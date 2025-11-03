#include <fstream>
#include <ios>
#include <memory>
#include <string>
#include <vector>
#include <filesystem>
#include <iostream>

#include "AlpakaCore/config.h"
#include "CondFormats/alpaka/SiPixelGainForHLTonGPU.h"
#include "CondFormats/SiPixelGainCalibrationForHLTHost.h"
#include "Framework/ESProducer.h"
#include "Framework/EventSetup.h"
#include "Framework/ESPluginFactory.h"

//#define GPU_DEBUG

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  class SiPixelGainCalibrationForHLTHostFromGPUBinESProducer : public edm::ESProducer {
  public:
    explicit SiPixelGainCalibrationForHLTHostFromGPUBinESProducer(std::filesystem::path const& datadir)
        : data_(datadir) {}
    void produce(edm::EventSetup& eventSetup);

  private:
    std::filesystem::path data_;
  };

  void SiPixelGainCalibrationForHLTHostFromGPUBinESProducer::produce(edm::EventSetup& eventSetup) {
    using GPUFormat = SiPixelGainForHLTonGPU;
    using DecodingStructure = GPUFormat::DecodingStructure;

    auto inputFile = data_ / "gain.bin";
    auto outputFile = data_ / "SiPixelGainCalibrationForHLTHostRun2.bin";

    std::ifstream in(inputFile, std::ios::binary);
    in.exceptions(std::ifstream::badbit | std::ifstream::failbit | std::ifstream::eofbit);

#ifdef GPU_DEBUG
    std::cout << "[GPU_DEBUG] Reading old GPU-format gain from: " << inputFile << "\n";
#endif

    // --- read GPU object header
    GPUFormat gainGPU;
    in.read(reinterpret_cast<char*>(&gainGPU), sizeof(GPUFormat));

    unsigned int nbytes = 0;
    in.read(reinterpret_cast<char*>(&nbytes), sizeof(unsigned int));
    std::vector<DecodingStructure> gainData(nbytes / sizeof(DecodingStructure));
    in.read(reinterpret_cast<char*>(gainData.data()), nbytes);
    in.close();

    unsigned int size = gainData.size();

#ifdef GPU_DEBUG
    std::cout << "[GPU_DEBUG] gainData entries: " << size
              << " (" << nbytes << " bytes)\n";
    std::cout << "[GPU_DEBUG] Calibration scalars:\n";
    std::cout << "  minPed=" << gainGPU.minPed_ << " maxPed=" << gainGPU.maxPed_
              << " minGain=" << gainGPU.minGain_ << " maxGain=" << gainGPU.maxGain_ << "\n";
    std::cout << "  pedPrecision=" << gainGPU.pedPrecision_
              << " gainPrecision=" << gainGPU.gainPrecision_ << "\n";
    std::cout << "  nRowsAvg=" << gainGPU.numberOfRowsAveragedOver_
              << " nBinsToUse=" << gainGPU.nBinsToUseForEncoding_
              << " deadFlag=" << gainGPU.deadFlag_
              << " noisyFlag=" << gainGPU.noisyFlag_ << "\n";
#endif

    // --- allocate the new host SoA
    auto gains = std::make_unique<SiPixelGainCalibrationForHLTHost>(size, cms::alpakatools::host());
    auto view = gains->view();

    // fill decoding data (manual copy since struct types differ)
    for (size_t i = 0; i < size; ++i) {
    view[i].v_pedestals().gain = gainData[i].gain;
    view[i].v_pedestals().ped  = gainData[i].ped;
    }

    // fill arrays (modStarts, modEnds, modCols) from the rangeAndCols_ in the GPU object
    const size_t nModules = phase1PixelTopology::numberOfModules;
    auto& modStarts = view.modStarts();
    auto& modEnds   = view.modEnds();
    auto& modCols   = view.modCols();

    for (size_t i = 0; i < nModules; ++i) {
      modStarts[i] = gainGPU.rangeAndCols_[i].first.first;
      modEnds[i]   = gainGPU.rangeAndCols_[i].first.second;
      modCols[i]   = gainGPU.rangeAndCols_[i].second;
    }

    // assign calibration scalars
    view.minPed() = gainGPU.minPed_;
    view.maxPed() = gainGPU.maxPed_;
    view.minGain() = gainGPU.minGain_;
    view.maxGain() = gainGPU.maxGain_;
    view.pedPrecision() = gainGPU.pedPrecision_;
    view.gainPrecision() = gainGPU.gainPrecision_;
    view.numberOfRowsAveragedOver() = gainGPU.numberOfRowsAveragedOver_;
    view.nBinsToUseForEncoding() = gainGPU.nBinsToUseForEncoding_;
    view.deadFlag() = gainGPU.deadFlag_;
    view.noisyFlag() = gainGPU.noisyFlag_;
    view.link() = 0.0f;  // not stored in old format, set to default

    // --- Write out the converted host format to a new binary file
    std::ofstream out(outputFile, std::ios::binary);
    out.exceptions(std::ofstream::badbit | std::ofstream::failbit);
    out.write(reinterpret_cast<const char*>(&size), sizeof(unsigned int));
    out.write(reinterpret_cast<const char*>(view.v_pedestals().data()),
              sizeof(siPixelGainsSoA::DecodingStructure) * size);
    out.write(reinterpret_cast<const char*>(modStarts.data()), sizeof(modStarts));
    out.write(reinterpret_cast<const char*>(modEnds.data()), sizeof(modEnds));
    out.write(reinterpret_cast<const char*>(modCols.data()), sizeof(modCols));

    out.write(reinterpret_cast<const char*>(&gainGPU.minPed_), sizeof(float));
    out.write(reinterpret_cast<const char*>(&gainGPU.maxPed_), sizeof(float));
    out.write(reinterpret_cast<const char*>(&gainGPU.minGain_), sizeof(float));
    out.write(reinterpret_cast<const char*>(&gainGPU.maxGain_), sizeof(float));
    out.write(reinterpret_cast<const char*>(&gainGPU.pedPrecision_), sizeof(float));
    out.write(reinterpret_cast<const char*>(&gainGPU.gainPrecision_), sizeof(float));
    out.write(reinterpret_cast<const char*>(&gainGPU.numberOfRowsAveragedOver_), sizeof(unsigned int));
    out.write(reinterpret_cast<const char*>(&gainGPU.nBinsToUseForEncoding_), sizeof(unsigned int));
    out.write(reinterpret_cast<const char*>(&gainGPU.deadFlag_), sizeof(unsigned int));
    out.write(reinterpret_cast<const char*>(&gainGPU.noisyFlag_), sizeof(unsigned int));
    float linkDummy = 0.0f;
    out.write(reinterpret_cast<const char*>(&linkDummy), sizeof(float));
    out.close();

#ifdef GPU_DEBUG
    std::cout << "[GPU_DEBUG] Wrote converted host-format gain to: " << outputFile << "\n";
#endif

    // --- put the new object into the EventSetup
    eventSetup.put(std::move(gains));
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

DEFINE_FWK_ALPAKA_EVENTSETUP_MODULE(SiPixelGainCalibrationForHLTHostFromGPUBinESProducer);
