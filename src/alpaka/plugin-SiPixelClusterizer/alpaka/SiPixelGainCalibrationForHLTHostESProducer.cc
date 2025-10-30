#include <fstream>
#include <ios>
#include <memory>
#include <string>
#include <vector>
#include <filesystem>

#include "AlpakaCore/config.h"
#include "CondFormats/SiPixelGainCalibrationForHLTHost.h"
#include "Framework/ESProducer.h"
#include "Framework/EventSetup.h"
#include "Framework/ESPluginFactory.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  class SiPixelGainCalibrationForHLTHostESProducer : public edm::ESProducer {
  public:
    explicit SiPixelGainCalibrationForHLTHostESProducer(std::filesystem::path const& datadir) : data_(datadir) {}
    void produce(edm::EventSetup& eventSetup);

  private:
    std::filesystem::path data_;
  };

  void SiPixelGainCalibrationForHLTHostESProducer::produce(edm::EventSetup& eventSetup) {
    std::ifstream in(data_ / "SiPixelGainCalibrationForHLTHost.bin", std::ios::binary);
    in.exceptions(std::ifstream::badbit | std::ifstream::failbit);

    unsigned int size;
    in.read(reinterpret_cast<char*>(&size), sizeof(unsigned int));

    auto gains = std::make_unique<SiPixelGainCalibrationForHLTHost>(size, cms::alpakatools::host());
    auto view = gains->view();

    in.read(reinterpret_cast<char*>(view.v_pedestals().data()),
            sizeof(siPixelGainsSoA::DecodingStructure) * size);

    siPixelGainsSoA::Ranges modStarts, modEnds;
    siPixelGainsSoA::Cols modCols;
    in.read(reinterpret_cast<char*>(&modStarts), sizeof(modStarts));
    in.read(reinterpret_cast<char*>(&modEnds), sizeof(modEnds));
    in.read(reinterpret_cast<char*>(&modCols), sizeof(modCols));

    float minPed, maxPed, minGain, maxGain, pedPrecision, gainPrecision;
    in.read(reinterpret_cast<char*>(&minPed), sizeof(float));
    in.read(reinterpret_cast<char*>(&maxPed), sizeof(float));
    in.read(reinterpret_cast<char*>(&minGain), sizeof(float));
    in.read(reinterpret_cast<char*>(&maxGain), sizeof(float));
    in.read(reinterpret_cast<char*>(&pedPrecision), sizeof(float));
    in.read(reinterpret_cast<char*>(&gainPrecision), sizeof(float));

    unsigned int numberOfRowsAveragedOver, nBinsToUseForEncoding, deadFlag, noisyFlag;
    in.read(reinterpret_cast<char*>(&numberOfRowsAveragedOver), sizeof(unsigned int));
    in.read(reinterpret_cast<char*>(&nBinsToUseForEncoding), sizeof(unsigned int));
    in.read(reinterpret_cast<char*>(&deadFlag), sizeof(unsigned int));
    in.read(reinterpret_cast<char*>(&noisyFlag), sizeof(unsigned int));

    float link;
    in.read(reinterpret_cast<char*>(&link), sizeof(float));
    in.close();

    // assign scalar fields
    view.modStarts() = modStarts;
    view.modEnds() = modEnds;
    view.modCols() = modCols;

    view.minPed() = minPed;
    view.maxPed() = maxPed;
    view.minGain() = minGain;
    view.maxGain() = maxGain;
    view.pedPrecision() = pedPrecision;
    view.gainPrecision() = gainPrecision;
    view.numberOfRowsAveragedOver() = numberOfRowsAveragedOver;
    view.nBinsToUseForEncoding() = nBinsToUseForEncoding;
    view.deadFlag() = deadFlag;
    view.noisyFlag() = noisyFlag;
    view.link() = link;

    eventSetup.put(std::move(gains));
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

DEFINE_FWK_ALPAKA_EVENTSETUP_MODULE(SiPixelGainCalibrationForHLTHostESProducer);
