#include <fstream>
#include <ios>
#include <memory>
#include <string>
#include <vector>
#include <filesystem>
#include <iostream>   // for debug output
#include <algorithm>  // for std::min

#include "AlpakaCore/config.h"
#include "CondFormats/SiPixelGainCalibrationForHLTHost.h"
#include "Framework/ESProducer.h"
#include "Framework/EventSetup.h"
#include "Framework/ESPluginFactory.h"

#define GPU_DEBUG 

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  class SiPixelGainCalibrationForHLTHostESProducer : public edm::ESProducer {
  public:
    explicit SiPixelGainCalibrationForHLTHostESProducer(std::filesystem::path const& datadir) : data_(datadir) {}
    void produce(edm::EventSetup& eventSetup);

  private:
    std::filesystem::path data_;
  };

  void SiPixelGainCalibrationForHLTHostESProducer::produce(edm::EventSetup& eventSetup) {
    auto filepath = data_ / "SiPixelGainCalibrationForHLTHostRun2.bin"; // "SiPixelGainCalibrationForHLTHost.bin";
    std::ifstream in(filepath, std::ios::binary);
    in.exceptions(std::ifstream::badbit | std::ifstream::failbit);

#ifdef GPU_DEBUG
    std::cout << "[GPU_DEBUG] Reading SiPixelGainCalibrationForHLTHost from: " << filepath << "\n";
#endif

    // 1. number of decoding structures
    unsigned int size = 0;
    in.read(reinterpret_cast<char*>(&size), sizeof(unsigned int));

#ifdef GPU_DEBUG
    std::cout << "[GPU_DEBUG] Decoding structures count: " << size << "\n";
#endif

    // 2. allocate SoA on host
    auto gains = std::make_unique<SiPixelGainCalibrationForHLTHost>(size, cms::alpakatools::host());
    auto view = gains->view();

    // 3. read column (ped/gain pairs)
    in.read(reinterpret_cast<char*>(view.v_pedestals().data()),
            sizeof(siPixelGainsSoA::DecodingStructure) * size);

    // 4. read module metadata (arrays, not pairs!)
    auto& modStarts = view.modStarts();  // std::array<uint32_t, N>
    auto& modEnds   = view.modEnds();    // std::array<uint32_t, N>
    auto& modCols   = view.modCols();    // std::array<int,      N>

    in.read(reinterpret_cast<char*>(modStarts.data()), sizeof(modStarts));
    in.read(reinterpret_cast<char*>(modEnds.data()),   sizeof(modEnds));
    in.read(reinterpret_cast<char*>(modCols.data()),   sizeof(modCols));

    // 5. read scalar calibration info
    float minPed, maxPed, minGain, maxGain, pedPrecision, gainPrecision;
    in.read(reinterpret_cast<char*>(&minPed),        sizeof(float));
    in.read(reinterpret_cast<char*>(&maxPed),        sizeof(float));
    in.read(reinterpret_cast<char*>(&minGain),       sizeof(float));
    in.read(reinterpret_cast<char*>(&maxGain),       sizeof(float));
    in.read(reinterpret_cast<char*>(&pedPrecision),  sizeof(float));
    in.read(reinterpret_cast<char*>(&gainPrecision), sizeof(float));

    unsigned int numberOfRowsAveragedOver, nBinsToUseForEncoding, deadFlag, noisyFlag;
    in.read(reinterpret_cast<char*>(&numberOfRowsAveragedOver), sizeof(unsigned int));
    in.read(reinterpret_cast<char*>(&nBinsToUseForEncoding),    sizeof(unsigned int));
    in.read(reinterpret_cast<char*>(&deadFlag),                 sizeof(unsigned int));
    in.read(reinterpret_cast<char*>(&noisyFlag),                sizeof(unsigned int));

    float link;
    in.read(reinterpret_cast<char*>(&link), sizeof(float));
    in.close();

#ifdef GPU_DEBUG
    std::cout << "[GPU_DEBUG] Calibration scalars:\n";
    std::cout << "  minPed=" << minPed << " maxPed=" << maxPed
              << " minGain=" << minGain << " maxGain=" << maxGain << "\n";
    std::cout << "  pedPrecision=" << pedPrecision << " gainPrecision=" << gainPrecision << "\n";
    std::cout << "  numberOfRowsAveragedOver=" << numberOfRowsAveragedOver
              << " nBinsToUseForEncoding=" << nBinsToUseForEncoding << "\n";
    std::cout << "  deadFlag=" << deadFlag << " noisyFlag=" << noisyFlag << "\n";
    std::cout << "  link=" << link << "\n";

    // print just first few modules from the arrays
    const std::size_t nModulesToPrint =
        std::min<std::size_t>(5, modStarts.size());  // phase1PixelTopology::numberOfModules is the full size

    std::cout << "[GPU_DEBUG] First " << nModulesToPrint << " module entries (modStarts / modEnds / modCols):\n";
    for (std::size_t i = 0; i < nModulesToPrint; ++i) {
      std::cout << "  mod " << i
                << ": start=" << modStarts[i]
                << " end="   << modEnds[i]
                << " cols="  << modCols[i]
                << "\n";
    }

    // sample few decoding structures
    unsigned int nDecPrint = std::min(size, 8u);
    std::cout << "[GPU_DEBUG] First " << nDecPrint << " DecodingStructure entries:\n";
    for (unsigned int i = 0; i < nDecPrint; ++i) {
      auto const& d = view.v_pedestals()[i];
      std::cout << "  [" << i << "] ped=" << static_cast<unsigned>(d.ped)
                << " gain=" << static_cast<unsigned>(d.gain) << "\n";
    }
#endif

    // 6. assign scalar fields into the view
    view.minPed()                  = minPed;
    view.maxPed()                  = maxPed;
    view.minGain()                 = minGain;
    view.maxGain()                 = maxGain;
    view.pedPrecision()            = pedPrecision;
    view.gainPrecision()           = gainPrecision;
    view.numberOfRowsAveragedOver()= numberOfRowsAveragedOver;
    view.nBinsToUseForEncoding()   = nBinsToUseForEncoding;
    view.deadFlag()                = deadFlag;
    view.noisyFlag()               = noisyFlag;
    view.link()                    = link;

    // the arrays (modStarts/modEnds/modCols) are already written in place since we read into view.*().data()

    eventSetup.put(std::move(gains));
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

DEFINE_FWK_ALPAKA_EVENTSETUP_MODULE(SiPixelGainCalibrationForHLTHostESProducer);
