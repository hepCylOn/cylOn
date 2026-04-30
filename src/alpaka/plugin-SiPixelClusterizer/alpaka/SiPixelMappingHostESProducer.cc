#include <fstream>
#include <ios>
#include <memory>
#include <string>
#include <vector>
#include <filesystem>
#include <iostream>  // for debug output

#include "AlpakaCore/config.h"
#include "CondFormats/SiPixelMappingHost.h"
#include "Framework/ESProducer.h"
#include "Framework/EventSetup.h"
#include "Framework/ESPluginFactory.h"
#include "Framework/ConfigRegistry.h"
#include "Framework/StreamFileUtils.h"

// #define GPU_DEBUG

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  class SiPixelMappingHostESProducer : public edm::ESProducer {
  public:
    explicit SiPixelMappingHostESProducer(edm::Config const& cfg) : data_(static_cast<std::string>(cfg.value("data", defaultPath_))) {
    }
    void produce(edm::EventSetup& eventSetup);

  private:
    std::filesystem::path data_;
    std::filesystem::path defaultPath_ = "data/SiPixelMappingHost.bin";
  };

  void SiPixelMappingHostESProducer::produce(edm::EventSetup& eventSetup) {
    auto in = edm::utils::openInputFile(data_);
    in.exceptions(std::ifstream::badbit | std::ifstream::failbit);

    unsigned int size;
    bool hasQuality;
    in.read(reinterpret_cast<char*>(&size), sizeof(unsigned int));
    in.read(reinterpret_cast<char*>(&hasQuality), sizeof(bool));

#ifdef GPU_DEBUG
    std::cout << "[GPU_DEBUG] Reading SiPixelMappingHost.bin\n";
    std::cout << "[GPU_DEBUG] Size: " << size << "\n";
    std::cout << "[GPU_DEBUG] hasQuality: " << std::boolalpha << hasQuality << "\n";
#endif

    // allocate host SoA
    auto mapping = std::make_unique<SiPixelMappingHost>(size, cms::alpakatools::host());
    auto view = mapping->view();

    // read all arrays
    in.read(reinterpret_cast<char*>(view.fed().data()), sizeof(unsigned int) * size);
    in.read(reinterpret_cast<char*>(view.link().data()), sizeof(unsigned int) * size);
    in.read(reinterpret_cast<char*>(view.roc().data()), sizeof(unsigned int) * size);
    in.read(reinterpret_cast<char*>(view.rawId().data()), sizeof(unsigned int) * size);
    in.read(reinterpret_cast<char*>(view.rocInDet().data()), sizeof(unsigned int) * size);
    in.read(reinterpret_cast<char*>(view.moduleId().data()), sizeof(unsigned int) * size);
    in.read(reinterpret_cast<char*>(view.badRocs().data()), sizeof(uint8_t) * size);
    in.read(reinterpret_cast<char*>(view.modToUnpDefault().data()), sizeof(unsigned char) * size);
    in.close();

    view.hasQuality() = hasQuality;

#ifdef GPU_DEBUG
    // print a few sample entries (up to first 5)
    unsigned int nPrint = std::min(size, 5u);
    std::cout << "[GPU_DEBUG] First " << nPrint << " entries:\n";
    for (unsigned int i = 0; i < nPrint; ++i) {
      std::cout << "  [" << i << "] fed=" << view.fed()[i]
                << " link=" << view.link()[i]
                << " roc=" << view.roc()[i]
                << " rawId=" << view.rawId()[i]
                << " rocInDet=" << view.rocInDet()[i]
                << " moduleId=" << view.moduleId()[i]
                << " badRocs=" << static_cast<unsigned>(view.badRocs()[i])
                << " modToUnpDefault=" << static_cast<unsigned>(view.modToUnpDefault()[i])
                << "\n";
    }
#endif

    eventSetup.put(std::move(mapping));
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

DEFINE_FWK_ALPAKA_EVENTSETUP_MODULE(SiPixelMappingHostESProducer);
