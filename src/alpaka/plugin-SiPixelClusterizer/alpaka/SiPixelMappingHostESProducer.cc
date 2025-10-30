#include <fstream>
#include <ios>
#include <memory>
#include <string>
#include <vector>
#include <filesystem>

#include "AlpakaCore/config.h"
#include "CondFormats/SiPixelMappingHost.h"
#include "Framework/ESProducer.h"
#include "Framework/EventSetup.h"
#include "Framework/ESPluginFactory.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  class SiPixelMappingHostESProducer : public edm::ESProducer {
  public:
    explicit SiPixelMappingHostESProducer(std::filesystem::path const& datadir) : data_(datadir) {}
    void produce(edm::EventSetup& eventSetup);

  private:
    std::filesystem::path data_;
  };

  void SiPixelMappingHostESProducer::produce(edm::EventSetup& eventSetup) {
    std::ifstream in(data_ / "SiPixelMappingHost.bin", std::ios::binary);
    in.exceptions(std::ifstream::badbit | std::ifstream::failbit);

    unsigned int size;
    bool hasQuality;
    in.read(reinterpret_cast<char*>(&size), sizeof(unsigned int));
    in.read(reinterpret_cast<char*>(&hasQuality), sizeof(bool));

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

    eventSetup.put(std::move(mapping));
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

DEFINE_FWK_ALPAKA_EVENTSETUP_MODULE(SiPixelMappingHostESProducer);
