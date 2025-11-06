#include <fstream>
#include <ios>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "AlpakaCore/config.h"
#include "CondFormats/SiPixelFedCablingMapGPU.h"
#include "CondFormats/alpaka/SiPixelFedCablingMapGPUWrapper.h"
#include "Framework/ESPluginFactory.h"
#include "Framework/ESProducer.h"
#include "Framework/EventSetup.h"
#include "Framework/ConfigRegistry.h"
#include "Framework/StreamFileUtils.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  class SiPixelFedCablingMapESProducer : public edm::ESProducer {
  public:
    explicit SiPixelFedCablingMapESProducer(edm::Config const& cfg) : data_(static_cast<std::string>(cfg.value("data", defaultPath_))) {}
    void produce(edm::EventSetup& eventSetup);

  private:
    std::filesystem::path data_;
    std::filesystem::path defaultPath_ = "data/cablingMap.bin";
  };

  void SiPixelFedCablingMapESProducer::produce(edm::EventSetup& eventSetup) {
    auto in = edm::utils::openInputFile(data_);

    in.exceptions(std::ifstream::badbit | std::ifstream::failbit | std::ifstream::eofbit);
    SiPixelFedCablingMapGPU obj;
    in.read(reinterpret_cast<char*>(&obj), sizeof(SiPixelFedCablingMapGPU));
    unsigned int modToUnpDefSize;
    in.read(reinterpret_cast<char*>(&modToUnpDefSize), sizeof(unsigned int));
    std::vector<unsigned char> modToUnpDefault(modToUnpDefSize);
    in.read(reinterpret_cast<char*>(modToUnpDefault.data()), modToUnpDefSize);
    eventSetup.put(std::make_unique<SiPixelFedCablingMapGPUWrapper>(obj, std::move(modToUnpDefault)));
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

DEFINE_FWK_ALPAKA_EVENTSETUP_MODULE(SiPixelFedCablingMapESProducer);
