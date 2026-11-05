#include <memory>
#include <string>
#include <iostream>  // for debug output

#include "AlpakaCore/config.h"
#include "CondFormats/alpaka/PixelCPEFast.h"
#include "Framework/ESPluginFactory.h"
#include "Framework/ESProducer.h"
#include "Framework/EventSetup.h"
#include "Geometry/SimplePixelTopology.h"
#include "Framework/ConfigRegistry.h"
#include "Framework/StreamFileUtils.h"

#define GPU_DEBUG 

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  template <typename TrackerTraits>
  class PixelCPEFastESProducer : public edm::ESProducer {
  public:
    explicit PixelCPEFastESProducer(edm::Config const& cfg) : data_(static_cast<std::string>(cfg.value("data", defaultPath_))) {}
    void produce(edm::EventSetup& eventSetup);

  private:
    std::string data_;
    std::string defaultPath_ = "data/PixelCPEFastRun2.bin";
  };

  template <typename TrackerTraits>
  void PixelCPEFastESProducer<TrackerTraits>::produce(edm::EventSetup& eventSetup) {

    try {
      auto cpeFast = std::make_unique<PixelCPEFast<TrackerTraits>>(data_.c_str());
      eventSetup.put(std::move(cpeFast));
    #ifdef GPU_DEBUG
      std::cout << "[GPU_DEBUG] Successfully constructed and stored PixelCPEFast<"
                << TrackerTraits::nameModifier << ">\n";
    #endif
    } catch (std::exception const& e) {
      std::cerr << "[ERROR] PixelCPEFast<" << TrackerTraits::nameModifier
                << "> construction failed: " << e.what() << std::endl;
      throw;
    }

#ifdef GPU_DEBUG
    std::cout << "[GPU_DEBUG] PixelCPEFast<" << TrackerTraits::nameModifier << "> added to EventSetup\n";
#endif
  }

  // using PixelCPEFastESProducerPhase1 = PixelCPEFastESProducer<pixelTopology::Phase1>;
  // using PixelCPEFastESProducerPhase2 = PixelCPEFastESProducer<pixelTopology::Phase2>;

  /// FIXME: These are needed to make these plugins visible when building the plugins.txt list
  /// see: src/alpaka/Makefile:204. This is a workaround but it works for the moment.
  class PixelCPEFastESProducerPhase1 : public PixelCPEFastESProducer<pixelTopology::Phase1> {
  public:
    using PixelCPEFastESProducer<pixelTopology::Phase1>::PixelCPEFastESProducer;
  };

  class PixelCPEFastESProducerPhase2 : public PixelCPEFastESProducer<pixelTopology::Phase2> {
  public:
    using PixelCPEFastESProducer<pixelTopology::Phase2>::PixelCPEFastESProducer;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

DEFINE_FWK_ALPAKA_EVENTSETUP_MODULE(PixelCPEFastESProducerPhase1);
DEFINE_FWK_ALPAKA_EVENTSETUP_MODULE(PixelCPEFastESProducerPhase2);
