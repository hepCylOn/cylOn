#include <memory>
#include <string>

#include "AlpakaCore/config.h"
#include "CondFormats/alpaka/PixelCPEFast.h"
#include "Framework/ESPluginFactory.h"
#include "Framework/ESProducer.h"
#include "Framework/EventSetup.h"

#include "Geometry/SimplePixelTopology.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  template <typename TrackerTraits>
  class PixelCPEFastESProducer : public edm::ESProducer {
  public:
    explicit PixelCPEFastESProducer(std::string const &datadir) : data_(datadir) {}
    void produce(edm::EventSetup &eventSetup);

  private:
    std::string data_;
  };

  template <typename TrackerTraits>
  void PixelCPEFastESProducer<TrackerTraits>::produce(edm::EventSetup &eventSetup) {
    eventSetup.put(std::make_unique<PixelCPEFast<TrackerTraits>>((data_ + "/cpefast.bin").c_str()));
  }
  
  using PixelCPEFastESProducerPhase1 = PixelCPEFastESProducer<pixelTopology::Phase1>;
  using PixelCPEFastESProducerPhase2 = PixelCPEFastESProducer<pixelTopology::Phase2>;

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE


DEFINE_FWK_ALPAKA_EVENTSETUP_MODULE(PixelCPEFastESProducerPhase1);
DEFINE_FWK_ALPAKA_EVENTSETUP_MODULE(PixelCPEFastESProducerPhase2);
