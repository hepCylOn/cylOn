#include <memory>
#include <string>

#include "AlpakaCore/config.h"
#include "CondFormats/alpaka/CAGeometry.h"
#include "Framework/ESPluginFactory.h"
#include "Framework/ESProducer.h"
#include "Framework/EventSetup.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  class CAGeometryESProducerFromHits : public edm::ESProducer {
  public:
    explicit CAGeometryESProducerFromHits(std::string const &datadir /*, std::string const &geometry*/) : data_(datadir)/*, geometry_(geometry)*/ {}
    void produce(edm::EventSetup &eventSetup);

  private:
    std::string data_;
    // std::string geometry_;
  };

  void CAGeometryESProducerFromHits::produce(edm::EventSetup &eventSetup) {
    eventSetup.put(std::make_unique<CAGeometry>((data_ + "/CMSPhase1GeometryFromHits.bin").c_str()));
  }
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

DEFINE_FWK_ALPAKA_EVENTSETUP_MODULE(CAGeometryESProducerFromHits);
