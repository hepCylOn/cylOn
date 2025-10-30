#include <utility>

#include "AlpakaCore/Product.h"
#include "AlpakaCore/ScopedContext.h"
#include "AlpakaCore/memory.h"
#include "AlpakaDataFormats/alpaka/BeamSpotSoACollection.h"
#include "Framework/EDProducer.h"
#include "Framework/Event.h"
#include "Framework/EventSetup.h"
#include "Framework/PluginFactory.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  class BeamSpotToSoA : public edm::EDProducer {
  public:
    explicit BeamSpotToSoA(edm::ProductRegistry& reg);
    ~BeamSpotToSoA() override = default;

    void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override;

  private:
    const edm::EDPutTokenT<cms::alpakatools::Product<Queue, BeamSpotSoACollection>> bsPutToken_;

    cms::alpakatools::host_buffer<BeamSpotPOD> bsHost_;
  };

  BeamSpotToSoA::BeamSpotToSoA(edm::ProductRegistry& reg)
      : bsPutToken_{reg.produces<cms::alpakatools::Product<Queue, BeamSpotSoACollection>>()},
        bsHost_{cms::alpakatools::make_host_buffer<BeamSpotPOD, Platform>()} {}

  void BeamSpotToSoA::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    *bsHost_ = iSetup.get<BeamSpotPOD>();

    cms::alpakatools::ScopedContextProduce<Queue> ctx{iEvent.streamID()};

    BeamSpotSoACollection bsDevice(ctx.stream());
    alpaka::memcpy(ctx.stream(), bsDevice.buffer(), bsHost_);

    ctx.emplace(iEvent, bsPutToken_, std::move(bsDevice));
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

DEFINE_FWK_ALPAKA_MODULE(BeamSpotToSoA);
