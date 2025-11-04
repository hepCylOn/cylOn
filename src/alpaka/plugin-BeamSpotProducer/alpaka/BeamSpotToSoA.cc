#include <utility>
#include <iostream>  // for std::cout, std::endl

#include "AlpakaCore/Product.h"
#include "AlpakaCore/ScopedContext.h"
#include "AlpakaCore/memory.h"
#include "AlpakaDataFormats/alpaka/BeamSpotSoACollection.h"
#include "Framework/EDProducer.h"
#include "Framework/Event.h"
#include "Framework/EventSetup.h"
#include "Framework/PluginFactory.h"

#define GPU_DEBUG

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
        bsHost_{cms::alpakatools::make_host_buffer<BeamSpotPOD, Platform>()} {
#ifdef GPU_DEBUG
    std::cout << "[GPU_DEBUG] BeamSpotToSoA: constructor called, host buffer allocated at "
              << static_cast<const void*>(bsHost_.data()) << std::endl;
#endif
  }

  void BeamSpotToSoA::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
#ifdef GPU_DEBUG
    std::cout << "[GPU_DEBUG] BeamSpotToSoA::produce start for stream "
              << iEvent.streamID() << std::endl;
#endif

    *bsHost_ = iSetup.get<BeamSpotPOD>();

#ifdef GPU_DEBUG
    std::cout << "[GPU_DEBUG] BeamSpotToSoA: BeamSpot read from EventSetup -> "
              << "x=" << bsHost_->x << "  y=" << bsHost_->y << "  z=" << bsHost_->z << " sigmaZ = " << bsHost_->sigmaZ
              << std::endl;
#endif

    cms::alpakatools::ScopedContextProduce<Queue> ctx{iEvent.streamID()};

#ifdef GPU_DEBUG
    std::cout << "[GPU_DEBUG] BeamSpotToSoA: ScopedContextProduce created for stream "
              << iEvent.streamID() << std::endl;
#endif

    // FIXME: have the copy only if needed
    BeamSpotSoACollection bsDevice(ctx.stream());
    alpaka::memcpy(ctx.stream(), bsDevice.buffer(), bsHost_);

#ifdef GPU_DEBUG
    std::cout << "[GPU_DEBUG] BeamSpotToSoA: memcpy to device done (device buffer "
              << static_cast<const void*>(bsDevice.buffer().data()) << ")"
              << std::endl;
#endif

    ctx.emplace(iEvent, bsPutToken_, std::move(bsDevice));

#ifdef GPU_DEBUG
    std::cout << "[GPU_DEBUG] BeamSpotToSoA::produce completed successfully" << std::endl;
#endif
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

DEFINE_FWK_ALPAKA_MODULE(BeamSpotToSoA);
