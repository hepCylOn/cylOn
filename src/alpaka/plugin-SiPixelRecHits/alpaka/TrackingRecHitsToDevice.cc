#include <utility>
#include <alpaka/alpaka.hpp>
#include <optional>

#include "AlpakaCore/ScopedContext.h"
#include "AlpakaCore/config.h"
#include "AlpakaCore/memory.h"
#include "AlpakaCore/Product.h"

#include "AlpakaDataFormats/TrackingRecHitsHost.h"
#include "AlpakaDataFormats/alpaka/TrackingRecHitsSoACollection.h"
#include "AlpakaDataFormats/TrackingRecHitsDevice.h"

#include "Framework/EDProducer.h"
#include "Framework/Event.h"
#include "Framework/EventSetup.h"
#include "Framework/PluginFactory.h"

#define GPU_DEBUG

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  class TrackingRecHitsToDevice : public edm::EDProducerExternalWork {
  public:
    explicit TrackingRecHitsToDevice(edm::ProductRegistry& reg);
    ~TrackingRecHitsToDevice() override = default;

    using HitsOnDevice = reco::TrackingRecHitsSoACollection;
    using HitsOnHost = ::reco::TrackingRecHitHost;

  private:
    void acquire(const edm::Event& iEvent,
                 const edm::EventSetup& iSetup,
                 edm::WaitingTaskWithArenaHolder waitingTaskHolder) override;
    void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override;
    
    cms::alpakatools::ContextState<Queue> ctxState_;
        
    const edm::EDGetTokenT<HitsOnHost> input_;
    const edm::EDPutTokenT<cms::alpakatools::Product<Queue, HitsOnDevice>> output_;

    std::optional<HitsOnDevice> dHits_;

  };

  TrackingRecHitsToDevice::TrackingRecHitsToDevice(edm::ProductRegistry& reg)
      : input_(reg.consumes<HitsOnHost>()),
        output_(reg.produces<cms::alpakatools::Product<Queue, HitsOnDevice>>()) {}

  
  void TrackingRecHitsToDevice::acquire(const edm::Event& iEvent,
                 const edm::EventSetup& iSetup,
                 edm::WaitingTaskWithArenaHolder waitingTaskHolder) {

#ifdef GPU_DEBUG
    std::cout << "[TrackingRecHitsToDevice::GPU_DEBUG] acquire started " << std::endl;
#endif
    
   auto const& hHits = iEvent.get(input_);
   cms::alpakatools::ScopedContextAcquire<Queue> ctx{iEvent.streamID(), std::move(waitingTaskHolder), ctxState_};

#if defined ALPAKA_ACC_CPU_B_SEQ_T_SEQ_ENABLED or defined ALPAKA_ACC_CPU_B_TBB_T_SEQ_ENABLED
  return; //dHits_.emplace(std::move(hHits));
#else
  dHits_.emplace(
      cms::alpakatools::CopyToDevice<HitsOnHost>::copyAsync(ctx.stream(), hHits)
  );
  alpaka::wait(ctx.stream());
#endif

    }
  void TrackingRecHitsToDevice::produce(edm::Event& iEvent, edm::EventSetup const& iSetup) {

    cms::alpakatools::ScopedContextProduce ctx{ctxState_};
    ctx.emplace(iEvent, output_, std::move(*dHits_));

#ifdef GPU_DEBUG
    std::cout << "[TrackingRecHitsToDevice::GPU_DEBUG] Finished successfully.\n" << std::endl;
#endif
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

DEFINE_FWK_ALPAKA_MODULE(TrackingRecHitsToDevice);
