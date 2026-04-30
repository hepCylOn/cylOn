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

// #define GPU_DEBUG

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  class TrackingRecHitsToDevice : public edm::EDProducerExternalWork {
  public:
    explicit TrackingRecHitsToDevice(edm::ProductRegistry& reg, edm::Config const& cfg);
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

  TrackingRecHitsToDevice::TrackingRecHitsToDevice(edm::ProductRegistry& reg, edm::Config const& cfg)
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
  //  std::cout << __LINE__ << " -- " << __FILE__ << " -- AAAAAAAAAAAAAAHHHHHHHHHHHHHH" << std::endl;
  //  std::cout << "hHits.view()[0].xLocal(): " << hHits.view()[0].xLocal() << std::endl;
  //  std::cout << "hHits.view()[0].yLocal(): " << hHits.view()[0].yLocal() << std::endl;
  //  std::cout << "hHits.view()[0].xerrLocal(): " << hHits.view()[0].xerrLocal() << std::endl;
  //  std::cout << "hHits.view()[0].yerrLocal(): " << hHits.view()[0].yerrLocal() << std::endl;
  //  std::cout << "hHits.view()[0].xGlobal(): " << hHits.view()[0].xGlobal() << std::endl;
  //  std::cout << "hHits.view()[0].yGlobal(): " << hHits.view()[0].yGlobal() << std::endl;
  //  std::cout << "hHits.view()[0].zGlobal(): " << hHits.view()[0].zGlobal() << std::endl;
  //  std::cout << "hHits.view()[0].rGlobal(): " << hHits.view()[0].rGlobal() << std::endl;
  //  std::cout << "hHits.view()[0].iphi(): " << hHits.view()[0].iphi() << std::endl;
  //  std::cout << "hHits.view()[0].chargeAndStatus().charge: " << hHits.view()[0].chargeAndStatus().charge << std::endl;
  //  std::cout << "hHits.view()[0].clusterSizeX(): " << hHits.view()[0].clusterSizeX() << std::endl;
  //  std::cout << "hHits.view()[0].clusterSizeY(): " << hHits.view()[0].clusterSizeY() << std::endl;
  //  std::cout << "hHits.view()[0].detectorIndex(): " << hHits.view()[0].detectorIndex() << std::endl;
  //  std::cout << "hHits.view().offsetBPIX2(): " << hHits.view().offsetBPIX2() << std::endl;
  //  std::cout << "hHits.view<::reco::HitModuleSoA>()[0].moduleStart(): " << hHits.view<::reco::HitModuleSoA>()[0].moduleStart() << std::endl;
  //  std::cout << "hHits.view<::reco::HitModuleSoA>()[1].moduleStart(): " << hHits.view<::reco::HitModuleSoA>()[1].moduleStart() << std::endl;
  //  std::cout << "hHits.view<::reco::HitModuleSoA>()[2].moduleStart(): " << hHits.view<::reco::HitModuleSoA>()[2].moduleStart() << std::endl;
  //  std::cout << "hHits.view<::reco::HitModuleSoA>()[3].moduleStart(): " << hHits.view<::reco::HitModuleSoA>()[3].moduleStart() << std::endl;
  //  std::cout << "hHits.view<::reco::HitModuleSoA>()[4].moduleStart(): " << hHits.view<::reco::HitModuleSoA>()[4].moduleStart() << std::endl;
  //  std::cout << "hHits.view<::reco::HitModuleSoA>()[5].moduleStart(): " << hHits.view<::reco::HitModuleSoA>()[5].moduleStart() << std::endl;
  //  std::cout << "hHits.view<::reco::HitModuleSoA>()[6].moduleStart(): " << hHits.view<::reco::HitModuleSoA>()[6].moduleStart() << std::endl;
  //  std::cout << "hHits.view<::reco::HitModuleSoA>()[7].moduleStart(): " << hHits.view<::reco::HitModuleSoA>()[7].moduleStart() << std::endl;
  //  std::cout << "hHits.view<::reco::HitModuleSoA>()[8].moduleStart(): " << hHits.view<::reco::HitModuleSoA>()[8].moduleStart() << std::endl;
  //  std::cout << "hHits.view<::reco::HitModuleSoA>()[9].moduleStart(): " << hHits.view<::reco::HitModuleSoA>()[9].moduleStart() << std::endl;
  //  std::cout << "hHits.view<::reco::HitModuleSoA>()[10].moduleStart(): " << hHits.view<::reco::HitModuleSoA>()[10].moduleStart() << std::endl;
  //  std::cout << __LINE__ << " -- " << __FILE__ << " -- AAAAAAAAAAAAAAHHHHHHHHHHHHHH" << std::endl;

#if defined ALPAKA_ACC_CPU_B_SEQ_T_SEQ_ENABLED or defined ALPAKA_ACC_CPU_B_TBB_T_SEQ_ENABLED
  // return; //dHits_.emplace(std::move(hHits));
  // dHits_.emplace(std::move(hHits));
  dHits_.emplace(
      cms::alpakatools::CopyToHost<HitsOnHost>::copyAsync(ctx.stream(), hHits)
  );
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
