// #include <utility>

// #include <alpaka/alpaka.hpp>

// #include "AlpakaCore/ScopedContext.h"
// #include "AlpakaCore/config.h"
// #include "AlpakaCore/memory.h"
// #include "AlpakaDataFormats/alpaka/SiPixelRecHitAlpaka.h"
// #include "AlpakaDataFormats/SiPixelRecHitHost.h"
// #include "Framework/EDProducer.h"
// #include "Framework/Event.h"
// #include "Framework/EventSetup.h"
// #include "Framework/PluginFactory.h"

// namespace ALPAKA_ACCELERATOR_NAMESPACE {

//   class SiPixelRecHitSoAFromAlpaka : public edm::EDProducerExternalWork {
//   public:
//     explicit SiPixelRecHitSoAFromAlpaka(edm::ProductRegistry& reg, edm::Config const& cfg);
//     ~SiPixelRecHitSoAFromAlpaka() override = default;

//   private:
//     void acquire(edm::Event const& iEvent,
//                  edm::EventSetup const& iSetup,
//                  edm::WaitingTaskWithArenaHolder waitingTaskHolder) override;
//     void produce(edm::Event& iEvent, edm::EventSetup const& iSetup) override;

//     edm::EDGetTokenT<cms::alpakatools::Product<Queue, SiPixelRecHitAlpaka>> tokenDevice_;
//     edm::EDPutTokenT<SiPixelRecHitHost> tokenHost_;

//     cms::alpakatools::host_buffer<SiPixelRecHit::TrackSoA> soa_;
//   };

//   SiPixelRecHitSoAFromAlpaka::SiPixelRecHitSoAFromAlpaka(edm::ProductRegistry& reg, edm::Config const& cfg)
//       : tokenDevice_(reg.consumes<cms::alpakatools::Product<Queue, SiPixelRecHitAlpaka>>()),
//         tokenHost_(reg.produces<SiPixelRecHitHost>()),
//         soa_{cms::alpakatools::make_host_buffer<SiPixelRecHit::TrackSoA, Platform>()} {}

//   void SiPixelRecHitSoAFromAlpaka::acquire(edm::Event const& iEvent,
//                                         edm::EventSetup const& iSetup,
//                                         edm::WaitingTaskWithArenaHolder waitingTaskHolder) {
//     cms::alpakatools::Product<Queue, SiPixelRecHitAlpaka> const& inputDataWrapped = iEvent.get(tokenDevice_);
//     cms::alpakatools::ScopedContextAcquire ctx{inputDataWrapped, std::move(waitingTaskHolder)};
//     auto const& inputData = ctx.get(inputDataWrapped);
//     soa_ = cms::alpakatools::make_host_buffer<SiPixelRecHit::TrackSoA>(ctx.stream());
//     alpaka::memcpy(ctx.stream(), soa_, inputData);
//   }

//   void SiPixelRecHitSoAFromAlpaka::produce(edm::Event& iEvent, edm::EventSetup const& iSetup) {
//     iEvent.emplace(tokenHost_, std::move(soa_));
//   }

// }  // namespace ALPAKA_ACCELERATOR_NAMESPACE

// DEFINE_FWK_ALPAKA_MODULE(SiPixelRecHitSoAFromAlpaka);
