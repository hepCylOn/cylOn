#include "AlpakaCore/Product.h"
#include "AlpakaCore/ScopedContext.h"
#include "AlpakaCore/config.h"
#include "AlpakaDataFormats/alpaka/PixelTrackAlpaka.h"
#include "AlpakaDataFormats/alpaka/TrackingRecHit2DAlpaka.h"
#include "CAHitNtupletGeneratorOnGPU.h"
#include "Framework/EDProducer.h"
#include "Framework/Event.h"
#include "Framework/EventSetup.h"
#include "Framework/PluginFactory.h"
#include "CondFormats/alpaka/CAGeometry.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  class CAHitNtupletAlpaka : public edm::EDProducer {
  public:
    explicit CAHitNtupletAlpaka(edm::ProductRegistry& reg);
    ~CAHitNtupletAlpaka() override = default;

  private:
    void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override;

    edm::EDGetTokenT<cms::alpakatools::Product<Queue, TrackingRecHit2DAlpaka>> tokenHitGPU_;
    edm::EDPutTokenT<cms::alpakatools::Product<Queue, PixelTrackAlpaka>> tokenTrackGPU_;

    CAHitNtupletGeneratorOnGPU gpuAlgo_;
  };

  CAHitNtupletAlpaka::CAHitNtupletAlpaka(edm::ProductRegistry& reg)
      : tokenHitGPU_{reg.consumes<cms::alpakatools::Product<Queue, TrackingRecHit2DAlpaka>>()},
        tokenTrackGPU_{reg.produces<cms::alpakatools::Product<Queue, PixelTrackAlpaka>>()},
        gpuAlgo_(reg) {}

  void CAHitNtupletAlpaka::produce(edm::Event& iEvent, const edm::EventSetup& es) {
    auto bf = 0.0114256972711507;  // 1/fieldInGeV
    // auto bf = 0.0166990960116818;  // 1/fieldInGeV

    auto const& geo = es.get<CAGeometry>();

    auto const& phits = iEvent.get(tokenHitGPU_);
    cms::alpakatools::ScopedContextProduce<Queue> ctx{phits};
    auto const& hits = ctx.get(phits);

    // // Lines below are used to write input files to be used with fromHits;
    // // Remember to also uncomment lines 28-31 in CAHitNtupletGeneratorKernels.cc
    // auto const& hits_h = hits.view_h();
    // std::cout << "hits:" << hits_h->nHits() << std::endl;
    // for(uint32_t i = 0; i < hits_h->nHits(); ++i){
    //   // std::cout << hits_h->xLocal(i) << "," << hits_h->yLocal(i) << "," << hits_h->xerrLocal(i) << "," << hits_h->yerrLocal(i) << "," << hits_h->xGlobal(i) << "," << hits_h->yGlobal(i) << "," << hits_h->zGlobal(i) << "," << hits_h->rGlobal(i) << "," << hits_h->iphi(i) << "," << hits_h->charge(i) << "," << hits_h->clusterSizeX(i) << "," << hits_h->clusterSizeY(i) << "," << hits_h->detectorIndex(i) << "," << hits_h->particleIndex(i) << std::endl;
    //   std::cout << hits_h->xLocal(i) << "," << hits_h->yLocal(i) << "," << hits_h->xerrLocal(i) << "," << hits_h->yerrLocal(i) << "," << hits_h->xGlobal(i) << "," << hits_h->yGlobal(i) << "," << hits_h->zGlobal(i) << "," << hits_h->rGlobal(i) << "," << hits_h->iphi(i) << "," << hits_h->charge(i) << "," << hits_h->clusterSizeX(i) << "," << hits_h->clusterSizeY(i) << "," << hits_h->detectorIndex(i) << "," << i << std::endl;
    // }
    // std::cout << "module:10" << std::endl;

    ctx.emplace(iEvent, tokenTrackGPU_, gpuAlgo_.makeTuplesAsync(hits, bf, geo.getGPUProductAsync(ctx.stream()), geo.sizes(), ctx.stream()));
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

DEFINE_FWK_ALPAKA_MODULE(CAHitNtupletAlpaka);
