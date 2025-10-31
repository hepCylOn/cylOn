#include "AlpakaCore/ScopedContext.h"
#include "AlpakaCore/config.h"
#include "AlpakaDataFormats/alpaka/PixelTrackAlpaka.h"
#include "AlpakaDataFormats/alpaka/ZVertexAlpaka.h"
#include "Framework/EDProducer.h"
#include "Framework/Event.h"
#include "Framework/EventSetup.h"
#include "Framework/PluginFactory.h"
#include "Framework/RunningAverage.h"

#include "gpuVertexFinder.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  class PixelVertexProducerAlpaka : public edm::EDProducer {
  public:
    explicit PixelVertexProducerAlpaka(edm::ProductRegistry& reg);
    ~PixelVertexProducerAlpaka() override = default;

  private:
    void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override;

    edm::EDGetTokenT<cms::alpakatools::Product<Queue, PixelTrackAlpaka>> tokenTrack_;
    edm::EDPutTokenT<cms::alpakatools::Product<Queue, ZVertexAlpaka>> tokenVertex_;

    const gpuVertexFinder::Producer m_gpuAlgo;

    // Tracking cuts before sending tracks to vertex algo
    const float m_ptMin;
  };

  PixelVertexProducerAlpaka::PixelVertexProducerAlpaka(edm::ProductRegistry& reg)
      : tokenTrack_(reg.consumes<cms::alpakatools::Product<Queue, PixelTrackAlpaka>>()),
        tokenVertex_(reg.produces<cms::alpakatools::Product<Queue, ZVertexAlpaka>>()),
        m_gpuAlgo(true,   // oneKernel
                  true,   // useDensity
                  false,  // useDBSCAN
                  false,  // useIterative
                  2,      // minT
                  0.07,   // eps
                  0.01,   // errmax
                  // 11.*0.01,   // errmax (for colliderML)
                  9       // chi2max
                  ),
        m_ptMin(0.5)  // 0.5 GeV
  {}

  void PixelVertexProducerAlpaka::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    cms::alpakatools::Product<Queue, PixelTrackAlpaka> const& tracksWrapped = iEvent.get(tokenTrack_);
    cms::alpakatools::ScopedContextProduce<Queue> ctx{tracksWrapped};
    auto const& tracks = ctx.get(tracksWrapped);

    // for(uint32_t i = 0; i < tracks.data()->partIndices.nbins(); ++i){
    //   if(tracks.data()->nHits(i) < 5) continue;
    //   if(tracks.data()->chi2(i) > 30) continue;
    //   if(tracks.data()->pt(i) < 0.3) continue;
    //   std::cout << "track index: " << i << std::endl;
    //   std::cout << "tracks.data()->nHits(" << i << "): " << tracks.data()->nHits(i) << std::endl;
    //   std::cout << "tracks.data()->chi2(" << i << "): " << tracks.data()->chi2(i) << std::endl;
    //   std::cout << "tracks.data()->pt(" << i << "): " << tracks.data()->pt(i) << std::endl;
    //   std::cout << "tracks.data()->eta(" << i << "): " << tracks.data()->eta(i) << std::endl;
    //   std::cout << "tracks.data()->phi(" << i << "): " << tracks.data()->phi(i) << std::endl;
    //   std::cout << "tracks.data()->charge(" << i << "): " << tracks.data()->charge(i) << std::endl;
    //   std::cout << "tracks.data()->tip(" << i << "): " << tracks.data()->tip(i) << std::endl;
    //   std::cout << "tracks.data()->zip(" << i << "): " << tracks.data()->zip(i) << std::endl;
    //   std::cout << "tracks.data()->detIndices->bins[" << i << "]: ";
    //   for(auto j = tracks.data()->detIndices.begin(i); j < tracks.data()->detIndices.end(i); ++j) std::cout << tracks.data()->detIndices.begin(i)[*j] << ", ";
    //   std::cout << std::endl;
    //   std::cout << "tracks.data()->partIndices->bins[" << i << "]: ";
    //   for(auto j = tracks.data()->partIndices.begin(i); j < tracks.data()->partIndices.end(i); ++j) std::cout << tracks.data()->partIndices.begin(i)[*j] << ", ";
    //   std::cout << std::endl;
    //   std::cout << "====================================" << std::endl;
    // }

    ctx.emplace(iEvent, tokenVertex_, m_gpuAlgo.makeAsync(tracks.data(), m_ptMin, ctx.stream()));
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

DEFINE_FWK_ALPAKA_MODULE(PixelVertexProducerAlpaka);
