#include <alpaka/alpaka.hpp>

#include "Geometry/SimplePixelTopology.h"
// #include "FWCore/Framework/interface/Frameworkfwd.h"
// #include "FWCore/Utilities/interface/StreamID.h"
// #include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
// #include "FWCore/ParameterSet/interface/ParameterSet.h"
// #include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
// #include "FWCore/Utilities/interface/InputTag.h"
// #include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDPutToken.h"
// #include "HeterogeneousCore/AlpakaCore/interface/alpaka/ESGetToken.h"
#include "AlpakaCore/config.h"
// #include "HeterogeneousCore/AlpakaCore/interface/alpaka/Event.h"
// #include "HeterogeneousCore/AlpakaCore/interface/alpaka/EventSetup.h"
// #include "HeterogeneousCore/AlpakaCore/interface/alpaka/global/EDProducer.h"

#include "AlpakaDataFormats/alpaka/TracksSoACollection.h"
#include "AlpakaDataFormats/TracksDevice.h"

#include "AlpakaDataFormats/ZVertexHost.h"
#include "AlpakaDataFormats/alpaka/ZVertexSoACollection.h"
#include "AlpakaDataFormats/ZVertexDevice.h"
// #include "AlpakaDataFormats/PixelVertexWorkSpaceLayout.h"
// #include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"

#include "AlpakaCore/ScopedContext.h"
#include "AlpakaCore/config.h"
// #include "AlpakaDataFormats/alpaka/PixelTrackAlpaka.h"
// #include "AlpakaDataFormats/alpaka/ZVertexAlpaka.h"
#include "Framework/EDProducer.h"
#include "Framework/Event.h"
#include "Framework/EventSetup.h"
#include "Framework/PluginFactory.h"
#include "Framework/RunningAverage.h"

#include "vertexFinder.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  using namespace cms::alpakatools;

  template <typename TrackerTraits>
  class PixelVertexProducerAlpaka : public edm::EDProducer {
    using TkSoADevice = reco::TracksSoACollection;
    using VtxSoADevice = ZVertexSoACollection;
    using Algo = vertexFinder::Producer<TrackerTraits>;

  public:
    explicit PixelVertexProducerAlpaka(edm::ProductRegistry& reg);
    ~PixelVertexProducerAlpaka() override = default;

    // static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
    void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override;

    const Algo algo_;

    // Maximum number of vertices that will be reconstructed
    const int maxVertices_;

    // Tracking cuts before sending tracks to vertex algo
    const float ptMin_;
    const float ptMax_;

    edm::EDGetTokenT<cms::alpakatools::Product<Queue, TkSoADevice>> tokenTrack_;
    edm::EDPutTokenT<cms::alpakatools::Product<Queue, VtxSoADevice>> tokenVertex_;

  };

  template <typename TrackerTraits>
  PixelVertexProducerAlpaka<TrackerTraits>::PixelVertexProducerAlpaka(edm::ProductRegistry& reg)
      : algo_(/* oneKernel   */ true,
              /* useDensity  */ true,
              /* useDBSCAN   */ false,
              /* useIterative*/ false,
              /* doSplitting */ true,
              /* minT        */ 2,
              /* eps         */ 0.07,
              /* errmax      */ 0.01,
              /* chi2max     */ 9.0),
        maxVertices_(256),
        ptMin_(0.5),
        ptMax_(75.0),
        tokenTrack_(reg.consumes<cms::alpakatools::Product<Queue, TkSoADevice>>()),
        tokenVertex_(reg.produces<cms::alpakatools::Product<Queue, VtxSoADevice>>()) {}

  // template <typename TrackerTraits>
  // void PixelVertexProducerAlpaka<TrackerTraits>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //   edm::ParameterSetDescription desc;

  //   // Only one of these three algos can be used at once.
  //   // Maybe this should become a Plugin Factory
  //   desc.add<bool>("oneKernel", true);
  //   desc.add<bool>("useDensity", true);
  //   desc.add<bool>("useDBSCAN", false);
  //   desc.add<bool>("useIterative", false);
  //   desc.add<bool>("doSplitting", true);

  //   desc.add<int>("minT", 2);          // min number of neighbours to be "core"
  //   desc.add<double>("eps", 0.07);     // max absolute distance to cluster
  //   desc.add<double>("errmax", 0.01);  // max error to be "seed"
  //   desc.add<double>("chi2max", 9.);   // max normalized distance to cluster

  //   desc.add<int>("maxVertices", 256);
  //   desc.add<double>("PtMin", 0.5);
  //   desc.add<double>("PtMax", 75.);
  //   desc.add<edm::InputTag>("pixelTrackSrc", edm::InputTag("pixelTracksAlpaka"));

  //   descriptions.addWithDefaultLabel(desc);
  // }

  template <typename TrackerTraits>
  void PixelVertexProducerAlpaka<TrackerTraits>::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    
    auto const& tracksWrapped = iEvent.get(tokenTrack_);
    cms::alpakatools::ScopedContextProduce<Queue> ctx{tracksWrapped};
    auto const& tracks = ctx.get(tracksWrapped);

    ctx.emplace(iEvent, tokenVertex_, algo_.makeAsync(ctx.stream(), tracks.view(), maxVertices_, ptMin_, ptMax_));
  }

  using PixelVertexProducerAlpakaPhase1 = PixelVertexProducerAlpaka<pixelTopology::Phase1>;
  using PixelVertexProducerAlpakaPhase2 = PixelVertexProducerAlpaka<pixelTopology::Phase2>;
  using PixelVertexProducerAlpakaHIonPhase1 = PixelVertexProducerAlpaka<pixelTopology::HIonPhase1>;

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

DEFINE_FWK_ALPAKA_MODULE(PixelVertexProducerAlpakaPhase1);
DEFINE_FWK_ALPAKA_MODULE(PixelVertexProducerAlpakaPhase2);
DEFINE_FWK_ALPAKA_MODULE(PixelVertexProducerAlpakaHIonPhase1);
