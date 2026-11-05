#include "AlpakaDataFormats/BeamSpotPOD.h"
#include "AlpakaDataFormats/alpaka/BeamSpotSoACollection.h"
#include "AlpakaDataFormats/SiPixelClustersDevice.h"
#include "AlpakaDataFormats/alpaka/SiPixelClustersSoACollection.h"
#include "AlpakaDataFormats/SiPixelDigisDevice.h"
#include "AlpakaDataFormats/alpaka/SiPixelDigisSoACollection.h"
#include "AlpakaDataFormats/TrackingRecHitsDevice.h"
#include "AlpakaDataFormats/alpaka/TrackingRecHitsSoACollection.h"
#include "Framework/EDProducer.h"
#include "Framework/Event.h"
#include "Framework/EventSetup.h"
#include "Framework/PluginFactory.h"
// #include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
// #include "FWCore/ParameterSet/interface/ParameterSet.h"
// #include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
// #include "FWCore/Utilities/interface/InputTag.h"
#include "Geometry/SimplePixelTopology.h"
// #include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
// #include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
// #include "HeterogeneousCore/AlpakaCore/interface/alpaka/Event.h"
// #include "HeterogeneousCore/AlpakaCore/interface/alpaka/EventSetup.h"
// #include "HeterogeneousCore/AlpakaCore/interface/alpaka/global/EDProducer.h"
#include "AlpakaCore/config.h"
#include "AlpakaCore/memory.h"
#include "AlpakaCore/Product.h"
#include "AlpakaCore/ScopedContext.h"

// #include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
// #include "FWCore/ParameterSet/interface/ParameterSet.h"
// #include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
// #include "FWCore/Utilities/interface/InputTag.h"
// #include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
// #include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
// #include "RecoLocalTracker/Records/interface/PixelCPEFastParamsRecord.h"

// #include "RecoLocalTracker/SiPixelRecHits/interface/PixelCPEBase.h"
#include "CondFormats/pixelCPEforDevice.h"
#include "CondFormats/alpaka/PixelCPEFast.h"

#include "PixelRecHitKernel.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  template <typename TrackerTraits>
  class SiPixelRecHit : public edm::EDProducer {
  public:
    explicit SiPixelRecHit(edm::ProductRegistry& reg, edm::Config const& cfg);
    ~SiPixelRecHit() override = default;

    // static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
     void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override;

    // const edm::ESGetToken<PixelCPEFastParams<TrackerTraits>, PixelCPEFastParamsRecord> cpeToken_;
    const edm::EDGetTokenT<cms::alpakatools::Product<Queue, BeamSpotSoACollection>> tBeamSpot;
    const edm::EDGetTokenT<cms::alpakatools::Product<Queue, SiPixelClustersSoACollection>> tokenClusters_;
    const edm::EDGetTokenT<cms::alpakatools::Product<Queue, SiPixelDigisSoACollection>> tokenDigi_;
    const edm::EDPutTokenT<cms::alpakatools::Product<Queue, reco::TrackingRecHitsSoACollection>> tokenHit_;

    pixelgpudetails::PixelRecHitKernel<TrackerTraits> Algo_;
  };

  template <typename TrackerTraits>
  SiPixelRecHit<TrackerTraits>::SiPixelRecHit(edm::ProductRegistry& reg, edm::Config const& cfg)
      : tBeamSpot(reg.consumes<cms::alpakatools::Product<Queue, BeamSpotSoACollection>>()),
        tokenClusters_(reg.consumes<cms::alpakatools::Product<Queue, SiPixelClustersSoACollection>>()),
        tokenDigi_(reg.consumes<cms::alpakatools::Product<Queue, SiPixelDigisSoACollection>>()),
        tokenHit_(reg.produces<cms::alpakatools::Product<Queue, reco::TrackingRecHitsSoACollection>>()) 
        {}

  // template <typename TrackerTraits>
  // void SiPixelRecHit<TrackerTraits>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //   edm::ParameterSetDescription desc;

  //   desc.add<edm::InputTag>("beamSpot", edm::InputTag("offlineBeamSpotSoACollection"));
  //   desc.add<edm::InputTag>("src", edm::InputTag("siPixelClustersPreSplittingAlpaka"));

  //   std::string cpe = "PixelCPEFastParams";
  //   cpe += TrackerTraits::nameModifier;
  //   desc.add<std::string>("CPE", cpe);

  //   descriptions.addWithDefaultLabel(desc);
  // }

  template <typename TrackerTraits>
  void SiPixelRecHit<TrackerTraits>::produce(edm::Event& iEvent, const edm::EventSetup& es){
                                                    
    auto const& fcpe = es.get<PixelCPEFast<TrackerTraits>>();

    auto const& pclusters = iEvent.get(tokenClusters_);
    cms::alpakatools::ScopedContextProduce<Queue> ctx{pclusters};

    auto const& clusters = ctx.get(pclusters);
    auto const& digis = ctx.get(iEvent, tokenDigi_);
    auto const& bs = ctx.get(iEvent, tBeamSpot);

#ifdef GPU_DEBUG
    std::cout << "[SiPixelRecHit::GPU_DEBUG<" << TrackerTraits::nameModifier
              << ">] Producing rechits with "
              << clusters.view().metadata().size() << " clusters and "
              << digis.view().metadata().size() << " digis" << std::endl;
#endif


    ctx.emplace(iEvent,
                tokenHit_,
                Algo_.makeHitsAsync(digis, clusters, bs.data(), fcpe.getGPUProductAsync(ctx.stream()), ctx.stream()));

#ifdef GPU_DEBUG
    std::cout << "[SiPixelRecHit::GPU_DEBUG<" << TrackerTraits::nameModifier
              << ">] Rechit production complete" << std::endl;
#endif

  }
  // using SiPixelRecHitPhase1 = SiPixelRecHit<pixelTopology::Phase1>;
  // using SiPixelRecHitHIonPhase1 = SiPixelRecHit<pixelTopology::HIonPhase1>;
  // using SiPixelRecHitPhase2 = SiPixelRecHit<pixelTopology::Phase2>;

  /// FIXME: These are needed to make these plugins visible when building the plugins.txt list
  /// see: src/alpaka/Makefile:204. This is a workaround but it works for the moment.
  class SiPixelRecHitPhase1 : public SiPixelRecHit<pixelTopology::Phase1> {
  public:
    using SiPixelRecHit<pixelTopology::Phase1>::SiPixelRecHit;
  };

  class SiPixelRecHitPhase2 : public SiPixelRecHit<pixelTopology::Phase2> {
  public:
    using SiPixelRecHit<pixelTopology::Phase2>::SiPixelRecHit;
  };

  class SiPixelRecHitHIonPhase1 : public SiPixelRecHit<pixelTopology::HIonPhase1> {
  public:
    using SiPixelRecHit<pixelTopology::HIonPhase1>::SiPixelRecHit;
  };
  
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

DEFINE_FWK_ALPAKA_MODULE(SiPixelRecHitPhase1);
DEFINE_FWK_ALPAKA_MODULE(SiPixelRecHitHIonPhase1);
DEFINE_FWK_ALPAKA_MODULE(SiPixelRecHitPhase2);
