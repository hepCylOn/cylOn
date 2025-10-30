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
  class SiPixelRecHitProducer : public edm::EDProducer {
  public:
    explicit SiPixelRecHitProducer(edm::ProductRegistry& reg);
    ~SiPixelRecHitProducer() override = default;

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
  SiPixelRecHitProducer<TrackerTraits>::SiPixelRecHitProducer(edm::ProductRegistry& reg)
      : tBeamSpot(reg.consumes<cms::alpakatools::Product<Queue, BeamSpotSoACollection>>()),
        tokenClusters_(reg.consumes<cms::alpakatools::Product<Queue, SiPixelClustersSoACollection>>()),
        tokenDigi_(reg.consumes<cms::alpakatools::Product<Queue, SiPixelDigisSoACollection>>()),
        tokenHit_(reg.produces<cms::alpakatools::Product<Queue, reco::TrackingRecHitsSoACollection>>()) 
        {}

  // template <typename TrackerTraits>
  // void SiPixelRecHitProducer<TrackerTraits>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //   edm::ParameterSetDescription desc;

  //   desc.add<edm::InputTag>("beamSpot", edm::InputTag("offlineBeamSpotSoACollection"));
  //   desc.add<edm::InputTag>("src", edm::InputTag("siPixelClustersPreSplittingAlpaka"));

  //   std::string cpe = "PixelCPEFastParams";
  //   cpe += TrackerTraits::nameModifier;
  //   desc.add<std::string>("CPE", cpe);

  //   descriptions.addWithDefaultLabel(desc);
  // }

  template <typename TrackerTraits>
  void SiPixelRecHitProducer<TrackerTraits>::produce(edm::Event& iEvent, const edm::EventSetup& es){
                                                    
    auto const& fcpe = es.get<PixelCPEFast<TrackerTraits>>();

    auto const& pclusters = iEvent.get(tokenClusters_);
    cms::alpakatools::ScopedContextProduce<Queue> ctx{pclusters};

    auto const& clusters = ctx.get(pclusters);
    auto const& digis = ctx.get(iEvent, tokenDigi_);
    auto const& bs = ctx.get(iEvent, tBeamSpot);

    ctx.emplace(iEvent,
                tokenHit_,
                Algo_.makeHitsAsync(digis, clusters, bs.data(), fcpe.getGPUProductAsync(ctx.stream()), ctx.stream()));

  }
  using SiPixelRecHitProducerPhase1 = SiPixelRecHitProducer<pixelTopology::Phase1>;
  using SiPixelRecHitProducerHIonPhase1 = SiPixelRecHitProducer<pixelTopology::HIonPhase1>;
  using SiPixelRecHitProducerPhase2 = SiPixelRecHitProducer<pixelTopology::Phase2>;
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

DEFINE_FWK_ALPAKA_MODULE(SiPixelRecHitProducerPhase1);
DEFINE_FWK_ALPAKA_MODULE(SiPixelRecHitProducerHIonPhase1);
DEFINE_FWK_ALPAKA_MODULE(SiPixelRecHitProducerPhase2);
