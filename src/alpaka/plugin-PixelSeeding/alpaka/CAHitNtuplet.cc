#include <alpaka/alpaka.hpp>

// #include <TFormula.h>
// #include "CommonTools/Utils/interface/FormulaEvaluator.h"

#include "AlpakaCore/Product.h"
#include "AlpakaCore/ScopedContext.h"
#include "AlpakaCore/config.h"

#include "AlpakaDataFormats/TracksHost.h"
#include "AlpakaDataFormats/alpaka/TracksSoACollection.h"
#include "AlpakaDataFormats/TracksDevice.h"

#include "AlpakaDataFormats/TrackingRecHitsHost.h"
#include "AlpakaDataFormats/alpaka/TrackingRecHitsSoACollection.h"
#include "AlpakaDataFormats/TrackingRecHitsDevice.h"

#include "Framework/EDProducer.h"
#include "Framework/Event.h"
#include "Framework/EventSetup.h"
#include "Framework/PluginFactory.h"

// #include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
// #include "RecoTracker/TkMSParametrization/interface/PixelRecoUtilities.h"

// #include "RecoTracker/Record/interface/TrackerRecoGeometryRecord.h"
// #include "plugin-PixelSeeding/alpaka/CAGeometrySoACollection.h"
// #include "plugin-PixelSeeding/CAGeometryHost.h"
#include "CAHitNtupletGenerator.h"

// #include "HeterogeneousCore/AlpakaCore/interface/MoveToDeviceCache.h"
// #include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
// #include "Geometry/Records/interface/TrackerTopologyRcd.h"
// #include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
// #include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "AlpakaDataFormats/CAGeometrySoA.h"
#include "AlpakaDataFormats/alpaka/CAGeometrySoACollection.h"


// #define GPU_DEBUG

// namespace reco {
//   struct CAGeometryParams {
//     //Constructor from ParameterSet
//     CAGeometryParams(edm::ParameterSet const& iConfig)
//         : caThetaCuts_(iConfig.getParameter<std::vector<double>>("caThetaCuts")),
//           caDCACuts_(iConfig.getParameter<std::vector<double>>("caDCACuts")),
//           pairGraph_(iConfig.getParameter<std::vector<unsigned int>>("pairGraph")),
//           startingPairs_(iConfig.getParameter<std::vector<unsigned int>>("startingPairs")),
//           phiCuts_(iConfig.getParameter<std::vector<int>>("phiCuts")),
//           minZ_(iConfig.getParameter<std::vector<double>>("minZ")),
//           maxZ_(iConfig.getParameter<std::vector<double>>("maxZ")),
//           maxR_(iConfig.getParameter<std::vector<double>>("maxR")) {}

//     // Layers params
//     const std::vector<double> caThetaCuts_;
//     const std::vector<double> caDCACuts_;

//     // Cells params
//     const std::vector<unsigned int> pairGraph_;
//     const std::vector<unsigned int> startingPairs_;
//     const std::vector<int> phiCuts_;
//     const std::vector<double> minZ_;
//     const std::vector<double> maxZ_;
//     const std::vector<double> maxR_;

//     mutable edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> tokenGeometry_;
//     mutable edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> tokenTopology_;
//   };

// }  // namespace reco

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  template <typename TrackerTraits>
  class CAHitNtupletAlpakaFromHits
      // : public stream::EDProducer<edm::GlobalCache<::reco::CAGeometryParams>,
      //                             edm::RunCache<cms::alpakatools::MoveToDeviceCache<Device, ::reco::CAGeometryHost>>> 
      : public edm::EDProducer{
    using HitsConstView = ::reco::TrackingRecHitConstView;
    using HitsOnDevice = reco::TrackingRecHitsSoACollection;
    using HitsOnHost = ::reco::TrackingRecHitHost;

    using CAGeometryOnDevice = reco::CAGeometrySoACollection;

    using TkSoAHost = ::reco::TracksHost;
    using TkSoADevice = reco::TracksSoACollection;

    using Algo = CAHitNtupletGenerator<TrackerTraits>;
    using Params = caHitNtupletGenerator::ParamsT<TrackerTraits>;

    // using CAGeometryCache = cms::alpakatools::MoveToDeviceCache<Device, ::reco::CAGeometryHost>;
    using Rotation = SOARotation<float>;
    using Frame = SOAFrame<float>;

  public:
    // explicit CAHitNtupletAlpakaFromHits(const edm::ParameterSet& iConfig, const ::reco::CAGeometryParams* iCache);
    explicit CAHitNtupletAlpakaFromHits(edm::ProductRegistry& reg);
    ~CAHitNtupletAlpakaFromHits() override = default;

    void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override;

  private:
    // const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> tokenField_;

    edm::EDGetTokenT<cms::alpakatools::Product<Queue, HitsOnDevice>> tokenHit_;
    edm::EDPutTokenT<cms::alpakatools::Product<Queue, TkSoADevice>> tokenTrack_;

    // const ::reco::FormulaEvaluator maxNumberOfDoublets_; //TODO TFormula
    // const ::reco::FormulaEvaluator maxNumberOfTuples_;
    const uint32_t maxNumberOfDoublets_ = 0;
    const uint32_t maxNumberOfTuples_ = 0;                                
    Algo deviceAlgo_;
  };

  template <typename TrackerTraits>
  CAHitNtupletAlpakaFromHits<TrackerTraits>::CAHitNtupletAlpakaFromHits(edm::ProductRegistry& reg)
  // (const edm::ParameterSet& iConfig,
                                                        // const ::reco::CAGeometryParams* iCache)
      :
      // : EDProducer(iConfig),
        // tokenField_(esConsumes()),
        tokenHit_{reg.consumes<cms::alpakatools::Product<Queue, HitsOnDevice>>()},
        tokenTrack_(reg.produces<cms::alpakatools::Product<Queue, TkSoADevice>>()),
        // maxNumberOfDoublets_(iConfig.getParameter<std::string>("maxNumberOfDoublets")),
        // maxNumberOfTuples_(iConfig.getParameter<std::string>("maxNumberOfTuples")),
        maxNumberOfDoublets_(500000),
        maxNumberOfTuples_(200000),
        deviceAlgo_(Params()) //default params
        // deviceAlgo_(iConfig) 
  {
    // iCache->tokenGeometry_ = esConsumes<edm::Transition::BeginRun>();
    // iCache->tokenTopology_ = esConsumes<edm::Transition::BeginRun>();
  }

  // template <typename TrackerTraits>
  // void CAHitNtupletAlpakaFromHits<TrackerTraits>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //   edm::ParameterSetDescription desc;

  //   desc.add<edm::InputTag>("pixelRecHitSrc", edm::InputTag("siPixelRecHitsPreSplittingAlpaka"));

  //   Algo::fillPSetDescription(desc);
  //   descriptions.addWithDefaultLabel(desc);
  // }

  template <typename TrackerTraits>
  void CAHitNtupletAlpakaFromHits<TrackerTraits>::produce(edm::Event& iEvent, const edm::EventSetup& es) {
    auto bf = 0.0114256972711507; //1. / es.getData(tokenField_).inverseBzAtOriginInGeV();

    auto const& geometry = es.get<CAGeometryOnDevice>();//runCache()->get(iEvent.queue());
    auto const& phits = iEvent.get(tokenHit_);
    cms::alpakatools::ScopedContextProduce<Queue> ctx{phits};
    auto const& hits = ctx.get(phits);
    // std::array<double, 1> nHitsV = {{double(hits.nHits())}};
    // std::array<double, 1> emptyV;

    uint32_t const maxTuples = maxNumberOfTuples_; //maxNumberOfTuples_.evaluate(nHitsV, emptyV);
    uint32_t const maxDoublets = maxNumberOfDoublets_; //maxNumberOfDoublets_.evaluate(nHitsV, emptyV);

    ctx.emplace(iEvent,tokenTrack_,
                deviceAlgo_.makeTuplesAsync(hits, geometry, bf, maxDoublets, maxTuples, ctx.stream()));
  }

  using CAHitNtupletAlpakaFromHitsPhase1 = CAHitNtupletAlpakaFromHits<pixelTopology::Phase1>;
  using CAHitNtupletAlpakaFromHitsHIonPhase1 = CAHitNtupletAlpakaFromHits<pixelTopology::HIonPhase1>;
  using CAHitNtupletAlpakaFromHitsPhase2 = CAHitNtupletAlpakaFromHits<pixelTopology::Phase2>;
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

// #include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"

DEFINE_FWK_ALPAKA_MODULE(CAHitNtupletAlpakaFromHitsPhase1);
DEFINE_FWK_ALPAKA_MODULE(CAHitNtupletAlpakaFromHitsHIonPhase1);
DEFINE_FWK_ALPAKA_MODULE(CAHitNtupletAlpakaFromHitsPhase2);
