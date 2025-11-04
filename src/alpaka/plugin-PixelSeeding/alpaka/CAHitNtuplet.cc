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
#include "AlpakaDataFormats/CAGeometryHost.h"

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

#define GPU_DEBUG

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  template <typename TrackerTraits>
  class CAHitNtuplet
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
    // explicit CAHitNtuplet(const edm::ParameterSet& iConfig, const ::reco::CAGeometryParams* iCache);
    explicit CAHitNtuplet(edm::ProductRegistry& reg);
    ~CAHitNtuplet() override = default;

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
  CAHitNtuplet<TrackerTraits>::CAHitNtuplet(edm::ProductRegistry& reg)
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
  // void CAHitNtuplet<TrackerTraits>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //   edm::ParameterSetDescription desc;

  //   desc.add<edm::InputTag>("pixelRecHitSrc", edm::InputTag("siPixelRecHitsPreSplittingAlpaka"));

  //   Algo::fillPSetDescription(desc);
  //   descriptions.addWithDefaultLabel(desc);
  // }

  template <typename TrackerTraits>
  void CAHitNtuplet<TrackerTraits>::produce(edm::Event& iEvent, const edm::EventSetup& es) {
    auto bf = 0.0114256972711507; //1. / es.getData(tokenField_).inverseBzAtOriginInGeV();

    auto const& hGeometry = es.get<::reco::CAGeometryHost>();//runCache()->get(iEvent.queue());
    auto const& phits = iEvent.get(tokenHit_);
    cms::alpakatools::ScopedContextProduce<Queue> ctx{phits};
    auto const& hits = ctx.get(phits);
    // std::array<double, 1> nHitsV = {{double(hits.nHits())}};
    // std::array<double, 1> emptyV;

/// TODO: make a helper for this, have the automatic mechamism for copy (later).
#if defined ALPAKA_ACC_CPU_B_SEQ_T_SEQ_ENABLED or defined ALPAKA_ACC_CPU_B_TBB_T_SEQ_ENABLED
    reco::CAGeometrySoACollection const& geometry = hGeometry;
#else
    reco::CAGeometrySoACollection geometry = cms::alpakatools::CopyToDevice<::reco::CAGeometryHost>::copyAsync(ctx.stream(), hGeometry);
    alpaka::wait(ctx.stream());
#endif

#ifdef GPU_DEBUG
  std::cout << "[CAHitNtuplet::GPU_DEBUG] Starting produce() for "
            << TrackerTraits::nameModifier << std::endl;
  std::cout << "[CAHitNtuplet::GPU_DEBUG] Number of hits: " << hits.nHits() << std::endl;
#endif

  uint32_t const maxTuples = maxNumberOfTuples_;
  uint32_t const maxDoublets = maxNumberOfDoublets_;

#ifdef GPU_DEBUG
  std::cout << "[CAHitNtuplet::GPU_DEBUG] maxTuples=" << maxTuples
            << " maxDoublets=" << maxDoublets << std::endl;
#endif

  ctx.emplace(iEvent,
              tokenTrack_,
              deviceAlgo_.makeTuplesAsync(hits, geometry, bf, maxDoublets, maxTuples, ctx.stream()));

#ifdef GPU_DEBUG
  std::cout << "[CAHitNtuplet::GPU_DEBUG] Finished produce() successfully." << std::endl;
#endif

  }

  // using CAHitNtupletPhase1 = CAHitNtuplet<pixelTopology::Phase1>;
  // using CAHitNtupletHIonPhase1 = CAHitNtuplet<pixelTopology::HIonPhase1>;
  // using CAHitNtupletPhase2 = CAHitNtuplet<pixelTopology::Phase2>;

  /// FIXME: These are needed to make these plugins visible when building the plugins.txt list
  /// see: src/alpaka/Makefile:204. This is a workaround but it works for the moment.
  class CAHitNtupletPhase1 : public CAHitNtuplet<pixelTopology::Phase1> {
  public:
    using CAHitNtuplet<pixelTopology::Phase1>::CAHitNtuplet;
  };

  class CAHitNtupletHIonPhase1 : public CAHitNtuplet<pixelTopology::HIonPhase1> {
  public:
    using CAHitNtuplet<pixelTopology::HIonPhase1>::CAHitNtuplet;
  };

  class CAHitNtupletPhase2 : public CAHitNtuplet<pixelTopology::Phase2> {
  public:
    using CAHitNtuplet<pixelTopology::Phase2>::CAHitNtuplet;
  };
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

// #include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"

DEFINE_FWK_ALPAKA_MODULE(CAHitNtupletPhase1);
DEFINE_FWK_ALPAKA_MODULE(CAHitNtupletHIonPhase1);
DEFINE_FWK_ALPAKA_MODULE(CAHitNtupletPhase2);
