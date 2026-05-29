#ifndef PixelSeeding_alpaka_CAHitNtupletGeneratorKernels_h
#define PixelSeeding_alpaka_CAHitNtupletGeneratorKernels_h

// #define GPU_DEBUG
// #define DUMP_GPU_TK_TUPLES

#include <cstdint>

#include <alpaka/alpaka.hpp>

#include "Framework/ConfigRegistry.h"

#include "AlpakaDataFormats/TrackDefinitions.h"
#include "AlpakaDataFormats/TracksHost.h"
#include "AlpakaDataFormats/alpaka/TrackUtilities.h"
#include "AlpakaDataFormats/TrackingRecHitsSoA.h"
#include "AlpakaCore/HistoContainer.h"
#include "AlpakaCore/config.h"
#include "AlpakaCore/memory.h"
#include "AlpakaDataFormats/CAGeometrySoA.h"
#include "AlpakaDataFormats/alpaka/CAPairSoACollection.h"

#include "CACell.h"
#include "CAPixelDoublets.h"
#include "CAStructures.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  using namespace ::caStructures;

  namespace caHitNtupletGenerator {

    //Counters
    struct Counters {
      unsigned long long nEvents;
      unsigned long long nHits;
      unsigned long long nCells;
      unsigned long long nTuples;
      unsigned long long nTrips;
      unsigned long long nCellTracks;
      unsigned long long nFitTracks;
      unsigned long long nLooseTracks;
      unsigned long long nGoodTracks;
      unsigned long long nUsedHits;
      unsigned long long nDupHits;
      unsigned long long nFishCells;
      unsigned long long nKilledCells;
      unsigned long long nEmptyCells;
      unsigned long long nZeroTrackCells;
    };

    //Full list of params = algo params + quality cuts
    //Generic template
    template <typename TrackerTraits, typename Enable = void>
    struct ParamsT {};

    template <typename TrackerTraits>
    struct ParamsT<TrackerTraits, pixelTopology::isPhase1Topology<TrackerTraits>> {
      using TT = TrackerTraits;
      using QualityCuts = ::pixelTrack::QualityCutsT<TT>;  //track quality cuts

      static constexpr AlgoParams defaultAlgoParams() {
      return {
          3.8f, // bField
          // Container sizes
          5.0f,   // avgHitsPerTrack_
          25.0f,  // avgCellsPerHit_
          2.0f,   // avgCellsPerCell_
          1.0f,   // avgTracksPerCell_

          // Algorithm parameters
          3,      // minHitsPerNtuplet_
          10,     // minHitsForSharingCut_
          0.9f,   // ptmin_
          1.0f / (0.35f * 87.0f),  // hardCurvCut_
          12.0f,  // cellZ0Cut_
          0.5f,   // cellPtCut_

          // Pixel cluster cut params
          8.0f * 0.0285f / 0.015f, // dzdrFact_
          36, //1, // minYsizeB1_
          28, //1, // minYsizeB2_
          28, //28, // maxDYsize12_
          20, //20, // maxDYsize_
          20, //20, // maxDYPred_

          // Flags
          false, // useRiemannFit_
          true, // fitNas4_
          true,  // earlyFishbone_
          false, // lateFishbone_
          false, // doStats_ (fillStatistics)
          false,  // doSharedHitCut_
          false, // dupPassThrough_
          true   // useSimpleTripletCleaner_
      };
    }

      AlgoParams makeAlgoParams(edm::Config const& cfg) const {
        return {
            static_cast<float>(cfg.value("BField", 3.8f)),
            // Container sizes
            static_cast<float>(cfg.value("avgHitsPerTrack", 5.0f)),
            static_cast<float>(cfg.value("avgCellsPerHit", 25.0f)),
            static_cast<float>(cfg.value("avgCellsPerCell", 2.0f)),
            static_cast<float>(cfg.value("avgTracksPerCell", 1.0f)),

            // Algorithm parameters
            static_cast<uint16_t>(cfg.value("minHitsPerNtuplet", 3)),
            static_cast<uint16_t>(cfg.value("minHitsForSharingCut", 10)),
            static_cast<float>(cfg.value("ptmin", 0.9f)),
            static_cast<float>(cfg.value("hardCurvCut", 1.0f / (0.35f * 87.0f))),
            static_cast<float>(cfg.value("cellZ0Cut", 12.0f)),
            static_cast<float>(cfg.value("cellPtCut", 0.5f)),

            // Pixel cluster cut params
            static_cast<float>(cfg.value("dzdrFact", 8.0f * 0.0285f / 0.015f)),
            static_cast<int16_t>(cfg.value("minYsizeB1", 36)),
            static_cast<int16_t>(cfg.value("minYsizeB2", 28)),
            static_cast<int16_t>(cfg.value("maxDYsize12", 28)),
            static_cast<int16_t>(cfg.value("maxDYsize", 20)),
            static_cast<int16_t>(cfg.value("maxDYPred", 20)),

            // Flags
            static_cast<bool>(cfg.value("useRiemannFit", false)),
            static_cast<bool>(cfg.value("fitNas4", true)),
            static_cast<bool>(cfg.value("earlyFishbone", true)),
            static_cast<bool>(cfg.value("lateFishbone", false)),
            static_cast<bool>(cfg.value("doStats", false)),
            static_cast<bool>(cfg.value("doSharedHitCut", false)),
            static_cast<bool>(cfg.value("dupPassThrough", false)),
            static_cast<bool>(cfg.value("useSimpleTripletCleaner", true))
        };
      }
    
    ///TODO: decommission these default methods
    static constexpr QualityCuts defaultQualityCuts() {
    return {
        // polynomial coefficients for pT-dependent chi2 cut
        {0.68177776, 0.74609577, -0.08035491, 0.00315399},
        // max pT used for chi2 cut
        10.,
        // chi2 scale factor
        30.,
        // triplet cuts
        {0.3, 0.5, 12.0}, // |tip|, pt, |zip| 
        // quadruplet cuts
        {0.5, 0.3, 12.0}};
  }
    QualityCuts makeQualityCuts(edm::Config const& cfg) const {
      return {
          // polynomial coefficients for pT-dependent chi2 cut
          {static_cast<float>(cfg.value("chi2Coeff0", 0.68177776f)),
          static_cast<float>(cfg.value("chi2Coeff1", 0.74609577f)),
          static_cast<float>(cfg.value("chi2Coeff2", -0.08035491f)),
          static_cast<float>(cfg.value("chi2Coeff3", 0.00315399f))},
          // max pT used for chi2 cut
          static_cast<float>(cfg.value("maxPtForChi2Cut", 10.f)),
          // chi2 scale factor
         static_cast<float>(cfg.value("chi2ScaleFactor", 30.f)),
          // triplet cuts
          {
              static_cast<float>(cfg.value("tripletTipCut", 0.3f)),
              static_cast<float>(cfg.value("tripletPtCut", 0.5f)),
              static_cast<float>(cfg.value("tripletZipCut", 12.0f))},
          // quadruplet cuts
          {
              static_cast<float>(cfg.value("quadrupletTipCut", 0.5f)),
              static_cast<float>(cfg.value("quadrupletPtCut", 0.3f)),
              static_cast<float>(cfg.value("quadrupletZipCut", 12.0f))
            }
            };
      }

      ParamsT() : algoParams_(defaultAlgoParams()), qualityCuts_(defaultQualityCuts()) {}

      ParamsT(edm::Config const& cfg) : algoParams_(makeAlgoParams(cfg)), qualityCuts_(makeQualityCuts(cfg)) {}

      ParamsT(AlgoParams const& commonCuts, QualityCuts const& qualityCuts)
          : algoParams_(commonCuts), qualityCuts_(qualityCuts) {}

      const AlgoParams algoParams_;
      const QualityCuts qualityCuts_{// polynomial coefficients for the pT-dependent chi2 cut
                                     {0.68177776, 0.74609577, -0.08035491, 0.00315399},
                                     // max pT used to determine the chi2 cut
                                     10.,
                                     // chi2 scale factor: 30 for broken line fit, 45 for Riemann fit
                                     30.,
                                     // regional cuts for triplets
                                     {
                                         0.3,  // |Tip| < 0.3 cm
                                         0.5,  // pT > 0.5 GeV
                                         12.0  // |Zip| < 12.0 cm
                                     },
                                     // regional cuts for quadruplets
                                     {
                                         0.5,  // |Tip| < 0.5 cm
                                         0.3,  // pT > 0.3 GeV
                                         12.0  // |Zip| < 12.0 cm
                                     }};

    };  // Params Phase1

    template <typename TrackerTraits>
    struct ParamsT<TrackerTraits, pixelTopology::isPhase2Topology<TrackerTraits>> : public AlgoParams {
      using TT = TrackerTraits;
      using QualityCuts = ::pixelTrack::QualityCutsT<TT>;

      static constexpr AlgoParams defaultAlgoParams() {
      return {
          3.8f, // bField
          // ---- Container sizes ----
          7.0f,   // avgHitsPerTrack_
          6.0f,   // avgCellsPerHit_
          0.151f, // avgCellsPerCell_
          0.040f, // avgTracksPerCell_

          // ---- Algorithm Parameters ----
          4,      // minHitsPerNtuplet_
          10,     // minHitsForSharingCut_
          0.9f,   // ptmin_  (kept same unless you want to change)
          1.0f / (0.35f * 87.0f),  // hardCurvCut_
          7.5f,   // cellZ0Cut_
          0.85f,  // cellPtCut_

          // ---- Pixel Cluster Cut Params ----
          8.0f * 0.0285f / 0.015f, // dzdrFact_
          25, // minYsizeB1_
          15, // minYsizeB2_
          12, // maxDYsize12_
          10, // maxDYsize_
          20, // maxDYPred_

          // ---- Flags ----
          false, // useRiemannFit_
          false, // fitNas4_
          true,  // earlyFishbone_
          false, // lateFishbone_
          false, // doStats_
          true,  // doSharedHitCut_
          false, // dupPassThrough_
          true   // useSimpleTripletCleaner_
      };
    }

      AlgoParams makeAlgoParams(edm::Config const& cfg) const {
        return {
            static_cast<float>(cfg.value("BField", 3.8f)),
            // ---- Container sizes ----
            static_cast<float>(cfg.value("avgHitsPerTrack", 7.0f)),
            static_cast<float>(cfg.value("avgCellsPerHit", 6.0f)),
            static_cast<float>(cfg.value("avgCellsPerCell", 0.151f)),
            static_cast<float>(cfg.value("avgTracksPerCell", 0.040f)),

            // ---- Algorithm Parameters ----
            static_cast<uint16_t>(cfg.value("minHitsPerNtuplet", 4)),
            static_cast<uint16_t>(cfg.value("minHitsForSharingCut", 10)),
            static_cast<float>(cfg.value("ptmin", 0.9f)),
            static_cast<float>(cfg.value("hardCurvCut", 1.0f / (0.35f * 87.0f))),
            static_cast<float>(cfg.value("cellZ0Cut", 7.5f)),
            static_cast<float>(cfg.value("cellPtCut", 0.85f)),

            // ---- Pixel Cluster Cut Params ----
            static_cast<float>(cfg.value("dzdrFact", 8.0f * 0.0285f / 0.015f)),
            static_cast<int16_t>(cfg.value("minYsizeB1", 25)),
            static_cast<int16_t>(cfg.value("minYsizeB2", 15)),
            static_cast<int16_t>(cfg.value("maxDYsize12", 12)),
            static_cast<int16_t>(cfg.value("maxDYsize", 10)),
            static_cast<int16_t>(cfg.value("maxDYPred", 20)),

            // ---- Flags ----
            static_cast<bool>(cfg.value("useRiemannFit", false)),
            static_cast<bool>(cfg.value("fitNas4", false)),
            static_cast<bool>(cfg.value("earlyFishbone", true)),
            static_cast<bool>(cfg.value("lateFishbone", false)),
            static_cast<bool>(cfg.value("doStats", false)),
            static_cast<bool>(cfg.value("doSharedHitCut", true)),
            static_cast<bool>(cfg.value("dupPassThrough", false)),
            static_cast<bool>(cfg.value("useSimpleTripletCleaner", true))
        };
      }

      QualityCuts makeQualityCuts(edm::Config const& cfg) const {
      return {
          static_cast<float>(cfg.value("maxChi2",5.0f)),
          static_cast<float>(cfg.value("minPtCut", 0.9f)),
          static_cast<float>(cfg.value("maxZip", 0.4f)),
          static_cast<float>(cfg.value("maxTip", 12.0f)),
      };
    }

      static constexpr QualityCuts defaultQualityCuts() {
      return {5.0f, /*chi2*/ 0.9f, /* pT in Gev*/ 0.4f, /*zip in cm*/ 12.0f /*tip in cm*/};
      }

      ParamsT() : algoParams_(defaultAlgoParams()), qualityCuts_(defaultQualityCuts()) {}
      
      ParamsT(edm::Config const& cfg) : algoParams_(makeAlgoParams(cfg)), qualityCuts_(makeQualityCuts(cfg)) {}

      ParamsT(AlgoParams const& commonCuts, QualityCuts const& qualityCuts)
          : algoParams_(commonCuts), qualityCuts_(qualityCuts) {}

      // quality cuts
      const AlgoParams algoParams_;
      const QualityCuts qualityCuts_{5.0f, /*chi2*/ 0.9f, /* pT in Gev*/ 0.4f, /*zip in cm*/ 12.0f /*tip in cm*/};

    };  // Params Phase2

    template <typename TrackerTraits>
    struct ParamsT<TrackerTraits, pixelTopology::isColliderMLPhase1Topology<TrackerTraits>> {
      using TT = TrackerTraits;
      using QualityCuts = ::pixelTrack::QualityCutsT<TT>;  //track quality cuts

      static constexpr AlgoParams defaultAlgoParams() {
      return {
          3.0f, // bField
          // Container sizes
          5.0f,   // avgHitsPerTrack_
          6.0f,   // avgCellsPerHit_
          0.151f, // avgCellsPerCell_
          0.040f, // avgTracksPerCell_

          // Algorithm parameters
          3,      // minHitsPerNtuplet_
          10,     // minHitsForSharingCut_
          0.9f,   // ptmin_
          0.015f,  // hardCurvCut_
          9.1f,  // cellZ0Cut_
          0.85f,   // cellPtCut_

          // Pixel cluster cut params
          8.0f * 0.0285f / 0.015f, // dzdrFact_
          0, //1, // minYsizeB1_
          0, //1, // minYsizeB2_
          0, //28, // maxDYsize12_
          0, //20, // maxDYsize_
          0, //20, // maxDYPred_

          // Flags
          false, // useRiemannFit_
          true, // fitNas4_
          true,  // earlyFishbone_
          false, // lateFishbone_
          false, // doStats_ (fillStatistics)
          false,  // doSharedHitCut_
          false, // dupPassThrough_
          true   // useSimpleTripletCleaner_
      };
    }

      AlgoParams makeAlgoParams(edm::Config const& cfg) const {
        return {
            static_cast<float>(cfg.value("BField", 3.0f)),
            // Container sizes
            static_cast<float>(cfg.value("avgHitsPerTrack", 5.0f)),
            static_cast<float>(cfg.value("avgCellsPerHit", 6.0f)),
            static_cast<float>(cfg.value("avgCellsPerCell", 0.151f)),
            static_cast<float>(cfg.value("avgTracksPerCell", 0.040f)),

            // Algorithm parameters
            static_cast<uint16_t>(cfg.value("minHitsPerNtuplet", 3)),
            static_cast<uint16_t>(cfg.value("minHitsForSharingCut", 10)),
            static_cast<float>(cfg.value("ptmin", 0.9f)),
            static_cast<float>(cfg.value("hardCurvCut", 0.015f)),
            static_cast<float>(cfg.value("cellZ0Cut", 9.1f)),
            static_cast<float>(cfg.value("cellPtCut", 0.85f)),

            // Pixel cluster cut params
            static_cast<float>(cfg.value("dzdrFact", 8.0f * 0.0285f / 0.015f)),
            static_cast<int16_t>(cfg.value("minYsizeB1", 0)),
            static_cast<int16_t>(cfg.value("minYsizeB2", 0)),
            static_cast<int16_t>(cfg.value("maxDYsize12", 0)),
            static_cast<int16_t>(cfg.value("maxDYsize", 0)),
            static_cast<int16_t>(cfg.value("maxDYPred", 0)),

            // Flags
            static_cast<bool>(cfg.value("useRiemannFit", false)),
            static_cast<bool>(cfg.value("fitNas4", true)),
            static_cast<bool>(cfg.value("earlyFishbone", true)),
            static_cast<bool>(cfg.value("lateFishbone", false)),
            static_cast<bool>(cfg.value("doStats", false)),
            static_cast<bool>(cfg.value("doSharedHitCut", false)),
            static_cast<bool>(cfg.value("dupPassThrough", false)),
            static_cast<bool>(cfg.value("useSimpleTripletCleaner", true))
        };
      }
    
    ///TODO: decommission these default methods
    static constexpr QualityCuts defaultQualityCuts() {
    return {
        // polynomial coefficients for pT-dependent chi2 cut
        {0.68177776, 0.74609577, -0.08035491, 0.00315399},
        // max pT used for chi2 cut
        10.,
        // chi2 scale factor
        30.,
        // triplet cuts
        {0.3, 0.5, 12.0}, // |tip|, pt, |zip| 
        // quadruplet cuts
        {0.5, 0.3, 12.0}};
  }
    QualityCuts makeQualityCuts(edm::Config const& cfg) const {
      return {
          // polynomial coefficients for pT-dependent chi2 cut
          {static_cast<float>(cfg.value("chi2Coeff0", 0.68177776f)),
          static_cast<float>(cfg.value("chi2Coeff1", 0.74609577f)),
          static_cast<float>(cfg.value("chi2Coeff2", -0.08035491f)),
          static_cast<float>(cfg.value("chi2Coeff3", 0.00315399f))},
          // max pT used for chi2 cut
          static_cast<float>(cfg.value("maxPtForChi2Cut", 10.f)),
          // chi2 scale factor
         static_cast<float>(cfg.value("chi2ScaleFactor", 30.f)),
          // triplet cuts
          {
              static_cast<float>(cfg.value("tripletTipCut", 0.3f)),
              static_cast<float>(cfg.value("tripletPtCut", 0.5f)),
              static_cast<float>(cfg.value("tripletZipCut", 12.0f))},
          // quadruplet cuts
          {
              static_cast<float>(cfg.value("quadrupletTipCut", 0.5f)),
              static_cast<float>(cfg.value("quadrupletPtCut", 0.3f)),
              static_cast<float>(cfg.value("quadrupletZipCut", 12.0f))
            }
            };
      }

      ParamsT() : algoParams_(defaultAlgoParams()), qualityCuts_(defaultQualityCuts()) {}

      ParamsT(edm::Config const& cfg) : algoParams_(makeAlgoParams(cfg)), qualityCuts_(makeQualityCuts(cfg)) {}

      ParamsT(AlgoParams const& commonCuts, QualityCuts const& qualityCuts)
          : algoParams_(commonCuts), qualityCuts_(qualityCuts) {}

      const AlgoParams algoParams_;
      const QualityCuts qualityCuts_{// polynomial coefficients for the pT-dependent chi2 cut
                                     {0.68177776, 0.74609577, -0.08035491, 0.00315399},
                                     // max pT used to determine the chi2 cut
                                     10.,
                                     // chi2 scale factor: 30 for broken line fit, 45 for Riemann fit
                                     30.,
                                     // regional cuts for triplets
                                     {
                                         0.3,  // |Tip| < 0.3 cm
                                         0.5,  // pT > 0.5 GeV
                                         12.0  // |Zip| < 12.0 cm
                                     },
                                     // regional cuts for quadruplets
                                     {
                                         0.5,  // |Tip| < 0.5 cm
                                         0.3,  // pT > 0.3 GeV
                                         12.0  // |Zip| < 12.0 cm
                                     }};

    };  // Params ColliderMLPhase1

    template <typename TrackerTraits>
    struct ParamsT<TrackerTraits, pixelTopology::isColliderMLPhase2Topology<TrackerTraits>> : public AlgoParams {
      using TT = TrackerTraits;
      using QualityCuts = ::pixelTrack::QualityCutsT<TT>;

      static constexpr AlgoParams defaultAlgoParams() {
      return {
          3.0f, // bField
          // ---- Container sizes ----
          7.0f,   // avgHitsPerTrack_
          6.0f,   // avgCellsPerHit_
          0.151f, // avgCellsPerCell_
          0.040f, // avgTracksPerCell_

          // ---- Algorithm Parameters ----
          4,      // minHitsPerNtuplet_
          10,     // minHitsForSharingCut_
          0.9f,   // ptmin_  (kept same unless you want to change)
          0.0125, // hardCurvCut_
          9.1f,   // cellZ0Cut_
          0.85f,  // cellPtCut_

          // ---- Pixel Cluster Cut Params ----
          8.0f * 0.0285f / 0.015f, // dzdrFact_
          0, // minYsizeB1_
          0, // minYsizeB2_
          0, // maxDYsize12_
          0, // maxDYsize_
          0, // maxDYPred_

          // ---- Flags ----
          false, // useRiemannFit_
          false, // fitNas4_
          true,  // earlyFishbone_
          false, // lateFishbone_
          false, // doStats_
          true,  // doSharedHitCut_
          false, // dupPassThrough_
          true   // useSimpleTripletCleaner_
      };
    }

      AlgoParams makeAlgoParams(edm::Config const& cfg) const {
        return {
            static_cast<float>(cfg.value("BField", 3.0f)),
            // ---- Container sizes ----
            static_cast<float>(cfg.value("avgHitsPerTrack", 7.0f)),
            static_cast<float>(cfg.value("avgCellsPerHit", 6.0f)),
            static_cast<float>(cfg.value("avgCellsPerCell", 0.151f)),
            static_cast<float>(cfg.value("avgTracksPerCell", 0.040f)),

            // ---- Algorithm Parameters ----
            static_cast<uint16_t>(cfg.value("minHitsPerNtuplet", 4)),
            static_cast<uint16_t>(cfg.value("minHitsForSharingCut", 10)),
            static_cast<float>(cfg.value("ptmin", 0.9f)),
            static_cast<float>(cfg.value("hardCurvCut", 0.0125)),
            static_cast<float>(cfg.value("cellZ0Cut", 9.1f)),
            static_cast<float>(cfg.value("cellPtCut", 0.85f)),

            // ---- Pixel Cluster Cut Params ----
            static_cast<float>(cfg.value("dzdrFact", 8.0f * 0.0285f / 0.015f)),
            static_cast<int16_t>(cfg.value("minYsizeB1", 0)),
            static_cast<int16_t>(cfg.value("minYsizeB2", 0)),
            static_cast<int16_t>(cfg.value("maxDYsize12", 0)),
            static_cast<int16_t>(cfg.value("maxDYsize", 0)),
            static_cast<int16_t>(cfg.value("maxDYPred", 0)),

            // ---- Flags ----
            static_cast<bool>(cfg.value("useRiemannFit", false)),
            static_cast<bool>(cfg.value("fitNas4", false)),
            static_cast<bool>(cfg.value("earlyFishbone", true)),
            static_cast<bool>(cfg.value("lateFishbone", false)),
            static_cast<bool>(cfg.value("doStats", false)),
            static_cast<bool>(cfg.value("doSharedHitCut", true)),
            static_cast<bool>(cfg.value("dupPassThrough", false)),
            static_cast<bool>(cfg.value("useSimpleTripletCleaner", true))
        };
      }

      QualityCuts makeQualityCuts(edm::Config const& cfg) const {
      return {
          static_cast<float>(cfg.value("maxChi2",5.0f)),
          static_cast<float>(cfg.value("minPtCut", 0.9f)),
          static_cast<float>(cfg.value("maxZip", 0.4f)),
          static_cast<float>(cfg.value("maxTip", 12.0f)),
      };
    }

      static constexpr QualityCuts defaultQualityCuts() {
      return {5.0f, /*chi2*/ 0.9f, /* pT in Gev*/ 0.4f, /*zip in cm*/ 12.0f /*tip in cm*/};
      }

      ParamsT() : algoParams_(defaultAlgoParams()), qualityCuts_(defaultQualityCuts()) {}
      
      ParamsT(edm::Config const& cfg) : algoParams_(makeAlgoParams(cfg)), qualityCuts_(makeQualityCuts(cfg)) {}

      ParamsT(AlgoParams const& commonCuts, QualityCuts const& qualityCuts)
          : algoParams_(commonCuts), qualityCuts_(qualityCuts) {}

      // quality cuts
      const AlgoParams algoParams_;
      const QualityCuts qualityCuts_{5.0f, /*chi2*/ 0.9f, /* pT in Gev*/ 0.4f, /*zip in cm*/ 12.0f /*tip in cm*/};

    };  // Params ColliderMLPhase2

  }  // namespace caHitNtupletGenerator
  template <typename TTTraits>
  class CAHitNtupletGeneratorKernels {
  public:
    using TrackerTraits = TTTraits;

    using SimpleCell = CACell<TrackerTraits>;
    using Params = caHitNtupletGenerator::ParamsT<TrackerTraits>;
    using Counters = caHitNtupletGenerator::Counters;
    // Track qualities
    using Quality = ::pixelTrack::Quality;
    using QualityCuts = ::pixelTrack::QualityCutsT<TrackerTraits>;

    // Histograms

    using PhiBinner = caStructures::PhiBinnerT<TrackerTraits>;  //the traits here define the number of layer/histograms
    using PhiBinnerStorageType = typename PhiBinner::index_type;
    using PhiBinnerView = typename PhiBinner::View;

    using HitToTuple = caStructures::GenericContainer;
    using HitContainer = caStructures::SequentialContainer;
    using TupleMultiplicity = caStructures::GenericContainer;
    using HitToCell = caStructures::GenericContainer;
    using CellToCell = caStructures::GenericContainer;
    using CellToTrack = caStructures::GenericContainer;

    using GenericContainer = caStructures::GenericContainer;
    using GenericContainerStorage = typename GenericContainer::index_type;
    using GenericContainerView = typename GenericContainer::View;
    using DeviceGenericContainerBuffer = std::optional<cms::alpakatools::device_buffer<Device, GenericContainer>>;
    using DeviceGenericStorageBuffer =
        std::optional<cms::alpakatools::device_buffer<Device, GenericContainerStorage[]>>;
    using DeviceGenericOffsetsBuffer =
        std::optional<cms::alpakatools::device_buffer<Device, GenericContainerOffsets[]>>;

    using SequentialContainer = caStructures::SequentialContainer;
    using SequentialContainerStorage = typename SequentialContainer::index_type;
    using SequentialContainerView = typename SequentialContainer::View;
    using DeviceSequentialContainerBuffer = std::optional<cms::alpakatools::device_buffer<Device, SequentialContainer>>;
    using DeviceSequentialStorageBuffer =
        std::optional<cms::alpakatools::device_buffer<Device, SequentialContainerStorage[]>>;
    using DeviceSequentialOffsetsBuffer =
        std::optional<cms::alpakatools::device_buffer<Device, SequentialContainerOffsets[]>>;

    CAHitNtupletGeneratorKernels(Params const& params,
                                 uint32_t nHits,
                                 uint32_t offsetBPIX2,
                                 uint32_t nDoublets,
                                 uint32_t nTracks,
                                 uint16_t nLayers,
                                 Queue& queue);
    ~CAHitNtupletGeneratorKernels() = default;

    TupleMultiplicity const* tupleMultiplicity() const { return device_tupleMultiplicity_->data(); }
    HitContainer const* hitContainer() const { return device_hitContainer_->data(); }
    HitToCell const* hitToCell() const { return device_hitToCell_->data(); }
    HitToTuple const* hitToTuple() const { return device_hitToTuple_->data(); }
    CellToCell const* cellToCell() const { return device_cellToNeighbors_->data(); }
    CellToTrack const* cellToTrack() const { return device_cellToTracks_->data(); }

    void prepareHits(const HitsConstView& hh,
                     const HitModulesConstView& mm,
                     const ::reco::CALayersSoAConstView& ll,
                     Queue& queue);

    void launchKernels(const HitsConstView& hh,
                       uint32_t offsetBPIX2,
                       uint16_t nLayers,
                       TkSoAView& track_view,
                       TkHitsSoAView& track_hits_view,
                       const ::reco::CALayersSoAConstView& ll,
                       const ::reco::CAGraphSoAConstView& cc,
                       Queue& queue);

    void classifyTuples(const HitsConstView& hh, TkSoAView& track_view, Queue& queue);

    void buildDoublets(const HitsConstView& hh,
                       const ::reco::CAGraphSoAConstView& cc,
                       const ::reco::CALayersSoAConstView& ll,
                       uint32_t offsetBPIX2,
                       Queue& queue);

    static void printCounters();

  private:
    // params
    Params const& m_params;
    std::optional<cms::alpakatools::device_buffer<Device, Counters>> counters_;

    // Hits->Track
    DeviceGenericContainerBuffer device_hitToTuple_;
    DeviceGenericStorageBuffer device_hitToTupleStorage_;
    DeviceGenericOffsetsBuffer device_hitToTupleOffsets_;
    GenericContainerView device_hitToTupleView_;

    // (Outer) Hits-> Cells
    DeviceGenericContainerBuffer device_hitToCell_;
    DeviceGenericStorageBuffer device_hitToCellStorage_;
    DeviceGenericOffsetsBuffer device_hitToCellOffsets_;
    GenericContainerView device_hitToCellView_;

    // Hits Phi Binner
    std::optional<cms::alpakatools::device_buffer<Device, PhiBinner>> device_hitPhiHist_;
    std::optional<cms::alpakatools::device_buffer<Device, PhiBinnerStorageType[]>> device_phiBinnerStorage_;
    PhiBinnerView device_hitPhiView_;
    std::optional<cms::alpakatools::device_buffer<Device, hindex_type[]>> device_layerStarts_;

    // Cells-> Neighbor Cells
    DeviceGenericContainerBuffer device_cellToNeighbors_;
    DeviceGenericStorageBuffer device_cellToNeighborsStorage_;
    DeviceGenericOffsetsBuffer device_cellToNeighborsOffsets_;
    GenericContainerView device_cellToNeighborsView_;

    // Cells-> Tracks
    DeviceGenericContainerBuffer device_cellToTracks_;
    DeviceGenericStorageBuffer device_cellToTracksStorage_;
    DeviceGenericOffsetsBuffer device_cellToTracksOffsets_;
    GenericContainerView device_cellToTracksView_;

    // Tracks->Hits
    DeviceSequentialContainerBuffer device_hitContainer_;
    DeviceGenericStorageBuffer device_hitContainerStorage_;
    DeviceSequentialOffsetsBuffer device_hitContainerOffsets_;
    SequentialContainerView device_hitContainerView_;

    // No.Hits -> Track (Multiplicity)
    DeviceGenericContainerBuffer device_tupleMultiplicity_;
    DeviceGenericStorageBuffer device_tupleMultiplicityStorage_;
    DeviceGenericOffsetsBuffer device_tupleMultiplicityOffsets_;
    GenericContainerView device_tupleMultiplicityView_;

    std::optional<cms::alpakatools::device_buffer<Device, SimpleCell[]>> device_simpleCells_;

    std::optional<cms::alpakatools::device_buffer<Device, cms::alpakatools::AtomicPairCounter::DoubleWord[]>>
        device_extraStorage_;
    cms::alpakatools::AtomicPairCounter* device_hitTuple_apc_;
    std::optional<cms::alpakatools::device_buffer<Device, uint32_t[]>> device_nCells_;
    std::optional<cms::alpakatools::device_buffer<Device, uint32_t[]>> device_nTriplets_;
    std::optional<cms::alpakatools::device_buffer<Device, uint32_t[]>> device_nCellTracks_;

    std::optional<CAPairSoACollection> deviceTriplets_;
    std::optional<CAPairSoACollection> deviceTracksCells_;

    // this could be inferred from the above buffers
    // but seems cleaner to have a dedicate variable
    uint32_t maxNumberOfDoublets_;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif  // PixelSeeding_alpaka_CAHitNtupletGeneratorKernels_h
