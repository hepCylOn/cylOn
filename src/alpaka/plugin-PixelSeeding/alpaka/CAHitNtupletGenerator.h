#ifndef PixelSeeding_alpaka_CAHitNtupletGenerator_h
#define PixelSeeding_alpaka_CAHitNtupletGenerator_h

#include <alpaka/alpaka.hpp>

#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "AlpakaDataFormats/TrackDefinitions.h"
#include "AlpakaDataFormats/alpaka/TracksSoACollection.h"
#include "AlpakaDataFormats/TracksHost.h"
#include "AlpakaDataFormats/TracksDevice.h"
#include "AlpakaDataFormats/TrackingRecHitsSoA.h"
#include "AlpakaDataFormats/alpaka/TrackingRecHitsSoACollection.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "AlpakaCore/config.h"
#include "plugin-PixelSeeding/alpaka/CAGeometrySoACollection.h"

#include "CACell.h"
#include "CAHitNtupletGeneratorKernels.h"
#include "HelixFit.h"

namespace edm {
  class ParameterSetDescription;
}  // namespace edm

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  template <typename TrackerTraits>
  class CAHitNtupletGenerator {
  public:
    using HitsView = ::reco::TrackingRecHitView;
    using HitsConstView = ::reco::TrackingRecHitConstView;
    using HitsOnDevice = reco::TrackingRecHitsSoACollection;
    using HitsOnHost = ::reco::TrackingRecHitHost;

    using TkSoADevice = reco::TracksSoACollection;
    using Quality = ::pixelTrack::Quality;

    using QualityCuts = ::pixelTrack::QualityCutsT<TrackerTraits>;
    using Params = caHitNtupletGenerator::ParamsT<TrackerTraits>;
    using Counters = caHitNtupletGenerator::Counters;

    using CAGeometryOnDevice = reco::CAGeometrySoACollection;

  public:
    CAHitNtupletGenerator(const edm::ParameterSet& cfg);

    static void fillPSetDescription(edm::ParameterSetDescription& desc);

    // NOTE: beginJob and endJob were meant to be used
    // to fill the statistics. This is still not implemented in Alpaka
    // since we are missing the begin/endJob functionality for the Alpaka
    // producers.
    //
    // void beginJob();
    // void endJob();

    TkSoADevice makeTuplesAsync(HitsOnDevice const& hits_d,
                                CAGeometryOnDevice const& params_d,
                                float bfield,
                                uint32_t maxDoublets,
                                uint32_t maxTuples,
                                Queue& queue) const;

  private:
    Params m_params;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif  // PixelSeeding_alpaka_CAHitNtupletGenerator_h
