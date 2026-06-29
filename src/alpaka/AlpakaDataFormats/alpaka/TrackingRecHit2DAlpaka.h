#ifndef AlpakaDataFormats_alpaka_TrackingRecHit2DAlpaka_h
#define AlpakaDataFormats_alpaka_TrackingRecHit2DAlpaka_h

#include "AlpakaCore/config.h"
#include "AlpakaCore/memory.h"
#include "AlpakaDataFormats/TrackingRecHitsSoA.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  using TrackingRecHit2DAlpaka = cms::alpakatools::device_buffer<Device, ::reco::TrackingRecHitSoA>;

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif  // AlpakaDataFormats_alpaka_TrackingRecHit2DAlpaka_h

