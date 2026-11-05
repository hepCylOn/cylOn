#ifndef RecoLocalTracker_SiPixelRecHits_PixelRecHitKernel_h
#define RecoLocalTracker_SiPixelRecHits_PixelRecHitKernel_h

#include <cstdint>

#include <alpaka/alpaka.hpp>

#include "AlpakaDataFormats/BeamSpotPOD.h"
#include "AlpakaDataFormats/alpaka/SiPixelClustersSoACollection.h"
#include "AlpakaDataFormats/SiPixelClustersDevice.h"
#include "AlpakaDataFormats/SiPixelDigisDevice.h"
#include "AlpakaDataFormats/alpaka/SiPixelDigisSoACollection.h"
#include "AlpakaDataFormats/TrackingRecHitsDevice.h"
#include "AlpakaDataFormats/alpaka/TrackingRecHitsSoACollection.h"
#include "AlpakaCore/config.h"
#include "Geometry/SimplePixelTopology.h"
#include "CondFormats/pixelCPEforDevice.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  namespace pixelgpudetails {
    using namespace cms::alpakatools;
    using namespace ALPAKA_ACCELERATOR_NAMESPACE::reco;
    template <typename TrackerTraits>
    class PixelRecHitKernel {
    public:
      PixelRecHitKernel() = default;
      ~PixelRecHitKernel() = default;

      PixelRecHitKernel(const PixelRecHitKernel&) = delete;
      PixelRecHitKernel(PixelRecHitKernel&&) = delete;
      PixelRecHitKernel& operator=(const PixelRecHitKernel&) = delete;
      PixelRecHitKernel& operator=(PixelRecHitKernel&&) = delete;

      using ParamsOnDevice = pixelCPEforDevice::ParamsOnDevice;

      reco::TrackingRecHitsSoACollection makeHitsAsync(SiPixelDigisSoACollection const& digis_d,
                                                       SiPixelClustersSoACollection const& clusters_d,
                                                       BeamSpotPOD const* bs_d,
                                                       ParamsOnDevice const* cpeParams,
                                                       Queue queue) const;
    };
  }  // namespace pixelgpudetails
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif  // RecoLocalTracker_SiPixelRecHits_PixelRecHitKernel_h
