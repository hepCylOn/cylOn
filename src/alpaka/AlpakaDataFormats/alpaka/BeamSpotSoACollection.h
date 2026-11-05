#ifndef AlpakaDataFormats_BeamSpotSoACollection_h
#define AlpakaDataFormats_BeamSpotSoACollection_h

#include <cstdint>

#include <alpaka/alpaka.hpp>

#include "Portable/PortableCollection.h"
#include "AlpakaDataFormats/BeamSpotDevice.h"
#include "AlpakaDataFormats/BeamSpotHost.h"
#include "AlpakaCore/CopyToHost.h"
#include "AlpakaCore/config.h"
#include "AlpakaCore/AssertDeviceMatchesHostCollection.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  using BeamSpotSoACollection =
      std::conditional_t<std::is_same_v<Device, alpaka::DevCpu>, BeamSpotHost, BeamSpotDevice<Device>>;

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

ASSERT_DEVICE_MATCHES_HOST_COLLECTION(BeamSpotSoACollection, BeamSpotHost);

#endif  // AlpakaDataFormats_BeamSpotSoACollection_h

