#ifndef AlpakaDataFormats_TracksSoACollection_h
#define AlpakaDataFormats_TracksSoACollection_h

#include <type_traits>

#include <alpaka/alpaka.hpp>

#include "Portable/PortableCollection.h"
#include "AlpakaDataFormats/TracksDevice.h"
#include "AlpakaDataFormats/TracksHost.h"
#include "Geometry/SimplePixelTopology.h"
#include "AlpakaCore/CopyToHost.h"
#include "AlpakaCore/config.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::reco {

  using ::reco::TracksDevice;
  using ::reco::TracksHost;
  using TracksSoACollection =
      std::conditional_t<std::is_same_v<Device, alpaka::DevCpu>, TracksHost, TracksDevice<Device>>;

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::reco

#include "AlpakaCore/AssertDeviceMatchesHostCollection.h"
ASSERT_DEVICE_MATCHES_HOST_COLLECTION(reco::TracksSoACollection, reco::TracksHost);

#endif  // AlpakaDataFormats_TracksSoACollection_h
