#ifndef AlpakaDataFormats_CAGeometrySoACollection_h
#define AlpakaDataFormats_CAGeometrySoACollection_h

#include <cstdint>

#include <alpaka/alpaka.hpp>

#include "Portable/PortableCollection.h"
#include "AlpakaDataFormats/CAGeometryDevice.h"
#include "AlpakaDataFormats/CAGeometryHost.h"
#include "AlpakaDataFormats/CAGeometrySoA.h"
#include "AlpakaCore/CopyToHost.h"
#include "AlpakaCore/config.h"
#include "AlpakaCore/AssertDeviceMatchesHostCollection.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::reco {

  using ::reco::CAGeometryDevice;
  using ::reco::CAGeometryHost;
  using CAGeometrySoACollection =
      std::conditional_t<std::is_same_v<Device, alpaka::DevCpu>, CAGeometryHost, CAGeometryDevice<Device>>;

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::reco

ASSERT_DEVICE_MATCHES_HOST_COLLECTION(reco::CAGeometrySoACollection, reco::CAGeometryHost);

#endif  // AlpakaDataFormats_CAGeometrySoACollection_h
