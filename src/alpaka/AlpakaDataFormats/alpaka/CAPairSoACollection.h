#ifndef AlpakaDataFormats_CAPairSoACollection_h
#define AlpakaDataFormats_CAPairSoACollection_h

#include <cstdint>

#include <alpaka/alpaka.hpp>

#include "Portable/PortableCollection.h"
#include "AlpakaDataFormats/CAPairDevice.h"
#include "AlpakaDataFormats/CAPairHost.h"
#include "AlpakaDataFormats/CAPairSoA.h"
#include "AlpakaCore/CopyToHost.h"
#include "AlpakaCore/config.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  using ::caStructures::CAPairDevice;
  using ::caStructures::CAPairHost;
  using CAPairSoACollection =
      std::conditional_t<std::is_same_v<Device, alpaka::DevCpu>, CAPairHost, CAPairDevice<Device>>;

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

ASSERT_DEVICE_MATCHES_HOST_COLLECTION(CAPairSoACollection, ::caStructures::CAPairHost);

#endif  // AlpakaDataFormats_CAPairSoACollection_h
