#ifndef AlpakaDataFormats_alpaka_VertexSoACollection_h
#define AlpakaDataFormats_alpaka_VertexSoACollection_h

#include <cstdint>

#include <alpaka/alpaka.hpp>

#include "Portable/PortableCollection.h"
#include "AlpakaDataFormats/VertexDevice.h"
#include "AlpakaDataFormats/VertexHost.h"
#include "AlpakaDataFormats/VertexSoA.h"
#include "AlpakaCore/CopyToHost.h"
#include "AlpakaCore/config.h"
#include "AlpakaCore/AssertDeviceMatchesHostCollection.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  
  using VertexSoACollection =
      std::conditional_t<std::is_same_v<Device, alpaka::DevCpu>, VertexHost, VertexDevice<Device>>;

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

ASSERT_DEVICE_MATCHES_HOST_COLLECTION(VertexSoACollection, VertexHost);

#endif  // AlpakaDataFormats_alpaka_VertexSoACollection_h
