#ifndef AlpakaDataFormats_alpaka_VertexSoACollection_h
#define AlpakaDataFormats_alpaka_VertexSoACollection_h

#include <cstdint>

#include <alpaka/alpaka.hpp>

#include "Portable/PortableCollection.h"
#include "AlpakaDataFormats/ZVertexDevice.h"
#include "AlpakaDataFormats/ZVertexHost.h"
#include "AlpakaDataFormats/ZVertexSoA.h"
#include "AlpakaCore/CopyToHost.h"
#include "AlpakaCore/config.h"
#include "AlpakaCore/AssertDeviceMatchesHostCollection.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  
  using ZVertexSoACollection =
      std::conditional_t<std::is_same_v<Device, alpaka::DevCpu>, ZVertexHost, ZVertexDevice<Device>>;

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

ASSERT_DEVICE_MATCHES_HOST_COLLECTION(ZVertexSoACollection, ZVertexHost);

#endif  // AlpakaDataFormats_alpaka_VertexSoACollection_h
