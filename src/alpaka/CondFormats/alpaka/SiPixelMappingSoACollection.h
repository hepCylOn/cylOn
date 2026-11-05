#ifndef CondFormats_alpaka_SiPixelMappingSoA_h
#define CondFormats_alpaka_SiPixelMappingSoA_h

#include "CondFormats/SiPixelMappingSoA.h"
#include "CondFormats/SiPixelMappingHost.h"
#include "CondFormats/SiPixelMappingDevice.h"
#include "Portable/PortableCollection.h"
#include "AlpakaCore/config.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  using SiPixelMappingSoACollection =
      std::conditional_t<std::is_same_v<Device, alpaka::DevCpu>, SiPixelMappingHost, SiPixelMappingDevice<Device>>;

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

ASSERT_DEVICE_MATCHES_HOST_COLLECTION(SiPixelMappingSoACollection, SiPixelMappingHost);

#endif  // CondFormats_alpaka_SiPixelMappingSoA_h
