#ifndef AlpakaDataFormats_alpaka_PixelVertexWorkSpaceSoACollection_h
#define AlpakaDataFormats_alpaka_PixelVertexWorkSpaceSoACollection_h

#include <alpaka/alpaka.hpp>

#include "Portable/PortableCollection.h"
#include "AlpakaCore/config.h"
#include "AlpakaDataFormats/PixelVertexWorkSpaceLayout.h"
#include "AlpakaDataFormats/PixelVertexWorkSpaceSoAHost.h"
#include "AlpakaDataFormats/PixelVertexWorkSpaceSoADevice.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::vertexFinder {
  
  using ::vertexFinder::PixelVertexWorkSpaceSoAHost;
  using ::vertexFinder::PixelVertexWorkSpaceSoADevice;

  using PixelVertexWorkSpaceSoACollection =
      std::conditional_t<std::is_same_v<Device, alpaka::DevCpu>, PixelVertexWorkSpaceSoAHost, PixelVertexWorkSpaceSoADevice<Device>>;

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif  // RecoVertex_PixelVertexFinding_plugins_alpaka_PixelVertexWorkSpaceSoADeviceAlpaka_h
