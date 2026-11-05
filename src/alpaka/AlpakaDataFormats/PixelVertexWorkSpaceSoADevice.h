#ifndef RecoVertex_PixelVertexFinding_plugins_alpaka_PixelVertexWorkSpaceSoADeviceAlpaka_h
#define RecoVertex_PixelVertexFinding_plugins_alpaka_PixelVertexWorkSpaceSoADeviceAlpaka_h

#include <alpaka/alpaka.hpp>

#include "Portable/PortableCollection.h"
#include "AlpakaCore/config.h"
#include "AlpakaDataFormats/PixelVertexWorkSpaceLayout.h"


namespace vertexFinder {
  
  template <typename TDev>
  using PixelVertexWorkSpaceSoADevice = PortableDeviceCollection<::vertexFinder::PixelVertexWSSoALayout<>, TDev>;

}  // namespace vertexFinder




#endif  // RecoVertex_PixelVertexFinding_plugins_alpaka_PixelVertexWorkSpaceSoADeviceAlpaka_h
