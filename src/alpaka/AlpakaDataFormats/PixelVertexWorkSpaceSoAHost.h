#ifndef RecoVertex_PixelVertexFinding_plugins_PixelVertexWorkSpaceSoAHostAlpaka_h
#define RecoVertex_PixelVertexFinding_plugins_PixelVertexWorkSpaceSoAHostAlpaka_h

#include <alpaka/alpaka.hpp>

#include "Portable/PortableHostCollection.h"
#include "AlpakaDataFormats/PixelVertexWorkSpaceLayout.h"

namespace vertexFinder {

  using PixelVertexWorkSpaceSoAHost = PortableHostCollection<PixelVertexWSSoALayout<>>;

}  // namespace vertexFinder

#endif  // RecoVertex_PixelVertexFinding_plugins_PixelVertexWorkSpaceSoAHostAlpaka_h
