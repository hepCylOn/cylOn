#ifndef DataFormats_VertexSoA_test_alpaka_VertexSoA_test_h
#define DataFormats_VertexSoA_test_alpaka_VertexSoA_test_h

#include "AlpakaDataFormats/ZVertexSoA.h"
#include "AlpakaCore/config.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::testVertexSoAT {

  void runKernels(reco::VertexSoAView vertex_view, reco::VertexTracksSoAView ztracks_view, Queue& queue);

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::testVertexSoAT

#endif  // DataFormats_VertexSoA_test_alpaka_VertexSoA_test_h
