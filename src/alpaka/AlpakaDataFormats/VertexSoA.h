#ifndef AlpakaDataFormats_VertexSoA_h
#define AlpakaDataFormats_VertexSoA_h

#include <alpaka/alpaka.hpp>

#include <Eigen/Core>

#include "SoATemplate/SoALayout.h"

namespace reco {

  GENERATE_SOA_LAYOUT(VertexLayout,
                      SOA_COLUMN(float, zv),          // output z-posistion of found vertices
                      SOA_COLUMN(float, wv),          // output weight (1/error^2) on the above
                      SOA_COLUMN(float, chi2),        // vertices chi2
                      SOA_COLUMN(float, ptv2),        // vertices pt^2
                      SOA_COLUMN(uint16_t, sortInd),  // sorted index (by pt2)  ascending
                      SOA_SCALAR(uint32_t, nvFinal))  // the number of vertices

  GENERATE_SOA_LAYOUT(VertexTracksLayout,
                      SOA_COLUMN(int16_t, idv),   // vertex index for each associated (original) track
                                                  // (-1 == not associate)
                      SOA_COLUMN(int32_t, ndof))  // vertices number of dof
                                                  // FIXME: reused as workspace for the number of nearest neighbours

  // Common types for both Host and Device code
  using VertexSoA = VertexLayout<>;
  using VertexSoAView = VertexSoA::View;
  using VertexSoAConstView = VertexSoA::ConstView;

  // Common types for both Host and Device code
  using VertexTracksSoA = VertexTracksLayout<>;
  using VertexTracksSoAView = VertexTracksSoA::View;
  using VertexTracksSoAConstView = VertexTracksSoA::ConstView;

  ALPAKA_FN_HOST_ACC ALPAKA_FN_INLINE void init(VertexSoAView &vertices) { vertices.nvFinal() = 0; }

}  // namespace reco

#endif  // DataFormats_VertexSoA_interface_VertexSoA_h
