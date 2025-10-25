// #include <alpaka/alpaka.hpp>

// #include "AlpakaDataFormats/VertexDevice.h"
// #include "AlpakaDataFormats/VertexHost.h"
// #include "AlpakaDataFormats/alpaka/VertexSoACollection.h"
// #include "AlpakaCore/config.h"
// #include "AlpakaCore/workdivisionAdvanced.h"

// namespace ALPAKA_ACCELERATOR_NAMESPACE::testVertexSoAT {

//   class TestFillKernel {
//   public:
//     ALPAKA_FN_ACC void operator()(Acc1D const& acc,
//                                   reco::VertexSoAView vertex_view,
//                                   reco::VertexTracksSoAView ztracks_view) const {
//       if (cms::alpakatools::once_per_grid(acc)) {
//         vertex_view.nvFinal() = 420;
//       }

//       for (int32_t j : cms::alpakatools::uniform_elements(acc, vertex_view.metadata().size())) {
//         vertex_view[j].zv() = (float)j;
//         vertex_view[j].wv() = (float)j;
//         vertex_view[j].chi2() = (float)j;
//         vertex_view[j].ptv2() = (float)j;
//         vertex_view[j].sortInd() = (uint16_t)j;
//       }
//       for (int32_t j : cms::alpakatools::uniform_elements(acc, ztracks_view.metadata().size())) {
//         ztracks_view[j].idv() = (int16_t)j;
//         ztracks_view[j].ndof() = (int32_t)j;
//       }
//     }
//   };

//   class TestVerifyKernel {
//   public:
//     ALPAKA_FN_ACC void operator()(Acc1D const& acc,
//                                   reco::VertexSoAView vertex_view,
//                                   reco::VertexTracksSoAView ztracks_view) const {
//       if (cms::alpakatools::once_per_grid(acc)) {
//         ALPAKA_ASSERT_ACC(vertex_view.nvFinal() == 420);
//       }

//       for (int32_t j : cms::alpakatools::uniform_elements(acc, vertex_view.nvFinal())) {
//         ALPAKA_ASSERT(vertex_view[j].zv() - (float)j < 0.0001);
//         ALPAKA_ASSERT(vertex_view[j].wv() - (float)j < 0.0001);
//         ALPAKA_ASSERT(vertex_view[j].chi2() - (float)j < 0.0001);
//         ALPAKA_ASSERT(vertex_view[j].ptv2() - (float)j < 0.0001);
//         ALPAKA_ASSERT(vertex_view[j].sortInd() == uint32_t(j));
//       }
//       for (int32_t j : cms::alpakatools::uniform_elements(acc, ztracks_view.metadata().size())) {
//         ALPAKA_ASSERT(ztracks_view[j].idv() == j);
//         ALPAKA_ASSERT(ztracks_view[j].ndof() == j);
//       }
//     }
//   };

//   void runKernels(reco::VertexSoAView vertex_view, reco::VertexTracksSoAView ztracks_view, Queue& queue) {
//     uint32_t items = 64;
//     uint32_t groups = cms::alpakatools::divide_up_by(vertex_view.metadata().size(), items);
//     auto workDiv = cms::alpakatools::make_workdiv<Acc1D>(groups, items);
//     alpaka::exec<Acc1D>(queue, workDiv, TestFillKernel{}, vertex_view, ztracks_view);
//     alpaka::exec<Acc1D>(queue, workDiv, TestVerifyKernel{}, vertex_view, ztracks_view);
//   }

// }  // namespace ALPAKA_ACCELERATOR_NAMESPACE::testVertexSoAT
