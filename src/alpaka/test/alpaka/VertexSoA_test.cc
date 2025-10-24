// /**
//    Simple test for the reco::VertexSoA data structure
//    which inherits from Portable{Host}Collection.

//    Creates an instance of the class (automatically allocates
//    memory on device), passes the view of the SoA data to
//    the kernels which:
//    - Fill the SoA with data.
//    - Verify that the data written is correct.

//    Then, the SoA data are copied back to Host, where
//    a temporary host-side view (tmp_view) is created using
//    the same Layout to access the data on host and print it.
//  */

// #include <cstdlib>
// #include <unistd.h>

// #include <alpaka/alpaka.hpp>

// #include "AlpakaDataFormats/VertexHost.h"
// #include "AlpakaDataFormats/alpaka/VertexSoACollection.h"
// #include "Framework/stringize.h"
// #include "AlpakaCore/config.h"
// #include "AlpakaCore/alpaka/devices.h"
// #include "AlpakaCore/memory.h"
// #include "AlpakaCore/workdivisionAdvanced.h"

// #include "Portable/PortableCollection.h"
// #include "VertexSoA_test.h"

// using namespace ALPAKA_ACCELERATOR_NAMESPACE;

// // Run 3 values, used for testing
// constexpr uint32_t maxTracks = 32 * 1024;
// constexpr uint32_t maxVertices = 1024;

// int main() {
//   // Get the list of devices on the current platform
//   auto const& devices = cms::alpakatools::devices<Platform>();
//   if (devices.empty()) {
//     std::cerr << "No devices available for the " EDM_STRINGIZE(ALPAKA_ACCELERATOR_NAMESPACE) " backend, "
//       "the test will be skipped.\n";
//     exit(EXIT_FAILURE);
//   }

//   // Run the test on each device
//   for (const auto& device : devices) {
//     Queue queue(device);

//     // Inner scope to deallocate memory before destroying the stream
//     {
//       // Instantiate vertices on device. PortableCollection allocates
//       // SoA on device automatically.
//       VertexSoACollection Vertex_d({{maxTracks, maxVertices}}, queue);
//       testVertexSoAT::runKernels(Vertex_d.view(), Vertex_d.view<reco::VertexTracksSoA>(), queue);

//       // If the device is actually the host, use the collection as-is.
//       // Otherwise, copy the data from the device to the host.
// #ifdef ALPAKA_ACC_CPU_B_SEQ_T_SEQ_ENABLED //TODO this was ALPAKA_ACC_CPU_B_SEQ_T_SEQ_ENABLED
//       VertexHost Vertex_h = std::move(Vertex_d);
// #else
//       VertexHost Vertex_h = cms::alpakatools::CopyToHost<VertexSoACollection>::copyAsync(queue, Vertex_d);
// #endif
//       alpaka::wait(queue);
//       std::cout << Vertex_h.view().metadata().size() << std::endl;

//       // Print results
//       std::cout << "idv\t"
//                 << "zv\t"
//                 << "wv\t"
//                 << "chi2\t"
//                 << "ptv2\t"
//                 << "ndof\t"
//                 << "sortInd\t"
//                 << "nvFinal\n";

//       auto vtx_v = Vertex_h.view<reco::VertexSoA>();
//       auto trk_v = Vertex_h.view<reco::VertexTracksSoA>();
//       for (int i = 0; i < 10; ++i) {
//         auto vi = vtx_v[i];
//         auto ti = trk_v[i];
//         std::cout << (int)ti.idv() << "\t" << vi.zv() << "\t" << vi.wv() << "\t" << vi.chi2() << "\t" << vi.ptv2()
//                   << "\t" << (int)ti.ndof() << "\t" << vi.sortInd() << "\t" << (int)vtx_v.nvFinal() << std::endl;
//       }
//     }
//   }

//   return EXIT_SUCCESS;
// }
