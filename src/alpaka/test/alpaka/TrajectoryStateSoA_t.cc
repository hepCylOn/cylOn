/* Simple test for the copyFromDense and copyToDense utilities from AlpakaDataFormats/alpaka/TrackUtilities.h .
 *
 * Creates an instance of TracksSoACollection<pixelTopology::Phase1> (automatically allocates memory on device),
 * passes the view of the SoA data to the kernel that:
 *   - fill the SoA with covariance data;
 *   - copy the covariance data to the dense representation, and back to the matrix representation;
 *   - verify that the data is copied back and forth correctly.
 */

#include <cstdlib>
#include <iostream>

#include <alpaka/alpaka.hpp>

#include "AlpakaDataFormats/TracksSoA.h"
#include "AlpakaDataFormats/alpaka/TracksSoACollection.h"
#include "Framework/stringize.h"
#include "Geometry/SimplePixelTopology.h"
#include "AlpakaCore/config.h"
#include "AlpakaCore/alpaka/devices.h"


/// TODO: move back to common header + dev.cc
using Vector5d = Eigen::Matrix<double, 5, 1>;
using Matrix5d = Eigen::Matrix<double, 5, 5>;

using namespace cms::alpakatools;

namespace ALPAKA_ACCELERATOR_NAMESPACE::test {

  namespace {

    ALPAKA_FN_ACC Matrix5d buildCovariance(Vector5d const& e) {
      Matrix5d cov;
      for (int i = 0; i < 5; ++i)
        cov(i, i) = e(i) * e(i);
      for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < i; ++j) {
          // this makes the matrix positive defined
          double v = 0.3 * std::sqrt(cov(i, i) * cov(j, j));
          cov(i, j) = (i + j) % 2 ? -0.4 * v : 0.1 * v;
          cov(j, i) = cov(i, j);
        }
      }
      return cov;
    }

    template <typename TrackerTraits>
    struct TestTrackSoA {
      ALPAKA_FN_ACC void operator()(Acc1D const& acc, ::reco::TrackSoAView tracks) const {
        Vector5d par0;
        par0 << 0.2, 0.1, 3.5, 0.8, 0.1;
        Vector5d e0;
        e0 << 0.01, 0.01, 0.035, -0.03, -0.01;
        Matrix5d cov0 = buildCovariance(e0);

        for (auto i : uniform_elements(acc, tracks.metadata().size())) {
          ::reco::copyFromDense(tracks, par0, cov0, i);
          Vector5d par1;
          Matrix5d cov1;
          ::reco::copyToDense(tracks, par1, cov1, i);
          Vector5d deltaV = par1 - par0;
          Matrix5d deltaM = cov1 - cov0;
          for (int j = 0; j < 5; ++j) {
            ALPAKA_ASSERT(std::abs(deltaV(j)) < 1.e-5);
            for (int k = j; k < 5; ++k) {
              ALPAKA_ASSERT(cov0(k, j) == cov0(j, k));
              ALPAKA_ASSERT(cov1(k, j) == cov1(j, k));
              ALPAKA_ASSERT(std::abs(deltaM(k, j)) < 1.e-5);
            }
          }
        }
      }
    };

  }  // namespace


  template <typename TrackerTraits>
  void testTrackSoA(Queue& queue, ::reco::TrackSoAView& tracks) {
    auto grid = make_workdiv<Acc1D>(1, 64);
    alpaka::exec<Acc1D>(queue, grid, TestTrackSoA<TrackerTraits>{}, tracks);
  }

  template void testTrackSoA<pixelTopology::Phase1>(Queue& queue, ::reco::TrackSoAView& tracks);
  template void testTrackSoA<pixelTopology::Phase2>(Queue& queue, ::reco::TrackSoAView& tracks);

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::test

////////////////////////////////////////////////////////////////////////////

// Each test binary is built for a single Alpaka backend.
using namespace ALPAKA_ACCELERATOR_NAMESPACE;
using namespace ALPAKA_ACCELERATOR_NAMESPACE::reco;
int main() {
  // Get the list of devices on the current platform.
  auto const& devices = cms::alpakatools::devices<Platform>();
  if (devices.empty()) {
    std::cerr << "No devices available for the " EDM_STRINGIZE(ALPAKA_ACCELERATOR_NAMESPACE) " backend, "
      "the test will be skipped.\n";
    exit(EXIT_FAILURE);
  }

  // Run the test on each device.
  for (const auto& device : devices) {
    Queue queue(device);

    // Inner scope to deallocate memory before destroying the stream.
    {
      TracksSoACollection tracks_d({{1000, 5000}}, queue);

      test::testTrackSoA<pixelTopology::Phase1>(queue, tracks_d.view());

      // Wait for the tests to complete.
      alpaka::wait(queue);
    }
  }

  std::cout << "All tests passed.\n";
  return EXIT_SUCCESS;
}
