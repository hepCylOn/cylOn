/**
   Simple test for the pixelTrack::TrackSoA data structure
   which inherits from PortableDeviceCollection.

   Creates an instance of the class (automatically allocates
   memory on device), passes the view of the SoA data to
   the kernels which:
   - Fill the SoA with data.
   - Verify that the data written is correct.

   Then, the SoA data are copied back to Host, where
   a temporary host-side view (tmp_view) is created using
   the same Layout to access the data on host and print it.
 */

#include <cstdlib>
#include <unistd.h>

#include <alpaka/alpaka.hpp>

#include "AlpakaDataFormats/TracksDevice.h"
#include "AlpakaDataFormats/TracksHost.h"
#include "AlpakaDataFormats/alpaka/TracksSoACollection.h"
#include "Framework/stringize.h"
#include "Geometry/SimplePixelTopology.h"
#include "AlpakaCore/config.h"
#include "AlpakaCore/alpaka/devices.h"
#include "AlpakaCore/memory.h"
#include "AlpakaCore/workdivisionAdvanced.h"

/// TODO: move back to common header + dev.cc
// #include "TrackSoAHeterogeneous_test.h"

using namespace reco;

using Quality = pixelTrack::Quality;
namespace ALPAKA_ACCELERATOR_NAMESPACE {
  using namespace cms::alpakatools;
  namespace testTrackSoA {

    // Kernel which fills the TrackSoAView with data
    // to test writing to it
    class TestFillKernel {
    public:
      ALPAKA_FN_ACC void operator()(Acc1D const& acc, TrackSoAView tracks_view, int32_t nTracks) const {
        if (cms::alpakatools::once_per_grid(acc)) {
          tracks_view.nTracks() = nTracks;
        }

        for (int32_t j : uniform_elements(acc, nTracks)) {
          tracks_view[j].pt() = (float)j;
          tracks_view[j].eta() = (float)j;
          tracks_view[j].chi2() = (float)j;
          tracks_view[j].quality() = (Quality)(j % 256);
          tracks_view[j].nLayers() = j % 128;
          tracks_view[j].hitOffsets() = j;
        }
      }
    };

    // Kernel which reads from the TrackSoAView to verify
    // that it was written correctly from the fill kernel
    class TestVerifyKernel {
    public:
      ALPAKA_FN_ACC void operator()(Acc1D const& acc, TrackSoAConstView tracks_view, int32_t nTracks) const {
        if (cms::alpakatools::once_per_grid(acc)) {
          ALPAKA_ASSERT(tracks_view.nTracks() == nTracks);
        }
        for (int32_t j : uniform_elements(acc, tracks_view.nTracks())) {
          ALPAKA_ASSERT(abs(tracks_view[j].pt() - (float)j) < .0001);
          ALPAKA_ASSERT(abs(tracks_view[j].eta() - (float)j) < .0001);
          ALPAKA_ASSERT(abs(tracks_view[j].chi2() - (float)j) < .0001);
          ALPAKA_ASSERT(tracks_view[j].quality() == (Quality)(j % 256));
          ALPAKA_ASSERT(tracks_view[j].nLayers() == j % 128);
          ALPAKA_ASSERT(tracks_view[j].hitOffsets() == uint32_t(j));
        }
      }
    };

    // Host function which invokes the two kernels above
    void runKernels(TrackSoAView tracks_view, Queue& queue) {
      int32_t tracks = 420;
      uint32_t items = 64;
      uint32_t groups = divide_up_by(tracks, items);
      auto workDiv = make_workdiv<Acc1D>(groups, items);
      alpaka::exec<Acc1D>(queue, workDiv, TestFillKernel{}, tracks_view, tracks);
      alpaka::exec<Acc1D>(queue, workDiv, TestVerifyKernel{}, tracks_view, tracks);
    }

  }  // namespace testTrackSoA

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

////////////////////////////////////////////////////////////////////////////
// Each test binary is built for a single Alpaka backend.
using namespace ALPAKA_ACCELERATOR_NAMESPACE;
using namespace ALPAKA_ACCELERATOR_NAMESPACE::reco;

int main() {
  // Get the list of devices on the current platform
  auto const& devices = cms::alpakatools::devices<Platform>();
  if (devices.empty()) {
    std::cerr << "No devices available for the " EDM_STRINGIZE(ALPAKA_ACCELERATOR_NAMESPACE) " backend, "
      "the test will be skipped.\n";
    exit(EXIT_FAILURE);
  }

  // Run the test on each device
  for (const auto& device : devices) {
    Queue queue(device);

    // Inner scope to deallocate memory before destroying the stream
    {
      // Instantiate tracks on device. PortableDeviceCollection allocates
      // SoA on device automatically.
      constexpr auto nTracks = 1000;
      constexpr auto nHits = nTracks * 5;

      TracksSoACollection tracks_d({{nTracks, nHits}}, queue);
      testTrackSoA::runKernels(tracks_d.view(), queue);

      // Instantate tracks on host. This is where the data will be
      // copied to from device.
      ::reco::TracksHost tracks_h({{nTracks, nHits}}, queue);

      std::cout << "no. of tracks = " << tracks_h.view().metadata().size() << std::endl;
      alpaka::memcpy(queue, tracks_h.buffer(), tracks_d.const_buffer());
      alpaka::wait(queue);

      // Print results
      std::cout << "pt"
                << "\t"
                << "eta"
                << "\t"
                << "chi2"
                << "\t"
                << "quality"
                << "\t"
                << "nLayers"
                << "\t"
                << "hitIndices off" << std::endl;

      for (int i = 0; i < 10; ++i) {
        std::cout << tracks_h.view()[i].pt() << "\t" << tracks_h.view()[i].eta() << "\t" << tracks_h.view()[i].chi2()
                  << "\t" << (int)tracks_h.view()[i].quality() << "\t" << (int)tracks_h.view()[i].nLayers() << "\t"
                  << tracks_h.view()[i].hitOffsets() << std::endl;
      }
    }
  }

  return EXIT_SUCCESS;
}
