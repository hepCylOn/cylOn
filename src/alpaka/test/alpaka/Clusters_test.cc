#include <cstdlib>

#include <alpaka/alpaka.hpp>

#include "AlpakaDataFormats/SiPixelClustersDevice.h"
#include "AlpakaDataFormats/SiPixelClustersHost.h"
#include "AlpakaDataFormats/SiPixelClustersSoA.h"
#include "AlpakaDataFormats/alpaka/SiPixelClustersSoACollection.h"

#include "Framework/stringize.h"
#include "AlpakaCore/config.h"
#include "AlpakaCore/alpaka/devices.h"
#include "AlpakaCore/memory.h"
#include "AlpakaCore/workdivisionAdvanced.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
/// TODO: This used to be a separte dev.cc + header  
// #include "Clusters_test.h"


namespace ALPAKA_ACCELERATOR_NAMESPACE::testClusterSoA {

  class TestFillKernel {
  public:
    ALPAKA_FN_ACC void operator()(Acc1D const& acc, SiPixelClustersSoAView clust_view) const {
      for (int32_t j : cms::alpakatools::uniform_elements(acc, clust_view.metadata().size())) {
        clust_view[j].moduleStart() = j;
        clust_view[j].clusInModule() = j * 2;
        clust_view[j].moduleId() = j * 3;
        clust_view[j].clusModuleStart() = j * 4;
      }
    }
  };

  class TestVerifyKernel {
  public:
    ALPAKA_FN_ACC void operator()(Acc1D const& acc, SiPixelClustersSoAConstView clust_view) const {
      for (uint32_t j : cms::alpakatools::uniform_elements(acc, clust_view.metadata().size())) {
        ALPAKA_ASSERT_ACC(clust_view[j].moduleStart() == j);
        ALPAKA_ASSERT_ACC(clust_view[j].clusInModule() == j * 2);
        ALPAKA_ASSERT_ACC(clust_view[j].moduleId() == j * 3);
        ALPAKA_ASSERT_ACC(clust_view[j].clusModuleStart() == j * 4);
      }
    }
  };

  void runKernels(SiPixelClustersSoAView clust_view, Queue& queue) {
    uint32_t items = 64;
    uint32_t groups = cms::alpakatools::divide_up_by(clust_view.metadata().size(), items);
    auto workDiv = cms::alpakatools::make_workdiv<Acc1D>(groups, items);
    alpaka::exec<Acc1D>(queue, workDiv, TestFillKernel{}, clust_view);
    alpaka::exec<Acc1D>(queue, workDiv, TestVerifyKernel{}, clust_view);
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::testClusterSoA

//////////////////////////////////////////////////////////////////////////////////////////////////

using namespace ALPAKA_ACCELERATOR_NAMESPACE;

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
      SiPixelClustersSoACollection clusters_d(100, queue);
      testClusterSoA::runKernels(clusters_d.view(), queue);

      // Instantate tracks on host. This is where the data will be
      // copied to from device.
      SiPixelClustersHost clusters_h(clusters_d.view().metadata().size(), queue);

      std::cout << clusters_h.view().metadata().size() << std::endl;
      alpaka::memcpy(queue, clusters_h.buffer(), clusters_d.const_buffer());
      alpaka::wait(queue);
    }
  }

  return EXIT_SUCCESS;
}
