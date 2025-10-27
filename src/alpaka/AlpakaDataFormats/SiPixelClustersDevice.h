#ifndef AlpakaDataFormats_SiPixelClustersDevice_h
#define AlpakaDataFormats_SiPixelClustersDevice_h

#include <cstdint>

#include <alpaka/alpaka.hpp>

#include "Framework/Uninitialized.h"
#include "Portable/PortableDeviceCollection.h"
#include "AlpakaDataFormats/SiPixelClustersHost.h"
#include "AlpakaDataFormats/SiPixelClustersSoA.h"
#include "AlpakaCore/CopyToHost.h"
#include "AlpakaCore/config.h"

template <typename TDev>
class SiPixelClustersDevice : public PortableDeviceCollection<SiPixelClustersSoA, TDev> {
public:
  SiPixelClustersDevice(edm::Uninitialized) : PortableDeviceCollection<SiPixelClustersSoA, TDev>{edm::kUninitialized} {}

  template <typename TQueue>
  explicit SiPixelClustersDevice(size_t maxModules, TQueue queue)
      : PortableDeviceCollection<SiPixelClustersSoA, TDev>(maxModules + 1, queue) {}

  // Constructor which specifies the SoA size
  explicit SiPixelClustersDevice(size_t maxModules, TDev const &device)
      : PortableDeviceCollection<SiPixelClustersSoA, TDev>(maxModules + 1, device) {}

  void setNClusters(uint32_t nClusters, int32_t offsetBPIX2) {
    nClusters_h = nClusters;
    offsetBPIX2_h = offsetBPIX2;
  }

  uint32_t nClusters() const { return nClusters_h; }
  int32_t offsetBPIX2() const { return offsetBPIX2_h; }

private:
  uint32_t nClusters_h = 0;
  int32_t offsetBPIX2_h = 0;
};

#endif  // AlpakaDataFormats_SiPixelClustersDevice_h
