#ifndef AlpakaDataFormats_alpaka_SiPixelDigisSoACollection_h
#define AlpakaDataFormats_alpaka_SiPixelDigisSoACollection_h

#include <cstdint>

#include <alpaka/alpaka.hpp>

#include "Portable/PortableCollection.h"
#include "AlpakaDataFormats/SiPixelDigisDevice.h"
#include "AlpakaDataFormats/SiPixelDigisHost.h"
#include "AlpakaCore/CopyToHost.h"
#include "AlpakaCore/config.h"
#include "AlpakaCore/AssertDeviceMatchesHostCollection.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  using SiPixelDigisSoACollection =
      std::conditional_t<std::is_same_v<Device, alpaka::DevCpu>, SiPixelDigisHost, SiPixelDigisDevice<Device>>;

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

namespace cms::alpakatools {
  template <typename TDevice>
  struct CopyToHost<SiPixelDigisDevice<TDevice>> {
    template <typename TQueue>
    static auto copyAsync(TQueue &queue, SiPixelDigisDevice<TDevice> const &srcData) {
      SiPixelDigisHost dstData(srcData.view().metadata().size() - 1, queue);
      alpaka::memcpy(queue, dstData.buffer(), srcData.buffer());
      dstData.setNModules(srcData.nModules());
      return dstData;
    }
  };
}  // namespace cms::alpakatools

ASSERT_DEVICE_MATCHES_HOST_COLLECTION(SiPixelDigisSoACollection, SiPixelDigisHost);

#endif  // AlpakaDataFormats_alpaka_SiPixelDigisSoACollection_h
