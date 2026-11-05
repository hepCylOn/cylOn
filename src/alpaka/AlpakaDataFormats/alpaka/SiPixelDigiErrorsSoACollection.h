#ifndef AlpakaDataFormats_alpaka_SiPixelDigiErrorsSoACollection_h
#define AlpakaDataFormats_alpaka_SiPixelDigiErrorsSoACollection_h

#include <cstdint>

#include <alpaka/alpaka.hpp>

#include "Portable/PortableCollection.h"
#include "AlpakaDataFormats/SiPixelDigiErrorsHost.h"
#include "AlpakaDataFormats/SiPixelDigiErrorsDevice.h"
#include "AlpakaCore/config.h"
#include "AlpakaCore/memory.h"
#include "AlpakaCore/CopyToHost.h"
#include "AlpakaCore/AssertDeviceMatchesHostCollection.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  using SiPixelDigiErrorsSoACollection =
      std::conditional_t<std::is_same_v<Device, alpaka::DevCpu>, SiPixelDigiErrorsHost, SiPixelDigiErrorsDevice<Device>>;

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

namespace cms::alpakatools {
  template <typename TDevice>
  struct CopyToHost<SiPixelDigiErrorsDevice<TDevice>> {
    template <typename TQueue>
    static auto copyAsync(TQueue& queue, SiPixelDigiErrorsDevice<TDevice> const& srcData) {
      SiPixelDigiErrorsHost dstData(srcData.maxFedWords(), queue);
      alpaka::memcpy(queue, dstData.buffer(), srcData.buffer());
#ifdef GPU_DEBUG
      printf("SiPixelDigiErrorsSoACollection: I'm copying to host.\n");
#endif
      return dstData;
    }
  };
}  // namespace cms::alpakatools

ASSERT_DEVICE_MATCHES_HOST_COLLECTION(SiPixelDigiErrorsSoACollection, SiPixelDigiErrorsHost);

#endif  // AlpakaDataFormats_alpaka_SiPixelDigiErrorsSoACollection_h