#ifndef CondFormats_alpaka_SiPixelGainCalibrationForHLTSoACollection_h
#define CondFormats_alpaka_SiPixelGainCalibrationForHLTSoACollection_h

#include "CondFormats/SiPixelGainCalibrationForHLTSoA.h"
#include "CondFormats/SiPixelGainCalibrationForHLTHost.h"
#include "CondFormats/SiPixelGainCalibrationForHLTDevice.h"
#include "Portable/PortableCollection.h"
#include "AlpakaCore/config.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  using SiPixelGainCalibrationForHLTSoACollection =
      std::conditional_t<std::is_same_v<Device, alpaka::DevCpu>, SiPixelGainCalibrationForHLTHost, SiPixelGainCalibrationForHLTDevice<Device>>;

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

ASSERT_DEVICE_MATCHES_HOST_COLLECTION(SiPixelGainCalibrationForHLTSoACollection, SiPixelGainCalibrationForHLTHost);

#endif  // CondFormats_alpaka_SiPixelGainCalibrationForHLTSoACollection_h
