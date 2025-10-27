#ifndef CondFormats_SiPixelGainCalibrationForHLTDevice_h
#define CondFormats_SiPixelGainCalibrationForHLTDevice_h

#include "CondFormats/SiPixelGainCalibrationForHLTSoA.h"
#include "Portable/PortableCollection.h"
#include "AlpakaCore/config.h"

template <typename TDev>
using SiPixelGainCalibrationForHLTDevice = PortableDeviceCollection<SiPixelGainCalibrationForHLTSoA, TDev>;

#endif  // CondFormats_alpaka_SiPixelGainCalibrationForHLTDevice_h
