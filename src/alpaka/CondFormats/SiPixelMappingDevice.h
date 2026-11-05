#ifndef CondFormats_SiPixelMappingDevice_h
#define CondFormats_SiPixelMappingDevice_h

#include <cstdint>

#include "Portable/PortableCollection.h"
#include "CondFormats/SiPixelMappingSoA.h"
#include "AlpakaCore/config.h"

template <typename TDev>
using SiPixelMappingDevice = PortableDeviceCollection<SiPixelMappingSoA, TDev>;

#endif  // CondFormats_SiPixelMappingDevice_h
