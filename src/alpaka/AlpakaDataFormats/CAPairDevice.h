#ifndef AlpakaDataFormats_CAPairDevice_H
#define AlpakaDataFormats_CAPairDevice_H

#include <cstdint>

#include <alpaka/alpaka.hpp>

#include "Portable/PortableDeviceCollection.h"
#include "AlpakaDataFormats/CAPairSoA.h" // move me somewhere else  
#include "AlpakaCore/config.h"

namespace caStructures {
  template <typename TDev>
  using CAPairDevice = PortableDeviceCollection<CAPairSoA, TDev>;
}

#endif  // AlpakaDataFormats_CAPairDevice_H
