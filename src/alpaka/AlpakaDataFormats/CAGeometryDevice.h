#ifndef AlpakaDataFormats_CAGeometryDevice_H
#define AlpakaDataFormats_CAGeometryDevice_H

#include <cstdint>

#include <alpaka/alpaka.hpp>

#include "Portable/PortableDeviceCollection.h"
#include "AlpakaDataFormats/CAGeometrySoA.h"
#include "AlpakaCore/config.h"

namespace reco {
  template <typename TDev>
  using CAGeometryDevice = PortableDeviceMultiCollection<TDev, CALayersSoA, CAGraphSoA, CAModulesSoA>;
}
#endif  // AlpakaDataFormats_CAGeometryDevice_H
