#ifndef AlpakaDataFormats_CAGeometryHost_H
#define AlpakaDataFormats_CAGeometryHost_H

#include <cstdint>

#include <alpaka/alpaka.hpp>

#include "Portable/PortableHostCollection.h"
#include "AlpakaDataFormats/CAGeometrySoA.h"
#include "AlpakaCore/config.h"

namespace reco {
  using CAGeometryHost = PortableHostMultiCollection<CALayersSoA, CAGraphSoA, CAModulesSoA>;
}
#endif  // AlpakaDataFormats_CAGeometryHost_H
