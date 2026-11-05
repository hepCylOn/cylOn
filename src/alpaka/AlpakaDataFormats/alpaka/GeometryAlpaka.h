#ifndef AlpakaDataFormats_alpaka_GeometryAlpaka_h
#define AlpakaDataFormats_alpaka_GeometryAlpaka_h

#include "AlpakaCore/config.h"
#include "AlpakaCore/memory.h"
#include "AlpakaDataFormats/GeometrySoA.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  using GeometryAlpaka = cms::alpakatools::device_buffer<Device, caGeometry::GeometrySoA>;

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif  // AlpakaDataFormats_alpaka_GeometryAlpaka_h
