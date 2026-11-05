#ifndef AlpakaDataFormats_GeometryHost_h
#define AlpakaDataFormats_GeometryHost_h

#include "AlpakaCore/config.h"
#include "AlpakaCore/memory.h"
#include "AlpakaDataFormats/GeometrySoA.h"

using GeometryHost = cms::alpakatools::host_buffer<caGeometry::GeometrySoA>;

#endif  // AlpakaDataFormats_GeometryHost_h
