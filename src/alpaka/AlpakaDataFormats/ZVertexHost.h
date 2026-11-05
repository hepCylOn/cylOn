#ifndef AlpakaDataFormats_VertexHost_H
#define AlpakaDataFormats_VertexHost_H

#include <cstdint>

#include <alpaka/alpaka.hpp>

#include "Portable/PortableHostCollection.h"
#include "AlpakaDataFormats/ZVertexSoA.h"
#include "AlpakaCore/config.h"

using ZVertexHost = PortableHostCollection2<reco::ZVertexSoA, reco::ZVertexTracksSoA>;

#endif  // AlpakaDataFormats_VertexHost_H
