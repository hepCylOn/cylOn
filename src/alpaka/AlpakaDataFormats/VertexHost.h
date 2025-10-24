#ifndef AlpakaDataFormats_VertexHost_H
#define AlpakaDataFormats_VertexHost_H

#include <cstdint>

#include <alpaka/alpaka.hpp>

#include "Portable/PortableHostCollection.h"
#include "AlpakaDataFormats/VertexSoA.h"
#include "AlpakaCore/config.h"

using VertexHost = PortableHostCollection2<reco::VertexSoA, reco::VertexTracksSoA>;

#endif  // AlpakaDataFormats_VertexHost_H
