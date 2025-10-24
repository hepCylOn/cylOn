#ifndef DataFormats_VertexSoA_interface_ZVertexDevice_h
#define DataFormats_VertexSoA_interface_ZVertexDevice_h

#include <cstdint>

#include <alpaka/alpaka.hpp>

#include "AlpakaDataFormats/VertexSoA.h"
#include "AlpakaDataFormats/VertexHost.h"
#include "Portable/PortableDeviceCollection.h"

template <typename TDev>
using VertexDevice = PortableDeviceMultiCollection<TDev, reco::VertexSoA, reco::VertexTracksSoA>;

#endif  // DataFormats_VertexSoA_interface_ZVertexDevice_h
