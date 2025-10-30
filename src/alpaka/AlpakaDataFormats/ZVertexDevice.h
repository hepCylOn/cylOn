#ifndef DataFormats_VertexSoA_interface_ZVertexDevice_h
#define DataFormats_VertexSoA_interface_ZVertexDevice_h

#include <cstdint>

#include <alpaka/alpaka.hpp>

#include "AlpakaDataFormats/ZVertexSoA.h"
#include "AlpakaDataFormats/ZVertexHost.h"
#include "Portable/PortableDeviceCollection.h"

template <typename TDev>
using ZVertexDevice = PortableDeviceMultiCollection<TDev, reco::ZVertexSoA, reco::ZVertexTracksSoA>;

#endif  // DataFormats_VertexSoA_interface_ZVertexDevice_h
