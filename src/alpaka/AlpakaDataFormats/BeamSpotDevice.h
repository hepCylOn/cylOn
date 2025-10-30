#ifndef AlpakaDataFormats_BeamSpotDevice_h
#define AlpakaDataFormats_BeamSpotDevice_h

#include "AlpakaDataFormats/BeamSpotHost.h"
#include "AlpakaDataFormats/BeamSpotPOD.h"
#include "Portable/PortableDeviceObject.h"
#include "AlpakaCore/config.h"

template <typename TDev>
using BeamSpotDevice = PortableDeviceObject<BeamSpotPOD, TDev>;

#endif  // AlpakaDataFormats_BeamSpotDevice_h