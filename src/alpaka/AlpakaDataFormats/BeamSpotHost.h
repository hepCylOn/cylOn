#ifndef DataFormats_BeamSpot_interface_BeamSpotHost_h
#define DataFormats_BeamSpot_interface_BeamSpotHost_h

#include "AlpakaDataFormats/BeamSpotPOD.h"
#include "Portable/PortableHostObject.h"

// simplified representation of the beamspot data, in host memory
using BeamSpotHost = PortableHostObject<BeamSpotPOD>;

#endif  // DataFormats_BeamSpot_interface_BeamSpotHost_h
