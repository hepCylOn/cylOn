#ifndef AlpakaDataFormats_TracksDevice_h
#define AlpakaDataFormats_TracksDevice_h

#include <cstdint>

#include <alpaka/alpaka.hpp>

#include "Framework/Uninitialized.h"
#include "Portable/PortableDeviceCollection.h"
#include "AlpakaDataFormats/TrackDefinitions.h"
#include "AlpakaDataFormats/TracksSoA.h"

namespace reco {
  template <typename TDev>
  using TracksDevice = PortableDeviceMultiCollection<TDev, TrackSoA, TrackHitSoA>;
}

#endif  // AlpakaDataFormats_TracksDevice_h
