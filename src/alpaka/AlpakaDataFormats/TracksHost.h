#ifndef AlpakaDataFormats_TracksHost_H
#define AlpakaDataFormats_TracksHost_H

#include "Portable/PortableHostCollection.h"
#include "AlpakaDataFormats/TracksSoA.h"

namespace reco {
  using TracksHost = PortableHostMultiCollection<TrackSoA, TrackHitSoA>;
}

#endif  // AlpakaDataFormats_TracksHost_H
