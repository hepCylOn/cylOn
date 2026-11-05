#ifndef AlpakaDataFormats_ParticleHost_H
#define AlpakaDataFormats_ParticleHost_H

#include "Portable/PortableHostCollection.h"
#include "AlpakaDataFormats/ParticleSoA.h"

namespace sim {
    using ParticleHost = PortableHostCollection<ParticleSoA>;
}

#endif  // AlpakaDataFormats_ParticleHost_H


