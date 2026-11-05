#ifndef AlpakaDataFormats_CAPairHost_h
#define AlpakaDataFormats_CAPairHost_h

#include <cstdint>

#include <alpaka/alpaka.hpp>

#include "Portable/PortableHostCollection.h"
#include "AlpakaDataFormats/CAPairSoA.h"
#include "AlpakaCore/config.h"

namespace caStructures {
  using CAPairHost = PortableHostCollection<CAPairSoA>;
}
#endif  // AlpakaDataFormats_CAPairHost_h
