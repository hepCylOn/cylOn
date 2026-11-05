#ifndef AlpakaDataFormats_SimpleMapHost_h
#define AlpakaDataFormats_SimpleMapHost_h

#include <cstdint>

#include <alpaka/alpaka.hpp>

#include "Portable/PortableHostCollection.h"
#include "AlpakaDataFormats/SimpleMapSoA.h"
#include "AlpakaCore/config.h"

namespace utils {
  using SimpleMapHost = PortableHostCollection<SimpleMapSoA>;
}
#endif  // AlpakaDataFormats_SimpleMapHost_h
