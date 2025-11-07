#ifndef AlpakaDataFormats_SimpleMapSoA_h
#define AlpakaDataFormats_SimpleMapSoA_h

#include <Eigen/Core>

#include <alpaka/alpaka.hpp>

#include "SoATemplate/SoALayout.h"

namespace utils {

  GENERATE_SOA_LAYOUT(SimpleMapLayout, SOA_COLUMN(uint32_t, id))

  using SimpleMapSoA = SimpleMapLayout<>;
  using SimpleMapSoAView = SimpleMapSoA::View;
  using SimpleMapSoAConstView = SimpleMapSoA::ConstView;

}  // namespace utils

#endif  // AlpakaDataFormats_SimpleMapSoA_h
