#ifndef AlpakaDataFormats_CAPairSoA_h
#define AlpakaDataFormats_CAPairSoA_h

#include <Eigen/Core>

#include <alpaka/alpaka.hpp>

#include "SoATemplate/SoALayout.h"

///TODO: since we have started to use this also outside of the CA maybe is worth a renaming
namespace caStructures {

  GENERATE_SOA_LAYOUT(CAPairLayout, SOA_COLUMN(uint32_t, inner), SOA_COLUMN(uint32_t, outer))

  using CAPairSoA = CAPairLayout<>;
  using CAPairSoAView = CAPairSoA::View;
  using CAPairSoAConstView = CAPairSoA::ConstView;

}  // namespace caStructures

#endif  // AlpakaDataFormats_CAPairSoA_h
