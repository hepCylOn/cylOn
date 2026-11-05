#ifndef AlpakaDataFormats_SiPixelDigiErrorsSoA_h
#define AlpakaDataFormats_SiPixelDigiErrorsSoA_h

#include "SoATemplate/SoALayout.h"
#include "DataFormats/PixelErrors.h"
#include "AlpakaCore/SimpleVector.h"

GENERATE_SOA_LAYOUT(SiPixelDigiErrorsLayout, SOA_COLUMN(PixelErrorCompact, pixelErrors), SOA_SCALAR(uint32_t, size))

using SiPixelDigiErrorsSoA = SiPixelDigiErrorsLayout<>;
using SiPixelDigiErrorsSoAView = SiPixelDigiErrorsSoA::View;
using SiPixelDigiErrorsSoAConstView = SiPixelDigiErrorsSoA::ConstView;

#endif  // AlpakaDataFormats_SiPixelDigiErrorsSoA_h