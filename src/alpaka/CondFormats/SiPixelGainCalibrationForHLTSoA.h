#ifndef CondFormats_interface_SiPixelGainCalibrationForHLTSoA_h
#define CondFormats_interface_SiPixelGainCalibrationForHLTSoA_h

#include <array>
#include "SoATemplate/SoALayout.h"
#include "Geometry/SimplePixelTopology.h"

namespace siPixelGainsSoA {
  struct DecodingStructure {
    uint8_t gain;
    uint8_t ped;
  };

  using Ranges = std::array<uint32_t, phase1PixelTopology::numberOfModules>;
  using Cols = std::array<int, phase1PixelTopology::numberOfModules>;
}  // namespace siPixelGainsSoA

GENERATE_SOA_LAYOUT(SiPixelGainCalibrationForHLTLayout,
                    SOA_COLUMN(siPixelGainsSoA::DecodingStructure, v_pedestals),

                    SOA_SCALAR(siPixelGainsSoA::Ranges, modStarts),
                    SOA_SCALAR(siPixelGainsSoA::Ranges, modEnds),
                    SOA_SCALAR(siPixelGainsSoA::Cols, modCols),

                    SOA_SCALAR(float, minPed),
                    SOA_SCALAR(float, maxPed),
                    SOA_SCALAR(float, minGain),
                    SOA_SCALAR(float, maxGain),
                    SOA_SCALAR(float, pedPrecision),
                    SOA_SCALAR(float, gainPrecision),

                    SOA_SCALAR(unsigned int, numberOfRowsAveragedOver),
                    SOA_SCALAR(unsigned int, nBinsToUseForEncoding),
                    SOA_SCALAR(unsigned int, deadFlag),
                    SOA_SCALAR(unsigned int, noisyFlag),
                    SOA_SCALAR(float, link))

using SiPixelGainCalibrationForHLTSoA = SiPixelGainCalibrationForHLTLayout<>;
using SiPixelGainCalibrationForHLTSoAView = SiPixelGainCalibrationForHLTSoA::View;
using SiPixelGainCalibrationForHLTSoAConstView = SiPixelGainCalibrationForHLTSoA::ConstView;

#endif  // CondFormats_interface_SiPixelGainCalibrationForHLTSoA_h
