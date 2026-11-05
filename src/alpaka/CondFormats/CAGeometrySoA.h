#ifndef CondFormats_CAGeometry_h
#define CondFormats_CAGeometrySoA_h

#include <utility>

#include <cstdint>
#include "DataFormats/SOARotation.h"

//TODO: make this a real SoA

namespace caGeometry
{
  using CAModule = SOAFrame<float>;

  struct CAPair
  {
    uint32_t innerLayer;
    uint32_t outerLayer;
    int16_t phiCut;
    float minz;
    float maxz;
    float maxr;
    uint8_t  startingPair;
  };

  struct CALayer
  {
      int32_t layerStarts;
      float caThetaCut;
      float caDCACut;
  };

  constexpr uint16_t MAXMODULES = 10000;
  constexpr uint16_t MAXLAYERS  = 100;
  constexpr uint16_t MAXPAIRS  = 200;
  
  using CAModules = CAModule[MAXMODULES];
  using CALayers = CALayer[MAXLAYERS];
  using CAPairs = CAPair[MAXPAIRS];

  struct CASizes
  {
    uint16_t nModules;
    uint16_t nLayers;
    uint16_t nPairs;
  };

  struct CAGeometrySoA
  {

    uint16_t m_nModules;
    uint16_t m_nLayers;
    uint16_t m_nPairs;

    CAModule  const* __restrict__ m_modules;
    CALayer   const* __restrict__ m_layers;
    CAPair    const* __restrict__ m_pairs;

  };

}

#endif