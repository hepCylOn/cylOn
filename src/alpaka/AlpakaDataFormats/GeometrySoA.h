#ifndef AlpakaDataFormats_GeometrySoA_h
#define AlpakaDataFormats_GeometrySoA_h

#include <cstdint>

#include "AlpakaCore/config.h"
#include "DataFormats/SOARotation.h"

namespace caGeometry
{
  using Frame = SOAFrame<float>;
  using Rotation = SOARotation<float>;
  using GraphNode = std::array<uint32_t, 2>;

  struct GeometrySoA {

    static constexpr uint16_t MAXMODULES = 10000;
    static constexpr uint16_t MAXLAYERS  = 100;
    static constexpr uint16_t MAXPAIRS  = 200;

    //CA Modules
    Frame detFrame[MAXMODULES];
    
    //CA Layers
    int32_t layerStarts[MAXLAYERS];
    float caThetaCut[MAXLAYERS];
    float caDCACut[MAXLAYERS];

    //CA Pairs
    GraphNode graph[MAXPAIRS];
    bool startingPair[MAXPAIRS];
    int16_t phiCuts[MAXPAIRS];
    float minz[MAXPAIRS];
    float maxz[MAXPAIRS];
    float maxr[MAXPAIRS];
};

}


#endif  // AlpakaDataFormats_GeometrySoA_h
