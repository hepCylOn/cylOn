#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>
#include <cassert>
#include <filesystem>

#include "DataFormats/SOARotation.h"
#include "CondFormats/CAGeometrySoA.h"

using namespace caGeometry;
using Frame = SOAFrame<float>;
namespace fs = std::filesystem;

struct DetParams {
    bool isBarrel;
    bool isPosZ;
    uint16_t layer;
    uint16_t index;
    uint32_t rawId;

    float shiftX;
    float shiftY;
    float chargeWidthX;
    float chargeWidthY;
    // CMSSW 11.2.x adds
    //uint16_t pixmx;  // max pix charge
    // which would break reading the binary dumps

    float x0, y0, z0;  // the vertex in the local coord of the detector

    float sx[3], sy[3];  // the errors...

    Frame frame;
  };

// all modules are identical!
struct CommonParams {
    float theThicknessB;
    float theThicknessE;
    float thePitchX;
    float thePitchY;
  };
  
class PixelCPEFast {
  public:
    PixelCPEFast(std::string const &inputPath, std::string const &outputPath){
      std::ifstream in(inputPath, std::ios::binary);
      in.exceptions(std::ifstream::badbit | std::ifstream::failbit | std::ifstream::eofbit);
      in.read(reinterpret_cast<char *>(&m_commonParamsGPU), sizeof(CommonParams));
      unsigned int ndetParams;
      in.read(reinterpret_cast<char *>(&ndetParams), sizeof(unsigned int));
      m_detParamsGPU.resize(ndetParams);
      in.read(reinterpret_cast<char *>(m_detParamsGPU.data()), ndetParams * sizeof(DetParams));

      std::ofstream out(outputPath + "/AdHocCMSPhase1Modules.bin", std::ios::binary);
      if (!out) throw std::runtime_error("Cannot open file for writing");

      std::cout << "=== Writing " <<  ndetParams << " modules" << std::endl;
      for(auto const &v: m_detParamsGPU)
      {
        auto f = v.frame;
        out.write(reinterpret_cast<const char*>(&f), sizeof(f));
      }

      out.close();
      std::cout << "AdHocCMSPhase1Modules.bin written.\n";

    }

    ~PixelCPEFast() = default;

      private:
    // allocate it with posix malloc to be compatible with cpu wf
    std::vector<DetParams> m_detParamsGPU;
    CommonParams m_commonParamsGPU;


};

// --- Constants ---
constexpr uint16_t nModules = 1856;
constexpr uint16_t nLayers  = 18;
constexpr uint16_t nPairs   = 41;

// phi cut constants
// constexpr int16_t phi0p01 = 106;
// constexpr int16_t phi0p02 = 210;
// constexpr int16_t phi0p03 = 314;
// constexpr int16_t phi0p04 = 418;
constexpr int16_t phi0p05 = 522;
constexpr int16_t phi0p06 = 626;
constexpr int16_t phi0p07 = 730;

// Phi cuts for layer pairs
constexpr int16_t phicuts[nPairs]{

    //0,  1    0,  4    0,  11
      phi0p05, phi0p05, phi0p05, 
    //1,  2    1,  4    1,  11
      phi0p06, phi0p07, phi0p07,
    //2,  3    2,  4    2,  11
      phi0p06, phi0p07, phi0p07, 
      
    //4,  5    5,  6    6,  7    7,  8    8,  9    9,  10
      phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, 
    //11, 12   12, 13   13, 14   14, 15   15, 16   16, 17
      phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05,
      
    //0,  2    0,  5    0,  12   0,  6    0,  13
      phi0p05, phi0p07, phi0p07, phi0p07, phi0p07,
    //1,  3    1,  5    1,  12   1,  6    1,  13 
      phi0p05, phi0p07, phi0p07, phi0p07, phi0p07, 
      
    //4,  6    5,  7    6,  8    7,  9    8,  10
      phi0p07, phi0p07, phi0p07, phi0p07, phi0p07, 
    //11, 13   12, 14   13, 15   14, 16   15, 17
      phi0p07, phi0p07, phi0p07, phi0p07, phi0p07

    // //0,  1    0,  4    0,  11
    //   phi0p03, phi0p03, phi0p03, 
    // //1,  2    1,  4    1,  11
    //   phi0p03, phi0p04, phi0p04,
    // //2,  3    2,  4    2,  11
    //   phi0p04, phi0p04, phi0p04, 
      
    // //4,  5    5,  6    6,  7    7,  8    8,  9    9,  10
    //   phi0p02, phi0p02, phi0p02, phi0p02, phi0p02, phi0p02, 
    // //11, 12   12, 13   13, 14   14, 15   15, 16   16, 17
    //   phi0p02, phi0p02, phi0p02, phi0p02, phi0p02, phi0p02,
      
    // //0,  2    0,  5    0,  12   0,  6    0,  13
    //   phi0p02, phi0p02, phi0p02, phi0p02, phi0p02,
    // //1,  3    1,  5    1,  12   1,  6    1,  13 
    //   phi0p02, phi0p02, phi0p02, phi0p02, phi0p02, 
      
    // //4,  6    5,  7    6,  8    7,  9    8,  10
    //   phi0p02, phi0p02, phi0p02, phi0p02, phi0p02, 
    // //11, 13   12, 14   13, 15   14, 16   15, 17
    //   phi0p02, phi0p02, phi0p02, phi0p02, phi0p02

};

// Min z for layer pairs
constexpr float minz_vals[nPairs] = {

    //0,  1     0,  4     0,  11
      -34.0,    10.0,      2.0*-48.0, 
    //1,  2     1,  4     1,  11      
      2.0*-44.0,    1.0*6.0,      2.0*-72.0, 
    //2,  3     2,  4     2,  11      
      2.0*-54.0,    1.0*11.0,     2.0*-96.0,
      
    //4,  5     5,  6     6,  7    7,  8    8,  9    9,  10      
      1.0*23.0,     1.0*30.0,     1.0*39.0,    1.0*50.0,    1.0*65.0,    1.0*82.0, 
    //11, 12    12, 13    13, 14   14, 15   15, 16   16, 17      
      2.0*-84.0,    2.0*-105.0,   2.0*-132.0,  2.0*-165.0,  2.0*-210.0,  2.0*-327.0, 
      
    //0,  2     0,  5   0,  12    0,  6    0,  13      
      -17.0,    7.0,    -24.0,    11.0,    -24.0,
    //1,  3     1,  5   1,  12    1,  6    1,  13       
      -17.0,    9.0,    -24.0,    13.0,    -24.0, 
      
    //4,  6     5,  7    6,  8    7,  9    8,  10      
      23.0,     30.0,    39.0,    50.0,    65.0, 
    //11, 13    12, 14   13, 15   14, 16   15, 17      
      -84.0,    -105.0,  -132.0,  -165.0,  -210.0

    // //0,  1     0,  4     0,  11
    //   -26.0,    20.0,      -49.0, 
    // //1,  2     1,  4     1,  11      
    //   -31.0,    28.0,      -52.0, 
    // //2,  3     2,  4     2,  11      
    //   -36.0,    38.0,      -52.0,
      
    // //4,  5     5,  6     6,  7    7,  8    8,  9    9,  10      
    //   59.0,     69.0,     80.0,    93.0,    109.0,   128.0, 
    // //11, 12    12, 13    13, 14   14, 15   15, 16   16, 17      
    //   -65.0,    -75.0,    -86.0,   -102.0,  -115.0,  -134.0,
      
    // //0,  2     0,  5   0,  12    0,  6    0,  13      
    //   -29.0,    17.0,    -52.0,    17.0,    -52.0,
    // //1,  3     1,  5   1,  12    1,  6    1,  13       
    //   -33.0,    25.0,    -55.0,    25.0,    -55.0, 
      
    // //4,  6     5,  7    6,  8    7,  9    8,  10      
    //   56.0,     66.0,    77.0,    90.0,    106.0
    // //11, 13    12, 14   13, 15   14, 16   15, 17      
    //   -68.0,    -78.0,   -89.0,   -105.0,  -118.0
    
};

// Max z for layer pairs
constexpr float maxz_vals[nPairs] = {

    
    //0,  1     0,  4     0,  11
      34.0,     50.0,     1.0*-4.0,
    //1,  2     1,  4     1,  11
      2.0*44.0,     2.0*72.0,     1.0*-6.0,
    //2,  3     2,  4     2,  11
      2.0*54.0,     2.0*96.0,     1.0*-11.0,
      
    //4,  5    5,  6    6,  7    7,  8    8,  9    9,  10
      2.0*84.0,    2.0*105.0,   2.0*132.0,   2.0*165.0,   2.0*210.0,   2.0*327.0, 
    //11, 12   12, 13   13, 14   14, 15   15, 16   16, 17
      1.0*-23.0,   1.0*-30.0,   1.0*-39.0,   1.0*-50.0,   1.0*-65.0,   1.0*-82.0,
      
    //0,  2    0,  5    0,  12   0,  6    0,  13 
      17.0,    24.0,    -7.0,    24.0,    -11.0,
    //1,  3    1,  5    1,  12   1,  6    1,  13
      17.0,    24.0,    -9.0,    24.0,    -13.0,
      
    //4,  6    5,  7    6,  8    7,  9    8,  10
      84.0,    105.0,   132.0,   165.0,   210.0,
    //11, 13   12, 14   13, 15   14, 16   15, 17
      -23.0,   -30.0,   -39.0,   -50.0,   -65.0

    // //0,  1     0,  4     0,  11
    //   26.0,     49.0,     -20.0,
    // //1,  2     1,  4     1,  11
    //   31.0,     52.0,     -28.0,
    // //2,  3     2,  4     2,  11
    //   36.0,     52.0,     -38.0,
      
    // //4,  5    5,  6    6,  7    7,  8    8,  9    9,  10
    //   65.0,    75.0,    86.0,    102.0,   115.0,   134.0, 
    // //11, 12   12, 13   13, 14   14, 15   15, 16   16, 17
    //   -59.0,   -69.0,   -80.0,   -93.0,   -109.0,  -128.0, 
      
    // //0,  2    0,  5    0,  12   0,  6    0,  13 
    //   29.0,    52.0,    -17.0,    52.0,    -17.0,
    // //1,  3    1,  5    1,  12   1,  6    1,  13
    //   33.0,    55.0,    -25.0,    55.0,    -25.0,
      
    // //4,  6    5,  7    6,  8    7,  9    8,  10
    //   68.0,    78.0,    89.0,    105.0,   118.0
    // //11, 13   12, 14   13, 15   14, 16   15, 17
    //   -56.0,   -66.0,   -77.0,   -90.0,   -106.0

};

// Max r for layer pairs
constexpr float maxr_vals[nPairs] = {
    
  // //0,  1   0,  4   0,  11
  //   5.0,    6.0,    6.0,
  // //1,  2   1,  4   1,  11
  //   7.0,    8.0,    8.0,
  // //2,  3   2,  4   2,  11
  //   7.0,    8.0,    8.0,

  // //4,  5   5,  6   6,  7   7,  8   8,  9   9,  10
  //   6.0,    6.0,    6.0,    6.0,    6.0,    6.0,
  // //11, 12  12, 13  13, 14  14, 15  15, 16  16, 17
  //   6.0,    6.0,    6.0,    6.0,    6.0,    6.0,

  // //0,  2    0,  5   0,  12  0,  6   0,  13 
  //   10.0,    5.0,    5.0,    5.0,    5.0,
  // //1,  3    1,  5   1,  12  1,  6   1,  13
  //   12.0,    8.0,    8.0,    8.0,    8.0,

  // //4,  6   5,  7   6,  8   7,  9   8,  10
  //   9.0,    9.0,    9.0,    8.0,    8.0,
  // //11, 13  12, 14  13, 15  14, 16  15, 17
  //   9.0,    9.0,    9.0,    8.0,    8.0

  //0,  1   0,  4   0,  11
    4.0,    5.0,    5.0,
  //1,  2   1,  4   1,  11
    // 5.0,    7.0,    7.0,
    5.0,    6.0,    6.0,
  //2,  3   2,  4   2,  11
    6.0,    6.0,    6.0,

  //4,  5   5,  6   6,  7   7,  8   8,  9   9,  10
    3.0,    3.0,    3.0,    3.0,    3.0,    3.0,
  //11, 12  12, 13  13, 14  14, 15  15, 16  16, 17
    3.0,    3.0,    3.0,    3.0,    3.0,    3.0,

  //0,  2    0,  5   0,  12  0,  6   0,  13 
    4.0,     5.0,    5.0,    5.0,    5.0,
  //1,  3    1,  5   1,  12  1,  6   1,  13
    5.0,     5.0,    5.0,    5.0,    5.0,

  //4,  6   5,  7   6,  8   7,  9   8,  10
    3.0,    3.0,    3.0,    3.0,    3.0,
  //11, 13  12, 14  13, 15  14, 16  15, 17
    3.0,    3.0,    3.0,    3.0,    3.0,

};

// Layer pairs
uint8_t layerPairs[2 * nPairs] = {

    0,  1,  0,  4,  0,  11,  // BPIX1 (3)
    1,  2,  1,  4,  1,  11,  // BPIX2 (6)
    2,  3,  2,  4,  2,  11,  // BPIX3 (9)

    4,  5,  5,  6,  6,  7,  7,  8,  8,  9,  9,  10,  // POS (15)
    11, 12, 12, 13, 13, 14, 14, 15, 15, 16, 16, 17,  // NEG (21)

    0,  2,  0,  5,  0,  12, 0,  6,  0,  13,  // BPIX1 Jump (26)
    1,  3,  1,  5,  1,  12, 1,  6,  1,  13,  // BPIX2 Jump (31)

    4,  6,  5,  7,  6,  8,  7,  9,  8,  10,  // POS Jump (36)
    11, 13, 12, 14, 13, 15, 14, 16, 15, 17,  // NEG Jump (41)

};

  // startingPairs: first three true, rest false
constexpr uint8_t startingPairs_flags[nPairs] = {
    1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};

// layer start indices (nLayers+1 elements)
constexpr uint32_t layerStart[nLayers + 1] = {
    0,
    1,
    2,
    3,  // barrel
    4,
    5,
    6,
    7,
    8,
    9,
    10, // positive endcap
    11,
    12,
    13,
    14,
    15,
    16,
    17,
    nLayers // negative endcap
};

// caThetaCuts and caDCACuts
// constexpr float caDCACuts_vals[nLayers] = {0.2*0.15,  //BPix1
//                                            0.2*0.25, 0.2*0.25, 0.2*0.25, 0.2*0.25, 0.2*0.25, 0.2*0.25, 0.2*0.25, 0.2*0.25, 0.2*0.25,
//                                            0.2*0.25, 0.2*0.25, 0.2*0.25, 0.2*0.25, 0.2*0.25, 0.2*0.25, 0.2*0.25, 0.2*0.25};
// constexpr float caThetaCuts_vals[nLayers] = {0.45*0.002, 0.45*0.002, 0.45*0.002, 0.45*0.002,  // BPix
//                                              0.45*0.003, 0.45*0.003, 0.45*0.003, 0.45*0.003, 0.45*0.003, 0.45*0.003, 0.45*0.003, 0.45*0.003,
//                                              0.45*0.003, 0.45*0.003, 0.45*0.003, 0.45*0.003, 0.45*0.003, 0.45*0.003};
// // caThetaCuts and caDCACuts
// constexpr float caDCACuts_vals[nLayers] = {0.15,  //BPix1
//                                            0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25,
//                                            0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25};
// constexpr float caThetaCuts_vals[nLayers] = {0.002, 0.002, 0.002, 0.002,  // BPix
//                                              0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003,
//                                              0.003, 0.003, 0.003, 0.003, 0.003, 0.003};
// caThetaCuts and caDCACuts
constexpr float caDCACuts_vals[nLayers] = {0.16*0.15,  //BPix1
                                           0.16*0.25, 0.16*0.25, 0.16*0.25, 0.18*0.25, 0.18*0.25, 0.18*0.25, 0.18*0.25, 0.18*0.25, 0.18*0.25,
                                           0.18*0.25, 0.18*0.25, 0.18*0.25, 0.18*0.25, 0.18*0.25, 0.18*0.25, 0.18*0.25, 0.18*0.25};
constexpr float caThetaCuts_vals[nLayers] = {0.47*0.002, 0.47*0.002, 0.47*0.002, 0.47*0.002,  // BPix
                                             0.48*0.003, 0.48*0.003, 0.48*0.003, 0.48*0.003, 0.48*0.003, 0.48*0.003, 0.48*0.003, 0.48*0.003,
                                             0.48*0.003, 0.48*0.003, 0.48*0.003, 0.48*0.003, 0.48*0.003, 0.48*0.003};

int writeModules(std::string dataDir) {
    PixelCPEFast dummyPixelCPEFast(dataDir + "/cpefast.bin", dataDir);
    return 1;
}

int write(std::string dataDir) {
    std::ifstream ifs(dataDir + "/AdHocCMSPhase1Modules.bin", std::ios::binary);
    if (!ifs) {
        std::cerr << "Error: cannot open AdHocCMSPhase1Modules.bin for reading.\n";
        return 1;
    }

    std::vector<CAModule> modules(nModules);
    ifs.read(reinterpret_cast<char*>(modules.data()), modules.size() * sizeof(CAModule));
    if (!ifs) {
        std::cerr << "Error: failed to read all modules from file.\n";
        return 1;
    }
    ifs.close();

    std::vector<CALayer> layers(nLayers + 1);
    for (uint16_t i = 0; i < nLayers + 1; i++) {
        layers[i].layerStarts = static_cast<int32_t>(layerStart[i]);
        if (i < nLayers) {
            layers[i].caThetaCut = caThetaCuts_vals[i];
            layers[i].caDCACut   = caDCACuts_vals[i];
        } else {
            layers[i].caThetaCut = 0.0f;
            layers[i].caDCACut   = 0.0f;
        }
    }

    std::vector<CAPair> pairs(nPairs);
    for (uint16_t i = 0; i < nPairs; i++) {
        pairs[i].innerLayer   = layerPairs[2 * i];
        pairs[i].outerLayer   = layerPairs[2 * i + 1];
        pairs[i].phiCut       = phicuts[i];
        pairs[i].minz         = minz_vals[i];
        pairs[i].maxz         = maxz_vals[i];
        pairs[i].maxr         = maxr_vals[i];
        pairs[i].startingPair = startingPairs_flags[i];
    }

    CAGeometrySoA geo;
    geo.m_nModules = nModules;
    geo.m_nLayers  = nLayers;  // detector layers
    geo.m_nPairs   = nPairs;
    geo.m_modules  = modules.data();
    geo.m_layers   = layers.data();
    geo.m_pairs    = pairs.data();

    std::ofstream ofs(dataDir + "/AdHocCMSPhase1Geometry.bin", std::ios::binary);
    ofs.write(reinterpret_cast<const char*>(&geo.m_nModules), sizeof(geo.m_nModules));
    ofs.write(reinterpret_cast<const char*>(&geo.m_nLayers),  sizeof(geo.m_nLayers));
    ofs.write(reinterpret_cast<const char*>(&geo.m_nPairs),   sizeof(geo.m_nPairs));

    ofs.write(reinterpret_cast<const char*>(modules.data()), modules.size() * sizeof(CAModule));
    ofs.write(reinterpret_cast<const char*>(layers.data()),  layers.size()  * sizeof(CALayer));
    ofs.write(reinterpret_cast<const char*>(pairs.data()),   pairs.size()   * sizeof(CAPair));

    ofs.close();
    std::cout << "AdHocCMSPhase1Geometry.bin written.\n";

    return 0;
}

int verify(std::string dataDir) {
    std::ifstream ifs(dataDir + "/AdHocCMSPhase1Geometry.bin", std::ios::binary);

    if (!ifs) {
        std::cerr << "Error opening AdHocCMSPhase1Geometry.bin\n";
        return 1;
    }

    uint16_t nMods, nLays, nP;
    ifs.read(reinterpret_cast<char*>(&nMods), sizeof(nMods));
    ifs.read(reinterpret_cast<char*>(&nLays), sizeof(nLays));
    ifs.read(reinterpret_cast<char*>(&nP), sizeof(nP));

    assert(nMods == nModules);
    assert(nLays == nLayers);
    assert(nP == nPairs);

    std::vector<CAModule> modules(nModules);
    std::vector<CALayer>  layers(nLayers + 1);
    std::vector<CAPair>   pairs(nPairs);

    ifs.read(reinterpret_cast<char*>(modules.data()), modules.size() * sizeof(CAModule));
    ifs.read(reinterpret_cast<char*>(layers.data()),  layers.size()  * sizeof(CALayer));
    ifs.read(reinterpret_cast<char*>(pairs.data()),   pairs.size()   * sizeof(CAPair));
    ifs.close();

    for (uint16_t i = 0; i < nLayers; ++i) {
        assert(layers[i].layerStarts == static_cast<int32_t>(layerStart[i]));
        assert(layers[i].caThetaCut == caThetaCuts_vals[i]);
        assert(layers[i].caDCACut   == caDCACuts_vals[i]);
    }

    for (uint16_t i = 0; i < nPairs; ++i) {
        assert(pairs[i].innerLayer   == layerPairs[2*i]);
        assert(pairs[i].outerLayer   == layerPairs[2*i+1]);
        assert(pairs[i].phiCut       == phicuts[i]);
        assert(pairs[i].minz         == minz_vals[i]);
        assert(pairs[i].maxz         == maxz_vals[i]);
        assert(pairs[i].maxr         == maxr_vals[i]);
        assert(pairs[i].startingPair == startingPairs_flags[i]);
    }

    std::cout << "All values verified successfully!\n";

    return 0;
}

int main(int argc, char* argv[]) {
    fs::path exePath = fs::absolute(argv[0]);
    fs::path exeDir = exePath.parent_path();
    fs::path dataDir = exeDir / "../../../../data/";
    dataDir = fs::canonical(dataDir);
    writeModules(dataDir.string());
    write(dataDir.string());
    verify(dataDir.string());
    return 0;
}
