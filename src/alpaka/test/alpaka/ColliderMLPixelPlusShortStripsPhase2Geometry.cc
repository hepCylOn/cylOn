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

      std::ofstream out(outputPath + "/ColliderMLPixelPlusShortStripsPhase2Modules.bin", std::ios::binary);
      if (!out) throw std::runtime_error("Cannot open file for writing");

      std::cout << "=== Writing " <<  ndetParams << " modules" << std::endl;
      for(auto const &v: m_detParamsGPU)
      {
        auto f = v.frame;
        out.write(reinterpret_cast<const char*>(&f), sizeof(f));
      }

      out.close();
      std::cout << "ColliderMLPixelPlusShortStripsPhase2Modules.bin written.\n";

    }

    ~PixelCPEFast() = default;

      private:
    // allocate it with posix malloc to be compatible with cpu wf
    std::vector<DetParams> m_detParamsGPU;
    CommonParams m_commonParamsGPU;


};

// --- Constants ---
constexpr uint16_t nModules = 34;
constexpr uint16_t nLayers  = 34;
constexpr uint16_t nPairs   = 5 + 14 + 22 + 3 + 4 + 26 + 5 + 14 + 18; // pixel barrel only + pixel barrel-endcap + pixel endcap only + pixel-shortStrips barrel only + pixel-shortStrips endcap-barrel + pixel-shortStrips endcap only + shortStrips barrel only + shortStrips barrel-endcap + shortStrips endcap-only

// Layer pairs
uint8_t layerPairs[2 * nPairs] = {

    0, 1, 1, 2, 2, 3, 0, 2, 1, 3,  // pixel barrel only (5)

    0, 4, 0, 5, 0, 6, 1, 4, 1, 5, 1, 6, 2, 4,        // pixel barrel-endcap pos (12)
    0, 11, 0, 12, 0, 13, 1, 11, 1, 12, 1, 13, 2, 11, // pixel barrel-endcap neg (19)

    4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 4, 6, 5, 7, 6, 8, 7, 9, 8, 10,                     // pixel endcap only pos (30)
    11, 12, 12, 13, 13, 14, 14, 15, 15, 16, 16, 17, 11, 13, 12, 14, 13, 15, 14, 16, 15, 17, // pixel endcap only neg (41)

    3, 18, 2, 18, 3, 19, // pixel-shortStrips barrel only (44)

    4, 18, 5, 18,   // pixel-shortStrips endcap-barrel pos (46)
    11, 18, 12, 18, // pixel-shortStrips endcap-barrel neg (48)

    4, 22, 5, 22, 6, 22, 5, 23, 6, 23, 7, 22, 7, 23, 8, 23, 8, 24, 9, 24, 9, 25, 10, 25, 10, 26,            // pixel-shortStrips endcap only pos (61)
    11, 28, 12, 28, 13, 28, 12, 29, 13, 29, 14, 28, 14, 29, 15, 29, 15, 30, 16, 30, 16, 31, 17, 31, 17, 32, // pixel-shortStrips endcap only neg (74)

    18, 19, 19, 20, 20, 21, 18, 20, 19, 21, // shortStrips barrel only (79)

    18, 22, 18, 23, 18, 24, 19, 22, 19, 23, 19, 24, 20, 22, // shortStrips barrel-endcap pos (86)
    18, 28, 18, 29, 18, 30, 19, 28, 19, 29, 19, 30, 20, 28, // shortStrips barrel-endcap neg (93)

    22, 23, 23, 24, 24, 25, 25, 26, 26, 27, 22, 24, 23, 25, 24, 26, 25, 27,  // shortStrips endcap only pos (102)
    28, 29, 29, 30, 30, 31, 31, 32, 32, 33, 28, 30, 29, 31, 30, 32, 31, 33,  // shortStrips endcap only neg (111)

};

// phi cut constants
constexpr int16_t phi0p05 = 522;
constexpr int16_t phi0p06 = 626;
constexpr int16_t phi0p07 = 730;

// Phi cuts for layer pairs
constexpr int16_t phicuts[nPairs]{

    phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, // pixel barrel only (5)

    phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, // pixel barrel-endcap pos (12)
    phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, // pixel barrel-endcap neg (19)

    phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, // pixel endcap only pos (30)
    phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, // pixel endcap only neg (41)

    phi0p05, phi0p05, phi0p05, // pixel-shortStrips barrel only (44)

    phi0p05, phi0p05, // pixel-shortStrips endcap-barrel pos (46)
    phi0p05, phi0p05, // pixel-shortStrips endcap-barrel neg (48)

    phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, // pixel-shortStrips endcap only pos (61)
    phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, // pixel-shortStrips endcap only neg (74)

    phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, // shortStrips barrel only (79)

    phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, // shortStrips barrel-endcap pos (86)
    phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, // shortStrips barrel-endcap neg (93)

    phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, // shortStrips endcap only pos (102)
    phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, // shortStrips endcap only neg (111)

};

constexpr int16_t pH = 0.0;

// Min z for layer pairs
constexpr float minz_vals[nPairs] = {

    pH, pH, pH, pH, pH, // pixel barrel only (5)

    pH, pH, pH, pH, pH, pH, pH, // pixel barrel-endcap pos (12)
    pH, pH, pH, pH, pH, pH, pH, // pixel barrel-endcap neg (19)

    pH, pH, pH, pH, pH, pH, pH, pH, pH, pH, pH, // pixel endcap only pos (30)
    pH, pH, pH, pH, pH, pH, pH, pH, pH, pH, pH, // pixel endcap only neg (41)

    pH, pH, pH, // pixel-shortStrips barrel only (44)

    pH, pH, // pixel-shortStrips endcap-barrel pos (46)
    pH, pH, // pixel-shortStrips endcap-barrel neg (48)

    pH, pH, pH, pH, pH, pH, pH, pH, pH, pH, pH, pH, pH, // pixel-shortStrips endcap only pos (61)
    pH, pH, pH, pH, pH, pH, pH, pH, pH, pH, pH, pH, pH, // pixel-shortStrips endcap only neg (74)

    pH, pH, pH, pH, pH, // shortStrips barrel only (79)

    pH, pH, pH, pH, pH, pH, pH, // shortStrips barrel-endcap pos (86)
    pH, pH, pH, pH, pH, pH, pH, // shortStrips barrel-endcap neg (93)

    pH, pH, pH, pH, pH, pH, pH, pH, pH, // shortStrips endcap only pos (102)
    pH, pH, pH, pH, pH, pH, pH, pH, pH, // shortStrips endcap only neg (111)

};

// Max z for layer pairs
constexpr float maxz_vals[nPairs] = {

    pH, pH, pH, pH, pH, // pixel barrel only (5)

    pH, pH, pH, pH, pH, pH, pH, // pixel barrel-endcap pos (12)
    pH, pH, pH, pH, pH, pH, pH, // pixel barrel-endcap neg (19)

    pH, pH, pH, pH, pH, pH, pH, pH, pH, pH, pH, // pixel endcap only pos (30)
    pH, pH, pH, pH, pH, pH, pH, pH, pH, pH, pH, // pixel endcap only neg (41)

    pH, pH, pH, // pixel-shortStrips barrel only (44)

    pH, pH, // pixel-shortStrips endcap-barrel pos (46)
    pH, pH, // pixel-shortStrips endcap-barrel neg (48)

    pH, pH, pH, pH, pH, pH, pH, pH, pH, pH, pH, pH, pH, // pixel-shortStrips endcap only pos (61)
    pH, pH, pH, pH, pH, pH, pH, pH, pH, pH, pH, pH, pH, // pixel-shortStrips endcap only neg (74)

    pH, pH, pH, pH, pH, // shortStrips barrel only (79)

    pH, pH, pH, pH, pH, pH, pH, // shortStrips barrel-endcap pos (86)
    pH, pH, pH, pH, pH, pH, pH, // shortStrips barrel-endcap neg (93)

    pH, pH, pH, pH, pH, pH, pH, pH, pH, // shortStrips endcap only pos (102)
    pH, pH, pH, pH, pH, pH, pH, pH, pH, // shortStrips endcap only neg (111)

};

// Max r for layer pairs
constexpr float maxr_vals[nPairs] = {
    
    pH, pH, pH, pH, pH, // pixel barrel only (5)

    pH, pH, pH, pH, pH, pH, pH, // pixel barrel-endcap pos (12)
    pH, pH, pH, pH, pH, pH, pH, // pixel barrel-endcap neg (19)

    pH, pH, pH, pH, pH, pH, pH, pH, pH, pH, pH, // pixel endcap only pos (30)
    pH, pH, pH, pH, pH, pH, pH, pH, pH, pH, pH, // pixel endcap only neg (41)

    pH, pH, pH, // pixel-shortStrips barrel only (44)

    pH, pH, // pixel-shortStrips endcap-barrel pos (46)
    pH, pH, // pixel-shortStrips endcap-barrel neg (48)

    pH, pH, pH, pH, pH, pH, pH, pH, pH, pH, pH, pH, pH, // pixel-shortStrips endcap only pos (61)
    pH, pH, pH, pH, pH, pH, pH, pH, pH, pH, pH, pH, pH, // pixel-shortStrips endcap only neg (74)

    pH, pH, pH, pH, pH, // shortStrips barrel only (79)

    pH, pH, pH, pH, pH, pH, pH, // shortStrips barrel-endcap pos (86)
    pH, pH, pH, pH, pH, pH, pH, // shortStrips barrel-endcap neg (93)

    pH, pH, pH, pH, pH, pH, pH, pH, pH, // shortStrips endcap only pos (102)
    pH, pH, pH, pH, pH, pH, pH, pH, pH, // shortStrips endcap only neg (111)

};

  // startingPairs: first three true, rest false
constexpr uint8_t startingPairs_flags[nPairs] = {
    1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};

// layer start indices (nLayers+1 elements)
constexpr uint32_t layerStart[nLayers + 1] = {
    0,
    1,
    2,
    3,  // pixel barrel
    4,
    5,
    6,
    7,
    8,
    9,
    10, // pixel positive endcap
    11,
    12,
    13,
    14,
    15,
    16,
    17, // pixel negative endcap
    18,
    19,
    20,
    21, // shortStrips barrel
    22,
    23,
    24,
    25,
    26,
    27, // shortStrips positive endcap
    28,
    29,
    30,
    31,
    32,
    33, // shortStrips negative endcap
    nLayers 
};

// caThetaCuts and caDCACuts
constexpr float caDCACuts_vals[nLayers] = {pH, pH, pH, pH, pH, pH, pH, pH, pH,
                                           pH, pH, pH, pH, pH, pH, pH, pH, pH,
                                           pH, pH, pH, pH, pH, pH, pH, pH, pH,
                                           pH, pH, pH, pH, pH, pH, pH};
constexpr float caThetaCuts_vals[nLayers] = {pH, pH, pH, pH, pH, pH, pH, pH, pH,
                                             pH, pH, pH, pH, pH, pH, pH, pH, pH,
                                             pH, pH, pH, pH, pH, pH, pH, pH, pH,
                                             pH, pH, pH, pH, pH, pH, pH};

int writeModules(std::string dataDir) {
    PixelCPEFast dummyPixelCPEFast(dataDir + "/cpefast.bin", dataDir);
    return 1;
}

int write(std::string dataDir) {
    std::ifstream ifs(dataDir + "/ColliderMLPixelPlusShortStripsPhase2Modules.bin", std::ios::binary);
    if (!ifs) {
        std::cerr << "Error: cannot open ColliderMLPixelPlusShortStripsPhase2Modules.bin for reading.\n";
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

    std::ofstream ofs(dataDir + "/ColliderMLPixelPlusShortStripsPhase2Geometry.bin", std::ios::binary);
    ofs.write(reinterpret_cast<const char*>(&geo.m_nModules), sizeof(geo.m_nModules));
    ofs.write(reinterpret_cast<const char*>(&geo.m_nLayers),  sizeof(geo.m_nLayers));
    ofs.write(reinterpret_cast<const char*>(&geo.m_nPairs),   sizeof(geo.m_nPairs));

    ofs.write(reinterpret_cast<const char*>(modules.data()), modules.size() * sizeof(CAModule));
    ofs.write(reinterpret_cast<const char*>(layers.data()),  layers.size()  * sizeof(CALayer));
    ofs.write(reinterpret_cast<const char*>(pairs.data()),   pairs.size()   * sizeof(CAPair));

    ofs.close();
    std::cout << "ColliderMLPixelPlusShortStripsPhase2Geometry.bin written.\n";

    return 0;
}

int verify(std::string dataDir) {
    std::ifstream ifs(dataDir + "/ColliderMLPixelPlusShortStripsPhase2Geometry.bin", std::ios::binary);

    if (!ifs) {
        std::cerr << "Error opening ColliderMLPixelPlusShortStripsPhase2Geometry.bin\n";
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