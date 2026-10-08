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
      m_detParamsGPU.resize(4000);

      std::cout << "Number of elements in m_detParamsGPU: " << m_detParamsGPU.size() << std::endl;
      std::cout << "Size of ndetParams: " << ndetParams << std::endl;

      std::ofstream out(outputPath + "/CMSPhase2Modules.bin", std::ios::binary);
      if (!out) throw std::runtime_error("Cannot open file for writing");

      std::cout << "=== Writing " <<  ndetParams << " modules" << std::endl;
      for(auto const &v: m_detParamsGPU)
      {
        auto f = v.frame;
        out.write(reinterpret_cast<const char*>(&f), sizeof(f));
      }

      out.close();
      std::cout << "CMSPhase2Modules.bin written.\n";

    }

    ~PixelCPEFast() = default;

      private:
    // allocate it with posix malloc to be compatible with cpu wf
    std::vector<DetParams> m_detParamsGPU;
    CommonParams m_commonParamsGPU;


};

// --- Constants ---
constexpr uint16_t nModules = 4000;
constexpr uint16_t nLayers  = 28;
constexpr uint16_t nPairs   = 23 + 6 + 14 + 8 + 4;

// phi cut constants
constexpr int16_t phi0p05 = 522;
constexpr int16_t phi0p06 = 626;
constexpr int16_t phi0p07 = 730;
// constexpr int16_t phi0p09 = 900;

constexpr int16_t phicuts[nPairs]{
    phi0p05, phi0p05, phi0p05, phi0p06, phi0p07, phi0p07, phi0p06, phi0p07, phi0p07, phi0p05, phi0p05,
      phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05,
      phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p07, phi0p07, phi0p07, phi0p07,
      phi0p07, phi0p07, phi0p07, phi0p07, phi0p07, phi0p07, phi0p07, phi0p07, phi0p07, phi0p07, phi0p07,
      phi0p07, phi0p07, phi0p07, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05
};

constexpr float minz_vals[nPairs]{
    -16.0, 4.0,   -22.0, -17.0, 6.0,   -22.0, -18.0, 11.0,  -22.0,  23.0,   30.0,   39.0,   50.0,   65.0,
      82.0,  109.0, -28.0, -35.0, -44.0, -55.0, -70.0, -87.0, -113.0, -16.,   7.0,    -22.0,  11.0,   -22.0,
      -17.0, 9.0,   -22.0, 13.0,  -22.0, 137.0, 173.0, 199.0, 229.0,  -142.0, -177.0, -203.0, -233.0, 23.0,
      30.0,  39.0,  50.0,  65.0,  82.0,  109.0, -28.0, -35.0, -44.0,  -55.0,  -70.0,  -87.0,  -113.0
};
constexpr float maxz_vals[nPairs]{
    17.0, 22.0,  -4.0,  17.0,  22.0,  -6.0,  18.0,  22.0,  -11.0,  28.0,   35.0,   44.0,   55.0,   70.0,
      87.0, 113.0, -23.0, -30.0, -39.0, -50.0, -65.0, -82.0, -109.0, 17.0,   22.0,   -7.0,   22.0,   -10.0,
      17.0, 22.0,  -9.0,  22.0,  -13.0, 142.0, 177.0, 203.0, 233.0,  -137.0, -173.0, -199.0, -229.0, 28.0,
      35.0, 44.0,  55.0,  70.0,  87.0,  113.0, -23.0, -30.0, -39.0,  -50.0,  -65.0,  -82.0,  -109.0
};
constexpr float maxr_vals[nPairs]{
    5.0, 5.0, 5.0, 7.0, 8.0, 8.0,  7.0, 7.0, 7.0, 6.0, 6.0, 6.0, 6.0, 5.0,
                                             6.0, 5.0, 6.0, 6.0, 6.0, 6.0,  5.0, 6.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0,
                                             5.0, 8.0, 8.0, 8.0, 8.0, 6.0,  5.0, 5.0, 5.0, 6.0, 5.0, 5.0, 5.0, 9.0,
                                             9.0, 9.0, 8.0, 8.0, 8.0, 11.0, 9.0, 9.0, 9.0, 8.0, 8.0, 8.0, 11.0
};

// Layer pairs
constexpr uint8_t layerPairs[2 * nPairs] = {
    0,  1,  0,  4,  0,  16,  // BPIX1 (3)
      1,  2,  1,  4,  1,  16,  // BPIX2 (6)
      2,  3,  2,  4,  2,  16,  // BPIX3 & Forward (9)

      4,  5,  5,  6,  6,  7,  7,  8,  8,  9,  9,  10, 10, 11,  // POS (16)
      16, 17, 17, 18, 18, 19, 19, 20, 20, 21, 21, 22, 22, 23,  // NEG (23)

      0,  2,  0,  5,  0,  17, 0,  6,  0,  18,  // BPIX1 Jump (28)
      1,  3,  1,  5,  1,  17, 1,  6,  1,  18,  // BPIX2 Jump (33)

      11, 12, 12, 13, 13, 14, 14, 15,  // Late POS (37)
      23, 24, 24, 25, 25, 26, 26, 27,  // Late NEG (41)

      4,  6,  5,  7,  6,  8,  7,  9,  8,  10, 9,  11, 10, 12,  // POS Jump (48)
      16, 18, 17, 19, 18, 20, 19, 21, 20, 22, 21, 23, 22, 24,  // NEG Jump (55)
};

// startingPairs: first three true, rest false
constexpr uint8_t startingPairs_flags[nPairs] = {
    1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};

// layer start indices (nLayers+1 elements)
constexpr uint32_t layerStart[nLayers + 1] = {
    0,
    216,
    432,
    612,  // Barrel
    864,
    972,
    1080,
    1188,
    1296,
    1404,
    1512,
    1620,
    1728,
    1904,
    2080,
    2256,  // Fp
    2432,
    2540,
    2648,
    2756,
    2864,
    2972,
    3080,
    3188,
    3296,
    3472,
    3648,
    3824,  // Np
    nModules
};

// caThetaCuts and caDCACuts
constexpr float caDCACuts_vals[nLayers] = {
    0.15,  //BPix1
                                                        0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25,
                                                        0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25,
                                                        0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25
};
constexpr float caThetaCuts_vals[nLayers] = {
    0.002, 0.002, 0.002, 0.002,  // BPix
                                                          0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003,
                                                          0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003,
                                                          0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003
};

int writeModules(std::string dataDir) {
    PixelCPEFast dummyPixelCPEFast(dataDir + "/cpefast.bin", dataDir);
    return 1;
}

int write(std::string dataDir) {
    std::ifstream ifs(dataDir + "/CMSPhase2Modules.bin", std::ios::binary);
    if (!ifs) {
        std::cerr << "Error: cannot open CMSPhase2Modules.bin for reading.\n";
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

    std::ofstream ofs(dataDir + "/CMSPhase2Geometry.bin", std::ios::binary);
    ofs.write(reinterpret_cast<const char*>(&geo.m_nModules), sizeof(geo.m_nModules));
    ofs.write(reinterpret_cast<const char*>(&geo.m_nLayers),  sizeof(geo.m_nLayers));
    ofs.write(reinterpret_cast<const char*>(&geo.m_nPairs),   sizeof(geo.m_nPairs));

    ofs.write(reinterpret_cast<const char*>(modules.data()), modules.size() * sizeof(CAModule));
    ofs.write(reinterpret_cast<const char*>(layers.data()),  layers.size()  * sizeof(CALayer));
    ofs.write(reinterpret_cast<const char*>(pairs.data()),   pairs.size()   * sizeof(CAPair));

    ofs.close();
    std::cout << "CMSPhase2Geometry.bin written.\n";

    return 0;
}

int verify(std::string dataDir) {
    std::ifstream ifs(dataDir + "/CMSPhase2Geometry.bin", std::ios::binary);

    if (!ifs) {
        std::cerr << "Error opening CMSPhase2Geometry.bin\n";
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
