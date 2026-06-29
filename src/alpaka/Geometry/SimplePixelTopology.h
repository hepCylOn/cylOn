#ifndef Geometry_SimplePixelTopology_h
#define Geometry_SimplePixelTopology_h

#include <array>
#include <cstdint>
#include <type_traits>
#include "Framework/HostDeviceConstant.h"

namespace pixelTopology {

  constexpr auto maxNumberOfLadders = 160;
  constexpr uint8_t maxLayers = 60;
  constexpr uint8_t maxPairs = 150;

  // TODO
  // Once CUDA is dropped this could be wrapped in #ifdef CA_TRIPLETS_HOLE
  // see DataFormats/TrackingRecHitSoa/interface/TrackingRecHitSoA.h

  template <typename TrackerTraits>
  struct AverageGeometryT {
    //
    float ladderZ[TrackerTraits::numberOfLaddersInBarrel];
    float ladderX[TrackerTraits::numberOfLaddersInBarrel];
    float ladderY[TrackerTraits::numberOfLaddersInBarrel];
    float ladderR[TrackerTraits::numberOfLaddersInBarrel];
    float ladderMinZ[TrackerTraits::numberOfLaddersInBarrel];
    float ladderMaxZ[TrackerTraits::numberOfLaddersInBarrel];
    float endCapZ[2];  // just for pos and neg Layer1
  };

  constexpr int16_t phi0p05 = 522;  // round(521.52189...) = phi2short(0.05);
  constexpr int16_t phi0p06 = 626;  // round(625.82270...) = phi2short(0.06);
  constexpr int16_t phi0p07 = 730;  // round(730.12648...) = phi2short(0.07);
  constexpr int16_t phi0p09 = 900;

  constexpr float pH = 0.0;

  constexpr uint16_t last_barrel_layer = 3;  // this is common between all the topologies

  template <class Function, std::size_t... Indices>
  constexpr auto map_to_array_helper(Function f, std::index_sequence<Indices...>)
      -> std::array<std::invoke_result_t<Function, std::size_t>, sizeof...(Indices)> {
    return {{f(Indices)...}};
  }

  template <int N, class Function>
  constexpr auto map_to_array(Function f) -> std::array<std::invoke_result_t<Function, std::size_t>, N> {
    return map_to_array_helper(f, std::make_index_sequence<N>{});
  }

  template <typename TrackerTraits>
  constexpr uint16_t findMaxModuleStride() {
    bool go = true;
    int n = 2;
    while (go) {
      for (uint8_t i = 1; i < TrackerTraits::numberOfLayers + 1; ++i) {
        if (TrackerTraits::layerStart[i] % n != 0) {
          go = false;
          break;
        }
      }
      if (!go)
        break;
      n *= 2;
    }
    return n / 2;
  }

  template <typename TrackerTraits>
  constexpr uint16_t maxModuleStride = findMaxModuleStride<TrackerTraits>();

  template <typename TrackerTraits>
  constexpr uint8_t findLayer(uint32_t detId, uint8_t sl = 0) {
    for (uint8_t i = sl; i < TrackerTraits::numberOfLayers + 1; ++i)
      if (detId < TrackerTraits::layerStart[i + 1])
        return i;
    return TrackerTraits::numberOfLayers + 1;
  }

  template <typename TrackerTraits>
  constexpr uint8_t findLayerFromCompact(uint32_t detId) {
    detId *= maxModuleStride<TrackerTraits>;
    for (uint8_t i = 0; i < TrackerTraits::numberOfLayers + 1; ++i)
      if (detId < TrackerTraits::layerStart[i + 1])
        return i;
    return TrackerTraits::numberOfLayers + 1;
  }

  template <typename TrackerTraits>
  constexpr uint32_t layerIndexSize = TrackerTraits::numberOfModules / maxModuleStride<TrackerTraits>;

  template <typename TrackerTraits>
#ifdef __CUDA_ARCH__
  __device__
#endif
      constexpr std::array<uint8_t, layerIndexSize<TrackerTraits>>
          layer = map_to_array<layerIndexSize<TrackerTraits>>(findLayerFromCompact<TrackerTraits>);

  template <typename TrackerTraits>
  constexpr uint8_t getLayer(uint32_t detId) {
    return layer<TrackerTraits>[detId / maxModuleStride<TrackerTraits>];
  }

  template <typename TrackerTraits>
  constexpr bool validateLayerIndex() {
    bool res = true;
    for (auto i = 0U; i < TrackerTraits::numberOfModules; ++i) {
      auto j = i / maxModuleStride<TrackerTraits>;
      res &= (layer<TrackerTraits>[j] < TrackerTraits::numberOfLayers);
      res &= (i >= TrackerTraits::layerStart[layer<TrackerTraits>[j]]);
      res &= (i < TrackerTraits::layerStart[layer<TrackerTraits>[j] + 1]);
    }
    return res;
  }

  template <typename TrackerTraits>
#ifdef __CUDA_ARCH__
  __device__
#endif
      constexpr inline uint32_t
      layerStart(uint32_t i) {
    return TrackerTraits::layerStart[i];
  }

  constexpr inline uint16_t divu52(uint16_t n) {
    n = n >> 2;
    uint16_t q = (n >> 1) + (n >> 4);
    q = q + (q >> 4) + (q >> 5);
    q = q >> 3;
    uint16_t r = n - q * 13;
    return q + ((r + 3) >> 4);
  }
}  // namespace pixelTopology

namespace phase1PixelTopology {

  using pixelTopology::phi0p05;
  using pixelTopology::phi0p06;
  using pixelTopology::phi0p07;

  constexpr uint32_t numberOfLayers = 10;
  constexpr int nPairs = 13 + 2 + 4;
  constexpr uint16_t numberOfModules = 1856;

  constexpr uint32_t maxNumClustersPerModules = 1024;

  constexpr uint32_t max_ladder_bpx0 = 12;
  constexpr uint32_t first_ladder_bpx0 = 0;
  constexpr float module_length_bpx0 = 6.7f;
  constexpr float module_tolerance_bpx0 = 0.4f;  // projection to cylinder is inaccurate on BPIX1
  constexpr uint32_t max_ladder_bpx4 = 64;
  constexpr uint32_t first_ladder_bpx4 = 84;
  constexpr float radius_even_ladder = 15.815f;
  constexpr float radius_odd_ladder = 16.146f;
  constexpr float module_length_bpx4 = 6.7f;
  constexpr float module_tolerance_bpx4 = 0.2f;
  constexpr float barrel_z_length = 26.f;
  constexpr float forward_z_begin = 32.f;

  HOST_DEVICE_CONSTANT uint8_t layerPairs[2 * nPairs] = {
      0, 1, 0, 4, 0, 7,              // BPIX1 (3)
      1, 2, 1, 4, 1, 7,              // BPIX2 (6)
      4, 5, 7, 8,                    // FPIX1 (8)
      2, 3, 2, 4, 2, 7, 5, 6, 8, 9,  // BPIX3 & FPIX2 (13)
      0, 2, 1, 3,                    // Jumping Barrel (15)
      0, 5, 0, 8,                    // Jumping Forward (BPIX1,FPIX2)
      4, 6, 7, 9                     // Jumping Forward (19)
  };

  HOST_DEVICE_CONSTANT int16_t phicuts[nPairs]{phi0p05,
                                               phi0p07,
                                               phi0p07,
                                               phi0p05,
                                               phi0p06,
                                               phi0p06,
                                               phi0p05,
                                               phi0p05,
                                               phi0p06,
                                               phi0p06,
                                               phi0p06,
                                               phi0p05,
                                               phi0p05,
                                               phi0p05,
                                               phi0p05,
                                               phi0p05,
                                               phi0p05,
                                               phi0p05,
                                               phi0p05};
  HOST_DEVICE_CONSTANT float minz[nPairs] = {
      -20., 0., -30., -22., 10., -30., -70., -70., -22., 15., -30, -70., -70., -20., -22., 0, -30., -70., -70.};
  HOST_DEVICE_CONSTANT float maxz[nPairs] = {
      20., 30., 0., 22., 30., -10., 70., 70., 22., 30., -15., 70., 70., 20., 22., 30., 0., 70., 70.};
  HOST_DEVICE_CONSTANT float maxr[nPairs] = {
      20., 9., 9., 20., 7., 7., 5., 5., 20., 6., 6., 5., 5., 20., 20., 9., 9., 9., 9.};

  HOST_DEVICE_CONSTANT float dcaCuts[numberOfLayers] = {0.15, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25};

  HOST_DEVICE_CONSTANT float thetaCuts[numberOfLayers] = {
      0.002, 0.002, 0.002, 0.002, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003};

  static constexpr uint32_t layerStart[numberOfLayers + 1] = {0,
                                                              96,
                                                              320,
                                                              672,  // barrel
                                                              1184,
                                                              1296,
                                                              1408,  // positive endcap
                                                              1520,
                                                              1632,
                                                              1744,  // negative endcap
                                                              numberOfModules};
}  // namespace phase1PixelTopology

namespace phase2PixelTopology {

  using pixelTopology::phi0p05;
  using pixelTopology::phi0p06;
  using pixelTopology::phi0p07;
  using pixelTopology::phi0p09;

  constexpr uint32_t numberOfLayers = 28;
  constexpr int nPairs = 23 + 6 + 14 + 8 + 4;  // include far forward layer pairs
  constexpr uint16_t numberOfModules = 4000;

  constexpr uint32_t maxNumClustersPerModules = 1024;

  HOST_DEVICE_CONSTANT uint8_t layerPairs[2 * nPairs] = {

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

  HOST_DEVICE_CONSTANT uint32_t layerStart[numberOfLayers + 1] = {0,
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
                                                                  numberOfModules};

  HOST_DEVICE_CONSTANT int16_t phicuts[nPairs]{
      phi0p05, phi0p05, phi0p05, phi0p06, phi0p07, phi0p07, phi0p06, phi0p07, phi0p07, phi0p05, phi0p05,
      phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05,
      phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p07, phi0p07, phi0p07, phi0p07,
      phi0p07, phi0p07, phi0p07, phi0p07, phi0p07, phi0p07, phi0p07, phi0p07, phi0p07, phi0p07, phi0p07,
      phi0p07, phi0p07, phi0p07, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05};

  HOST_DEVICE_CONSTANT float minz[nPairs] = {
      -16.0, 4.0,   -22.0, -17.0, 6.0,   -22.0, -18.0, 11.0,  -22.0,  23.0,   30.0,   39.0,   50.0,   65.0,
      82.0,  109.0, -28.0, -35.0, -44.0, -55.0, -70.0, -87.0, -113.0, -16.,   7.0,    -22.0,  11.0,   -22.0,
      -17.0, 9.0,   -22.0, 13.0,  -22.0, 137.0, 173.0, 199.0, 229.0,  -142.0, -177.0, -203.0, -233.0, 23.0,
      30.0,  39.0,  50.0,  65.0,  82.0,  109.0, -28.0, -35.0, -44.0,  -55.0,  -70.0,  -87.0,  -113.0};

  HOST_DEVICE_CONSTANT float maxz[nPairs] = {

      17.0, 22.0,  -4.0,  17.0,  22.0,  -6.0,  18.0,  22.0,  -11.0,  28.0,   35.0,   44.0,   55.0,   70.0,
      87.0, 113.0, -23.0, -30.0, -39.0, -50.0, -65.0, -82.0, -109.0, 17.0,   22.0,   -7.0,   22.0,   -10.0,
      17.0, 22.0,  -9.0,  22.0,  -13.0, 142.0, 177.0, 203.0, 233.0,  -137.0, -173.0, -199.0, -229.0, 28.0,
      35.0, 44.0,  55.0,  70.0,  87.0,  113.0, -23.0, -30.0, -39.0,  -50.0,  -65.0,  -82.0,  -109.0};

  HOST_DEVICE_CONSTANT float maxr[nPairs] = {5.0, 5.0, 5.0, 7.0, 8.0, 8.0,  7.0, 7.0, 7.0, 6.0, 6.0, 6.0, 6.0, 5.0,
                                             6.0, 5.0, 6.0, 6.0, 6.0, 6.0,  5.0, 6.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0,
                                             5.0, 8.0, 8.0, 8.0, 8.0, 6.0,  5.0, 5.0, 5.0, 6.0, 5.0, 5.0, 5.0, 9.0,
                                             9.0, 9.0, 8.0, 8.0, 8.0, 11.0, 9.0, 9.0, 9.0, 8.0, 8.0, 8.0, 11.0};

  HOST_DEVICE_CONSTANT float dcaCuts[numberOfLayers] = {0.15,  //BPix1
                                                        0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25,
                                                        0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25,
                                                        0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25};

  HOST_DEVICE_CONSTANT float thetaCuts[numberOfLayers] = {0.002, 0.002, 0.002, 0.002,  // BPix
                                                          0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003,
                                                          0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003,
                                                          0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003};

}  // namespace phase2PixelTopology

namespace phase1HIonPixelTopology {
  // Storing here the needed constants different w.r.t. pp Phase1 topology.
  // All the other defined by inheritance in the HIon topology struct.
  using pixelTopology::phi0p09;

  constexpr uint32_t maxNumClustersPerModules = 2048;

  HOST_DEVICE_CONSTANT int16_t phicuts[phase1PixelTopology::nPairs]{phi0p09,
                                                                    phi0p09,
                                                                    phi0p09,
                                                                    phi0p09,
                                                                    phi0p09,
                                                                    phi0p09,
                                                                    phi0p09,
                                                                    phi0p09,
                                                                    phi0p09,
                                                                    phi0p09,
                                                                    phi0p09,
                                                                    phi0p09,
                                                                    phi0p09,
                                                                    phi0p09,
                                                                    phi0p09,
                                                                    phi0p09,
                                                                    phi0p09,
                                                                    phi0p09,
                                                                    phi0p09};

  HOST_DEVICE_CONSTANT float dcaCuts[phase1PixelTopology::numberOfLayers] = {
      0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};

  HOST_DEVICE_CONSTANT float thetaCuts[phase1PixelTopology::numberOfLayers] = {
      0.001, 0.001, 0.001, 0.001, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002};

}  // namespace phase1HIonPixelTopology

namespace colliderMLPhase1PixelTopology {

  using pixelTopology::phi0p05;
  using pixelTopology::phi0p06;
  using pixelTopology::phi0p07;

  constexpr uint32_t numberOfLayers = 10;
  constexpr int nPairs = 13 + 2 + 4;
  constexpr uint16_t numberOfModules = 10;

  constexpr uint32_t maxNumClustersPerModules = 1024;

  constexpr uint32_t max_ladder_bpx0 = 12;
  constexpr uint32_t first_ladder_bpx0 = 0;
  constexpr float module_length_bpx0 = 6.7f;
  constexpr float module_tolerance_bpx0 = 0.4f;  // projection to cylinder is inaccurate on BPIX1
  constexpr uint32_t max_ladder_bpx4 = 64;
  constexpr uint32_t first_ladder_bpx4 = 84;
  constexpr float radius_even_ladder = 15.815f;
  constexpr float radius_odd_ladder = 16.146f;
  constexpr float module_length_bpx4 = 6.7f;
  constexpr float module_tolerance_bpx4 = 0.2f;
  constexpr float barrel_z_length = 26.f;
  constexpr float forward_z_begin = 32.f;

  HOST_DEVICE_CONSTANT uint8_t layerPairs[2 * nPairs] = {
      0, 1, 0, 4, 0, 7,              // BPIX1 (3)
      1, 2, 1, 4, 1, 7,              // BPIX2 (6)
      4, 5, 7, 8,                    // FPIX1 (8)
      2, 3, 2, 4, 2, 7, 5, 6, 8, 9,  // BPIX3 & FPIX2 (13)
      0, 2, 1, 3,                    // Jumping Barrel (15)
      0, 5, 0, 8,                    // Jumping Forward (BPIX1,FPIX2)
      4, 6, 7, 9                     // Jumping Forward (19)
  };

  // ------ Begin for SimPixelTracks geometry ------
  HOST_DEVICE_CONSTANT uint8_t startingPairs[nPairs] = {
    1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
  };

  HOST_DEVICE_CONSTANT bool isBarrel[numberOfLayers] = {
    1, 1, 1, 1, 0, 0, 0, 0, 0, 0
  };

  HOST_DEVICE_CONSTANT float ptCuts[nPairs] = {
    0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
    0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5
  };
  // ------ End for SimPixelTracks geometry ------

  // // Using values from simPixelTracksAnalyser 99.5%
  // HOST_DEVICE_CONSTANT int16_t phicuts[nPairs]{250,
  //                                              330,
  //                                              330,
  //                                              270,
  //                                              370,
  //                                              370,
  //                                              130,
  //                                              130,
  //                                              310,
  //                                              290,
  //                                              290,
  //                                              130,
  //                                              130,
  //                                              450,
  //                                              550,
  //                                              190,
  //                                              190,
  //                                              110,
  //                                              110};
  // Using values from simPixelTracksAnalyser 90%
  HOST_DEVICE_CONSTANT int16_t phicuts[nPairs]{170,
                                               190,
                                               190,
                                               230,
                                               250,
                                               250,
                                               90,
                                               90,
                                               270,
                                               230,
                                               230,
                                               90,
                                               90,
                                               0,
                                               0,
                                               0,
                                               0,
                                               0,
                                               0};

  // // Using values from simPixelTracksAnalyser 99.5%
  // HOST_DEVICE_CONSTANT float minz[nPairs] = {
  //     -25., 20., -50., -32., 28., -50., 60., -70., -34., 38., -50, 70., -72., -20., -25., 45, -52., -72., -72.};
  // HOST_DEVICE_CONSTANT float maxz[nPairs] = {
  //     25., 50., -20., 32., 50., -28., 70., -60., 34., 50., -38., 72., -70., 20., 25., 52., -45., -70., -70.};
  HOST_DEVICE_CONSTANT float maxr[nPairs] = {
      4., 5., 5., 5., 7., 7., 3., 3., 6., 6., 6., 3., 3., 9., 11., 2., 2., 2., 2.};
  // Using values from simPixelTracksAnalyser 90%
  HOST_DEVICE_CONSTANT float minz[nPairs] = {
      -19.5427, 23.8853, -44.1519, -22.4379, 31.1233, -49.9423, 61.5229, -62.9707, -28.2283, 39.8089, -49.9423, 71.6561, -71.6563, 10.8569, -49.9423, -71.6563, -71.6563, -71.6563, -71.6563};
  HOST_DEVICE_CONSTANT float maxz[nPairs] = {
      19.5425, 44.1517, -23.8855, 22.4377, 49.9421, -31.1235, 62.9705, -61.5231, 28.2281, 49.9421, -39.8091, 71.6561, -71.6563, 12.3045, 13.7521, -71.6563, -71.6563, -71.6563, -71.6563};
  // HOST_DEVICE_CONSTANT float maxr[nPairs] = {
  //     3.6492, 4.5743, 4.5743, 4.6771, 6.4246, 6.4246, 2.1073, 2.1073, 5.7051, 5.4995, 5.4995, 2.2101, 2.2101, 8.2749, 10.2280, 0.0514, 0.0514, 0.0514, 0.0514};

  // // Using values from simPixelTracksAnalyser 99.5%
  // HOST_DEVICE_CONSTANT float dcaCuts[numberOfLayers] = {0.11, 0.14, 0.23, 0.0, 0.25, 0.0, 0.0, 0.25, 0.0, 0.0};
  // Using values from simPixelTracksAnalyser 90%
  HOST_DEVICE_CONSTANT float dcaCuts[numberOfLayers] = {0.0225, 0.0405, 0.1305, 0.0, 0.1515, 0.0, 0.0, 0.1515, 0.0, 0.0};

  // Using values from simPixelTracksAnalyser 99.5%
  HOST_DEVICE_CONSTANT float thetaCuts[numberOfLayers] = {
      0.0, 0.002, 0.002, 0.0, 0.002, 0.002, 0.0, 0.002, 0.002, 0.0};
  // // Using values from simPixelTracksAnalyser 90%
  // HOST_DEVICE_CONSTANT float thetaCuts[numberOfLayers] = {
  //     0.0, 0.0007, 0.0007, 0.0, 0.0009, 0.0013, 0.0, 0.0009, 0.0013, 0.0};

  static constexpr uint32_t layerStart[numberOfLayers + 1] = {0,
                                                              1,
                                                              2,
                                                              3,
                                                              4,
                                                              5,
                                                              6,
                                                              7,
                                                              8,
                                                              9,  // negative endcap
                                                              10};
}  // namespace colliderMLPhase1PixelTopology

namespace colliderMLPhase2PixelTopology {

  using pixelTopology::phi0p05;
  using pixelTopology::phi0p06;
  using pixelTopology::phi0p07;
  using pixelTopology::phi0p09;

  constexpr uint32_t numberOfLayers = 18;
  constexpr int nPairs = 19 + 12 + 10;  // include far forward layer pairs
  constexpr uint16_t numberOfModules = 18;

  constexpr uint32_t maxNumClustersPerModules = 1024;

  HOST_DEVICE_CONSTANT uint8_t layerPairs[2 * nPairs] = {

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

  // ------ Begin for SimPixelTracks geometry ------
  HOST_DEVICE_CONSTANT uint8_t startingPairs[nPairs] = {
    1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
  };

  HOST_DEVICE_CONSTANT bool isBarrel[numberOfLayers] = {
    1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
  };

  HOST_DEVICE_CONSTANT float ptCuts[nPairs] = {
    0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
    0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
    0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
    0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
    0.5
  };
  // ------ End for SimPixelTracks geometry ------

  HOST_DEVICE_CONSTANT uint32_t layerStart[numberOfLayers + 1] = {0,
                                                                  1,
                                                                  2,
                                                                  3, // Barrel
                                                                  4,
                                                                  5,
                                                                  6,
                                                                  7,
                                                                  8,
                                                                  9,
                                                                  10, // Fp
                                                                  11,
                                                                  12,
                                                                  13,
                                                                  14,
                                                                  15,
                                                                  16,
                                                                  17, // Np
                                                                  numberOfModules};

  // 99% simPixelTracks
  HOST_DEVICE_CONSTANT int16_t phicuts[nPairs]{
      //0,  1    0,  4    0,  11
      250, 330, 330, 
    //1,  2    1,  4    1,  11
      250, 370, 370,
    //2,  3    2,  4    2,  11
      310, 290, 290, 
      
    //4,  5    5,  6    6,  7    7,  8    8,  9    9,  10
      130, 130, 130, 110, 130, 130, 
    //11, 12   12, 13   13, 14   14, 15   15, 16   16, 17
      130, 130, 130, 110, 130, 130,
      
    //0,  2    0,  5    0,  12   0,  6    0,  13
      0, 0, 0, 0, 0,
    //1,  3    1,  5    1,  12   1,  6    1,  13 
      0, 0, 0, 0, 0, 
      
    //4,  6    5,  7    6,  8    7,  9    8,  10
      0, 0, 0, 0, 0, 
    //11, 13   12, 14   13, 15   14, 16   15, 17
      0, 0, 0, 0, 0};

  // // 90% simPixelTracks
  // HOST_DEVICE_CONSTANT int16_t phicuts[nPairs]{
  //     //0,  1    0,  4    0,  11
  //     170, 190, 190, 
  //   //1,  2    1,  4    1,  11
  //     230, 250, 250,
  //   //2,  3    2,  4    2,  11
  //     270, 230, 230, 
      
  //   //4,  5    5,  6    6,  7    7,  8    8,  9    9,  10
  //     90, 90, 90, 90, 110, 90, 
  //   //11, 12   12, 13   13, 14   14, 15   15, 16   16, 17
  //     90, 90, 90, 90, 110, 90,
      
  //   //0,  2    0,  5    0,  12   0,  6    0,  13
  //     0, 0, 0, 0, 0,
  //   //1,  3    1,  5    1,  12   1,  6    1,  13 
  //     0, 0, 0, 0, 0, 
      
  //   //4,  6    5,  7    6,  8    7,  9    8,  10
  //     0, 0, 0, 0, 0, 
  //   //11, 13   12, 14   13, 15   14, 16   15, 17
  //     0, 0, 0, 0, 0};

  // 99% simPixelTracks
  HOST_DEVICE_CONSTANT float minz[nPairs] = {
    //0,  1    0,  4    0,  11
      -25.1522,    19.8570,      -48.9806, 
    //1,  2     1,  4     1,  11      
      -30.4474,    27.7998,      -51.6282, 
    //2,  3     2,  4     2,  11      
      -35.7426,    38.3902,     -51.6282,
      
    //4,  5     5,  6     6,  7    7,  8    8,  9    9,  10      
      61.2186,     71.8090,     82.3994,    96.6374,    111.5230,    130.0560, 
    //11, 12    12, 13    13, 14   14, 15   15, 16   16, 17      
      -63.2186,    -73.8090,   -84.3994,  -99.2850,  -113.5230,  -132.0560, 
      
    //0,  2     0,  5   0,  12    0,  6    0,  13      
      -0.0,    0.0,    -0.0,    0.0,    -0.0,
    //1,  3     1,  5   1,  12    1,  6    1,  13       
      -0.0,    0.0,    -0.0,    0.0,    -0.0, 
      
    //4,  6     5,  7    6,  8    7,  9    8,  10      
      0.0,     0.0,    0.0,    0.0,    0.0, 
    //11, 13    12, 14   13, 15   14, 16   15, 17      
      -0.0,    -0.0,  -0.0,  -0.0,  -0.0};
  HOST_DEVICE_CONSTANT float maxz[nPairs] = {

    //0,  1    0,  4    0,  11
      25.1522,     48.9806,     -19.8570,
    //1,  2     1,  4     1,  11
      30.4474,     51.6282,     -27.7998,
    //2,  3     2,  4     2,  11
      35.7426,     51.6282,     -38.3902,
      
    //4,  5    5,  6    6,  7    7,  8    8,  9    9,  10
      63.2186,    73.8090,   84.3994,   99.2850,   113.5230,   132.0560, 
    //11, 12   12, 13   13, 14   14, 15   15, 16   16, 17
      -61.2186,   -71.8090,   -82.3994,   -96.6374,   -111.5230,   -130.0560,
      
    //0,  2    0,  5    0,  12   0,  6    0,  13 
      0.0,    0.0,    -0.0,    0.0,    -0.0,
    //1,  3    1,  5    1,  12   1,  6    1,  13
      0.0,    0.0,    -0.0,    0.0,    -0.0,
      
    //4,  6    5,  7    6,  8    7,  9    8,  10
      0.0,    0.0,   0.0,   0.0,   0.0,
    //11, 13   12, 14   13, 15   14, 16   15, 17
      -0.0,   -0.0,   -0.0,   -0.0,   -0.0};

  // // 90% simPixelTracks
  // HOST_DEVICE_CONSTANT float minz[nPairs] = {
  //   //0,  1    0,  4    0,  11
  //     -19.8570,    22.5046,      -43.6854, 
  //   //1,  2     1,  4     1,  11      
  //     -22.5046,    30.4474,      -48.9806, 
  //   //2,  3     2,  4     2,  11      
  //     -27.7998,    41.0378,     -48.9806,
      
  //   //4,  5     5,  6     6,  7    7,  8    8,  9    9,  10      
  //     61.2186,     71.8090,     82.3994,    96.6374,    111.5230,    130.0560, 
  //   //11, 12    12, 13    13, 14   14, 15   15, 16   16, 17      
  //     -63.2186,    -73.8090,   -84.3994,  -99.2850,  -113.5230,  -132.0560, 
      
  //   //0,  2     0,  5   0,  12    0,  6    0,  13      
  //     -0.0,    0.0,    -0.0,    0.0,    -0.0,
  //   //1,  3     1,  5   1,  12    1,  6    1,  13       
  //     -0.0,    0.0,    -0.0,    0.0,    -0.0, 
      
  //   //4,  6     5,  7    6,  8    7,  9    8,  10      
  //     0.0,     0.0,    0.0,    0.0,    0.0, 
  //   //11, 13    12, 14   13, 15   14, 16   15, 17      
  //     -0.0,    -0.0,  -0.0,  -0.0,  -0.0};
  // HOST_DEVICE_CONSTANT float maxz[nPairs] = {

  //   //0,  1    0,  4    0,  11
  //     19.8570,     43.6854,     -22.5046,
  //   //1,  2     1,  4     1,  11
  //     22.5046,     48.9806,     -30.4474,
  //   //2,  3     2,  4     2,  11
  //     27.7998,     48.9806,     -41.0378,
      
  //   //4,  5    5,  6    6,  7    7,  8    8,  9    9,  10
  //     63.2186,    73.8090,   84.3994,   99.2850,   113.5230,   132.0560, 
  //   //11, 12   12, 13   13, 14   14, 15   15, 16   16, 17
  //     -61.2186,   -71.8090,   -82.3994,   -96.6374,   -111.5230,   -130.0560,
      
  //   //0,  2    0,  5    0,  12   0,  6    0,  13 
  //     0.0,    0.0,    -0.0,    0.0,    -0.0,
  //   //1,  3    1,  5    1,  12   1,  6    1,  13
  //     0.0,    0.0,    -0.0,    0.0,    -0.0,
      
  //   //4,  6    5,  7    6,  8    7,  9    8,  10
  //     0.0,    0.0,   0.0,   0.0,   0.0,
  //   //11, 13   12, 14   13, 15   14, 16   15, 17
  //     -0.0,   -0.0,   -0.0,   -0.0,   -0.0};

  // 99% simPixelTracks
  HOST_DEVICE_CONSTANT float maxr[nPairs] = {
    //0,  1    0,  4    0,  11
    3.7237,    5.0873,    5.0873,
  //1,  2   1,  4   1,  11
    // 5.0,    7.0,    7.0,
    4.7727,    7.2901,    7.2901,
  //2,  3   2,  4   2,  11
    5.8216,    5.9265,    5.9265,

  //4,  5   5,  6   6,  7   7,  8   8,  9   9,  10
    2.5699,    2.5699,    2.5699,    2.2552,    2.6748,    2.3601,
  //11, 12  12, 13  13, 14  14, 15  15, 16  16, 17
    2.5699,    2.5699,    2.5699,    2.2552,    2.6748,    2.3601,

  //0,  2    0,  5   0,  12  0,  6   0,  13 
    8.2342,     0.0524,    0.0524,    0.0524,    0.0524,
  //1,  3    1,  5   1,  12  1,  6   1,  13
    10.2271,     0.0524,    0.0524,    0.0524,    2.9895,

  //4,  6   5,  7   6,  8   7,  9   8,  10
    0.0524,    0.0524,    0.0524,    0.0524,    0.0524,
  //11, 13  12, 14  13, 15  14, 16  15, 17
    0.0524,    1.7308,    0.0524,    0.0524,    1.5210};

  // // 90% simPixelTracks
  // HOST_DEVICE_CONSTANT float maxr[nPairs] = {
  //   //0,  1    0,  4    0,  11
  //   3.7237,    4.5629,    4.5629,
  // //1,  2   1,  4   1,  11
  //   4.6678,    6.4510,    6.4510,
  // //2,  3   2,  4   2,  11
  //   5.7167,    5.5069,    5.5069,

  // //4,  5   5,  6   6,  7   7,  8   8,  9   9,  10
  //   2.1503,    2.2552,    2.2552,    2.0454,    2.4650,    2.2552,
  // //11, 12  12, 13  13, 14  14, 15  15, 16  16, 17
  //   2.1503,    2.2552,    2.2552,    2.0454,    2.4650,    2.2552,

  // //0,  2    0,  5   0,  12  0,  6   0,  13 
  //   0.0524,     0.0524,    0.0524,    0.0524,    0.0524,
  // //1,  3    1,  5   1,  12  1,  6   1,  13
  //   0.0524,     0.0524,    0.0524,    0.0524,    0.0524,

  // //4,  6   5,  7   6,  8   7,  9   8,  10
  //   0.0524,    0.0524,    0.0524,    0.0524,    0.0524,
  // //11, 13  12, 14  13, 15  14, 16  15, 17
  //   0.0524,    0.0524,    0.0524,    0.0524,    0.0524};


  // // 99% simPixelTracks
  // HOST_DEVICE_CONSTANT float dcaCuts[numberOfLayers] = {0.1275, 0.1675, 0.2775, 0.0, 0.2925, 0.2975, 0.3175, 0.3175, 0.3175, 0.0, 0.0, 0.2925, 0.2925, 0.3225, 0.3275, 0.3125, 0.0, 0.0};
  // 90% simPixelTracks
  HOST_DEVICE_CONSTANT float dcaCuts[numberOfLayers] = {0.0225, 0.0425, 0.1325, 0.0, 0.1575, 0.1575, 0.1775, 0.1725, 0.1675, 0.0, 0.0, 0.1575, 0.1575, 0.1775, 0.1725, 0.1675, 0.0, 0.0};

  // 99% simPixelTracks
  HOST_DEVICE_CONSTANT float thetaCuts[numberOfLayers] = {0.0, 0.0022, 0.0018, 0.0, 0.0019, 0.0022, 0.0021, 0.0019, 0.0018, 0.0018, 0.0, 0.0019, 0.0022, 0.0021, 0.0019, 0.0018, 0.0018, 0.0};
  // // 90% simPixelTracks
  // HOST_DEVICE_CONSTANT float thetaCuts[numberOfLayers] = {0.0, 0.0008, 0.0006, 0.0, 0.0009, 0.0013, 0.0011, 0.0011, 0.0010, 0.0009, 0.0, 0.0009, 0.0013, 0.0011, 0.0011, 0.0010, 0.0009, 0.0};

}  // namespace colliderMLPhase2PixelTopology

namespace colliderMLPhase2PixelPlusShortStripsTopology {

  using pixelTopology::phi0p05;
  using pixelTopology::phi0p06;
  using pixelTopology::phi0p07;
  using pixelTopology::phi0p09;
  using pixelTopology::pH;

  constexpr uint32_t numberOfLayers = 34;
  constexpr int nPairs = 5 + 14 + 22 + 3 + 4 + 26 + 5 + 14 + 18;  // include far forward layer pairs
  constexpr uint16_t numberOfModules = 34;

  constexpr uint32_t maxNumClustersPerModules = 1024;

  HOST_DEVICE_CONSTANT uint8_t layerPairs[2 * nPairs] = {

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

  // ------ Begin for SimPixelTracks geometry ------
  HOST_DEVICE_CONSTANT uint8_t startingPairs[nPairs] = {
    1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
  };

  HOST_DEVICE_CONSTANT bool isBarrel[numberOfLayers] = {
    1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  };



  HOST_DEVICE_CONSTANT float ptCuts[nPairs] = {
    0.5, 0.5, 0.5, 0.5, 0.5, // pixel barrel only (5)

    0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, // pixel barrel-endcap pos (12)
    0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, // pixel barrel-endcap neg (19)

    0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, // pixel endcap only pos (30)
    0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, // pixel endcap only neg (41)

    0.5, 0.5, 0.5, // pixel-shortStrips barrel only (44)

    0.5, 0.5, // pixel-shortStrips endcap-barrel pos (46)
    0.5, 0.5, // pixel-shortStrips endcap-barrel neg (48)

    0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, // pixel-shortStrips endcap only pos (61)
    0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, // pixel-shortStrips endcap only neg (74)

    0.5, 0.5, 0.5, 0.5, 0.5, // shortStrips barrel only (79)

    0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, // shortStrips barrel-endcap pos (86)
    0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, // shortStrips barrel-endcap neg (93)

    0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, // shortStrips endcap only pos (102)
    0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, // shortStrips endcap only neg (111)
  };
  // ------ End for SimPixelTracks geometry ------

  HOST_DEVICE_CONSTANT uint32_t layerStart[numberOfLayers + 1] = {0,
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
                                                                  numberOfModules };

  // 99% simPixelTracks
  HOST_DEVICE_CONSTANT int16_t phicuts[nPairs]{
    // 99% simPixelTracks
    // 250, 250, 310, 550, 550, 310, 110, 10, 370, 10, 10, 310, 310, 90, 10, 350, 10, 10, 290, 130, 130, 130, 110, 130, 110, 10, 10, 10, 10, 10, 130, 130, 130, 110, 130, 130, 10, 10, 10, 10, 10, 490, 810, 1010, 610, 550, 590, 650, 630, 710, 530, 10, 710, 10, 550, 530, 670, 430, 630, 490, 610, 730, 690, 550, 10, 710, 10, 550, 690, 590, 410, 630, 450, 630, 550, 770, 890, 1290, 1630, 730, 10, 10, 1010, 1610, 10, 1010, 730, 10, 10, 1010, 1590, 10, 1010, 550, 550, 550, 470, 470, 850, 1010, 710, 410, 550, 550, 550, 470, 470, 1010, 630, 690, 10
    // 90% simPixelTracks
    170, 230, 270, 410, 510, 190, 110, 0, 250, 0, 0, 230, 190, 90, 0, 250, 0, 0, 230, 90, 90, 90, 90, 110, 90, 0, 0, 0, 0, 0, 90, 90, 90, 90, 110, 90, 0, 0, 0, 0, 0, 450, 730, 930, 510, 470, 490, 450, 630, 610, 450, 0, 630, 0, 470, 350, 530, 350, 530, 370, 530, 730, 610, 450, 0, 630, 0, 470, 350, 510, 370, 530, 370, 530, 490, 710, 830, 1230, 1490, 570, 0, 0, 810, 1470, 0, 790, 550, 0, 0, 790, 1490, 0, 790, 410, 430, 410, 350, 350, 850, 1010, 710, 410, 410, 410, 410, 370, 350, 950, 630, 690, 10
  };

  // 99% simPixelTracks
  HOST_DEVICE_CONSTANT float minz[nPairs] = {
    // 99% simPixelTracks
    -28.1325, -33.2475, -33.2475, -17.9025, -23.0175, 17.9025, 48.5925, -253.193, 28.1325, -253.193, -253.193, 38.3625, -48.5925, -48.5925, -253.193, -48.5925, -253.193, -253.193, -48.5925, 63.9375, 69.0525, 84.3975, 99.7425, 109.972, 130.432, -253.193, -253.193, -253.193, -253.193, -253.193, -63.9375, -74.1675, -84.3975, -99.7425, -109.973, -130.433, -253.193, -253.193, -253.193, -253.193, -253.193, -48.5925, -43.4775, -33.2475, 63.9375, 74.1675, -63.9375, -74.1675, 63.9375, 74.1675, 84.3975, -253.193, 84.3975, -253.193, 99.7425, 109.972, 109.972, 130.432, 130.432, 150.892, 150.892, -63.9375, -74.1675, -84.3975, -253.193, -84.3975, -253.193, -99.7425, -109.973, -109.973, -130.433, -130.433, -150.893, -150.893, -79.2825, -79.2825, -84.3975, -53.7075, -58.8225, 79.2825, -253.193, -253.193, 79.2825, 84.3975, -253.193, 89.5125, -115.088, -253.193, -253.193, -115.088, -89.5125, -253.193, -115.088, 130.432, 156.007, 186.697, 217.387, 253.192, 130.432, 156.007, 186.697, 222.502, -130.433, -156.008, -186.698, -222.503, -253.193, -130.433, -156.008, -186.698, -253.193

    // 90% simPixelTracks
    // -17.9025, -23.0175, -28.1325, -12.7875, -12.7875, 23.0175, 48.5925, -253.193, 28.1325, -253.193, -253.193, 38.3625, -43.4775, -48.5925, -253.193, -48.5925, -253.193, -253.193, -48.5925, 63.9375, 69.0525, 84.3975, 99.7425, 109.972, 130.432, -253.193, -253.193, -253.193, -253.193, -253.193, -63.9375, -74.1675, -84.3975, -99.7425, -109.973, -130.433, -253.193, -253.193, -253.193, -253.193, -253.193, -43.4775, -43.4775, -28.1325, 63.9375, 74.1675, -63.9375, -74.1675, 63.9375, 74.1675, 84.3975, -253.193, 84.3975, -253.193, 99.7425, 109.972, 109.972, 130.432, 130.432, 150.892, 150.892, -63.9375, -74.1675, -84.3975, -253.193, -84.3975, -253.193, -99.7425, -109.973, -109.973, -130.433, -130.433, -150.893, -150.893, -69.0525, -69.0525, -74.1675, -43.4775, -38.3625, 84.3975, -253.193, -253.193, 84.3975, 84.3975, -253.193, 94.6275, -109.973, -253.193, -253.193, -109.973, -89.5125, -253.193, -115.088, 130.432, 156.007, 186.697, 217.387, 253.192, 130.432, 156.007, 186.697, 222.502, -130.433, -156.008, -186.698, -222.503, -253.193, -130.433, -156.008, -186.698, -253.193
  };
  HOST_DEVICE_CONSTANT float maxz[nPairs] = {
    // 99% simPixelTracks
    28.1325, 33.2475, 33.2475, 17.9025, 23.0175, 48.5925, 48.5925, -253.193, 48.5925, -253.193, -253.193, 48.5925, -17.9025, -48.5925, -253.193, -28.1325, -253.193, -253.193, -38.3625, 63.9375, 74.1675, 84.3975, 99.7425, 109.972, 130.432, -253.193, -253.193, -253.193, -253.193, -253.193, -63.9375, -69.0525, -84.3975, -99.7425, -109.973, -130.433, -253.193, -253.193, -253.193, -253.193, -253.193, 48.5925, 43.4775, 38.3625, 63.9375, 74.1675, -63.9375, -74.1675, 63.9375, 74.1675, 84.3975, -253.193, 84.3975, -253.193, 99.7425, 109.972, 109.972, 130.432, 130.432, 150.892, 150.892, -63.9375, -74.1675, -84.3975, -253.193, -84.3975, -253.193, -99.7425, -109.973, -109.973, -130.433, -130.433, -150.893, -150.893, 79.2825, 79.2825, 84.3975, 63.9375, 79.2825, 115.087, -253.193, -253.193, 115.087, 94.6275, -253.193, 115.087, -79.2825, -253.193, -253.193, -79.2825, -79.2825, -253.193, -89.5125, 130.432, 156.007, 186.697, 222.502, 253.192, 130.432, 156.007, 186.697, 222.502, -130.433, -156.008, -186.698, -217.388, -253.193, -130.433, -156.008, -186.698, -253.193

    // 90% simPixelTracks
    // 17.9025, 23.0175, 28.1325, 12.7875, 17.9025, 43.4775, 48.5925, -253.193, 48.5925, -253.193, -253.193, 48.5925, -23.0175, -48.5925, -253.193, -28.1325, -253.193, -253.193, -38.3625, 63.9375, 74.1675, 84.3975, 99.7425, 109.972, 130.432, -253.193, -253.193, -253.193, -253.193, -253.193, -63.9375, -69.0525, -84.3975, -99.7425, -109.973, -130.433, -253.193, -253.193, -253.193, -253.193, -253.193, 43.4775, 43.4775, 28.1325, 63.9375, 74.1675, -63.9375, -74.1675, 63.9375, 74.1675, 84.3975, -253.193, 84.3975, -253.193, 99.7425, 109.972, 109.972, 130.432, 130.432, 150.892, 150.892, -63.9375, -74.1675, -84.3975, -253.193, -84.3975, -253.193, -99.7425, -109.973, -109.973, -130.433, -130.433, -150.893, -150.893, 63.9375, 69.0525, 74.1675, 43.4775, 53.7075, 109.972, -253.193, -253.193, 109.972, 94.6275, -253.193, 115.087, -84.3975, -253.193, -253.193, -84.3975, -79.2825, -253.193, -94.6275, 130.432, 156.007, 186.697, 222.502, 253.192, 130.432, 156.007, 186.697, 222.502, -130.433, -156.008, -186.698, -217.388, -253.193, -130.433, -156.008, -186.698, -253.193
  };

  HOST_DEVICE_CONSTANT float maxr[nPairs] = {
    // // 99% simPixelTracks
    // 3.87207, 4.80137, 5.73067, 8.2088, 10.3772, 5.11114, 2.01348, 0.154883, 7.2795, 0.154883, 0.154883, 6.04044, 5.11114, 2.01348, 0.154883, 7.2795, 0.154883, 0.154883, 6.04044, 2.63301, 2.63301, 2.63301, 2.32324, 2.63301, 2.32324, 0.154883, 0.154883, 0.154883, 0.154883, 0.154883, 2.63301, 2.63301, 2.63301, 2.32324, 2.63301, 2.32324, 0.154883, 0.154883, 0.154883, 0.154883, 0.154883, 9.44786, 15.0236, 19.3604, 11.6162, 10.6869, 11.3065, 10.6869, 14.0944, 14.0944, 10.3772, 0.154883, 13.7846, 0.154883, 10.9967, 8.51856, 11.926, 8.82833, 11.6162, 8.2088, 11.6162, 14.0944, 14.0944, 10.6869, 0.154883, 14.0944, 0.154883, 10.9967, 7.58927, 11.926, 7.58927, 11.926, 8.2088, 11.6162, 10.6869, 14.7139, 16.5725, 24.6264, 30.5119, 15.0236, 0.154883, 0.154883, 20.5994, 30.8217, 0.154883, 19.9799, 15.0236, 0.154883, 0.154883, 20.5994, 29.5826, 0.154883, 19.9799, 11.3065, 11.3065, 10.9967, 9.44786, 9.44786, 19.3604, 19.3604, 17.5018, 15.6432, 11.3065, 11.3065, 10.9967, 9.44786, 9.44786, 19.6701, 19.6701, 20.2897, 0.154883

    // 90% simPixelTracks
    3.56231, 4.80137, 5.73067, 8.2088, 10.3772, 4.49161, 2.01348, 0.154883, 6.3502, 0.154883, 0.154883, 5.4209, 4.49161, 2.01348, 0.154883, 6.3502, 0.154883, 0.154883, 5.4209, 2.01348, 2.32324, 2.32324, 2.01348, 2.63301, 2.32324, 0.154883, 0.154883, 0.154883, 0.154883, 0.154883, 2.01348, 2.32324, 2.32324, 2.01348, 2.63301, 2.32324, 0.154883, 0.154883, 0.154883, 0.154883, 0.154883, 9.44786, 15.0236, 19.3604, 10.9967, 9.75763, 10.9967, 10.0674, 14.0944, 13.1651, 9.75763, 0.154883, 13.1651, 0.154883, 10.0674, 7.58927, 11.3065, 7.2795, 11.3065, 7.89903, 10.9967, 14.0944, 13.4748, 9.75763, 0.154883, 13.4748, 0.154883, 10.0674, 7.2795, 11.3065, 7.2795, 11.3065, 7.89903, 10.9967, 10.3772, 14.4041, 16.2627, 24.3166, 30.2022, 13.7846, 0.154883, 0.154883, 19.0506, 30.2022, 0.154883, 18.7408, 13.7846, 0.154883, 0.154883, 19.0506, 29.5826, 0.154883, 18.7408, 10.0674, 10.0674, 9.75763, 8.51856, 8.51856, 19.3604, 19.3604, 17.5018, 15.6432, 10.0674, 10.0674, 9.75763, 8.51856, 8.51856, 19.6701, 19.6701, 20.2897, 0.154883
  };

  // 99% simPixelTracks
  // HOST_DEVICE_CONSTANT float dcaCuts[numberOfLayers] = {0.1275, 0.1625, 0.2175, 0.2475, 0.2875, 0.2825, 0.3025, 0.3225, 0.3225, 0.3825, 0.3775, 0.2925, 0.2825, 0.3025, 0.2925, 0.3175, 0.3825, 0.4025, 0.3625, 0.4125, 0.4825, 0.0025, 0.4625, 0.4575, 0.4675, 0.4775, 0.0025, 0.0025, 0.4525, 0.4625, 0.4675, 0.4775, 0.0025, 0.0025};
  // 90% simPixelTracks
  HOST_DEVICE_CONSTANT float dcaCuts[numberOfLayers] = {0.0175, 0.0425, 0.0625, 0.0825, 0.1425, 0.1425, 0.1625, 0.1575, 0.1575, 0.1625, 0.1725, 0.1425, 0.1425, 0.1575, 0.1575, 0.1525, 0.1675, 0.1675, 0.1575, 0.2125, 0.3775, 0.0025, 0.2925, 0.2925, 0.3275, 0.3525, 0.0025, 0.0025, 0.2925, 0.2975, 0.3225, 0.3525, 0.0025, 0.0025};

  // 99% simPixelTracks
  // HOST_DEVICE_CONSTANT float thetaCuts[numberOfLayers] = {5e-05, 0.00225, 0.00205, 0.00215, 0.00195, 0.00175, 0.00175, 0.00165, 0.00185, 0.00225, 0.00375, 0.00185, 0.00175, 0.00175, 0.00165, 0.00175, 0.00225, 0.00475, 0.00305, 0.00335, 0.00465, 5e-05, 0.00425, 0.00405, 0.00415, 0.00425, 0.00435, 5e-05, 0.00415, 0.00405, 0.00415, 0.00435, 0.00445, 5e-05};
  // 90% simPixelTracks
  HOST_DEVICE_CONSTANT float thetaCuts[numberOfLayers] = {5e-05, 0.00075, 0.00065, 0.00095, 0.00075, 0.00095, 0.00095, 0.00095, 0.00085, 0.00085, 0.00125, 0.00075, 0.00095, 0.00095, 0.00085, 0.00085, 0.00085, 0.00135, 0.00155, 0.00205, 0.00315, 5e-05, 0.00265, 0.00265, 0.00275, 0.00275, 0.00295, 5e-05, 0.00265, 0.00265, 0.00275, 0.00285, 0.00295, 5e-05};

}  // namespace colliderMLPhase2PixelPlusShortStripsTopology

namespace colliderMLPhase2AllTrackerTopology {

  constexpr uint32_t numberOfLayers = 48;
  constexpr int nPairs = 171;  // include far forward layer pairs
  constexpr uint16_t numberOfModules = 48;

  constexpr uint32_t maxNumClustersPerModules = 1024;

  HOST_DEVICE_CONSTANT uint8_t layerPairs[2 * nPairs] = {

      0,1,0,2,0,4,0,5,0,11,0,12,1,2,1,3,1,4,1,5,1,11,1,12,1,18,2,3,2,4,2,5,2,11,2,12,2,18,2,19,2,22,2,28,3,18,3,19,4,5,4,6,4,18,4,22,5,6,5,7,5,18,5,22,5,23,6,7,6,8,6,22,6,23,7,8,7,9,7,23,7,24,8,9,8,10,8,24,8,25,9,10,9,25,9,26,9,27,10,26,10,27,11,12,11,13,11,18,11,28,12,13,12,14,12,18,12,28,12,29,13,14,13,15,13,28,13,29,14,15,14,16,14,29,14,30,15,16,15,17,15,30,15,31,16,17,16,31,16,32,16,33,17,32,17,33,18,19,18,20,18,22,18,23,18,28,18,29,19,20,19,21,19,22,19,23,19,28,19,29,19,37,19,43,20,21,20,22,20,23,20,28,20,29,20,34,20,36,20,37,20,38,20,42,20,43,20,44,21,34,21,35,21,36,21,37,21,42,21,43,22,23,22,24,22,37,22,38,22,39,23,24,23,25,23,38,23,39,23,40,24,25,24,26,24,39,24,40,24,41,25,26,25,27,25,40,25,41,26,27,26,41,28,29,28,30,28,43,28,44,28,45,29,30,29,31,29,44,29,45,29,46,30,31,30,32,30,45,30,46,30,47,31,32,31,33,31,46,31,47,32,33,32,47,34,35,34,36,34,42,36,37,37,38,37,39,38,39,38,40,39,40,39,41,40,41,42,43,43,44,43,45,44,45,44,46,45,46,45,47,46,47,
      
  };

  // ------ Begin for SimPixelTracks geometry ------
  HOST_DEVICE_CONSTANT uint8_t startingPairs[nPairs] = {
    1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  };

  HOST_DEVICE_CONSTANT bool isBarrel[numberOfLayers] = {
    1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,
  };

  HOST_DEVICE_CONSTANT float ptCuts[nPairs] = {
    0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,
  };
  // ------ End for SimPixelTracks geometry ------

  HOST_DEVICE_CONSTANT uint32_t layerStart[numberOfLayers + 1] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48};

  HOST_DEVICE_CONSTANT int16_t phicuts[nPairs]{
    250, 550, 310, 110, 310, 90, 250, 550, 370, 10, 350, 10, 1070, 310, 310, 10, 290, 10, 810, 10, 10, 10, 490, 1010, 130, 10, 610, 630, 130, 10, 550, 710, 10, 130, 10, 530, 710, 110, 10, 550, 810, 130, 10, 670, 730, 110, 630, 10, 10, 610, 730, 130, 10, 590, 730, 130, 10, 650, 690, 10, 130, 10, 550, 710, 110, 10, 550, 1410, 130, 10, 590, 1110, 130, 630, 10, 10, 630, 750, 550, 1290, 730, 10, 730, 10, 770, 1630, 1010, 1610, 1010, 1590, 1950, 10, 890, 1010, 10, 1010, 10, 1810, 1490, 1990, 10, 1350, 1990, 10, 930, 1990, 1450, 1430, 1470, 1530, 550, 850, 850, 1510, 10, 550, 1010, 810, 1510, 10, 550, 710, 790, 1330, 10, 470, 410, 670, 1230, 470, 650, 550, 1010, 850, 1510, 10, 550, 630, 830, 1490, 10, 550, 690, 790, 1330, 10, 470, 10, 670, 1230, 470, 650, 1230, 1430, 1370, 1230, 1030, 1470, 970, 1430, 830, 1510, 810, 1210, 1030, 1250, 990, 1490, 830, 1550, 810
  };

  HOST_DEVICE_CONSTANT float minz[nPairs] = {
    -24.6296, -19.3786, 17.3784, 48.8844, -50.8846, -50.8846, -29.8806, -24.6296, 27.8804, -260.925, -50.8846, -260.925, -24.6296, -35.1316, 38.3824, -260.925, -50.8846, -260.925, -45.6336, -260.925, -260.925, -260.925, -50.8846, -35.1316, 59.3864, -260.925, 59.3864, 59.3864, 69.8884, -260.925, 69.8884, 69.8884, -260.925, 80.3904, -260.925, 85.6414, 85.6414, 96.1434, -260.925, 96.1434, 96.1434, 111.896, -260.925, 111.896, 111.896, 132.9, 132.9, -260.925, -260.925, 148.653, 148.653, -61.3866, -260.925, -61.3866, -61.3866, -71.8886, -260.925, -71.8886, -71.8886, -260.925, -87.6416, -260.925, -87.6416, -87.6416, -98.1436, -260.925, -98.1436, -98.1436, -113.897, -260.925, -113.897, -113.897, -134.901, -134.901, -260.925, -260.925, -155.905, -155.905, -82.3906, -56.1356, 80.3904, -260.925, -113.897, -260.925, -82.3906, -56.1356, 80.3904, 80.3904, -113.897, -87.6416, 64.6374, -260.925, -87.6416, 85.6414, -260.925, -113.897, -260.925, -82.3906, 75.1394, 80.3904, 85.6414, -92.8926, -98.1436, -87.6416, -92.8926, -92.8926, 90.8924, 96.1434, -113.897, -113.897, 127.649, 127.649, 127.649, 127.649, -260.925, 153.904, 153.904, 153.904, 153.904, -260.925, 185.41, 185.41, 185.41, 185.41, -260.925, 216.916, 222.167, 216.916, 216.916, 253.673, 253.673, -129.65, -129.65, -129.65, -129.65, -260.925, -155.905, -155.905, -155.905, -155.905, -260.925, -187.411, -187.411, -187.411, -187.411, -260.925, -224.168, -260.925, -218.917, -218.917, -255.674, -255.674, -87.6416, 64.6374, -108.646, 127.649, 153.904, 159.155, 185.41, 185.41, 222.167, 222.167, 258.924, -129.65, -161.156, -161.156, -192.662, -192.662, -229.419, -224.168, -260.925
  };

  HOST_DEVICE_CONSTANT float maxz[nPairs] = {
    24.6294, 19.3784, 50.8844, 50.8844, -17.3786, -48.8846, 29.8804, 24.6294, 50.8844, -258.925, -27.8806, -258.925, 29.8804, 35.1314, 50.8844, -258.925, -38.3826, -258.925, 45.6334, -258.925, -258.925, -258.925, 50.8844, 40.3824, 61.3864, -258.925, 61.3864, 61.3864, 71.8884, -258.925, 71.8884, 71.8884, -258.925, 87.6414, -258.925, 87.6414, 87.6414, 98.1434, -258.925, 98.1434, 98.1434, 113.896, -258.925, 113.896, 113.896, 134.9, 134.9, -258.925, -258.925, 155.904, 155.904, -59.3866, -258.925, -59.3866, -59.3866, -69.8886, -258.925, -69.8886, -69.8886, -258.925, -80.3906, -258.925, -85.6416, -85.6416, -96.1436, -258.925, -96.1436, -96.1436, -111.897, -258.925, -111.897, -111.897, -132.901, -132.901, -258.925, -258.925, -148.654, -148.654, 82.3904, 66.6374, 113.896, -258.925, -80.3906, -258.925, 82.3904, 61.3864, 113.896, 92.8924, -80.3906, -80.3906, 71.8884, -258.925, 87.6414, 113.896, -258.925, -90.8926, -258.925, 71.8884, 92.8924, 98.1434, 92.8924, -64.6376, -80.3906, -85.6416, 92.8924, 8.87644, 113.896, 113.896, -90.8926, -106.646, 129.649, 129.649, 129.649, 129.649, -258.925, 155.904, 155.904, 155.904, 155.904, -258.925, 187.41, 187.41, 187.41, 187.41, -258.925, 224.167, 224.167, 218.916, 218.916, 255.673, 255.673, -127.65, -127.65, -127.65, -127.65, -258.925, -153.905, -153.905, -153.905, -153.905, -258.925, -185.411, -185.411, -185.411, -185.411, -258.925, -216.917, -258.925, -216.917, -216.917, -253.674, -253.674, 87.6414, 108.645, -96.1436, 129.649, 161.155, 161.155, 192.661, 192.661, 229.418, 229.418, 260.924, -127.65, -153.905, -159.156, -185.411, -185.411, -222.168, -222.168, -258.925
  };

  HOST_DEVICE_CONSTANT float maxr[nPairs] = {
    3.47759, 8.53591, 5.37446, 2.21301, 5.37446, 1.58072, 4.74217, 10.4328, 7.27133, 0.316145, 7.27133, 0.316145, 19.2848, 6.00675, 6.00675, 0.316145, 6.00675, 0.316145, 14.8588, 0.316145, 0.316145, 0.316145, 9.1682, 19.2848, 2.21301, 0.316145, 11.6974, 14.2265, 2.8453, 0.316145, 10.4328, 14.2265, 0.316145, 2.8453, 0.316145, 10.4328, 13.5942, 2.21301, 0.316145, 11.0651, 13.5942, 2.8453, 0.316145, 11.6974, 13.5942, 2.21301, 11.6974, 0.316145, 0.316145, 11.6974, 13.5942, 2.8453, 0.316145, 11.6974, 14.2265, 2.8453, 0.316145, 11.0651, 14.2265, 0.316145, 2.8453, 0.316145, 10.4328, 14.2265, 2.21301, 0.316145, 11.0651, 13.5942, 2.8453, 0.316145, 11.6974, 13.5942, 2.21301, 11.6974, 0.316145, 0.316145, 11.6974, 13.5942, 10.4328, 24.3432, 14.8588, 0.316145, 14.8588, 0.316145, 14.8588, 30.6661, 20.5494, 30.6661, 20.5494, 29.4015, 62.9129, 0.316145, 16.7557, 19.9171, 0.316145, 19.9171, 0.316145, 33.1952, 32.5629, 49.6348, 49.0025, 49.0025, 49.6348, 49.0025, 17.388, 37.6213, 33.1952, 33.8275, 33.1952, 33.1952, 11.0651, 19.2848, 21.814, 23.7109, 0.316145, 11.0651, 19.2848, 21.1817, 23.7109, 0.316145, 11.0651, 17.388, 20.5494, 23.0786, 0.316145, 9.1682, 15.4911, 19.2848, 21.814, 9.1682, 18.6526, 11.0651, 19.9171, 21.814, 23.7109, 0.316145, 11.0651, 19.9171, 21.1817, 23.7109, 0.316145, 11.0651, 20.5494, 20.5494, 23.0786, 0.316145, 9.80049, 0.316145, 19.2848, 21.814, 9.1682, 18.6526, 21.814, 18.0203, 18.0203, 17.388, 17.388, 16.7557, 17.388, 17.388, 17.388, 17.388, 17.388, 17.388, 17.388, 16.7557, 17.388, 17.388, 17.388, 17.388, 17.388
  };

  HOST_DEVICE_CONSTANT float dcaCuts[numberOfLayers] = {0.1275, 0.1625, 0.2175, 0.2475, 0.2875, 0.2825, 0.3025, 0.3225, 0.3225, 0.3575, 0.3525, 0.2925, 0.2825, 0.3025, 0.3025, 0.3175, 0.3825, 0.3925, 0.3625, 0.4175, 0.4875, 0.4925, 0.4725, 0.4675, 0.4725, 0.4775, 0.0025, 0.0025, 0.4675, 0.4725, 0.4725, 0.4775, 0.0025, 0.0025, 0.2775, 0.0025, 0.0025, 0.0025, 0.0025, 0.0025, 0.0025, 0.0025, 0.0025, 0.0025, 0.0025, 0.0025, 0.0025, 0.0025};

  HOST_DEVICE_CONSTANT float thetaCuts[numberOfLayers] = {5e-05, 0.00225, 0.00205, 0.00215, 0.00195, 0.00175, 0.00175, 0.00165, 0.00185, 0.00225, 0.00305, 0.00185, 0.00175, 0.00175, 0.00165, 0.00175, 0.00225, 0.00375, 0.00305, 0.00335, 0.00475, 0.00995, 0.00715, 0.00705, 0.00825, 0.00805, 0.00505, 5e-05, 0.00735, 0.00635, 0.00825, 0.00775, 0.00515, 5e-05, 0.00985, 5e-05, 0.00995, 0.00985, 0.00995, 0.00985, 0.00995, 5e-05, 0.00995, 0.00995, 0.00985, 0.00995, 0.00975, 5e-05};

}  // namespace colliderMLPhase2AllTrackerTopology

namespace pixelTopology {

  struct Phase2 {
    // types
    using hindex_type = uint32_t;  // FIXME from siPixelRecHitsHeterogeneousProduct
    using tindex_type = uint32_t;  // for tuples
    using cindex_type = uint32_t;  // for cells

    static constexpr uint32_t maxCellNeighbors = 256;
    static constexpr uint32_t maxCellTracks = 302;
    static constexpr uint32_t maxHitsOnTrack = 20;
    static constexpr uint32_t maxHitsOnTrackForFullFit = 6;
    static constexpr uint32_t avgHitsPerTrack = 7;
    static constexpr uint32_t maxCellsPerHit = 256;
    static constexpr uint32_t avgTracksPerHit = 10;
    static constexpr uint32_t maxNumberOfTuples = 256 * 1024;
    // this is well above thanks to maxNumberOfTuples
    static constexpr uint32_t maxHitsForContainers = avgHitsPerTrack * maxNumberOfTuples;
    static constexpr uint32_t maxNumberOfDoublets = 5 * 512 * 1024;
    static constexpr uint32_t maxNumOfActiveDoublets = maxNumberOfDoublets / 2;
    static constexpr uint32_t maxNumberOfQuadruplets = maxNumberOfTuples;
    static constexpr uint32_t maxDepth = 12;
    static constexpr uint32_t numberOfLayers = 28;

    static constexpr uint32_t maxSizeCluster = 2047;

    static constexpr uint32_t getDoubletsFromHistoMaxBlockSize = 128;  // for both x and y
    static constexpr uint32_t getDoubletsFromHistoMinBlocksPerMP = 16;

    static constexpr uint16_t last_bpix1_detIndex = 216;
    static constexpr uint16_t last_bpix2_detIndex = 432;
    static constexpr uint16_t last_barrel_detIndex = 864;

    static constexpr uint32_t maxPixInModule = 6000;
    static constexpr uint32_t maxPixInModuleForMorphing = 0;
    static constexpr uint32_t maxIterClustering = 16;

    static constexpr uint32_t maxNumClustersPerModules = phase2PixelTopology::maxNumClustersPerModules;
    static constexpr uint32_t maxHitsInModule = phase2PixelTopology::maxNumClustersPerModules;

    static constexpr float moduleLength = 4.345f;
    static constexpr float endcapCorrection = 0.0f;

    static constexpr float xerr_barrel_l1_def = 0.00035f;
    static constexpr float yerr_barrel_l1_def = 0.00125f;
    static constexpr float xerr_barrel_ln_def = 0.00035f;
    static constexpr float yerr_barrel_ln_def = 0.00125f;
    static constexpr float xerr_endcap_def = 0.00060f;
    static constexpr float yerr_endcap_def = 0.00180f;

    static constexpr float bigPixXCorrection = 0.0f;
    static constexpr float bigPixYCorrection = 0.0f;

    static constexpr float dzdrFact = 8 * 0.0285 / 0.015;  // from dz/dr to "DY"

    static constexpr int nPairsMinimal = 33;
    static constexpr int nPairsFarForwards = nPairsMinimal + 8;  // include barrel "jumping" layer pairs
    static constexpr int nPairs = phase2PixelTopology::nPairs;   // include far forward layer pairs

    static constexpr int maxDYsize12 = 12;
    static constexpr int maxDYsize = 10;
    static constexpr int maxDYPred = 20;

    static constexpr uint16_t numberOfModules = phase2PixelTopology::numberOfModules;

    // 1000 bins < 1024 bins (10 bits) must be:
    // - < 32*32 (warpSize*warpSize for block prefix scan for CUDA)
    // - > number of columns (y) in any module. This is due to the fact
    //     that in pixel clustering we give for granted that in each
    //     bin we only have the pixel belonging to the same column.
    //     See RecoLocalTracker/SiPixelClusterizer/plugins/alpaka/PixelClustering.h#L325-L347
    static constexpr uint16_t clusterBinning = 1000;
    static constexpr uint16_t clusterBits = 10;

    static constexpr uint16_t numberOfModulesInBarrel = 756;
    static constexpr uint16_t numberOfModulesInLadder = 9;
    static constexpr uint16_t numberOfLaddersInBarrel = numberOfModulesInBarrel / numberOfModulesInLadder;

    static constexpr uint16_t firstEndcapPos = 4;
    static constexpr uint16_t firstEndcapNeg = 16;

    static constexpr int16_t xOffset = -1e4;  // not used actually, to suppress static analyzer warnings

    static constexpr char const *nameModifier = "Phase2";
    static constexpr char const *cpeModules = "PixelCPEFastParamsPhase2";

    static constexpr uint32_t const *layerStart = phase2PixelTopology::layerStart;
    static constexpr float const *minz = phase2PixelTopology::minz;
    static constexpr float const *maxz = phase2PixelTopology::maxz;
    static constexpr float const *maxr = phase2PixelTopology::maxr;

    static constexpr uint8_t const *layerPairs = phase2PixelTopology::layerPairs;
    static constexpr int16_t const *phicuts = phase2PixelTopology::phicuts;
    static constexpr float const *thetaCuts = phase2PixelTopology::thetaCuts;
    static constexpr float const *dcaCuts = phase2PixelTopology::dcaCuts;

    static constexpr inline bool isBigPixX(uint16_t px) { return false; }
    static constexpr inline bool isBigPixY(uint16_t py) { return false; }

    static constexpr inline uint16_t localX(uint16_t px) { return px; }
    static constexpr inline uint16_t localY(uint16_t py) { return py; }
  };

  struct Phase1 {
    // types
    using hindex_type = uint32_t;  // FIXME from siPixelRecHitsHeterogeneousProduct
    using tindex_type = uint16_t;  // for tuples
    using cindex_type = uint32_t;  // for cells

    static constexpr uint32_t maxCellNeighbors = 36;
    static constexpr uint32_t maxCellTracks = 48;
    static constexpr uint32_t maxHitsOnTrack = 10;
    static constexpr uint32_t maxHitsOnTrackForFullFit = 6;
    static constexpr uint32_t avgHitsPerTrack = 5;
    static constexpr uint32_t maxCellsPerHit = 256;
    static constexpr uint32_t avgTracksPerHit = 6;
    static constexpr uint32_t maxNumberOfTuples = 32 * 1024;
    static constexpr uint32_t maxHitsForContainers = avgHitsPerTrack * maxNumberOfTuples;
    static constexpr uint32_t maxNumberOfDoublets = 512 * 1024;
    static constexpr uint32_t maxNumOfActiveDoublets = maxNumberOfDoublets / 8;
    static constexpr uint32_t maxNumberOfQuadruplets = maxNumberOfTuples;
    static constexpr uint32_t maxDepth = 6;
    static constexpr uint32_t numberOfLayers = 10;

    static constexpr uint32_t maxSizeCluster = 1023;

    static constexpr uint32_t getDoubletsFromHistoMaxBlockSize = 64;  // for both x and y
    static constexpr uint32_t getDoubletsFromHistoMinBlocksPerMP = 16;

    static constexpr uint16_t last_bpix1_detIndex = 96;
    static constexpr uint16_t last_bpix2_detIndex = 320;
    static constexpr uint16_t last_barrel_detIndex = 1184;

    static constexpr uint32_t maxPixInModule = 6000;
    static constexpr uint32_t maxPixInModuleForMorphing = maxPixInModule * 2 / 5;
    static constexpr uint32_t maxIterClustering = 24;

    static constexpr uint32_t maxNumClustersPerModules = phase1PixelTopology::maxNumClustersPerModules;
    static constexpr uint32_t maxHitsInModule = phase1PixelTopology::maxNumClustersPerModules;

    static constexpr float moduleLength = 6.7f;
    static constexpr float endcapCorrection = 1.5f;

    static constexpr float xerr_barrel_l1_def = 0.00200f;
    static constexpr float yerr_barrel_l1_def = 0.00210f;
    static constexpr float xerr_barrel_ln_def = 0.00200f;
    static constexpr float yerr_barrel_ln_def = 0.00210f;
    static constexpr float xerr_endcap_def = 0.0020f;
    static constexpr float yerr_endcap_def = 0.00210f;

    static constexpr float bigPixXCorrection = 1.0f;
    static constexpr float bigPixYCorrection = 8.0f;

    static constexpr float dzdrFact = 8 * 0.0285 / 0.015;  // from dz/dr to "DY"

    static constexpr int nPairsForQuadruplets = 13;                     // quadruplets require hits in all layers
    static constexpr int nPairsForTriplets = nPairsForQuadruplets + 2;  // include barrel "jumping" layer pairs
    static constexpr int nPairs = nPairsForTriplets + 4;                // include forward "jumping" layer pairs

    static constexpr int maxDYsize12 = 28;
    static constexpr int maxDYsize = 20;
    static constexpr int maxDYPred = 20;

    static constexpr uint16_t numberOfModules = phase1PixelTopology::numberOfModules;

    static constexpr uint16_t numRowsInRoc = 80;
    static constexpr uint16_t numColsInRoc = 52;
    static constexpr uint16_t lastRowInRoc = numRowsInRoc - 1;
    static constexpr uint16_t lastColInRoc = numColsInRoc - 1;

    static constexpr uint16_t numRowsInModule = 2 * numRowsInRoc;
    static constexpr uint16_t numColsInModule = 8 * numColsInRoc;
    static constexpr uint16_t lastRowInModule = numRowsInModule - 1;
    static constexpr uint16_t lastColInModule = numColsInModule - 1;

    // 418 bins < 512, 9 bits are enough
    static constexpr uint16_t clusterBinning = numColsInModule + 2;
    static constexpr uint16_t clusterBits = 9;

    static constexpr uint16_t numberOfModulesInBarrel = 1184;
    static constexpr uint16_t numberOfModulesInLadder = 8;
    static constexpr uint16_t numberOfLaddersInBarrel = numberOfModulesInBarrel / numberOfModulesInLadder;

    static constexpr uint16_t firstEndcapPos = 4;
    static constexpr uint16_t firstEndcapNeg = 7;

    static constexpr int16_t xOffset = -81;

    static constexpr char const *nameModifier = "Phase1";
    static constexpr char const *cpeModules = "PixelCPEFastParams";

    static constexpr uint32_t const *layerStart = phase1PixelTopology::layerStart;
    static constexpr float const *minz = phase1PixelTopology::minz;
    static constexpr float const *maxz = phase1PixelTopology::maxz;
    static constexpr float const *maxr = phase1PixelTopology::maxr;

    static constexpr uint8_t const *layerPairs = phase1PixelTopology::layerPairs;
    static constexpr int16_t const *phicuts = phase1PixelTopology::phicuts;
    static constexpr float const *thetaCuts = phase1PixelTopology::thetaCuts;
    static constexpr float const *dcaCuts = phase1PixelTopology::dcaCuts;

    static constexpr inline bool isEdgeX(uint16_t px) { return (px == 0) | (px == lastRowInModule); }

    static constexpr inline bool isEdgeY(uint16_t py) { return (py == 0) | (py == lastColInModule); }

    static constexpr inline uint16_t toRocX(uint16_t px) { return (px < numRowsInRoc) ? px : px - numRowsInRoc; }

    static constexpr inline uint16_t toRocY(uint16_t py) {
      auto roc = divu52(py);
      return py - 52 * roc;
    }

    static constexpr inline bool isBigPixX(uint16_t px) { return (px == 79) | (px == 80); }
    static constexpr inline bool isBigPixY(uint16_t py) {
      auto ly = toRocY(py);
      return (ly == 0) | (ly == lastColInRoc);
    }

    static constexpr inline uint16_t localX(uint16_t px) {
      auto shift = 0;
      if (px > lastRowInRoc)
        shift += 1;
      if (px > numRowsInRoc)
        shift += 1;
      return px + shift;
    }

    static constexpr inline uint16_t localY(uint16_t py) {
      auto roc = divu52(py);
      auto shift = 2 * roc;
      auto yInRoc = py - 52 * roc;
      if (yInRoc > 0)
        shift += 1;
      return py + shift;
    }
  };

  struct HIonPhase1 : public Phase1 {
    // Storing here the needed constants different w.r.t. pp Phase1 topology.
    // All the other defined by inheritance in the HIon topology struct.

    using tindex_type = uint32_t;  // for tuples

    static constexpr uint32_t maxCellNeighbors = 90;
    static constexpr uint32_t maxCellTracks = 90;
    static constexpr uint32_t maxNumberOfTuples = 256 * 1024;
    static constexpr uint32_t maxNumberOfDoublets = 6 * 512 * 1024;
    static constexpr uint32_t maxHitsForContainers = avgHitsPerTrack * maxNumberOfTuples;
    static constexpr uint32_t maxNumberOfQuadruplets = maxNumberOfTuples;

    static constexpr uint32_t maxPixInModule = 10000;
    static constexpr uint32_t maxPixInModuleForMorphing = maxPixInModule * 1 / 10;
    static constexpr uint32_t maxIterClustering = 32;

    static constexpr uint32_t maxNumOfActiveDoublets =
        maxNumberOfDoublets / 4;  // TODO need to think a better way to avoid this duplication
    static constexpr uint32_t maxCellsPerHit = 256;

    static constexpr uint32_t maxNumClustersPerModules = phase1HIonPixelTopology::maxNumClustersPerModules;
    static constexpr uint32_t maxHitsInModule = phase1HIonPixelTopology::maxNumClustersPerModules;

    static constexpr char const *nameModifier = "HIonPhase1";
  };

  struct GenericUpgrade : public Phase2
  {
    static constexpr uint32_t numberOfModules = 20000; // this is just a maximum when running from hits
    static constexpr uint32_t numberOfLayers = 50;
    static constexpr uint32_t maxDepth = 30;
    static constexpr uint32_t maxHitsOnTrack = 30;
    static constexpr uint32_t maxHitsOnTrackForFullFit = 15;

    static constexpr char const *nameModifier = "Upgrade";

    static constexpr uint32_t const *layerStart = phase2PixelTopology::layerStart;
    static constexpr float const *minz = phase2PixelTopology::minz;
    static constexpr float const *maxz = phase2PixelTopology::maxz;
    static constexpr float const *maxr = phase2PixelTopology::maxr;

    static constexpr uint8_t const *layerPairs = phase2PixelTopology::layerPairs;
    static constexpr int16_t const *phicuts = phase2PixelTopology::phicuts;
    static constexpr float const *thetaCuts = phase2PixelTopology::thetaCuts;
    static constexpr float const *dcaCuts = phase2PixelTopology::dcaCuts;
  };

  struct ColliderMLPhase1 {
    // types
    using hindex_type = uint32_t;  // FIXME from siPixelRecHitsHeterogeneousProduct
    using tindex_type = uint32_t;  // for tuples
    using cindex_type = uint32_t;  // for cells

    static constexpr uint32_t maxCellNeighbors = 256;
    static constexpr uint32_t maxCellTracks = 302;
    static constexpr uint32_t maxHitsOnTrack = 20;
    static constexpr uint32_t maxHitsOnTrackForFullFit = 6;
    static constexpr uint32_t avgHitsPerTrack = 7;
    static constexpr uint32_t maxCellsPerHit = 256;
    static constexpr uint32_t avgTracksPerHit = 10;
    static constexpr uint32_t maxNumberOfTuples = 256 * 1024;
    // this is well above thanks to maxNumberOfTuples
    static constexpr uint32_t maxHitsForContainers = avgHitsPerTrack * maxNumberOfTuples;
    static constexpr uint32_t maxNumberOfDoublets = 5 * 512 * 1024;
    static constexpr uint32_t maxNumOfActiveDoublets = maxNumberOfDoublets / 2;
    static constexpr uint32_t maxNumberOfQuadruplets = maxNumberOfTuples;
    static constexpr uint32_t maxDepth = 6;
    static constexpr uint32_t numberOfLayers = 10;

    static constexpr uint32_t maxSizeCluster = 2047;

    static constexpr uint32_t getDoubletsFromHistoMaxBlockSize = 128;  // for both x and y
    static constexpr uint32_t getDoubletsFromHistoMinBlocksPerMP = 16;

    static constexpr uint16_t last_bpix1_detIndex = 216;
    static constexpr uint16_t last_bpix2_detIndex = 432;
    static constexpr uint16_t last_barrel_detIndex = 864;

    static constexpr uint32_t maxPixInModule = 6000;
    static constexpr uint32_t maxPixInModuleForMorphing = 0;
    static constexpr uint32_t maxIterClustering = 16;

    static constexpr uint32_t maxNumClustersPerModules = colliderMLPhase1PixelTopology::maxNumClustersPerModules;
    static constexpr uint32_t maxHitsInModule = 20000;

    static constexpr float moduleLength = 4.345f;
    static constexpr float endcapCorrection = 0.0f;

    static constexpr float xerr_barrel_l1_def = 0.00035f;
    static constexpr float yerr_barrel_l1_def = 0.00125f;
    static constexpr float xerr_barrel_ln_def = 0.00035f;
    static constexpr float yerr_barrel_ln_def = 0.00125f;
    static constexpr float xerr_endcap_def = 0.00060f;
    static constexpr float yerr_endcap_def = 0.00180f;

    static constexpr float bigPixXCorrection = 0.0f;
    static constexpr float bigPixYCorrection = 0.0f;

    static constexpr float dzdrFact = 8 * 0.0285 / 0.015;  // from dz/dr to "DY"

    static constexpr int nPairsMinimal = 33;
    static constexpr int nPairsFarForwards = nPairsMinimal + 8;  // include barrel "jumping" layer pairs
    static constexpr int nPairs = colliderMLPhase1PixelTopology::nPairs;   // include far forward layer pairs

    static constexpr int maxDYsize12 = 12;
    static constexpr int maxDYsize = 10;
    static constexpr int maxDYPred = 20;

    static constexpr uint16_t numberOfModules = colliderMLPhase1PixelTopology::numberOfModules;

    // 1000 bins < 1024 bins (10 bits) must be:
    // - < 32*32 (warpSize*warpSize for block prefix scan for CUDA)
    // - > number of columns (y) in any module. This is due to the fact
    //     that in pixel clustering we give for granted that in each
    //     bin we only have the pixel belonging to the same column.
    //     See RecoLocalTracker/SiPixelClusterizer/plugins/alpaka/PixelClustering.h#L325-L347
    static constexpr uint16_t clusterBinning = 1000;
    static constexpr uint16_t clusterBits = 10;

    static constexpr uint16_t numberOfModulesInBarrel = 756;
    static constexpr uint16_t numberOfModulesInLadder = 9;
    static constexpr uint16_t numberOfLaddersInBarrel = numberOfModulesInBarrel / numberOfModulesInLadder;

    static constexpr uint16_t firstPixelBarrelPos = 0;
    static constexpr uint16_t lastPixelBarrelPos = 3;
    static constexpr uint16_t firstPixelEndcapPos = 4;
    static constexpr uint16_t lastPixelEndcapPos = 6;
    static constexpr uint16_t firstPixelEndcapNeg = 7;
    static constexpr uint16_t lastPixelEndcapNeg = 9;

    static constexpr uint16_t firstShortStripsBarrelPos = 0;
    static constexpr uint16_t lastShortStripsBarrelPos = 0;
    static constexpr uint16_t firstShortStripsEndcapPos = 0;
    static constexpr uint16_t lastShortStripsEndcapPos = 0;
    static constexpr uint16_t firstShortStripsEndcapNeg = 0;
    static constexpr uint16_t lastShortStripsEndcapNeg = 0;

    static constexpr uint16_t firstLongStripsBarrelPos = 0;
    static constexpr uint16_t lastLongStripsBarrelPos = 0;
    static constexpr uint16_t firstLongStripsEndcapPos = 0;
    static constexpr uint16_t lastLongStripsEndcapPos = 0;
    static constexpr uint16_t firstLongStripsEndcapNeg = 0;
    static constexpr uint16_t lastLongStripsEndcapNeg = 0;

    static constexpr int16_t xOffset = -1e4;  // not used actually, to suppress static analyzer warnings

    static constexpr char const *nameModifier = "ColliderMLPhase1";
    static constexpr char const *cpeModules = "PixelCPEFastParams";

    static constexpr uint32_t const *layerStart = colliderMLPhase1PixelTopology::layerStart;
    static constexpr float const *minz = colliderMLPhase1PixelTopology::minz;
    static constexpr float const *maxz = colliderMLPhase1PixelTopology::maxz;
    static constexpr float const *maxr = colliderMLPhase1PixelTopology::maxr;

    static constexpr uint8_t const *layerPairs = colliderMLPhase1PixelTopology::layerPairs;
    static constexpr int16_t const *phicuts = colliderMLPhase1PixelTopology::phicuts;
    static constexpr float const *thetaCuts = colliderMLPhase1PixelTopology::thetaCuts;
    static constexpr float const *dcaCuts = colliderMLPhase1PixelTopology::dcaCuts;

    static constexpr uint8_t const *startingPairs = colliderMLPhase1PixelTopology::startingPairs;
    static constexpr bool const *isBarrel = colliderMLPhase1PixelTopology::isBarrel;
    static constexpr float const *ptCuts = colliderMLPhase1PixelTopology::ptCuts;

    static constexpr float const cellZ0Cut = 12.0;
    static constexpr float const cellPtCut = 0.5;
    static constexpr float const hardCurvCut = 0.015;

    static constexpr inline bool isBigPixX(uint16_t px) { return false; }
    static constexpr inline bool isBigPixY(uint16_t py) { return false; }

    static constexpr inline uint16_t localX(uint16_t px) { return px; }
    static constexpr inline uint16_t localY(uint16_t py) { return py; }
  };

  struct ColliderMLPhase2 {
    // types
    using hindex_type = uint32_t;  // FIXME from siPixelRecHitsHeterogeneousProduct
    using tindex_type = uint32_t;  // for tuples
    using cindex_type = uint32_t;  // for cells

    static constexpr uint32_t maxCellNeighbors = 256;
    static constexpr uint32_t maxCellTracks = 302;
    static constexpr uint32_t maxHitsOnTrack = 20;
    static constexpr uint32_t maxHitsOnTrackForFullFit = 6;
    static constexpr uint32_t avgHitsPerTrack = 7;
    static constexpr uint32_t maxCellsPerHit = 256;
    static constexpr uint32_t avgTracksPerHit = 10;
    static constexpr uint32_t maxNumberOfTuples = 256 * 1024;
    // this is well above thanks to maxNumberOfTuples
    static constexpr uint32_t maxHitsForContainers = avgHitsPerTrack * maxNumberOfTuples;
    static constexpr uint32_t maxNumberOfDoublets = 5 * 512 * 1024;
    static constexpr uint32_t maxNumOfActiveDoublets = maxNumberOfDoublets / 2;
    static constexpr uint32_t maxNumberOfQuadruplets = maxNumberOfTuples;
    static constexpr uint32_t maxDepth = 10;
    static constexpr uint32_t numberOfLayers = 18;

    static constexpr uint32_t maxSizeCluster = 2047;

    static constexpr uint32_t getDoubletsFromHistoMaxBlockSize = 128;  // for both x and y
    static constexpr uint32_t getDoubletsFromHistoMinBlocksPerMP = 16;

    static constexpr uint16_t last_bpix1_detIndex = 216;
    static constexpr uint16_t last_bpix2_detIndex = 432;
    static constexpr uint16_t last_barrel_detIndex = 864;

    static constexpr uint32_t maxPixInModule = 6000;
    static constexpr uint32_t maxPixInModuleForMorphing = 0;
    static constexpr uint32_t maxIterClustering = 16;

    static constexpr uint32_t maxNumClustersPerModules = colliderMLPhase2PixelTopology::maxNumClustersPerModules;
    static constexpr uint32_t maxHitsInModule = colliderMLPhase2PixelTopology::maxNumClustersPerModules;

    static constexpr float moduleLength = 4.345f;
    static constexpr float endcapCorrection = 0.0f;

    static constexpr float xerr_barrel_l1_def = 0.00035f;
    static constexpr float yerr_barrel_l1_def = 0.00125f;
    static constexpr float xerr_barrel_ln_def = 0.00035f;
    static constexpr float yerr_barrel_ln_def = 0.00125f;
    static constexpr float xerr_endcap_def = 0.00060f;
    static constexpr float yerr_endcap_def = 0.00180f;

    static constexpr float bigPixXCorrection = 0.0f;
    static constexpr float bigPixYCorrection = 0.0f;

    static constexpr float dzdrFact = 8 * 0.0285 / 0.015;  // from dz/dr to "DY"

    static constexpr int nPairsMinimal = 33;
    static constexpr int nPairsFarForwards = nPairsMinimal + 8;  // include barrel "jumping" layer pairs
    static constexpr int nPairs = colliderMLPhase2PixelTopology::nPairs;   // include far forward layer pairs

    static constexpr int maxDYsize12 = 12;
    static constexpr int maxDYsize = 10;
    static constexpr int maxDYPred = 20;

    static constexpr uint16_t numberOfModules = colliderMLPhase2PixelTopology::numberOfModules;

    // 1000 bins < 1024 bins (10 bits) must be:
    // - < 32*32 (warpSize*warpSize for block prefix scan for CUDA)
    // - > number of columns (y) in any module. This is due to the fact
    //     that in pixel clustering we give for granted that in each
    //     bin we only have the pixel belonging to the same column.
    //     See RecoLocalTracker/SiPixelClusterizer/plugins/alpaka/PixelClustering.h#L325-L347
    static constexpr uint16_t clusterBinning = 1000;
    static constexpr uint16_t clusterBits = 10;

    static constexpr uint16_t numberOfModulesInBarrel = 756;
    static constexpr uint16_t numberOfModulesInLadder = 9;
    static constexpr uint16_t numberOfLaddersInBarrel = numberOfModulesInBarrel / numberOfModulesInLadder;

    static constexpr uint16_t firstPixelBarrelPos = 0;
    static constexpr uint16_t lastPixelBarrelPos = 3;
    static constexpr uint16_t firstPixelEndcapPos = 4;
    static constexpr uint16_t lastPixelEndcapPos = 10;
    static constexpr uint16_t firstPixelEndcapNeg = 11;
    static constexpr uint16_t lastPixelEndcapNeg = 17;

    static constexpr uint16_t firstShortStripsBarrelPos = 0;
    static constexpr uint16_t lastShortStripsBarrelPos = 0;
    static constexpr uint16_t firstShortStripsEndcapPos = 0;
    static constexpr uint16_t lastShortStripsEndcapPos = 0;
    static constexpr uint16_t firstShortStripsEndcapNeg = 0;
    static constexpr uint16_t lastShortStripsEndcapNeg = 0;

    static constexpr uint16_t firstLongStripsBarrelPos = 0;
    static constexpr uint16_t lastLongStripsBarrelPos = 0;
    static constexpr uint16_t firstLongStripsEndcapPos = 0;
    static constexpr uint16_t lastLongStripsEndcapPos = 0;
    static constexpr uint16_t firstLongStripsEndcapNeg = 0;
    static constexpr uint16_t lastLongStripsEndcapNeg = 0;

    static constexpr int16_t xOffset = -1e4;  // not used actually, to suppress static analyzer warnings

    static constexpr char const *nameModifier = "ColliderMLPhase2";
    static constexpr char const *cpeModules = "PixelCPEFastParamsPhase2";

    static constexpr uint32_t const *layerStart = colliderMLPhase2PixelTopology::layerStart;
    static constexpr float const *minz = colliderMLPhase2PixelTopology::minz;
    static constexpr float const *maxz = colliderMLPhase2PixelTopology::maxz;
    static constexpr float const *maxr = colliderMLPhase2PixelTopology::maxr;

    static constexpr uint8_t const *layerPairs = colliderMLPhase2PixelTopology::layerPairs;
    static constexpr int16_t const *phicuts = colliderMLPhase2PixelTopology::phicuts;
    static constexpr float const *thetaCuts = colliderMLPhase2PixelTopology::thetaCuts;
    static constexpr float const *dcaCuts = colliderMLPhase2PixelTopology::dcaCuts;

    static constexpr uint8_t const *startingPairs = colliderMLPhase2PixelTopology::startingPairs;
    static constexpr bool const *isBarrel = colliderMLPhase2PixelTopology::isBarrel;
    static constexpr float const *ptCuts = colliderMLPhase2PixelTopology::ptCuts;

    static constexpr float const cellZ0Cut = 14.5;
    static constexpr float const cellPtCut = 0.85;
    static constexpr float const hardCurvCut = 0.0125;

    static constexpr inline bool isBigPixX(uint16_t px) { return false; }
    static constexpr inline bool isBigPixY(uint16_t py) { return false; }

    static constexpr inline uint16_t localX(uint16_t px) { return px; }
    static constexpr inline uint16_t localY(uint16_t py) { return py; }
  };

  struct ColliderMLPixelPlusShortStripsPhase2 {
    static constexpr uint32_t numberOfLayers = 34;
    static constexpr uint32_t maxDepth = 15;
    static constexpr uint32_t maxHitsOnTrack = 30;
    static constexpr uint32_t maxHitsOnTrackForFullFit = 15;

    static constexpr char const *nameModifier = "ColliderMLPixelPlusShortStripsPhase2";

    static constexpr uint32_t maxNumClustersPerModules = colliderMLPhase2PixelPlusShortStripsTopology::maxNumClustersPerModules;
    static constexpr uint32_t maxHitsInModule = colliderMLPhase2PixelPlusShortStripsTopology::maxNumClustersPerModules;

    static constexpr int nPairsMinimal = 95;
    static constexpr int nPairsFarForwards = nPairsMinimal + 16;  // include barrel "jumping" layer pairs
    static constexpr int nPairs = colliderMLPhase2PixelPlusShortStripsTopology::nPairs;   // include far forward layer pairs

    static constexpr uint16_t numberOfModules = colliderMLPhase2PixelPlusShortStripsTopology::numberOfModules;

    static constexpr uint32_t const *layerStart = colliderMLPhase2PixelPlusShortStripsTopology::layerStart;

    static constexpr float const *minz = colliderMLPhase2PixelPlusShortStripsTopology::minz;
    static constexpr float const *maxz = colliderMLPhase2PixelPlusShortStripsTopology::maxz;
    static constexpr float const *maxr = colliderMLPhase2PixelPlusShortStripsTopology::maxr;

    static constexpr uint8_t const *layerPairs = colliderMLPhase2PixelPlusShortStripsTopology::layerPairs;
    static constexpr int16_t const *phicuts = colliderMLPhase2PixelPlusShortStripsTopology::phicuts;
    static constexpr float const *thetaCuts = colliderMLPhase2PixelPlusShortStripsTopology::thetaCuts;
    static constexpr float const *dcaCuts = colliderMLPhase2PixelPlusShortStripsTopology::dcaCuts;

    static constexpr uint8_t const *startingPairs = colliderMLPhase2PixelPlusShortStripsTopology::startingPairs;
    static constexpr bool const *isBarrel = colliderMLPhase2PixelPlusShortStripsTopology::isBarrel;
    static constexpr float const *ptCuts = colliderMLPhase2PixelPlusShortStripsTopology::ptCuts;

    static constexpr uint16_t firstPixelBarrelPos = 0;
    static constexpr uint16_t lastPixelBarrelPos = 3;
    static constexpr uint16_t firstPixelEndcapPos = 4;
    static constexpr uint16_t lastPixelEndcapPos = 10;
    static constexpr uint16_t firstPixelEndcapNeg = 11;
    static constexpr uint16_t lastPixelEndcapNeg = 17;

    static constexpr uint16_t firstShortStripsBarrelPos = 18;
    static constexpr uint16_t lastShortStripsBarrelPos = 21;
    static constexpr uint16_t firstShortStripsEndcapPos = 22;
    static constexpr uint16_t lastShortStripsEndcapPos = 27;
    static constexpr uint16_t firstShortStripsEndcapNeg = 28;
    static constexpr uint16_t lastShortStripsEndcapNeg = 33;

    static constexpr uint16_t firstLongStripsBarrelPos = 0;
    static constexpr uint16_t lastLongStripsBarrelPos = 0;
    static constexpr uint16_t firstLongStripsEndcapPos = 0;
    static constexpr uint16_t lastLongStripsEndcapPos = 0;
    static constexpr uint16_t firstLongStripsEndcapNeg = 0;
    static constexpr uint16_t lastLongStripsEndcapNeg = 0;

    // types
    using hindex_type = uint32_t;  // FIXME from siPixelRecHitsHeterogeneousProduct
    using tindex_type = uint32_t;  // for tuples
    using cindex_type = uint32_t;  // for cells

    static constexpr uint32_t maxCellNeighbors = 256;
    static constexpr uint32_t maxCellTracks = 302;
    static constexpr uint32_t avgHitsPerTrack = 20;
    static constexpr uint32_t maxCellsPerHit = 256;
    static constexpr uint32_t avgTracksPerHit = 20;
    static constexpr uint32_t maxNumberOfTuples = 256 * 1024;
    // this is well above thanks to maxNumberOfTuples
    static constexpr uint32_t maxHitsForContainers = avgHitsPerTrack * maxNumberOfTuples;
    static constexpr uint32_t maxNumberOfDoublets = 5 * 512 * 1024;
    static constexpr uint32_t maxNumOfActiveDoublets = maxNumberOfDoublets / 2;
    static constexpr uint32_t maxNumberOfQuadruplets = maxNumberOfTuples;

    static constexpr uint32_t maxSizeCluster = 2047;

    static constexpr uint32_t getDoubletsFromHistoMaxBlockSize = 128;  // for both x and y
    static constexpr uint32_t getDoubletsFromHistoMinBlocksPerMP = 16;

    static constexpr uint16_t last_bpix1_detIndex = 216;
    static constexpr uint16_t last_bpix2_detIndex = 432;
    static constexpr uint16_t last_barrel_detIndex = 864;

    static constexpr uint32_t maxPixInModule = 6000;
    static constexpr uint32_t maxPixInModuleForMorphing = 0;
    static constexpr uint32_t maxIterClustering = 16;

    static constexpr float moduleLength = 4.345f;
    static constexpr float endcapCorrection = 0.0f;

    static constexpr float xerr_barrel_l1_def = 0.00035f;
    static constexpr float yerr_barrel_l1_def = 0.00125f;
    static constexpr float xerr_barrel_ln_def = 0.00035f;
    static constexpr float yerr_barrel_ln_def = 0.00125f;
    static constexpr float xerr_endcap_def = 0.00060f;
    static constexpr float yerr_endcap_def = 0.00180f;

    static constexpr float bigPixXCorrection = 0.0f;
    static constexpr float bigPixYCorrection = 0.0f;

    static constexpr float dzdrFact = 8 * 0.0285 / 0.015;  // from dz/dr to "DY"

    static constexpr int maxDYsize12 = 12;
    static constexpr int maxDYsize = 10;
    static constexpr int maxDYPred = 20;

    // 1000 bins < 1024 bins (10 bits) must be:
    // - < 32*32 (warpSize*warpSize for block prefix scan for CUDA)
    // - > number of columns (y) in any module. This is due to the fact
    //     that in pixel clustering we give for granted that in each
    //     bin we only have the pixel belonging to the same column.
    //     See RecoLocalTracker/SiPixelClusterizer/plugins/alpaka/PixelClustering.h#L325-L347
    static constexpr uint16_t clusterBinning = 1000;
    static constexpr uint16_t clusterBits = 10;

    static constexpr uint16_t numberOfModulesInBarrel = 756;
    static constexpr uint16_t numberOfModulesInLadder = 9;
    static constexpr uint16_t numberOfLaddersInBarrel = numberOfModulesInBarrel / numberOfModulesInLadder;

    static constexpr int16_t xOffset = -1e4;  // not used actually, to suppress static analyzer warnings

    static constexpr char const *cpeModules = "PixelCPEFastParamsPhase2";

    static constexpr float const cellZ0Cut = 14.5;
    static constexpr float const cellPtCut = 0.85;
    static constexpr float const hardCurvCut = 0.0115;

    static constexpr inline bool isBigPixX(uint16_t px) { return false; }
    static constexpr inline bool isBigPixY(uint16_t py) { return false; }

    static constexpr inline uint16_t localX(uint16_t px) { return px; }
    static constexpr inline uint16_t localY(uint16_t py) { return py; }
  };

  struct ColliderMLAllTrackerPhase2 {
    static constexpr uint32_t numberOfLayers = 48;
    static constexpr uint32_t maxDepth = 20;
    static constexpr uint32_t maxHitsOnTrack = 40;
    static constexpr uint32_t maxHitsOnTrackForFullFit = 20;

    static constexpr char const *nameModifier = "ColliderMLAllTrackerPhase2";

    static constexpr uint32_t maxNumClustersPerModules = colliderMLPhase2AllTrackerTopology::maxNumClustersPerModules;
    static constexpr uint32_t maxHitsInModule = colliderMLPhase2AllTrackerTopology::maxNumClustersPerModules;

    static constexpr int nPairsMinimal = colliderMLPhase2AllTrackerTopology::nPairs - 16;
    static constexpr int nPairsFarForwards = nPairsMinimal + 16;  // include barrel "jumping" layer pairs
    static constexpr int nPairs = colliderMLPhase2AllTrackerTopology::nPairs;   // include far forward layer pairs

    static constexpr uint16_t numberOfModules = colliderMLPhase2AllTrackerTopology::numberOfModules;

    static constexpr uint32_t const *layerStart = colliderMLPhase2AllTrackerTopology::layerStart;

    static constexpr float const *minz = colliderMLPhase2AllTrackerTopology::minz;
    static constexpr float const *maxz = colliderMLPhase2AllTrackerTopology::maxz;
    static constexpr float const *maxr = colliderMLPhase2AllTrackerTopology::maxr;

    static constexpr uint8_t const *layerPairs = colliderMLPhase2AllTrackerTopology::layerPairs;
    static constexpr int16_t const *phicuts = colliderMLPhase2AllTrackerTopology::phicuts;
    static constexpr float const *thetaCuts = colliderMLPhase2AllTrackerTopology::thetaCuts;
    static constexpr float const *dcaCuts = colliderMLPhase2AllTrackerTopology::dcaCuts;

    static constexpr uint8_t const *startingPairs = colliderMLPhase2AllTrackerTopology::startingPairs;
    static constexpr bool const *isBarrel = colliderMLPhase2AllTrackerTopology::isBarrel;
    static constexpr float const *ptCuts = colliderMLPhase2AllTrackerTopology::ptCuts;

    static constexpr uint16_t firstPixelBarrelPos = 0;
    static constexpr uint16_t lastPixelBarrelPos = 3;
    static constexpr uint16_t firstPixelEndcapPos = 4;
    static constexpr uint16_t lastPixelEndcapPos = 10;
    static constexpr uint16_t firstPixelEndcapNeg = 11;
    static constexpr uint16_t lastPixelEndcapNeg = 17;

    static constexpr uint16_t firstShortStripsBarrelPos = 18;
    static constexpr uint16_t lastShortStripsBarrelPos = 21;
    static constexpr uint16_t firstShortStripsEndcapPos = 22;
    static constexpr uint16_t lastShortStripsEndcapPos = 27;
    static constexpr uint16_t firstShortStripsEndcapNeg = 28;
    static constexpr uint16_t lastShortStripsEndcapNeg = 33;

    static constexpr uint16_t firstLongStripsBarrelPos = 34;
    static constexpr uint16_t lastLongStripsBarrelPos = 35;
    static constexpr uint16_t firstLongStripsEndcapPos = 36;
    static constexpr uint16_t lastLongStripsEndcapPos = 41;
    static constexpr uint16_t firstLongStripsEndcapNeg = 42;
    static constexpr uint16_t lastLongStripsEndcapNeg = 47;

    // types
    using hindex_type = uint32_t;  // FIXME from siPixelRecHitsHeterogeneousProduct
    using tindex_type = uint32_t;  // for tuples
    using cindex_type = uint32_t;  // for cells

    static constexpr uint32_t maxCellNeighbors = 256;
    static constexpr uint32_t maxCellTracks = 302;
    static constexpr uint32_t avgHitsPerTrack = 20;
    static constexpr uint32_t maxCellsPerHit = 256;
    static constexpr uint32_t avgTracksPerHit = 20;
    static constexpr uint32_t maxNumberOfTuples = 256 * 1024;
    // this is well above thanks to maxNumberOfTuples
    static constexpr uint32_t maxHitsForContainers = avgHitsPerTrack * maxNumberOfTuples;
    static constexpr uint32_t maxNumberOfDoublets = 5 * 512 * 1024;
    static constexpr uint32_t maxNumOfActiveDoublets = maxNumberOfDoublets / 2;
    static constexpr uint32_t maxNumberOfQuadruplets = maxNumberOfTuples;

    static constexpr uint32_t maxSizeCluster = 2047;

    static constexpr uint32_t getDoubletsFromHistoMaxBlockSize = 128;  // for both x and y
    static constexpr uint32_t getDoubletsFromHistoMinBlocksPerMP = 16;

    static constexpr uint16_t last_bpix1_detIndex = 216;
    static constexpr uint16_t last_bpix2_detIndex = 432;
    static constexpr uint16_t last_barrel_detIndex = 864;

    static constexpr uint32_t maxPixInModule = 6000;
    static constexpr uint32_t maxPixInModuleForMorphing = 0;
    static constexpr uint32_t maxIterClustering = 16;

    static constexpr float moduleLength = 4.345f;
    static constexpr float endcapCorrection = 0.0f;

    static constexpr float xerr_barrel_l1_def = 0.00035f;
    static constexpr float yerr_barrel_l1_def = 0.00125f;
    static constexpr float xerr_barrel_ln_def = 0.00035f;
    static constexpr float yerr_barrel_ln_def = 0.00125f;
    static constexpr float xerr_endcap_def = 0.00060f;
    static constexpr float yerr_endcap_def = 0.00180f;

    static constexpr float bigPixXCorrection = 0.0f;
    static constexpr float bigPixYCorrection = 0.0f;

    static constexpr float dzdrFact = 8 * 0.0285 / 0.015;  // from dz/dr to "DY"

    static constexpr int maxDYsize12 = 12;
    static constexpr int maxDYsize = 10;
    static constexpr int maxDYPred = 20;

    // 1000 bins < 1024 bins (10 bits) must be:
    // - < 32*32 (warpSize*warpSize for block prefix scan for CUDA)
    // - > number of columns (y) in any module. This is due to the fact
    //     that in pixel clustering we give for granted that in each
    //     bin we only have the pixel belonging to the same column.
    //     See RecoLocalTracker/SiPixelClusterizer/plugins/alpaka/PixelClustering.h#L325-L347
    static constexpr uint16_t clusterBinning = 1000;
    static constexpr uint16_t clusterBits = 10;

    static constexpr uint16_t numberOfModulesInBarrel = 756;
    static constexpr uint16_t numberOfModulesInLadder = 9;
    static constexpr uint16_t numberOfLaddersInBarrel = numberOfModulesInBarrel / numberOfModulesInLadder;

    static constexpr int16_t xOffset = -1e4;  // not used actually, to suppress static analyzer warnings

    static constexpr char const *cpeModules = "PixelCPEFastParamsPhase2";

    static constexpr float const cellZ0Cut = 14.5;
    static constexpr float const cellPtCut = 0.85;
    static constexpr float const hardCurvCut = 0.0915;

    static constexpr inline bool isBigPixX(uint16_t px) { return false; }
    static constexpr inline bool isBigPixY(uint16_t py) { return false; }

    static constexpr inline uint16_t localX(uint16_t px) { return px; }
    static constexpr inline uint16_t localY(uint16_t py) { return py; }
  };

  template <typename T>
  using isPhase1Topology = typename std::enable_if<std::is_base_of<Phase1, T>::value>::type;

  template <typename T>
  using isPhase2Topology = typename std::enable_if<std::is_base_of<Phase2, T>::value>::type;

  template <typename T>
  using isColliderMLPhase1Topology = typename std::enable_if<std::is_base_of<ColliderMLPhase1, T>::value>::type;

  template <typename T>
  using isColliderMLPhase2Topology = typename std::enable_if<std::is_base_of<ColliderMLPhase2, T>::value>::type;

  template <typename T>
  using isColliderMLPixelPlusShortStripsPhase2Topology = typename std::enable_if<std::is_base_of<ColliderMLPixelPlusShortStripsPhase2, T>::value>::type;

  template <typename T>
  using isColliderMLAllTrackerPhase2Topology = typename std::enable_if<std::is_base_of<ColliderMLAllTrackerPhase2, T>::value>::type;

}  // namespace pixelTopology

#endif  // Geometry_SimplePixelTopology_h
