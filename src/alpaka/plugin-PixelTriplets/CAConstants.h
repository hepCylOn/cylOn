#ifndef plugin_PixelTriplets_CAConstants_h
#define plugin_PixelTriplets_CAConstants_h

#include <cstdint>

#include "AlpakaCore/HistoContainer.h"
#include "AlpakaCore/SimpleVector.h"
#include "AlpakaCore/VecArray.h"
#include "AlpakaDataFormats/gpuClusteringConstants.h"

#define ONLY_PHICUT

namespace CAConstants {

  // constants
#ifndef ONLY_PHICUT
#ifdef GPU_SMALL_EVENTS
  constexpr uint32_t maxNumberOfTuples() { return 3 * 1024; }
#else
  constexpr uint32_t maxNumberOfTuples() { return 96 * 1024; }
#endif
#else
  constexpr uint32_t maxNumberOfTuples() { return 256 * 1024; }
#endif
  constexpr uint32_t maxNumberOfQuadruplets() { return maxNumberOfTuples(); }
#ifndef ONLY_PHICUT
#ifndef GPU_SMALL_EVENTS
  constexpr uint32_t maxNumberOfDoublets() { return 512 * 1024; }
  constexpr uint32_t maxCellsPerHit() { return 128; }
#else
  constexpr uint32_t maxNumberOfDoublets() { return 256 * 1024; }
  constexpr uint32_t maxCellsPerHit() { return 128 / 2; }
#endif
#else
  constexpr uint32_t maxNumberOfDoublets() { return 2 * 1024 * 1024; }
  constexpr uint32_t maxCellsPerHit() { return 8 * 128; }
#endif
  constexpr uint32_t maxNumOfActiveDoublets() { return maxNumberOfDoublets() / 8; }

  constexpr uint32_t maxNumberOfLayerPairs() { return 50; }
  constexpr uint32_t maxNumberOfLayers() { return 10; }
  constexpr uint32_t maxTuples() { return maxNumberOfTuples(); }

  // types
  using hindex_type = uint32_t;  // FIXME from siPixelRecHitsHeterogeneousProduct
  using tindex_type = uint32_t;  //  for tuples

#ifndef ONLY_PHICUT
  using CellNeighbors = cms::alpakatools::VecArray<uint32_t, 36>;
  using CellTracks = cms::alpakatools::VecArray<tindex_type, 48>;
#else
  using CellNeighbors = cms::alpakatools::VecArray<uint32_t, 128>;
  using CellTracks = cms::alpakatools::VecArray<tindex_type, 128>;
#endif

  using PhiHist =
      cms::alpakatools::HistoContainer<int16_t, 128, gpuClustering::MaxNumClusters, 8 * sizeof(int16_t), uint32_t, 50>; //TODO make this templated

  using CellNeighborsVector = cms::alpakatools::SimpleVector<CellNeighbors>;
  using CellTracksVector = cms::alpakatools::SimpleVector<CellTracks>;

  using OuterHitOfCell = cms::alpakatools::VecArray<uint32_t, maxCellsPerHit()>;

  using TuplesContainer = cms::alpakatools::OneToManyAssoc<hindex_type, maxTuples(), 5 * maxTuples()>;
  using HitToTuple = cms::alpakatools::
      OneToManyAssoc<tindex_type, pixelGPUConstants::maxNumberOfHits, 4 * maxTuples()>;  // 3.5 should be enough
  using TupleMultiplicity = cms::alpakatools::OneToManyAssoc<tindex_type, 8, maxTuples()>;

}  // namespace CAConstants

#endif  // plugin_PixelTriplets_CAConstants_h
