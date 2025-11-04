#include "DataFormats/SimDoublets.h"

#include "DataFormats/TrackingRecHitSimpleSoA.h"
#include "DataFormats/ParticleSimpleSoA.h"

#include <numeric>

namespace simdoublets {

  // Function that gets the global position of a RecHit with respect to a reference point.
  std::vector<double> getGlobalHitPosition(uint32_t const& recHitID, TrackingRecHitSimpleSoA recHits, std::vector<double> const& referencePosition) {
    return {recHits.xg(recHitID) - referencePosition[0],recHits.yg(recHitID) - referencePosition[1],recHits.zg(recHitID) - referencePosition[2]};
  }

  // Function that determines the number of skipped layers for a given pair of RecHits.
  int getNumSkippedLayers(std::pair<int16_t, int16_t> const& layerIds) {
    // Possibility 0: invalid case (outer layer is not the outer one), set to -1 immediately
    if (layerIds.first >= layerIds.second) {
      return -1;
    }

    // determine where the RecHits are
    bool innerInBarrel = (layerIds.first < 4);
    bool outerInBarrel = (layerIds.second < 4);
    bool innerInBackward = (!innerInBarrel) && (layerIds.first >= 11);
    bool outerInBackward = (!outerInBarrel) && (layerIds.second >= 11);
    bool innerInForward = (!innerInBarrel) && (!innerInBackward);
    bool outerInForward = (!outerInBarrel) && (!outerInBackward);

    // Possibility 1: both RecHits lie in the same detector part (barrel, forward or backward)
    if ((innerInBarrel && outerInBarrel) || (innerInForward && outerInForward) ||
        (innerInBackward && outerInBackward)) {
      return (layerIds.second - layerIds.first - 1);
    }
    // Possibility 2: the inner RecHit is in the barrel while the outer is in either forward or backward
    else if (innerInBarrel) {
      if (outerInBackward) return ((layerIds.second - 10) - 1);
      if (outerInForward) return ((layerIds.second - 3) - 1);
      return -1;
    }
    // Possibility 3: invalid case (one is forward and the other in backward), set to -1
    else {
      return -1;
    }
  }

  // Function that, for a pair of two layers, gives a unique pair Id (innerLayerId * 100 + outerLayerId).
  int getLayerPairId(std::pair<int16_t, int16_t> const& layerIds) {
    // calculate the unique layer pair Id as (innerLayerId * 100 + outerLayerId)
    return (layerIds.first * 100 + layerIds.second);
  }
}  // end namespace simdoublets

// ------------------------------------------------------------------------------------------------------
// SimDoublets::Doublet class member functions
// ------------------------------------------------------------------------------------------------------

// constructor
SimDoublets::Doublet::Doublet(SimDoublets const& simDoublets,
                              size_t const innerIndex,
                              size_t const outerIndex)
    : simParticleID_(simDoublets.simParticle()), beamSpotPosition_(simDoublets.beamSpotPosition()) {
  // fill recHits and layers
  recHitRefs_ = std::make_pair(simDoublets.recHits(innerIndex), simDoublets.recHits(outerIndex));
  layerIds_ = std::make_pair(simDoublets.layerIds(innerIndex), simDoublets.layerIds(outerIndex));

  // determine number of skipped layers
  numSkippedLayers_ = simdoublets::getNumSkippedLayers(layerIds_);

  // determine Id of the layer pair
  layerPairId_ = simdoublets::getLayerPairId(layerIds_);
}

std::vector<double> SimDoublets::Doublet::innerGlobalPos(TrackingRecHitSimpleSoA recHits) const {
  // get the inner RecHit's global position
  return simdoublets::getGlobalHitPosition(recHitRefs_.first, recHits, beamSpotPosition_);
}

std::vector<double> SimDoublets::Doublet::outerGlobalPos(TrackingRecHitSimpleSoA recHits) const {
  // get the outer RecHit's global position
  return simdoublets::getGlobalHitPosition(recHitRefs_.second, recHits, beamSpotPosition_);
}

// ------------------------------------------------------------------------------------------------------
// SimDoublets class member functions
// ------------------------------------------------------------------------------------------------------

// method to sort the RecHits according to the position; needs particle position
void SimDoublets::sortRecHits(const ParticleSimpleSoA& particles, const TrackingRecHitSimpleSoA& recHits) {
  // get the production vertex of the TrackingParticle
  const std::vector<double> vertex{particles.vx(simParticleID_), particles.vy(simParticleID_), particles.vz(simParticleID_)};

  // get the vector of squared magnitudes of the global RecHit positions
  std::vector<double> recHitMag2;
  recHitMag2.reserve(layerIdVector_.size());
  for (const auto& recHit : recHitIdVector_) {
    // global RecHit position with respect to the production vertex
    std::vector<double> globalPosition = simdoublets::getGlobalHitPosition(recHit, recHits, vertex);
    recHitMag2.push_back((globalPosition[0]*globalPosition[0]) + (globalPosition[1]*globalPosition[1]) + (globalPosition[2]*globalPosition[2]));
  }

  // find the permutation vector that sort the magnitudes
  std::vector<std::size_t> sortedPerm(recHitMag2.size());
  std::iota(sortedPerm.begin(), sortedPerm.end(), 0);
  std::sort(sortedPerm.begin(), sortedPerm.end(), [&](std::size_t i, std::size_t j) {
    return (recHitMag2[i] < recHitMag2[j]);
  });

  // create the sorted recHitRefVector and the sorted layerIdVector accordingly
  std::vector<uint32_t> sorted_recHitIdVector;
  std::vector<int16_t> sorted_layerIdVector;
  sorted_recHitIdVector.reserve(sortedPerm.size());
  sorted_layerIdVector.reserve(sortedPerm.size());
  for (size_t i : sortedPerm) {
    sorted_recHitIdVector.push_back(recHitIdVector_[i]);
    sorted_layerIdVector.push_back(layerIdVector_[i]);
  }

  // swap them with the class member
  recHitIdVector_.swap(sorted_recHitIdVector);
  layerIdVector_.swap(sorted_layerIdVector);

  // set sorted bool to true
  recHitsAreSorted_ = true;
}

// method to produce the true doublets on the fly
std::vector<SimDoublets::Doublet> SimDoublets::getSimDoublets() const {
  // create output vector for the doublets
  std::vector<SimDoublets::Doublet> doubletVector;

  // confirm that the RecHits are sorted
  assert(recHitsAreSorted_);

  // check if there are at least two hits
  if (numRecHits() < 2) {
    return doubletVector;
  }

  // loop over the RecHits/layer Ids
  for (size_t i = 0; i < layerIdVector_.size(); i++) {
    int16_t innerLayerId = layerIdVector_[i];
    int16_t outerLayerId{};
    size_t outerLayerStart{layerIdVector_.size()};

    // find the next layer Id + at which hit this layer starts
    for (size_t j = i + 1; j < layerIdVector_.size(); j++) {
      if (innerLayerId != layerIdVector_[j]) {
        outerLayerId = layerIdVector_[j];
        outerLayerStart = j;
        break;
      }
    }

    // std::cout << innerLayerId << " -- " << outerLayerId << std::endl;

    // build the doublets of the inner hit i with all outer hits in the layer outerLayerId
    for (size_t j = outerLayerStart; j < layerIdVector_.size(); j++) {
      // break if the hit doesn't belong to the outer layer anymore
      if (outerLayerId != layerIdVector_[j]) {
        break;
      }

      doubletVector.push_back(SimDoublets::Doublet(*this, i, j));
    }
  }  // end loop over the RecHits/layer Ids

  return doubletVector;
}
