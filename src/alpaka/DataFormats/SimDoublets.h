#ifndef SimDataFormats_TrackingAnalysis_SimDoublets_h
#define SimDataFormats_TrackingAnalysis_SimDoublets_h

// #include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
// #include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
// #include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitFwd.h"

#include "DataFormats/TrackingRecHitSimpleSoA.h"
#include "DataFormats/ParticleSimpleSoA.h"
#include "DataFormats/BeamSpotPOD.h"

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>
#include <algorithm>

/** @brief Semi-Monte Carlo truth information used for pixel-tracking opimization.
 *
 * SimDoublets hold references to all pixel RecHits of a simulated TrackingParticle.
 * Ones those RecHits are sorted according to their position relative to the particle vertex
 * by the method sortRecHits(), you can create the true doublets of RecHits that the 
 * TrackingParticle left in the detector. These SimDoublets::Doublet objects can be used  
 * to optimize the doublet creation in the reconstruction.
 *
 * The Doublets are generated as the RecHit pairs between two consecutively hit layers.
 * I.e., if a TrackingParticle produces
 *  - 1 hit (A) in 1st layer
 *  - 2 hits (B, C) in 3rd layer
 *  - 1 hit (D) in 4th layer
 * then, the true Doublets are:
 *  (A-B), (A-C), (B-D) and (C-D).
 * So, neither does it matter that the 2nd layer got "skipped" as there are no hits,
 * nor is the Doublet of (A-D) formed since there is a layer with hits in between.
 * Doublets are not created between hits within the same layer.
 *
 * @author Jan Schulz (jan.gerrit.schulz@cern.ch)
 * @date January 2025
 */
class SimDoublets {
public:
  /**
    * Sub-class for true doublets of RecHits
    *  - first hit = inner RecHit
    *  - second hit = outer RecHit
    */
  class Doublet {
  public:
    // default constructor
    Doublet() = default;

    // constructor
    Doublet(SimDoublets const&, size_t const, size_t const);

    // method to get the layer pair
    std::pair<int16_t, int16_t> layerIds() const { return layerIds_; }

    // method to get the RecHit pair
    std::pair<uint32_t, uint32_t> recHits() const { return recHitRefs_; }

    // method to get the number of skipped layers
    int8_t numSkippedLayers() const { return numSkippedLayers_; }

    // method to get the layer pair ID
    int16_t layerPairId() const { return layerPairId_; }

    // method to get the inner layerId
    int16_t innerLayerId() const { return layerIds_.first; }

    // method to get the outer layerId
    int16_t outerLayerId() const { return layerIds_.second; }

    // method to get a reference to the inner RecHit
    uint32_t innerRecHit() const { return recHitRefs_.first; }

    // method to get a reference to the outer RecHit
    uint32_t outerRecHit() const { return recHitRefs_.second; }

    // method to get the global position of the inner RecHit
    std::vector<double> innerGlobalPos(TrackingRecHitSimpleSoA) const;

    // method to get the global position of the outer RecHit
    std::vector<double> outerGlobalPos(TrackingRecHitSimpleSoA) const;

  private:
    uint32_t simParticleID_;                   // simulated particle ID
    std::pair<uint32_t, uint32_t> recHitRefs_;       // reference pair to RecHits of the Doublet
    std::pair<int16_t, int16_t> layerIds_;           // pair of layer IDs corresponding to the RecHits
    uint8_t numSkippedLayers_;                       // number of layers skipped by the Doublet
    uint16_t layerPairId_;                           // ID of the layer pair as defined in the reconstruction for the doublets
    std::vector<double> beamSpotPosition_{0.,0.,0.};  // global position of the beam spot (needed to correct the global RecHit position)
  };

  // default contructor
  SimDoublets() = default;

  // constructor
  SimDoublets(uint32_t const simParticleID, BeamSpotPOD const& beamSpot)
      : simParticleID_(simParticleID), beamSpotPosition_{beamSpot.x, beamSpot.y, beamSpot.z} {}

  // method to add a recHitId with its layer
  void addRecHit(uint32_t const recHitId, int16_t const layerId) {
    recHitsAreSorted_ = false;  // set sorted-bool to false again

    // check if the layerId is not present in the layerIdVector yet
    if (std::find(layerIdVector_.begin(), layerIdVector_.end(), layerId) == layerIdVector_.end()) {
      // if it does not exist, increment number of layers
      numLayers_++;
    }

    // add recHit and layerId to the vectors
    recHitIdVector_.push_back(recHitId);
    layerIdVector_.push_back(layerId);
  }

  // method to get the reference to the sim particle
  uint32_t simParticle() const { return simParticleID_; }

  // method to get the reference vector to the RecHits
  std::vector<uint32_t> recHits() const { return recHitIdVector_; }

  // method to get a reference to the RecHit at index i
  uint32_t recHits(size_t i) const { return recHitIdVector_[i]; }

  // method to get the layer id vector
  std::vector<int16_t> layerIds() const { return layerIdVector_; }

  // method to get the layer id at index i
  int16_t layerIds(size_t i) const { return layerIdVector_[i]; }

  // method to get the beam spot position
  std::vector<double> beamSpotPosition() const { return beamSpotPosition_; }

  // method to get the number of layers
  int numLayers() const { return numLayers_; }

  // method to get number of RecHits in the SimDoublets
  int numRecHits() const { return layerIdVector_.size(); }

  // method to sort the RecHits according to the position
  void sortRecHits(const ParticleSimpleSoA& particles, const TrackingRecHitSimpleSoA& hits);

  // method to produce the SimDoublets from the RecHits
  std::vector<Doublet> getSimDoublets() const;

private:
  uint32_t simParticleID_;                   // simulated particle ID
  std::vector<uint32_t> recHitIdVector_;           // reference vector to rec hits ID associated to the simulated particle (sorted afer building)
  std::vector<int16_t> layerIdVector_;             // vector of layer IDs corresponding to the RecHits
  bool recHitsAreSorted_{false};                   // true if RecHits were sorted
  std::vector<double> beamSpotPosition_{0.,0.,0.};  // beam spot position (defined as 0,0,0 by default)
  int numLayers_{0};                               // number of layers hit by the simulated particle
};

// collection of SimDoublets
typedef std::vector<SimDoublets> SimDoubletsCollection;

#endif
