// -*- C++ -*-
//
// Package:    SimTracker/TrackerHitAssociation
// Class:      SimPixelTrackProducer
//

// user include files
#include "Framework/EDProducer.h"
#include "Framework/Event.h"
#include "Framework/EventSetup.h"
#include "Framework/PluginFactory.h"

#include "DataFormats/TrackingRecHitSimpleSoA.h"
#include "DataFormats/ParticleSimpleSoA.h"
#include "DataFormats/BeamSpotPOD.h"
#include "DataFormats/SimPixelTrack.h"

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>
#include <memory>
#include <typeinfo>
#include <algorithm>

/** Class: SimPixelTrackProducer
 * 
 * @brief Produces SimPixelTracks (MC-info based PixelRecHit doublets) for selected SimParticles.
 *
 * SimDoublets represent the true doublets of RecHits that a simulated particle (TrackingParticle) 
 * created in the pixel detector. They can be used to analyze cuts which are applied in the reconstruction
 * when producing doublets as the first part of patatrack pixel tracking.
 *
 * The SimPixelTrack are produced in the following way:
 * 1. We select reasonable SimParticles according to the criteria given in the config file as 
 *    "TrackingParticleSelectionConfig".
 * 2. For each selected particle, we create and append a new SimPixelTrack object to the SimPixelTrackCollection.
 * 3. We loop over all RecHits in the pixel tracker and check if the given RecHit is associated to one of
 *    the selected particles (association via TP to cluster association). If it is, we add a RecHit reference
 *    to the respective SimDoublet.
 * 4. In the end, we sort the RecHits in each SimPixelTrack object according to their global position.
 *
 * @author Jan Schulz (jan.gerrit.schulz@cern.ch)
 * @date January 2025
 */
class SimPixelTrackProducer : public edm::EDProducer {
public:
  explicit SimPixelTrackProducer(edm::ProductRegistry& reg);

  void produce(edm::Event&, const edm::EventSetup&) override;

private:
  bool dropEvenLayerRecHits_;  // if true, no RecHits from even layers are considered
  bool dropOddLayerRecHits_;   // if true, no RecHits from odd layers are considered

  BeamSpotPOD bs_;
  edm::EDGetTokenT<TrackingRecHitSimpleSoA> tSimpleHits_;
  edm::EDGetTokenT<ParticleSimpleSoA> tSimplePart_;

  const edm::EDPutTokenT<SimPixelTrackCollection> simPixelTracks_putToken_;

};

// constructor
SimPixelTrackProducer::SimPixelTrackProducer(edm::ProductRegistry& reg)
      : bs_{BeamSpotPOD()},
        tSimpleHits_(reg.consumes<TrackingRecHitSimpleSoA>()),
        tSimplePart_(reg.consumes<ParticleSimpleSoA>()),
        simPixelTracks_putToken_(reg.produces<SimPixelTrackCollection>()) {}

void SimPixelTrackProducer::produce(edm::Event& event, const edm::EventSetup& eventSetup) {

  // get information from the event
  bs_ = eventSetup.get<BeamSpotPOD>();
  auto const& hits = event.get(tSimpleHits_);
  auto const& particles = event.get(tSimplePart_);

  // create collection of SimPixelTrack
  // each element will correspond to one selected TrackingParticle
  SimPixelTrackCollection simPixelTrackCollection;

  // loop over SimParticles
  for (size_t i = 0; i < particles.nParticles(); ++i) {

    // select reasonable SimParticles for the study (e.g., only signal)
    simPixelTrackCollection.push_back(SimPixelTrack(i, bs_));
  }

//   // create a set of the keys of the selected SimParticles
//   edm::IndexSet selectedTrackingParticleKeys;
//   selectedTrackingParticleKeys.reserve(simPixelTrackCollection.size());
//   for (const auto& simPixelTrack : simPixelTrackCollection) {
//     TrackingParticleRef trackingParticleRef = simPixelTrack.trackingParticle();
//     selectedTrackingParticleKeys.insert(trackingParticleRef.key());
//   }

  // initialize a couple of counters
  int count_associatedRecHits{0}, count_RecHitsInSimPixelTrack{0};

  // initialize a couple of variables used in the following loop
  unsigned int layerId;

  dropEvenLayerRecHits_ = false;
  dropOddLayerRecHits_ = false;

  // loop over pixel RecHit collections of the different pixel modules
  for (uint32_t hitId = 0; hitId < hits.nHits(); ++hitId) {

    // determine layer Id from detector Id
    layerId = hits.detInd(hitId);

    // check if we would like to skip
    if (dropEvenLayerRecHits_ && (layerId % 2 == 0)) {
      continue;
    }
    if (dropOddLayerRecHits_ && (layerId % 2 == 1)) {
      continue;
    }

    auto it = std::find(particles.partIndVector().begin(), particles.partIndVector().end(), hits.partInd(hitId));
    uint32_t assocSimParticle = 0;
    if (it != particles.partIndVector().end()) {
      size_t index = std::distance(particles.partIndVector().begin(), it);
      assocSimParticle = index;
    //   std::cout << assocSimParticle << std::endl;
    }
    if (assocSimParticle != 0 && assocSimParticle < 2000000){
      count_associatedRecHits++;
      // loop over collection of SimDoublets and find the one of the associated TrackingParticle
      for (auto& simPixelTrack : simPixelTrackCollection) {
        uint32_t simParticle = simPixelTrack.simParticle();
        if (assocSimParticle == simParticle) {
          simPixelTrack.addRecHit(hitId, hits, layerId);
          count_RecHitsInSimPixelTrack++;
        }
      }
    }
  }  // end loop over pixel RecHit collections of the different pixel modules

  // loop over collection of SimPixelTrack and sort the RecHits according to their position
  for (auto& simPixelTrack : simPixelTrackCollection) {
    simPixelTrack.sortRecHits(simPixelTrack.simParticle(),particles);
  }

//   std::cout << "Size of SiPixelRecHitCollection : " << hits.nHits() << std::endl;
//   std::cout << count_associatedRecHits << " of " << hits.nHits()
//             << " RecHits are associated to selected SimParticles ("
//             << count_RecHitsInSimPixelTrack - count_associatedRecHits
//             << " of them were associated multiple times)." << std::endl;
//   std::cout << "Number of selected SimParticles : " << simPixelTrackCollection.size()
//                                     << std::endl;
//   std::cout << "Size of TrackingParticle Collection  : " << particles.nParticles()
//                                     << std::endl;

  // put the produced SimPixelTrack collection in the event
  event.emplace(simPixelTracks_putToken_, std::move(simPixelTrackCollection));
}

DEFINE_FWK_MODULE(SimPixelTrackProducer);