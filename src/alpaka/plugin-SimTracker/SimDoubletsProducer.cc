// -*- C++ -*-
//
// Package:    SimTracker/TrackerHitAssociation
// Class:      SimDoubletsProducer
//

// user include files
#include "Framework/EDProducer.h"
#include "Framework/Event.h"
#include "Framework/EventSetup.h"
#include "Framework/PluginFactory.h"

#include "DataFormats/TrackingRecHitSimpleSoA.h"
#include "DataFormats/ParticleSimpleSoA.h"
#include "DataFormats/BeamSpotPOD.h"
#include "DataFormats/SimDoublets.h"

#include <cstddef>
#include <utility>
#include <vector>
#include <memory>
#include <typeinfo>

/** Class: SimDoubletsProducer
 * 
 * @brief Produces SimDoublets (MC-info based PixelRecHit doublets) for selected TrackingParticles.
 *
 * SimDoublets represent the true doublets of RecHits that a simulated particle (TrackingParticle) 
 * created in the pixel detector. They can be used to analyze cuts which are applied in the reconstruction
 * when producing doublets as the first part of patatrack pixel tracking.
 *
 * The SimDoublets are produced in the following way:
 * 1. We select reasonable TrackingParticles according to the criteria given in the config file as 
 *    "TrackingParticleSelectionConfig".
 * 2. For each selected particle, we create and append a new SimDoublets object to the SimDoubletsCollection.
 * 3. We loop over all RecHits in the pixel tracker and check if the given RecHit is associated to one of
 *    the selected particles (association via TP to cluster association). If it is, we add a RecHit reference
 *    to the respective SimDoublet.
 * 4. In the end, we sort the RecHits in each SimDoublets object according to their global position.
 *
 * @author Jan Schulz (jan.gerrit.schulz@cern.ch)
 * @date January 2025
 */

class SimDoubletsProducer : public edm::EDProducer {
public:
  explicit SimDoubletsProducer(edm::ProductRegistry& reg);

  void produce(edm::Event&, const edm::EventSetup&) override;

private:
  BeamSpotPOD bs_;
  edm::EDGetTokenT<TrackingRecHitSimpleSoA> tSimpleHits_;
  edm::EDGetTokenT<ParticleSimpleSoA> tSimplePart_;

  const edm::EDPutTokenT<SimDoubletsCollection> simDoublets_putToken_;

  bool verbose_{false};
};

// constructor

SimDoubletsProducer::SimDoubletsProducer(edm::ProductRegistry& reg)
      : bs_{BeamSpotPOD()},
        tSimpleHits_(reg.consumes<TrackingRecHitSimpleSoA>()),
        tSimplePart_(reg.consumes<ParticleSimpleSoA>()),
        simDoublets_putToken_(reg.produces<SimDoubletsCollection>()) {}

void SimDoubletsProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // get information from the event
  bs_ = iSetup.get<BeamSpotPOD>();
  auto const& hits = iEvent.get(tSimpleHits_);
  auto const& particles = iEvent.get(tSimplePart_);

  // create collection of SimDoublets
  // each element will correspond to one selected TrackingParticle
  SimDoubletsCollection simDoubletsCollection;

  // loop over TrackingParticles
  for (uint32_t i = 0; i < particles.nParticles(); ++i) {
    simDoubletsCollection.push_back(SimDoublets(i, bs_));
  }

  // initialize a couple of counters
  int count_associatedRecHits{0}, count_RecHitsInSimDoublets{0};

  // loop over pixel RecHit collections of the different pixel modules
  for (uint32_t hitId = 0; hitId < hits.nHits(); ++hitId) {
    // determine layer Id
    unsigned int layerId = hits.detInd(hitId);

    auto it = std::find(particles.partIndVector().begin(), particles.partIndVector().end(), hits.partInd(hitId));
    uint32_t assocSimParticle = 0;
    if (it != particles.partIndVector().end()) {
      size_t index = std::distance(particles.partIndVector().begin(), it);
      assocSimParticle = index;
    }
    if (assocSimParticle != 0 && assocSimParticle < 2000000){
      count_associatedRecHits++;
      // loop over collection of SimDoublets and find the one of the associated TrackingParticle
      for (auto& simDoublets : simDoubletsCollection) {
        uint32_t simParticle = simDoublets.simParticle();
        if (assocSimParticle == simParticle) {
          simDoublets.addRecHit(hitId, layerId);
          count_RecHitsInSimDoublets++;
        }
      }
    }
  }  // end loop over pixel RecHit collections of the different pixel modules

  // loop over collection of SimDoublets and sort the RecHits according to their position
  for (auto& simDoublets : simDoubletsCollection) {
    simDoublets.sortRecHits(particles, hits);
  }

  if (verbose_){
    std::cout << "=====================================" << std::endl;
    std::cout << "Size of SiPixelRecHitCollection : " << hits.nHits() << std::endl;
    std::cout << count_associatedRecHits << " of " << hits.nHits()
                                    << " RecHits are associated to selected sim particles ("
                                    << count_RecHitsInSimDoublets - count_associatedRecHits
                                    << " of them were associated multiple times)." << std::endl;
    std::cout << "Number of selected sim particles : " << simDoubletsCollection.size()
                                    << std::endl;
    std::cout << "Size of sim particles Collection  : " << particles.nParticles()
                                    << std::endl;
  }

  // put the produced SimDoublets collection in the event
  iEvent.emplace(simDoublets_putToken_, std::move(simDoubletsCollection));
}

DEFINE_FWK_MODULE(SimDoubletsProducer);
