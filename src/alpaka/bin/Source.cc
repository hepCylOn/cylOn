#include <cassert>
#include <chrono>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <ios>
#include <memory>
#include <mutex>

#include "Source.h"

namespace {
  FEDRawDataCollection readRaw(std::ifstream &is, unsigned int nfeds) {
    FEDRawDataCollection rawCollection;
    for (unsigned int ifed = 0; ifed < nfeds; ++ifed) {
      unsigned int fedId;
      is.read(reinterpret_cast<char *>(&fedId), sizeof(unsigned int));
      unsigned int fedSize;
      is.read(reinterpret_cast<char *>(&fedSize), sizeof(unsigned int));
      FEDRawData &rawData = rawCollection.FEDData(fedId);
      rawData.resize(fedSize);
      is.read(reinterpret_cast<char *>(rawData.data()), fedSize);
    }
    return rawCollection;
  }

//   TrackingRecHitSimpleSoA readHits(std::ifstream &if) {

//     TrackingRecHitSimpleSoA hitSoA(1);
//     hitSoA.readText(if);
// // #ifdef HITSBIN
// //     hitSoA.readBinary(is);
// // #else
    
// // #endif
//     return hitSoA;
//   }

}  // namespace

namespace edm {
  Source::Source(
      int maxEvents, int runForMinutes, ProductRegistry &reg, std::filesystem::path const &datadir, bool validation, bool fromHits, bool isPhase2)
      : maxEvents_(maxEvents),
        runForMinutes_(runForMinutes),
        validation_(validation),
        fromHits_(fromHits),
        isPhase2_(isPhase2) {
    
    std::ifstream in_raw;
    std::ifstream in_hits;
      
    in_raw.open(datadir / "raw.bin", std::ios::binary);
    rawToken_ = reg.produces<FEDRawDataCollection>();

    if (not isPhase2_) in_hits.open(datadir / "hitsCMSPhase1Test.txt");
    else in_hits.open(datadir / "hits.txt");
    // TODO: remember to set this back to something more general
    // in_raw.open(datadir / "hitsTest.txt", std::ios::binary);
    hitToken_ = reg.produces<TrackingRecHitSimpleSoA>();
    
    std::ifstream in_digiclusters;
    std::ifstream in_tracks;
    std::ifstream in_vertices;

    if (validation_) {
      if (not fromHits_) digiClusterToken_ = reg.produces<DigiClusterCount>();
      trackToken_ = reg.produces<TrackCount>();
      vertexToken_ = reg.produces<VertexCount>();

      if (not fromHits_) in_digiclusters = std::ifstream(datadir / "digicluster.bin", std::ios::binary);
      in_tracks = std::ifstream(datadir / "tracks.bin", std::ios::binary);
      in_vertices = std::ifstream(datadir / "vertices.bin", std::ios::binary);
      if (not fromHits_) in_digiclusters.exceptions(std::ifstream::badbit | std::ifstream::failbit | std::ifstream::eofbit);
      in_tracks.exceptions(std::ifstream::badbit | std::ifstream::failbit | std::ifstream::eofbit);
      in_vertices.exceptions(std::ifstream::badbit | std::ifstream::failbit | std::ifstream::eofbit);
    }

    unsigned int nfeds;
    in_raw.exceptions(std::ifstream::badbit);
    in_raw.read(reinterpret_cast<char *>(&nfeds), sizeof(unsigned int));
    while (not in_raw.eof()) {
      in_raw.exceptions(std::ifstream::badbit | std::ifstream::failbit | std::ifstream::eofbit);

      raw_.emplace_back(readRaw(in_raw, nfeds));

      if (validation_) {
        unsigned int nm, nd, nc, nt, nv;
        if (not fromHits_) in_digiclusters.read(reinterpret_cast<char *>(&nm), sizeof(unsigned int));
        if (not fromHits_) in_digiclusters.read(reinterpret_cast<char *>(&nd), sizeof(unsigned int));
        if (not fromHits_) in_digiclusters.read(reinterpret_cast<char *>(&nc), sizeof(unsigned int));
        in_tracks.read(reinterpret_cast<char *>(&nt), sizeof(unsigned int));
        in_vertices.read(reinterpret_cast<char *>(&nv), sizeof(unsigned int));
        if (not fromHits_) digiclusters_.emplace_back(nm, nd, nc);
        tracks_.emplace_back(nt);
        vertices_.emplace_back(nv);
      }

      // next event
      in_raw.exceptions(std::ifstream::badbit);
      in_raw.read(reinterpret_cast<char *>(&nfeds), sizeof(unsigned int));
    }

    if (fromHits_)
    {
      while (true) {
        TrackingRecHitSimpleSoA soa;
        if (!soa.readText(in_hits)) break;
        hits_.push_back(std::move(soa));
      }
    }

    if (validation_ and not fromHits_) { //TODO allow for fromHits validation
      if (not fromHits_) assert(raw_.size() == digiclusters_.size());
      assert(raw_.size() == tracks_.size());
      assert(raw_.size() == vertices_.size());
    }

    if (runForMinutes_ < 0 and maxEvents_ < 0) {
      if (not fromHits_) maxEvents_ = raw_.size();
      else maxEvents_ = hits_.size();
    }
  }

  void Source::reconfigure(int maxEvents, int runForMinutes) {
    std::scoped_lock lock(timeMutex_);
    maxEvents_ = maxEvents;
    runForMinutes_ = runForMinutes;
    numEventsTimeLastCheck_ = 0;
    shouldStop_ = false;
    numEvents_ = 0;
  }

  void Source::startProcessing() {
    if (runForMinutes_ >= 0) {
      startTime_ = std::chrono::steady_clock::now();
    }
  }

  std::unique_ptr<Event> Source::produce(int streamId, ProductRegistry const &reg) {
    if (shouldStop_) {
      return nullptr;
    }

    const int old = numEvents_.fetch_add(1);
    const int iev = old + 1;
    if (runForMinutes_ < 0) {
      if (old >= maxEvents_) {
        shouldStop_ = true;
        --numEvents_;
        return nullptr;
      }
    } else {
      if (numEvents_ - numEventsTimeLastCheck_ > static_cast<int>(raw_.size())) {
        std::scoped_lock lock(timeMutex_);
        // if some other thread beat us, no need to do anything
        if (numEvents_ - numEventsTimeLastCheck_ > static_cast<int>(raw_.size())) {
          auto processingTime = std::chrono::steady_clock::now() - startTime_;
          if (std::chrono::duration_cast<std::chrono::minutes>(processingTime).count() >= runForMinutes_) {
            shouldStop_ = true;
          }
          numEventsTimeLastCheck_ = (numEvents_ / raw_.size()) * raw_.size();
        }
        if (shouldStop_) {
          --numEvents_;
          return nullptr;
        }
      }
    }
    auto ev = std::make_unique<Event>(streamId, iev, reg);
    // This was const. Is it really needed? Can it stay as not const because of the distinct running modes?
    int index = 0;
    if (not fromHits_) index = old % raw_.size();
    else index = old % hits_.size();

    if (not fromHits_)
      ev->emplace(rawToken_, raw_[index]);
    else 
      ev->emplace(hitToken_, hits_[index]);
    if (validation_) {
      if (not fromHits_) ev->emplace(digiClusterToken_, digiclusters_[index]);
      ev->emplace(trackToken_, tracks_[index]);
      ev->emplace(vertexToken_, vertices_[index]);
    }

    return ev;
  }
}  // namespace edm
