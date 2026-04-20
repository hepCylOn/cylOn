#include <cassert>
#include <chrono>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <ios>
#include <memory>
#include <mutex>

#include "Source.h"
#include "ParticleReader.h"
#include "HitReader.h"

#define INPUT_DEBUG

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

}  // namespace

namespace edm {
  Source::Source(
      int maxEvents, int runForMinutes, ProductRegistry &reg, std::filesystem::path const &datadir, bool validation, bool fromHits)
      : maxEvents_(maxEvents),
        runForMinutes_(runForMinutes),
        validation_(validation),
        fromHits_(fromHits) {

    // if(fromHits_ and validation_)
    //  throw std::runtime_error("--fromHits and --validation can't work together (yet)");
    
    std::ifstream in_file;

    if (not fromHits_)
    {
      in_file.open(datadir / "raw.bin", std::ios::binary);
      rawToken_ = reg.produces<FEDRawDataCollection>();
    }
    else
    {
      in_file.open(datadir / "hits.bin");
      // in_file.open("/data/user/borzari/cmssw/pixeltrack-standalone_withData/data/hits.txt");
      // TODO: remember to set this back to something more general
      // in_file.open(datadir / "hitsTest.txt", std::ios::binary);
      hitToken_ = reg.produces<reco::TrackingRecHitHost>();
    }
    std::ifstream in_digiclusters;
    std::ifstream in_tracks;
    std::ifstream in_vertices;
    std::ifstream in_particles;
    std::ifstream in_map;

    if (validation_) {
      digiClusterToken_ = reg.produces<DigiClusterCount>();
      trackToken_ = reg.produces<TrackCount>();
      vertexToken_ = reg.produces<VertexCount>();

      in_digiclusters = std::ifstream(datadir / "digicluster.bin", std::ios::binary);
      in_tracks = std::ifstream(datadir / "tracks.bin", std::ios::binary);
      in_vertices = std::ifstream(datadir / "vertices.bin", std::ios::binary);
      // in_particles    = std::ifstream(datadir / "particles.bin", std::ios::binary);
      // in_map = std::ifstream(datadir / "map.bin", std::ios::binary);

      in_digiclusters.exceptions(std::ifstream::badbit | std::ifstream::failbit | std::ifstream::eofbit);
      in_tracks.exceptions(std::ifstream::badbit | std::ifstream::failbit | std::ifstream::eofbit);
      in_vertices.exceptions(std::ifstream::badbit | std::ifstream::failbit | std::ifstream::eofbit);
      // in_particles.exceptions(std::ifstream::badbit | std::ifstream::failbit | std::ifstream::eofbit);
      
      // particleToken_ = reg.produces<sim::ParticleHost>();
      // mapToken_ = reg.produces<utils::SimpleMapHost>();
    }

    if(not fromHits_)
    {
      unsigned int nfeds;
      in_file.exceptions(std::ifstream::badbit);
      in_file.read(reinterpret_cast<char *>(&nfeds), sizeof(unsigned int));
      while (not in_file.eof()) {
        in_file.exceptions(std::ifstream::badbit | std::ifstream::failbit | std::ifstream::eofbit);

        raw_.emplace_back(readRaw(in_file, nfeds));

        if (validation_) {
          unsigned int nm, nd, nc, nt, nv;
          in_digiclusters.read(reinterpret_cast<char *>(&nm), sizeof(unsigned int));
          in_digiclusters.read(reinterpret_cast<char *>(&nd), sizeof(unsigned int));
          in_digiclusters.read(reinterpret_cast<char *>(&nc), sizeof(unsigned int));
          in_tracks.read(reinterpret_cast<char *>(&nt), sizeof(unsigned int));
          in_vertices.read(reinterpret_cast<char *>(&nv), sizeof(unsigned int));
          digiclusters_.emplace_back(nm, nd, nc);
          tracks_.emplace_back(nt);
          vertices_.emplace_back(nv);
        }

        // next event
        in_file.exceptions(std::ifstream::badbit);
        in_file.read(reinterpret_cast<char *>(&nfeds), sizeof(unsigned int));
      }
    }
    else
    {
      std::cout << "Reading hits from " << (datadir / "hits.bin") << std::endl;
      // std::cout << "Reading hits from /data/user/borzari/cmssw/pixeltrack-standalone_withData/data/hits.txt" << std::endl;
      hitReader::check_header(in_file);
      if (validation)
      {
        particleReader::check_header(in_particles);
        mapReader::check_header(in_map);
      }

      int32_t nEventsP, nEventsH, nEventsM;
      in_file.read(reinterpret_cast<char*>(&nEventsH), sizeof(nEventsH));
      nEventsP = nEventsM = nEventsH;

      if (validation)
      {
        in_particles.read(reinterpret_cast<char*>(&nEventsP), sizeof(nEventsP));
        in_map.read(reinterpret_cast<char*>(&nEventsM), sizeof(nEventsM));
      }

      // if(nEventsH != nEventsP or nEventsH != nEventsM)
      //   throw std::runtime_error("Error nEvents differs in the hits file and the particles file!");
      maxEvents_ = (maxEvents_ < 0) ? nEventsH : std::min(maxEvents_, nEventsH);
        
#ifdef INPUT_DEBUG
      std::cout << "File contains " << nEventsH << " events\n";
#endif
      for (int32_t ev = 0; ev < nEventsH && ev < maxEvents; ++ev) {
        hits_.emplace_back(hitReader::read_single_event(in_file));
        if (validation)
        {
          particles_.emplace_back(particleReader::read_single_event(in_particles));
          maps_.emplace_back(mapReader::read_single_event(in_map));
        }
          
#ifdef INPUT_DEBUG
        if (validation)
          std::cout << "Event " << ev << ": " << hits_[ev].nHits() << " hits, " << hits_[ev].nModules() << " modules - n. particles = " << particles_[ev].view().metadata().size() << std::endl;
        else
          std::cout << "Event " << ev << ": " << hits_[ev].nHits() << " hits, " << hits_[ev].nModules() << " modules\n";
#endif
      }
      if (!in_file.good() && !in_file.eof()) {
        throw std::runtime_error("I/O error while reading input file");
    }

      if (validation){
        if (!in_particles.good() && !in_particles.eof()) 
          throw std::runtime_error("I/O error while reading particles file");
        if (!in_map.good() && !in_map.eof()) 
          throw std::runtime_error("I/O error while reading hit-map file");
      }

    // std::cout << "Successfully read all events from " << filename << std::endl;
    }

    if (validation_ and not fromHits_) { //TODO allow for fromHits validation
      assert(raw_.size() == digiclusters_.size());
      assert(raw_.size() == tracks_.size());
      assert(raw_.size() == vertices_.size());
    } 
    else if (validation_)
    {
      assert(hits_.size() == particles_.size());
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
      ev->emplace(hitToken_, std::move(hits_[index]));
    if (validation_ and not fromHits_) {
      ev->emplace(digiClusterToken_, digiclusters_[index]);
      ev->emplace(trackToken_, tracks_[index]);
      ev->emplace(vertexToken_, vertices_[index]);
    }
    else if (validation_)
    {
      ev->emplace(particleToken_, std::move(particles_[index]));
      ev->emplace(mapToken_, std::move(maps_[index]));
    }

    return ev;
  }
}  // namespace edm
