#include <cassert>
#include <chrono>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <ios>
#include <memory>
#include <mutex>

#include "Source.h"

#define INPUT_DEBUG

///TODO: move this outside in an-hoc helper header
namespace hitReader {

  constexpr uint32_t kExpectedEndianness = 0x01020304;
  constexpr uint32_t kFormatVersion = 1;
  constexpr char kMagic[4] = {'T', 'R', 'H', '1'};

  inline void check_header(std::ifstream& in) {
    char magic[4];
    in.read(magic, 4);
    if (std::memcmp(magic, kMagic, 4) != 0)
      throw std::runtime_error("Invalid file magic (not a TRH1 binary)");

    uint32_t version;
    in.read(reinterpret_cast<char*>(&version), sizeof(version));
    if (version != kFormatVersion)
      throw std::runtime_error("Unsupported TRH version");

    uint32_t endian_marker;
    in.read(reinterpret_cast<char*>(&endian_marker), sizeof(endian_marker));
    if (endian_marker != kExpectedEndianness)
      throw std::runtime_error("Endianness mismatch — file not native endian");
#ifdef INPUT_DEBUG
    std::cout << "Input file is good!" << std::endl;
#endif
  }

  template <typename Span>
  void read_column(std::ifstream& in, Span view, uint32_t nHits) {
    using Elem = typename Span::value_type;
    if (view.size() < nHits) {
      throw std::runtime_error("Span too small for requested number of hits");
    }
    in.read(reinterpret_cast<char*>(view.data()), nHits * sizeof(Elem));
    if (!in)
      throw std::runtime_error("Error reading column data");
  }

  inline reco::TrackingRecHitHost read_single_event(std::ifstream& in)
  {
    uint32_t nHits, nModules;
    in.read(reinterpret_cast<char*>(&nHits), sizeof(nHits));
    in.read(reinterpret_cast<char*>(&nModules), sizeof(nModules));
    if (!in) throw std::runtime_error("Error reading event header");

    std::vector<uint32_t> moduleStart(nModules + 1);
    in.read(reinterpret_cast<char*>(moduleStart.data()), (nModules + 1) * sizeof(uint32_t));

    reco::TrackingRecHitHost recHitHost(cms::alpakatools::host(), nHits, nModules);

    auto hitView = recHitHost.view<reco::TrackingRecHitSoA>();
    auto modView = recHitHost.view<reco::HitModuleSoA>();

    // copy module starts
    std::memcpy(modView.moduleStart().data(), moduleStart.data(),
                (nModules + 1) * sizeof(uint32_t));

    read_column(in, hitView.xLocal(),       nHits);
    read_column(in, hitView.yLocal(),       nHits);
    read_column(in, hitView.xerrLocal(),    nHits);
    read_column(in, hitView.yerrLocal(),    nHits);
    read_column(in, hitView.xGlobal(),      nHits);
    read_column(in, hitView.yGlobal(),      nHits);
    read_column(in, hitView.zGlobal(),      nHits);
    read_column(in, hitView.rGlobal(),      nHits);
    read_column(in, hitView.iphi(),         nHits);
    read_column(in, hitView.chargeAndStatus(), nHits); 
    read_column(in, hitView.clusterSizeX(), nHits);
    read_column(in, hitView.clusterSizeY(), nHits);
    read_column(in, hitView.detectorIndex(), nHits);

#ifdef INPUT_DEBUG    
    std::cout << "  First hit global: ("
              << hitView.xGlobal()[0] << ", "
              << hitView.yGlobal()[0] << ", "
              << hitView.zGlobal()[0] << "), "
              << "r=" << hitView.rGlobal()[0]
              << ", detIdx=" << hitView.detectorIndex()[0] << '\n';
#endif

    return recHitHost;
  }


}
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
    
    
    if(fromHits_ and validation_)
     throw std::runtime_error("--fromHits and --validation can't work together (yet)");
    
    std::ifstream in_file;
      
    if (not fromHits_)
    {
      in_file.open(datadir / "raw.bin", std::ios::binary);
      rawToken_ = reg.produces<FEDRawDataCollection>();
    }
    else
    {
      in_file.open(datadir / "hits.bin");
      // TODO: remember to set this back to something more general
      // in_file.open(datadir / "hitsTest.txt", std::ios::binary);
      hitToken_ = reg.produces<reco::TrackingRecHitHost>();
    }
    std::ifstream in_digiclusters;
    std::ifstream in_tracks;
    std::ifstream in_vertices;

    if (validation_) {
      digiClusterToken_ = reg.produces<DigiClusterCount>();
      trackToken_ = reg.produces<TrackCount>();
      vertexToken_ = reg.produces<VertexCount>();

      in_digiclusters = std::ifstream(datadir / "digicluster.bin", std::ios::binary);
      in_tracks = std::ifstream(datadir / "tracks.bin", std::ios::binary);
      in_vertices = std::ifstream(datadir / "vertices.bin", std::ios::binary);
      in_digiclusters.exceptions(std::ifstream::badbit | std::ifstream::failbit | std::ifstream::eofbit);
      in_tracks.exceptions(std::ifstream::badbit | std::ifstream::failbit | std::ifstream::eofbit);
      in_vertices.exceptions(std::ifstream::badbit | std::ifstream::failbit | std::ifstream::eofbit);
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
      hitReader::check_header(in_file);
      int32_t nEvents;
      in_file.read(reinterpret_cast<char*>(&nEvents), sizeof(nEvents));
#ifdef INPUT_DEBUG
      std::cout << "File contains " << nEvents << " events\n";
#endif
      for (int32_t ev = 0; ev < nEvents && ev < maxEvents; ++ev) {
        hits_.emplace_back(hitReader::read_single_event(in_file));
#ifdef INPUT_DEBUG
        std::cout << "Event " << ev << ": " << hits_[ev].nHits() << " hits, " << hits_[ev].nModules() << " modules\n";
#endif
      }
      if (!in_file.good() && !in_file.eof()) {
        throw std::runtime_error("I/O error while reading file");
    }

    // std::cout << "Successfully read all events from " << filename << std::endl;
    }

    if (validation_ and not fromHits_) { //TODO allow for fromHits validation
      assert(raw_.size() == digiclusters_.size());
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
      ev->emplace(hitToken_, std::move(hits_[index]));
    if (validation_) {
      ev->emplace(digiClusterToken_, digiclusters_[index]);
      ev->emplace(trackToken_, tracks_[index]);
      ev->emplace(vertexToken_, vertices_[index]);
    }

    return ev;
  }
}  // namespace edm
