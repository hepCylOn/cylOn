#include <cassert>
#include <chrono>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <ios>
#include <memory>
#include <mutex>
#include <istream>
#include <streambuf>
#include <stdexcept>
#include <zlib.h>

#include "Source.h"
#include "ParticleReader.h"
#include "HitReader.h"
#include "HitReaderTest.h"

// #define INPUT_DEBUG

namespace {
  FEDRawDataCollection readRaw(std::istream &is, unsigned int nfeds) {
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

  class GzipInputBuffer : public std::streambuf {
  public:
    explicit GzipInputBuffer(std::filesystem::path const& path)
        : path_(path), file_(gzopen(path.c_str(), "rb")) {
      if (file_ == nullptr) {
        throw std::runtime_error("Cannot open gzip input file: " + path.string());
      }
      setg(buffer_, buffer_, buffer_);
    }

    ~GzipInputBuffer() override {
      if (file_ != nullptr) {
        gzclose(file_);
      }
    }

  protected:
    int_type underflow() override {
      if (gptr() < egptr()) {
        return traits_type::to_int_type(*gptr());
      }

      int count = gzread(file_, buffer_, static_cast<unsigned>(sizeof(buffer_)));
      if (count < 0) {
        int errorCode;
        char const* message = gzerror(file_, &errorCode);
        throw std::runtime_error(
            "Error reading gzip input file " + path_.string() + ": " + message);
      }
      if (count == 0) {
        return traits_type::eof();
      }

      setg(buffer_, buffer_, buffer_ + count);
      return traits_type::to_int_type(*gptr());
    }

  private:
    std::filesystem::path path_;
    gzFile file_;
    char buffer_[8192];
  };

  class GzipInputStream : public std::istream {
  public:
    explicit GzipInputStream(std::filesystem::path const& path)
        : std::istream(nullptr), buffer_(path) {
      rdbuf(&buffer_);
    }

  private:
    GzipInputBuffer buffer_;
  };

  std::unique_ptr<std::istream> openInputFile(
    std::filesystem::path const& path,
    bool binary = false) {
    if (path.extension() == ".gz") {
      return std::make_unique<GzipInputStream>(path);
    }

    auto mode = std::ios::in;
    if (binary) {
      mode |= std::ios::binary;
    }

    auto file = std::make_unique<std::ifstream>(path, mode);
    if (!*file) {
      throw std::runtime_error("Cannot open input file: " + path.string());
    }
    return file;
  }

}  // namespace

namespace edm {
  Source::Source(
      int maxEvents, int runForMinutes, ProductRegistry &reg, std::filesystem::path const &datadir, bool validation, bool fromHits, bool isPhase2, bool runSimTracks, bool isColliderML, std::string const& trackerType, std::string const& inputFileName, std::string const& mapFileName, std::string const& partFileName) // Change to work with Phase-2
      : maxEvents_(maxEvents),
        runForMinutes_(runForMinutes),
        validation_(validation),
        fromHits_(fromHits),
        isPhase2_(isPhase2),
        runSimTracks_(runSimTracks),
        isColliderML_(isColliderML),
        trackerType_(trackerType),
        inputFileName_(inputFileName),
        mapFileName_(mapFileName),
        partFileName_(partFileName) {

    std::unique_ptr<std::istream> in_file;

    std::string in_fileName;

    if (not fromHits_)
    {
      in_fileName = "raw.bin";
      if (!inputFileName.empty()) {
        in_fileName = inputFileName;
      }
      auto inputPath = datadir / in_fileName;
      in_file = openInputFile(inputPath, !fromHits_);
      std::cout << "Reading bin from " << datadir / in_fileName << std::endl;
      rawToken_ = reg.produces<FEDRawDataCollection>();
    }
    else
    {
      if (!inputFileName.empty()) {
        in_fileName = inputFileName;
      }
      auto inputPath = datadir / in_fileName;
      in_file = openInputFile(inputPath, fromHits_);
      std::cout << "Reading hits from " << datadir / in_fileName  << std::endl;
      hitToken_ = reg.produces<reco::TrackingRecHitHost>();
    }

    std::ifstream in_digiclusters;
    std::ifstream in_tracks;
    std::ifstream in_vertices;
    std::unique_ptr<std::istream> in_particles;
    std::unique_ptr<std::istream> in_map;

    if (validation_) {
      if (not fromHits_){

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

      else{

        if (not isColliderML) {
          trackToken_ = reg.produces<TrackCount>();
          vertexToken_ = reg.produces<VertexCount>();

          in_tracks = std::ifstream(datadir / "tracks.bin", std::ios::binary);
          in_vertices = std::ifstream(datadir / "vertices.bin", std::ios::binary);
          
          in_tracks.exceptions(std::ifstream::badbit | std::ifstream::failbit | std::ifstream::eofbit);
          in_vertices.exceptions(std::ifstream::badbit | std::ifstream::failbit | std::ifstream::eofbit);
        }
        
        else {
          if (!mapFileName.empty()) {
            auto mapPath = datadir / mapFileName;
            in_map = openInputFile(mapPath, !fromHits_);
          }
          if (!partFileName.empty()) {
            auto partPath = datadir / partFileName;
            in_particles = openInputFile(partPath, !fromHits_);
          }
        }

        particleToken_ = reg.produces<sim::ParticleHost>();
        mapToken_ = reg.produces<utils::SimpleMapHost>();

      }
      
    }

    if(not fromHits_)
    {
      unsigned int nfeds;
      in_file->exceptions(std::ifstream::badbit);
      in_file->read(reinterpret_cast<char *>(&nfeds), sizeof(unsigned int));
      while (not in_file->eof()) {
        in_file->exceptions(std::ifstream::badbit | std::ifstream::failbit | std::ifstream::eofbit);

        raw_.emplace_back(readRaw(*in_file, nfeds));

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
        in_file->exceptions(std::ifstream::badbit);
        in_file->read(reinterpret_cast<char *>(&nfeds), sizeof(unsigned int));
      }
    }
    else
    {
      int32_t nEventsP, nEventsH, nEventsM;

      std::string line_hits;
      while (std::getline(*in_file, line_hits)) {
        if (!line_hits.empty()) break;
      }
      if (line_hits.rfind("events:", 0) != 0) {
        std::cerr << "Hits check -- Expected 'events:NEVENTS', got: " << line_hits << "\n";
      }

      nEventsH = std::stoul(line_hits.substr(7)); // after "events:"

      if (isColliderML and validation) {
        std::string line_map;
        while (std::getline(*in_map, line_map)) {
          if (!line_map.empty()) break;
        }
        if (line_map.rfind("events:", 0) != 0) {
          std::cerr << "Map check -- Expected 'events:NEVENTS', got: " << line_map << "\n";
        }

        nEventsM = std::stoul(line_map.substr(7)); // after "events:"
        nEventsP = nEventsH;

        if(nEventsH != nEventsP or nEventsH != nEventsM)
          throw std::runtime_error("Error nEvents differs in the hits file and the particles file!");
      }
      maxEvents_ = (maxEvents_ < 0) ? nEventsH : std::min(maxEvents_, nEventsH);

#ifdef INPUT_DEBUG
      std::cout << "File contains " << nEventsH << " events\n";
#endif
      auto beginReadingTime = std::chrono::steady_clock::now();
      for (int32_t ev = 0; ev < nEventsH && ev < maxEvents_; ++ev) {
        // hits_.emplace_back(hitReader::read_single_event(in_file));

        // hits_.emplace_back(hitReader::read_single_event_fromText(in_file));
        hits_.emplace_back(hitReaderTest::readEvent(*in_file));
        if (validation and isColliderML)
        {
          particles_.emplace_back(particleReader::read_single_event_fromText(*in_particles));
          // map_.emplace_back(mapReader::read_single_event_fromText(in_map));
          map_.emplace_back(mapReaderTest::readEvent(*in_map));
        }
          
#ifdef INPUT_DEBUG
        if (validation and isColliderML)
          std::cout << "Event " << ev << ": " << hits_[ev].nHits() << " hits, " << hits_[ev].nModules() << " modules - n. particles = " << particles_[ev].view().metadata().size() << std::endl;
        else
          std::cout << "Event " << ev << ": " << hits_[ev].nHits() << " hits, " << hits_[ev].nModules() << " modules\n";
#endif
      }
      auto endReadingTime = std::chrono::steady_clock::now();
      auto readingTime = std::chrono::duration_cast<std::chrono::seconds>(endReadingTime - beginReadingTime);
      std::cout << "It took " << readingTime.count() << " seconds to read the input files" << std::endl;
      
#ifdef INPUT_DEBUG
      std::cout << "Successfully read all events from " << in_fileName << std::endl;
#endif
    }

    if (validation_ and not fromHits_) { //TODO allow for fromHits validation
      assert(raw_.size() == digiclusters_.size());
      assert(raw_.size() == tracks_.size());
      assert(raw_.size() == vertices_.size());
    } 
    else if (validation_ and isColliderML_)
    {
      assert(hits_.size() == particles_.size());
      assert(hits_.size() == map_.size());
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
    if (validation_ and not isColliderML_) {
      if (not fromHits_) ev->emplace(digiClusterToken_, digiclusters_[index]);
      ev->emplace(trackToken_, tracks_[index]);
      ev->emplace(vertexToken_, vertices_[index]);
    }
    else if (validation_ and isColliderML_)
    {
      ev->emplace(particleToken_, std::move(particles_[index]));
      ev->emplace(mapToken_, std::move(map_[index]));
    }

    return ev;
  }
}  // namespace edm
