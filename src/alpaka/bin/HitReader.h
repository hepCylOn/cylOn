#ifndef Framework_HitReader_h
#define Framework_HitReader_h

#include <cassert>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "AlpakaDataFormats/TrackingRecHitsHost.h"
#include "AlpakaDataFormats/SimpleMapSoA.h"

#define INPUT_DEBUG

namespace mapReader {

  constexpr uint32_t kExpectedEndianness = 0x01020304;
  constexpr uint32_t kFormatVersion      = 1;
  constexpr char     kMagic[4]           = {'M', 'A', 'P', '1'};  // "MAP1" magic

  inline void check_header(std::ifstream& in) {
    char magic[4];
    in.read(magic, 4);
    if (std::memcmp(magic, kMagic, 4) != 0)
      throw std::runtime_error("Invalid file magic (not a MAP1 binary)");

    uint32_t version;
    in.read(reinterpret_cast<char*>(&version), sizeof(version));
    if (version != kFormatVersion)
      throw std::runtime_error("Unsupported MAP version");

    uint32_t endian;
    in.read(reinterpret_cast<char*>(&endian), sizeof(endian));
    if (endian != kExpectedEndianness)
      throw std::runtime_error("Endianness mismatch — file not native endian");

#ifdef INPUT_DEBUG
    std::cout << "[mapReader] Input file header OK" << std::endl;
#endif
  }

  // --- read one event into a SimpleMapHost ---
  inline utils::SimpleMapHost read_single_event(std::ifstream& in) {
    uint32_t nEntries;
    in.read(reinterpret_cast<char*>(&nEntries), sizeof(nEntries));
    if (!in)
      throw std::runtime_error("Error reading map event header (nEntries)");

#ifdef INPUT_DEBUG
    std::cout << "[mapReader] Reading event with " << nEntries << " entries" << std::endl;
#endif

    utils::SimpleMapHost host(static_cast<int>(nEntries), cms::alpakatools::host());
    auto view = host.view();

    // Read the two columns (key, value)
    in.read(reinterpret_cast<char*>(view.id().data()), nEntries * sizeof(uint32_t));

    if (!in)
      throw std::runtime_error("Error reading map column data");

#ifdef INPUT_DEBUG
    if (nEntries > 0) {
      std::cout << "  First entry: particleIdx=" << view.id()[0] << std::endl;
    }
#endif

    return host;
  }

}

namespace hitReader {

  // --- file header constants ---
  constexpr uint32_t kExpectedEndianness = 0x01020304;
  constexpr uint32_t kFormatVersion      = 1;
  constexpr char     kMagic[4]           = {'T', 'R', 'H', '1'};  // "TRH1" = TrackingRecHits v1

  // --- header check ---
  inline void check_header(std::ifstream& in) {
    char magic[4];
    in.read(magic, 4);
	std::cout << magic << std::endl;
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
    std::cout << "[hitReader] Input file header OK" << std::endl;
#endif
  }

  // --- helper to read a column directly into the SoA view ---
  template <typename T>
  inline void read_column(std::ifstream& in, T* data, uint32_t n) {
    in.read(reinterpret_cast<char*>(data), n * sizeof(T));
    if (!in)
      throw std::runtime_error("Error reading hit column data");
  }

  // --- read a single event into a TrackingRecHitHost ---
  inline reco::TrackingRecHitHost read_single_event(std::ifstream& in) {
    uint32_t nHits, nModules;
    in.read(reinterpret_cast<char*>(&nHits), sizeof(nHits));
    in.read(reinterpret_cast<char*>(&nModules), sizeof(nModules));
    if (!in)
      throw std::runtime_error("Error reading event header (nHits/nModules)");

#ifdef INPUT_DEBUG
    std::cout << "[hitReader] Reading event with "
              << nHits << " hits, " << nModules << " modules" << std::endl;
#endif

    // Read moduleStart (nModules + 1 entries)
    std::vector<uint32_t> moduleStart(nModules + 1);
    in.read(reinterpret_cast<char*>(moduleStart.data()), (nModules + 1) * sizeof(uint32_t));
    if (!in)
      throw std::runtime_error("Error reading moduleStart array");

    // Construct the host collection (nHits, nModules)
    reco::TrackingRecHitHost recHitHost(cms::alpakatools::host(), nHits, nModules);

    auto hitsView = recHitHost.view<reco::TrackingRecHitSoA>();
    auto modsView = recHitHost.view<reco::HitModuleSoA>();

    // Copy moduleStart array
    std::memcpy(modsView.moduleStart().data(), moduleStart.data(),
                (nModules + 1) * sizeof(uint32_t));

    // Read hit columns sequentially (in the same order as written)
    std::cout << "Reading hit columns for " << nHits << " hits" << std::endl;
    read_column(in, hitsView.xLocal().data(), nHits);
    std::cout << "Read xLocal" << std::endl;
    read_column(in, hitsView.yLocal().data(), nHits);
    std::cout << "Read yLocal" << std::endl;
    read_column(in, hitsView.xerrLocal().data(), nHits);
    std::cout << "Read xerrLocal" << std::endl;
    read_column(in, hitsView.yerrLocal().data(), nHits);
    std::cout << "Read yerrLocal" << std::endl;
    read_column(in, hitsView.xGlobal().data(), nHits);
    std::cout << "Read xGlobal" << std::endl;
    read_column(in, hitsView.yGlobal().data(), nHits);
    std::cout << "Read yGlobal" << std::endl;
    read_column(in, hitsView.zGlobal().data(), nHits);
    std::cout << "Read zGlobal" << std::endl;
    read_column(in, hitsView.rGlobal().data(), nHits);
    std::cout << "Read rGlobal" << std::endl;
    read_column(in, hitsView.iphi().data(), nHits);
    std::cout << "Read iphi" << std::endl;
    read_column(in, hitsView.chargeAndStatus().data(), nHits);
    std::cout << "Read chargeAndStatus" << std::endl;
    read_column(in, hitsView.clusterSizeX().data(), nHits);
    std::cout << "Read clusterSizeX" << std::endl;
    read_column(in, hitsView.clusterSizeY().data(), nHits);
    std::cout << "Read clusterSizeY" << std::endl;
    read_column(in, hitsView.detectorIndex().data(), nHits);
    std::cout << "Read detectorIndex" << std::endl;

#ifdef INPUT_DEBUG
    if (nHits > 0) {
      std::cout << "  First hit global: ("
                << hitsView.xGlobal()[0] << ", "
                << hitsView.yGlobal()[0] << ", "
                << hitsView.zGlobal()[0]
                << ")  r=" << hitsView.rGlobal()[0]
                << " detIdx=" << hitsView.detectorIndex()[0]
                << std::endl;
    }
#endif

    return recHitHost;
  }

}  //

#endif
