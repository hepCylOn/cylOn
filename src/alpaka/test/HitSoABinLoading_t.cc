#include <fstream>
#include <iostream>
#include <vector>
#include <cstdint>
#include <stdexcept>
#include <cstring>
#include <limits>
#include <cmath>

#include "AlpakaDataFormats/TrackingRecHitsHost.h"
#include "AlpakaCore/config.h"
#include "AlpakaCore/host.h"

using namespace reco;

#define GPU_DEBUG

namespace {

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

  // Skip number of events will be read by caller
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

}  // namespace

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " hits.bin [maxEvents]\n";
    return 1;
  }

  const char* filename = argv[1];
  uint32_t maxEvents = (argc > 2) ? std::atoi(argv[2]) : std::numeric_limits<uint32_t>::max();

  std::ifstream in(filename, std::ios::binary);
  if (!in)
    throw std::runtime_error(std::string("Cannot open file ") + filename);

  // ---- Read and validate header ----
  check_header(in);

  uint32_t nEvents;
  in.read(reinterpret_cast<char*>(&nEvents), sizeof(nEvents));
  std::cout << "File contains " << nEvents << " events\n";

  // ---- Loop over events ----
  for (uint32_t ev = 0; ev < nEvents && ev < maxEvents; ++ev) {
    uint32_t nHits, nModules;
    in.read(reinterpret_cast<char*>(&nHits), sizeof(nHits));
    in.read(reinterpret_cast<char*>(&nModules), sizeof(nModules));
    if (!in) throw std::runtime_error("Error reading event header");

    std::vector<uint32_t> moduleStart(nModules + 1);
    in.read(reinterpret_cast<char*>(moduleStart.data()), (nModules + 1) * sizeof(uint32_t));

    TrackingRecHitHost recHitHost(cms::alpakatools::host(), nHits, nModules);

    auto hitView = recHitHost.view<TrackingRecHitSoA>();
    auto modView = recHitHost.view<HitModuleSoA>();

    // copy module starts
    std::memcpy(modView.moduleStart().data(), moduleStart.data(),
                (nModules + 1) * sizeof(uint32_t));

    // ---- read all columns ----
    // ---- read all columns (SoA) ----
    read_column(in, hitView.xLocal(),       nHits);
    read_column(in, hitView.yLocal(),       nHits);
    read_column(in, hitView.xerrLocal(),    nHits);
    read_column(in, hitView.yerrLocal(),    nHits);
    read_column(in, hitView.xGlobal(),      nHits);
    read_column(in, hitView.yGlobal(),      nHits);
    read_column(in, hitView.zGlobal(),      nHits);
    read_column(in, hitView.rGlobal(),      nHits);
    read_column(in, hitView.iphi(),         nHits);
    read_column(in, hitView.chargeAndStatus(), nHits);   // <-- FIXED: deduce type
    read_column(in, hitView.clusterSizeX(), nHits);
    read_column(in, hitView.clusterSizeY(), nHits);
    read_column(in, hitView.detectorIndex(), nHits);

#ifdef GPU_DEBUG
    std::cout << "Event " << ev << ": " << nHits << " hits, " << nModules << " modules\n";
    std::cout << "  First hit global: ("
              << hitView.xGlobal()[0] << ", "
              << hitView.yGlobal()[0] << ", "
              << hitView.zGlobal()[0] << "), "
              << "r=" << hitView.rGlobal()[0]
              << ", detIdx=" << hitView.detectorIndex()[0] << '\n';
#endif
  }

  if (!in.good() && !in.eof()) {
    throw std::runtime_error("I/O error while reading file");
  }

  std::cout << "Successfully read all events from " << filename << std::endl;
  return 0;
}
