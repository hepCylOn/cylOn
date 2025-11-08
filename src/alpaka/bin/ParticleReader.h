#ifndef bin_particleReader_h
#define bin_particleReader_h

#include <cassert>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdexcept>

#include "AlpakaDataFormats/ParticleHost.h"  

#define INPUT_DEBUG

namespace particleReader {

  // --- file header constants ---
  constexpr uint32_t kExpectedEndianness = 0x01020304;
  constexpr uint32_t kFormatVersion      = 1;
  constexpr char     kMagic[4]           = {'P', 'A', 'R', '1'};  // "PAR1" = Particle binary format v1

  // --- validate the header ---
  inline void check_header(std::ifstream& in) {
    char magic[4];
    in.read(magic, 4);
    if (std::memcmp(magic, kMagic, 4) != 0)
      throw std::runtime_error("Invalid file magic (not a PAR1 binary)");

    uint32_t version;
    in.read(reinterpret_cast<char*>(&version), sizeof(version));
    if (version != kFormatVersion)
      throw std::runtime_error("Unsupported PAR1 version");

    uint32_t endian_marker;
    in.read(reinterpret_cast<char*>(&endian_marker), sizeof(endian_marker));
    if (endian_marker != kExpectedEndianness)
      throw std::runtime_error("Endianness mismatch — file not native endian");

#ifdef INPUT_DEBUG
    std::cout << "[particleReader] Input file header OK" << std::endl;
#endif
  }

  // --- read one column directly into the host view ---
  template <typename T>
  inline void read_column(std::ifstream& in, T* data, uint32_t n) {
    in.read(reinterpret_cast<char*>(data), n * sizeof(T));
    if (!in)
      throw std::runtime_error("Error reading particle column data");
  }

  // --- read a single event into a ParticleHost ---
  inline sim::ParticleHost read_single_event(std::ifstream& in) {
    uint32_t nParticles;
    in.read(reinterpret_cast<char*>(&nParticles), sizeof(nParticles));
    if (!in)
      throw std::runtime_error("Error reading particle event header");

#ifdef INPUT_DEBUG
    std::cout << "[particleReader] Reading event with " << nParticles << " particles" << std::endl;
#endif

    // Create host collection (on CPU queue)
    sim::ParticleHost particleHost(int(nParticles), cms::alpakatools::host());

    auto view = particleHost.view();

    read_column(in, view.vx().data(),      nParticles);
    read_column(in, view.vy().data(),      nParticles);
    read_column(in, view.vz().data(),      nParticles);

    read_column(in, view.px().data(),      nParticles);
    read_column(in, view.py().data(),      nParticles);
    read_column(in, view.pz().data(),      nParticles);
    read_column(in, view.energy().data(),  nParticles);

    read_column(in, view.pt().data(),      nParticles);
    read_column(in, view.eta().data(),     nParticles);
    read_column(in, view.phi().data(),     nParticles);
    read_column(in, view.mass().data(),    nParticles);

    read_column(in, view.charge().data(),  nParticles);
    read_column(in, view.pdgID().data(),   nParticles);
    read_column(in, view.partInd().data(), nParticles);

#ifdef INPUT_DEBUG
    if (nParticles > 0) {
      std::cout << "  First particle: "
                << "pdgID=" << view.pdgID()[0]
                << " pt=" << view.pt()[0]
                << " eta=" << view.eta()[0]
                << " phi=" << view.phi()[0]
                << " charge=" << view.charge()[0]
                << std::endl;
    }
#endif

    return particleHost;
  }

}  // namespace particleReader

#endif  // bin_particleReader_h
