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

  // Split a line by delimiter (default = comma)
  std::vector<std::string> split(const std::string &line, char delimiter = ',') {
      std::vector<std::string> tokens;
      std::stringstream ss(line);
      std::string item;
      while (std::getline(ss, item, delimiter)) {
          if (!item.empty()) tokens.push_back(item);
      }
      return tokens;
  }
inline sim::ParticleHost read_single_event_fromText(std::ifstream& file) {

    std::string line;

    // Used to indicate something broke when reading file
    sim::ParticleHost auxParticleHost(0, cms::alpakatools::host());

    // skip empty lines
    while (std::getline(file, line)) {
        if (!line.empty()) break;
    }
    if (file.eof()) return auxParticleHost;

     if (line.rfind("particles:", 0) != 0) {
        std::cerr << "Expected 'particles:NPARTICLES', got: " << line << "\n";
        return auxParticleHost;
    }

    uint32_t nParticles;

    nParticles = std::stoul(line.substr(10)); // after "particles:"

    // Construct the host collection (nParticles, nModules)
    sim::ParticleHost particlesHost(static_cast<int>(nParticles), cms::alpakatools::host());

    auto particlesView = particlesHost.view();

    for (size_t i = 0; i < nParticles; ++i) {
        if (!std::getline(file, line)) {
            std::cerr << "Unexpected end of file while reading particles.\n";
            return auxParticleHost;
        }
        auto tokens = split(line);
        if (tokens.size() != 14) {
            std::cerr << "Hit row " << i << " has " << tokens.size()
                      << " columns, expected 14.\n";
            return auxParticleHost;
        }

        particlesView[i].vx() = std::stof(tokens[0]);
        particlesView[i].vy() = std::stof(tokens[1]);
        particlesView[i].vz() = std::stof(tokens[2]);
        particlesView[i].px() = std::stof(tokens[3]);
        particlesView[i].py() = std::stof(tokens[4]);
        particlesView[i].pz() = std::stof(tokens[5]);
        particlesView[i].energy() = std::stof(tokens[6]);
        particlesView[i].pt() = std::stof(tokens[7]);
        particlesView[i].eta() = std::stof(tokens[8]);
        particlesView[i].phi() = std::stof(tokens[9]);
        particlesView[i].mass() = std::stof(tokens[10]);
        particlesView[i].charge() = static_cast<int16_t>(std::stoi(tokens[11]));
        particlesView[i].pdgID() = static_cast<int32_t>(std::stoi(tokens[12]));
        particlesView[i].partInd() = static_cast<uint32_t>(std::stoi(tokens[13]));

    }

    return particlesHost;
  }

}  // namespace particleReader

#endif  // bin_particleReader_h
