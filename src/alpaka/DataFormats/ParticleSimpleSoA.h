#ifndef DataFormats_ParticleSimpleSoA_h
#define DataFormats_ParticleSimpleSoA_h

#define NDEBUG
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <sstream>

class ParticleSimpleSoA {
public:
  ParticleSimpleSoA() = default;
  explicit ParticleSimpleSoA(unsigned int nParticles)
  {
    m_vx.resize(nParticles);
    m_vy.resize(nParticles);
    m_vz.resize(nParticles);

    m_px.resize(nParticles);
    m_py.resize(nParticles);
    m_pz.resize(nParticles);
    m_energy.resize(nParticles);

    m_pt.resize(nParticles);
    m_eta.resize(nParticles);
    m_phi.resize(nParticles);
    m_mass.resize(nParticles);

    m_charge.resize(nParticles);
    m_pdgID.resize(nParticles);

    m_partInd.resize(nParticles);
  }

  explicit ParticleSimpleSoA(
      size_t nParticles,
      const float* vx, const float* vy, const float* vz, const float* px,
      const float* py, const float* pz, const float* energy, const float* pt, const float* eta,
      const float* phi, const float* mass, const int16_t* charge, const int32_t* pdgID,
      const uint32_t* partInd)
      : m_vx(vx, vx + nParticles),
        m_vy(vy, vy + nParticles),
        m_vz(vz, vz + nParticles),
        m_px(px, px + nParticles),
        m_py(py, py + nParticles),
        m_pz(pz, pz + nParticles),
        m_energy(energy, energy + nParticles),
        m_pt(pt, pt + nParticles),
        m_eta(eta, eta + nParticles),
        m_phi(phi, phi + nParticles),
        m_mass(mass, mass + nParticles),
        m_charge(charge, charge + nParticles),
        m_pdgID(pdgID, pdgID + nParticles),
        m_partInd(partInd, partInd + nParticles)
  {
    assert(m_vx.size() == nParticles);
  }
  ~ParticleSimpleSoA() = default;

  auto nParticles() const { return m_vx.size(); }

  float vx(size_t i) const { return m_vx[i]; }
  float vy(size_t i) const { return m_vy[i]; }
  float vz(size_t i) const { return m_vz[i]; }

  float px(size_t i) const { return m_px[i]; }
  float py(size_t i) const { return m_py[i]; }
  float pz(size_t i) const { return m_pz[i]; }
  float energy(size_t i) const { return m_energy[i]; }

  float pt(size_t i) const { return m_pt[i]; }
  float eta(size_t i) const { return m_eta[i]; }
  float phi(size_t i) const { return m_phi[i]; }
  float mass(size_t i) const { return m_mass[i]; }

  int16_t charge(size_t i) const { return m_charge[i]; }
  int32_t pdgID(size_t i) const { return m_pdgID[i]; }

  uint32_t partInd(size_t i) const { return m_partInd[i]; }

  std::vector<float>& vxVector() { return m_vx; }
  const std::vector<float>& vxVector() const { return m_vx; }

  std::vector<float>& vyVector() { return m_vy; }
  const std::vector<float>& vyVector() const { return m_vy; }

  std::vector<float>& vzVector() { return m_vz; }
  const std::vector<float>& vzVector() const { return m_vz; }

  std::vector<float>& pxVector() { return m_px; }
  const std::vector<float>& pxVector() const { return m_px; }

  std::vector<float>& pyVector() { return m_py; }
  const std::vector<float>& pyVector() const { return m_py; }

  std::vector<float>& pzVector() { return m_pz; }
  const std::vector<float>& pzVector() const { return m_pz; }

  std::vector<float>& energyVector() { return m_energy; }
  const std::vector<float>& energyVector() const { return m_energy; }

  std::vector<float>& ptVector() { return m_pt; }
  const std::vector<float>& ptVector() const { return m_pt; }

  std::vector<float>& etaVector() { return m_eta; }
  const std::vector<float>& etaVector() const { return m_eta; }

  std::vector<float>& phiVector() { return m_phi; }
  const std::vector<float>& phiVector() const { return m_phi; }

  std::vector<float>& massVector() { return m_mass; }
  const std::vector<float>& massVector() const { return m_mass; }

  std::vector<int16_t>& chargeVector() { return m_charge; }
  const std::vector<int16_t>& chargeVector() const { return m_charge; }

  std::vector<int32_t>& pdgIDVector() { return m_pdgID; }
  const std::vector<int32_t>& pdgIDVector() const { return m_pdgID; }

  std::vector<uint32_t>& partIndVector() { return m_partInd; }
  const std::vector<uint32_t>& partIndVector() const { return m_partInd; }

  void setParticles(uint32_t nParticles)
  {
    m_vx.resize(nParticles);
    m_vy.resize(nParticles);
    m_vz.resize(nParticles);

    m_px.resize(nParticles);
    m_py.resize(nParticles);
    m_pz.resize(nParticles);
    m_energy.resize(nParticles);

    m_pt.resize(nParticles);
    m_eta.resize(nParticles);
    m_phi.resize(nParticles);
    m_mass.resize(nParticles);

    m_charge.resize(nParticles);
    m_pdgID.resize(nParticles);

    m_partInd.resize(nParticles);
  }

  void readBinary(std::istream& is) {

    auto nParticles = this->nParticles();
    auto readVector = [&is, nParticles](auto& vec) {
      if (nParticles > 0) {
        is.read(reinterpret_cast<char*>(vec.data()),
                vec.size() * sizeof(typename std::decay<decltype(vec)>::type::value_type));
      }
    };

    readVector(m_vx);
    readVector(m_vy);
    readVector(m_vz);

    readVector(m_px);
    readVector(m_py);
    readVector(m_pz);
    readVector(m_energy);

    readVector(m_pt);
    readVector(m_eta);
    readVector(m_phi);
    readVector(m_mass);

    readVector(m_charge);
    readVector(m_pdgID);

    readVector(m_partInd);
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

  bool readText(std::istream& file) {
    
    std::string line;

    // skip empty lines
    while (std::getline(file, line)) {
        if (!line.empty()) break;
    }
    if (file.eof()) return false;

     if (line.rfind("particles:", 0) != 0) {
        std::cerr << "Expected 'particles:nParticles', got: " << line << "\n";
        return false;
    }

    size_t nParticles = std::stoul(line.substr(10)); // after "particles:"

    // --- read nParticles lines ---
    m_vx.reserve(nParticles);
    m_vy.reserve(nParticles);
    m_vz.reserve(nParticles);
    m_px.reserve(nParticles);
    m_py.reserve(nParticles);
    m_pz.reserve(nParticles);
    m_energy.reserve(nParticles);
    m_pt.reserve(nParticles);
    m_eta.reserve(nParticles);
    m_phi.reserve(nParticles);
    m_mass.reserve(nParticles);
    m_charge.reserve(nParticles);
    m_pdgID.reserve(nParticles);
    m_partInd.reserve(nParticles);

    for (size_t i = 0; i < nParticles; ++i) {
        if (!std::getline(file, line)) {
            std::cerr << "Unexpected end of file while reading particles.\n";
            return false;
        }
        auto tokens = split(line);
        if (tokens.size() != 14) {
            std::cerr << "Hit row " << i << " has " << tokens.size()
                      << " columns, expected 14.\n";
            return false;
        }

        m_vx.push_back(std::stof(tokens[0]));
        m_vy.push_back(std::stof(tokens[1]));
        m_vz.push_back(std::stof(tokens[2]));
        m_px.push_back(std::stof(tokens[3]));
        m_py.push_back(std::stof(tokens[4]));
        m_pz.push_back(std::stof(tokens[5]));
        m_energy.push_back(std::stof(tokens[6]));
        m_pt.push_back(std::stof(tokens[7]));
        m_eta.push_back(std::stof(tokens[8]));
        m_phi.push_back(std::stof(tokens[9]));
        m_mass.push_back(std::stof(tokens[10]));
        m_charge.push_back(static_cast<int16_t>(std::stoi(tokens[11])));
        m_pdgID.push_back(static_cast<int32_t>(std::stoi(tokens[12])));
        m_partInd.push_back(static_cast<uint32_t>(std::stoi(tokens[13])));
    }

    return true;
  }

private:

  // vertex coord
  std::vector<float> m_vx;
  std::vector<float> m_vy;
  std::vector<float> m_vz;

  // four-momentum
  std::vector<float> m_px;
  std::vector<float> m_py;
  std::vector<float> m_pz;
  std::vector<float> m_energy;

  // HEP coord and mass
  std::vector<float> m_pt;
  std::vector<float> m_eta;
  std::vector<float> m_phi;
  std::vector<float> m_mass;

  // charge & PDG ID
  std::vector<int16_t> m_charge;
  std::vector<int32_t> m_pdgID;

  // particle ID
  std::vector<uint32_t> m_partInd;

};

#endif  // DataFormats_ParticleSimpleSoA_h
