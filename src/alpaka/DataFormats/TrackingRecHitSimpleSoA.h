#ifndef DataFormats_TrackingRecHitSimpleSoA_h
#define DataFormats_TrackingRecHitSimpleSoA_h

#include <cstdint>
#include <cstdlib>
#include <vector>
#include <iostream>

class TrackingRecHitSimpleSoA {
public:
  TrackingRecHitSimpleSoA() = default;
  explicit TrackingRecHitSimpleSoA(unsigned int nHits)
  {
    m_xl.resize(nHits);
    m_yl.resize(nHits);
    m_xerr.resize(nHits);
    m_yerr.resize(nHits);

    m_xg.resize(nHits);
    m_yg.resize(nHits);
    m_zg.resize(nHits);
    m_rg.resize(nHits);
    m_iphi.resize(nHits);

    m_charge.resize(nHits);
    m_xsize.resize(nHits);
    m_ysize.resize(nHits);
    m_detInd.resize(nHits);
  }

  explicit TrackingRecHitSimpleSoA(
      size_t nHits,
      const float* xl, const float* yl, const float* xerr, const float* yerr,
      const float* xg, const float* yg, const float* zg, const float* rg, const int16_t* iphi,
      const int32_t* charge, const int16_t* xsize, const int16_t* ysize, const int16_t* detInd,
      const uint32_t* modStart)
      : m_xl(xl, xl + nHits),
        m_yl(yl, yl + nHits),
        m_xerr(xerr, xerr + nHits),
        m_yerr(yerr, yerr + nHits),
        m_xg(xg, xg + nHits),
        m_yg(yg, yg + nHits),
        m_zg(zg, zg + nHits),
        m_rg(rg, rg + nHits),
        m_iphi(iphi, iphi + nHits),
        m_charge(charge, charge + nHits),
        m_xsize(xsize, xsize + nHits),
        m_ysize(ysize, ysize + nHits),
        m_detInd(detInd, detInd + nHits),
        m_moduleStart(modStart, modStart + nHits)
  {
    assert(m_xl.size() == nHits);
    //TODO add assert for modStart[-1] == nHits;
  }
  ~TrackingRecHitSimpleSoA() = default;

  auto nHits() const { return m_xl.size(); }

  float xl(size_t i) const { return m_xl[i]; }
  float yl(size_t i) const { return m_yl[i]; }
  float xerr(size_t i) const { return m_xerr[i]; }
  float yerr(size_t i) const { return m_yerr[i]; }

  float xg(size_t i) const { return m_xg[i]; }
  float yg(size_t i) const { return m_yg[i]; }
  float zg(size_t i) const { return m_zg[i]; }
  float rg(size_t i) const { return m_rg[i]; }
  int16_t iphi(size_t i) const { return m_iphi[i]; }

  int32_t charge(size_t i) const { return m_charge[i]; }
  int16_t xsize(size_t i) const { return m_xsize[i]; }
  int16_t ysize(size_t i) const { return m_ysize[i]; }
  int16_t detInd(size_t i) const { return m_detInd[i]; }

  std::vector<float>& xlVector() { return m_xl; }
  const std::vector<float>& xlVector() const { return m_xl; }

  std::vector<float>& ylVector() { return m_yl; }
  const std::vector<float>& ylVector() const { return m_yl; }

  std::vector<float>& xerrVector() { return m_xerr; }
  const std::vector<float>& xerrVector() const { return m_xerr; }

  std::vector<float>& yerrVector() { return m_yerr; }
  const std::vector<float>& yerrVector() const { return m_yerr; }

  std::vector<float>& xgVector() { return m_xg; }
  const std::vector<float>& xgVector() const { return m_xg; }

  std::vector<float>& ygVector() { return m_yg; }
  const std::vector<float>& ygVector() const { return m_yg; }

  std::vector<float>& zgVector() { return m_zg; }
  const std::vector<float>& zgVector() const { return m_zg; }

  std::vector<float>& rgVector() { return m_rg; }
  const std::vector<float>& rgVector() const { return m_rg; }

  std::vector<int16_t>& iphiVector() { return m_iphi; }
  const std::vector<int16_t>& iphiVector() const { return m_iphi; }

  std::vector<int32_t>& chargeVector() { return m_charge; }
  const std::vector<int32_t>& chargeVector() const { return m_charge; }

  std::vector<int16_t>& xsizeVector() { return m_xsize; }
  const std::vector<int16_t>& xsizeVector() const { return m_xsize; }

  std::vector<int16_t>& ysizeVector() { return m_ysize; }
  const std::vector<int16_t>& ysizeVector() const { return m_ysize; }

  std::vector<int16_t>& detIndVector() { return m_detInd; }
  const std::vector<int16_t>& detIndVector() const { return m_detInd; }

  std::vector<uint32_t>& moduleStartVec() { return m_moduleStart; }
  const std::vector<uint32_t>& moduleStartVec() const { return m_moduleStart; }

  void setHits(uint32_t nHits)
  {
    m_xl.resize(nHits);
    m_yl.resize(nHits);
    m_xerr.resize(nHits);
    m_yerr.resize(nHits);

    m_xg.resize(nHits);
    m_yg.resize(nHits);
    m_zg.resize(nHits);
    m_rg.resize(nHits);
    m_iphi.resize(nHits);

    m_charge.resize(nHits);
    m_xsize.resize(nHits);
    m_ysize.resize(nHits);
    m_detInd.resize(nHits);
  }

  void readBinary(std::istream& is) {

    auto nHits = this->nHits();
    auto readVector = [&is, nHits](auto& vec) {
      if (nHits > 0) {
        is.read(reinterpret_cast<char*>(vec.data()),
                vec.size() * sizeof(typename std::decay<decltype(vec)>::type::value_type));
      }
    };

    readVector(m_xl);
    readVector(m_yl);
    readVector(m_xerr);
    readVector(m_yerr);

    readVector(m_xg);
    readVector(m_yg);
    readVector(m_zg);
    readVector(m_rg);
    readVector(m_iphi);

    readVector(m_charge);
    readVector(m_xsize);
    readVector(m_ysize);
    readVector(m_detInd);
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

     if (line.rfind("hits:", 0) != 0) {
        std::cerr << "Expected 'hits:NHITS', got: " << line << "\n";
        return false;
    }

    size_t nHits = std::stoul(line.substr(5)); // after "hits:"

    // --- read NHITS lines ---
    m_xl.reserve(nHits);
    m_yl.reserve(nHits);
    m_xerr.reserve(nHits);
    m_yerr.reserve(nHits);
    m_xg.reserve(nHits);
    m_yg.reserve(nHits);
    m_zg.reserve(nHits);
    m_rg.reserve(nHits);
    m_iphi.reserve(nHits);
    m_charge.reserve(nHits);
    m_xsize.reserve(nHits);
    m_ysize.reserve(nHits);
    m_detInd.reserve(nHits);

    for (size_t i = 0; i < nHits; ++i) {
        if (!std::getline(file, line)) {
            std::cerr << "Unexpected end of file while reading hits.\n";
            return false;
        }
        auto tokens = split(line);
        if (tokens.size() != 13) {
            std::cerr << "Hit row " << i << " has " << tokens.size()
                      << " columns, expected 13.\n";
            return false;
        }

        m_xl.push_back(std::stof(tokens[0]));
        m_yl.push_back(std::stof(tokens[1]));
        m_xerr.push_back(std::stof(tokens[2]));
        m_yerr.push_back(std::stof(tokens[3]));
        m_xg.push_back(std::stof(tokens[4]));
        m_yg.push_back(std::stof(tokens[5]));
        m_zg.push_back(std::stof(tokens[6]));
        m_rg.push_back(std::stof(tokens[7]));
        m_iphi.push_back(static_cast<int16_t>(std::stoi(tokens[8])));
        m_charge.push_back(std::stoi(tokens[9]));
        m_xsize.push_back(static_cast<int16_t>(std::stoi(tokens[10])));
        m_ysize.push_back(static_cast<int16_t>(std::stoi(tokens[11])));
        m_detInd.push_back(static_cast<int16_t>(std::stoi(tokens[12])));
    }

    // --- read module header ---
    if (!std::getline(file, line)) {
        std::cerr << "Missing module header.\n";
        return false;
    }
    if (line.rfind("module:", 0) != 0) {
        std::cerr << "Expected 'module:NMODULES' header, got: " << line << "\n";
        return false;
    }
    size_t nModules = std::stoul(line.substr(7)); // after "module:"

    // --- read one line with NMODULES+1 entries ---
    if (!std::getline(file, line)) {
        std::cerr << "Missing module start line.\n";
        return false;
    }
    auto tokens = split(line);
    if (tokens.size() != nModules + 1) {
        std::cerr << "Module start line has " << tokens.size()
                  << " entries, expected " << (nModules + 1) << ".\n";
        return false;
    }

    m_moduleStart.reserve(nModules + 1);
    for (auto &tok : tokens) {
        m_moduleStart.push_back(static_cast<uint32_t>(std::stoul(tok)));
    }

    return true;
  }

private:

    // local coord
  std::vector<float> m_xl;
  std::vector<float> m_yl;
  std::vector<float> m_xerr;
  std::vector<float> m_yerr;

  // global coord
  std::vector<float> m_xg;
  std::vector<float> m_yg;
  std::vector<float> m_zg;
  std::vector<float> m_rg;
  std::vector<int16_t> m_iphi;

  // cluster properties
  std::vector<int32_t> m_charge;
  std::vector<int16_t> m_xsize;
  std::vector<int16_t> m_ysize;
  std::vector<int16_t> m_detInd;

  std::vector<uint32_t> m_moduleStart;

};

#endif  // DataFormats_TrackingRecHitSimpleSoA_h
