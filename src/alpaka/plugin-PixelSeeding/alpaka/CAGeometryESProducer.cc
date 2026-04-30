#include <memory>
#include <string>
#include <vector>
#include <algorithm>

#include "AlpakaCore/config.h"
#include "Framework/ESProducer.h"
#include "Framework/EventSetup.h"
#include "Framework/ESPluginFactory.h"
#include "Framework/ConfigRegistry.h"

#include "AlpakaDataFormats/CAGeometryHost.h"
#include "AlpakaDataFormats/CAGeometrySoA.h"
#include "AlpakaDataFormats/alpaka/CAGeometrySoACollection.h"
#include "AlpakaDataFormats/SOARotation.h"
#include "Geometry/SimplePixelTopology.h"
#include "Framework/StreamFileUtils.h"

// #define GPU_DEBUG
namespace ALPAKA_ACCELERATOR_NAMESPACE {

  template <typename TrackerTraits>
  class CAGeometryHostESProducer : public edm::ESProducer {
  public:
    explicit CAGeometryHostESProducer(edm::Config const& cfg) : 
      data_(static_cast<std::string>(cfg.value("data", defaultPath_))),
      nLayers_(static_cast<int>(cfg.value("nLayers", TrackerTraits::numberOfLayers))),
      nPairs_(static_cast<int>(cfg.value("nPairs", TrackerTraits::nPairs))),
      nModules_(static_cast<int>(cfg.value("nModules", TrackerTraits::numberOfModules)))
      {

      auto getVectorOrDefault = [&](std::string const& key, auto const* defaults, size_t expected) {
      using T = std::remove_cvref_t<decltype(defaults[0])>;
      std::vector<T> vec;

      if (cfg.contains(key)) {
        try {
          vec = cfg.at(key).get<std::vector<T>>();
          if (vec.size() != expected) {
            std::cerr << "[CAGeometryHostESProducer ERROR] Size mismatch in " << key << ":\n"
                      << "  expected = " << expected << ", actual = " << vec.size() << std::endl;
            std::abort();
          }
        } catch (std::exception const& e) {
          std::cerr << "[CAGeometryHostESProducer WARNING] Failed to read key '" << key
                    << "': " << e.what() << ". Using default values." << std::endl;
          vec.assign(defaults, defaults + expected);
        }
      } else {
        vec.assign(defaults, defaults + expected);
      }

      return vec;
    };

  thetaCuts_   = getVectorOrDefault("thetaCuts", TrackerTraits::thetaCuts, nLayers_);
  dcaCuts_     = getVectorOrDefault("dcaCuts", TrackerTraits::dcaCuts, nLayers_);
  layerStarts_ = getVectorOrDefault("layerStarts", TrackerTraits::layerStart, nLayers_ + 1);
  phiCuts_     = getVectorOrDefault("phiCuts", TrackerTraits::phicuts, nPairs_);
  minZ_        = getVectorOrDefault("minZ", TrackerTraits::minz, nPairs_);
  maxZ_        = getVectorOrDefault("maxZ", TrackerTraits::maxz, nPairs_);
  maxR_        = getVectorOrDefault("maxR", TrackerTraits::maxr, nPairs_);
  pairGraph_   = getVectorOrDefault("pairGraph", TrackerTraits::layerPairs, nPairs_ * 2);
  startingPairs_ = cfg.contains("startingPairs")
                   ? cfg.at("startingPairs").get<std::vector<uint8_t>>()
                   : std::vector<uint8_t>{0u, 1u, 2u};

auto maxVal = std::ranges::max(startingPairs_);
  if (maxVal >= static_cast<unsigned int>(nPairs_)) {
    std::cerr << "[CAGeometryHostESProducer ERROR] Invalid 'startingPairs' values: "
              << "maximum value " << maxVal << " exceeds nPairs = " << nPairs_ << std::endl;
    std::abort();
  }


  #ifdef GPU_DEBUG

    auto vecToString = [](auto const& v) {
    std::ostringstream os;
    os << "{";
    for (size_t i = 0; i < v.size(); ++i) {
      os << v[i];
      if (i + 1 < v.size()) os << ", ";
    }
    os << "}";
    return os.str();
  };

  std::cout << "[GPU_DEBUG] Constructing CAGeometryHostESProducer for "
            << TrackerTraits::nameModifier << "\n"
            << "  - Data directory: " << data_ << "\n"
            << "  - Layers: " << nLayers_
            << ", Pairs: " << nPairs_
            << ", Starting pairs: " << vecToString(startingPairs_) << std::endl;
#endif
}

    void produce(edm::EventSetup& eventSetup);

  private:

    std::filesystem::path data_;
    int nLayers_, nPairs_, nModules_; //int(s) just because the PortableCollection wants so
    
    std::filesystem::path defaultPath_ = std::string("data/CAGeometryHostModules") + std::string(TrackerTraits::nameModifier) + ".bin";;
    // Geometry parameter data members
    std::vector<float> thetaCuts_;
    std::vector<float> dcaCuts_;
    std::vector<unsigned int> layerStarts_;
    
    std::vector<short int> phiCuts_;
    std::vector<float> minZ_;
    std::vector<float> maxZ_;
    std::vector<float> maxR_;
    std::vector<uint8_t> pairGraph_;
    std::vector<uint8_t> startingPairs_;
  };

  // ---------- Produce ----------
  template <typename TrackerTraits>
  void CAGeometryHostESProducer<TrackerTraits>::produce(edm::EventSetup& eventSetup) {

    #ifdef GPU_DEBUG
    std::cout << "[GPU_DEBUG] Producing CAGeometryHost for " << TrackerTraits::nameModifier << "\n"
              << "  - Reading modules from: " << data_ << std::endl;
    #endif
        
    // Construct full host geometry (layers + graph + modules)
    auto caGeometryHost = std::unique_ptr<reco::CAGeometryHost>(
      new reco::CAGeometryHost{{{nLayers_ + 1, nPairs_, nModules_}}, cms::alpakatools::host()});

    auto layers = caGeometryHost->template view<::reco::CALayersSoA>();
    auto graph = caGeometryHost->template view<::reco::CAGraphSoA>();
    auto modules = caGeometryHost->template view<::reco::CAModulesSoA>();

    // Fill layer-level data
    for (uint32_t i = 0; i < thetaCuts_.size(); ++i) {
      layers.layerStarts(i) = layerStarts_[i];  // placeholder
      layers.caThetaCut(i) = static_cast<float>(thetaCuts_[i]);
      layers.caDCACut(i) = static_cast<float>(dcaCuts_[i]);
    }

    layers.layerStarts(thetaCuts_.size()) = layerStarts_[thetaCuts_.size()];

#ifdef GPU_DEBUG
    std::cout << "[GPU_DEBUG] Filled " << nLayers_ << " layers with theta/DCACuts." << std::endl;
#endif


    // Fill graph-level data
    for (uint32_t i = 0; i < phiCuts_.size(); ++i) {
      graph.graph(i) = {{pairGraph_[2 * i], pairGraph_[2 * i + 1]}};
      graph.startingPair(i) =
          std::find(startingPairs_.begin(), startingPairs_.end(), i) != startingPairs_.end();
      graph.phiCuts(i) = static_cast<int16_t>(phiCuts_[i]);
      graph.minz(i) = static_cast<float>(minZ_[i]);
      graph.maxz(i) = static_cast<float>(maxZ_[i]);
      graph.maxr(i) = static_cast<float>(maxR_[i]);
    }

#ifdef GPU_DEBUG
    std::cout << "[GPU_DEBUG] Filled " << nPairs_ << " CA graph pairs." << std::endl;
#endif

    // --- Read CAModulesLayout from binary file ---

    using Rotation = SOARotation<float>;
    using Frame = SOAFrame<float>;

    auto in = edm::utils::openInputFile(data_);
    if (!in.is_open()) {
      throw std::runtime_error("CAGeometryHostESProducer: cannot open " + data_.string());
    }

    // read number of modules (and check consistency)
    int nModulesInFile = 0;
    in.read(reinterpret_cast<char*>(&nModulesInFile), sizeof(uint16_t));


    if (nModulesInFile < nModules_) {
      std::ostringstream msg;
      msg << "[CAGeometryHostESProducer ERROR] Module count mismatch when reading file:\n"
          << "  File: " << data_ << "\n"
          << "  Expected (TrackerTraits::numberOfModules) = " << nModules_ << "\n"
          << "  Found in file = " << nModulesInFile << " (too few!)\n";
    #ifdef GPU_DEBUG
      std::cerr << msg.str();
    #endif
      throw std::runtime_error(msg.str());
    }

    if (nModulesInFile > nModules_) {
      std::ostringstream msg;
      msg << "[CAGeometryHostESProducer WARNING] File contains more modules than expected.\n"
          << "  File: " << data_ << "\n"
          << "  Expected = " << nModules_ << ", Found = " << nModulesInFile << "\n"
          << "  Will load only the first " << nModules_ << " modules.\n";
    #ifdef GPU_DEBUG
      std::cerr << msg.str();
    #else
      // Print at least once even without GPU_DEBUG
      std::cerr << msg.str();
    #endif
      // Continue loading, but only up to nModules_
      nModulesInFile = nModules_;
    }

    #ifdef GPU_DEBUG
    std::cout << "[GPU_DEBUG] Loading " << nModulesInFile << " module frames from file..." << std::endl;
    #endif

    for (int i = 0; i < nModulesInFile; ++i) {
      Frame frame;
      in.read(reinterpret_cast<char*>(&frame), sizeof(frame));
      modules.detFrame(i) = frame;
    }

    // If file had extra entries, skip to end
    if (in.peek() != EOF) {
      in.ignore(std::numeric_limits<std::streamsize>::max());
    }

    in.close();

    #ifdef GPU_DEBUG
    std::cout << "[GPU_DEBUG] Finished reading " << nModulesInFile
              << " module frames from " << data_.filename() << std::endl;
    #endif

    /// TODO: allow for a queue to be here (porcoddue). And after, have the automatic mechamism for ESProducers (later).
    eventSetup.put(std::move(caGeometryHost));

#ifdef GPU_DEBUG
    std::cout << "[GPU_DEBUG] CAGeometryHost put into EventSetup successfully.\n";
#endif

  }

    // using CAGeometryHostESProducerPhase1 = CAGeometryHostESProducer<pixelTopology::Phase1>;

  /// FIXME: These are needed to make these plugins visible when building the plugins.txt list
  /// see: src/alpaka/Makefile:204. This is a workaround but it works for the moment.
  class CAGeometryHostESProducerPhase1 : public CAGeometryHostESProducer<pixelTopology::Phase1> {
  public:
    using CAGeometryHostESProducer<pixelTopology::Phase1>::CAGeometryHostESProducer;
  };

  class CAGeometryHostESProducerPhase1FromHits : public CAGeometryHostESProducer<pixelTopology::Phase1FromHits> {
  public:
    using CAGeometryHostESProducer<pixelTopology::Phase1FromHits>::CAGeometryHostESProducer;
  };

  class CAGeometryHostESProducerGenericUpgrade : public CAGeometryHostESProducer<pixelTopology::GenericUpgrade> {
  public:
    using CAGeometryHostESProducer<pixelTopology::GenericUpgrade>::CAGeometryHostESProducer;
  };

  class CAGeometryHostESProducerColliderMLPhase1 : public CAGeometryHostESProducer<pixelTopology::ColliderMLPhase1> {
  public:
    using CAGeometryHostESProducer<pixelTopology::ColliderMLPhase1>::CAGeometryHostESProducer;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

  // ---------- Explicit instantiation and registration ----------
  DEFINE_FWK_ALPAKA_EVENTSETUP_MODULE(CAGeometryHostESProducerPhase1);
  DEFINE_FWK_ALPAKA_EVENTSETUP_MODULE(CAGeometryHostESProducerPhase1FromHits);
  DEFINE_FWK_ALPAKA_EVENTSETUP_MODULE(CAGeometryHostESProducerGenericUpgrade);
  DEFINE_FWK_ALPAKA_EVENTSETUP_MODULE(CAGeometryHostESProducerColliderMLPhase1);