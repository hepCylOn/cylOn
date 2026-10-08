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
      data_(static_cast<std::string>(cfg.value("data", defaultPath_)))
      {

  auto getRequiredVector = [&]<typename T>(std::string const& key, size_t expected) {
    if (!cfg.contains(key)) {
      throw std::runtime_error("[CAGeometryHostESProducer - getRequiredVector ERROR] Missing required config key '" + key + "'");
    }

    std::vector<T> vec;
    try {
      vec = cfg.at(key).get<std::vector<T>>();
    } catch (std::exception const& e) {
      throw std::runtime_error("[CAGeometryHostESProducer - getRequiredVector ERROR] Invalid config key '" + key + "': " + e.what());
    }

    if (vec.size() != expected) {
      throw std::runtime_error("[CAGeometryHostESProducer - getRequiredVector ERROR] Size mismatch in config key '" + key +
                               "': expected " + std::to_string(expected) +
                               ", got " + std::to_string(vec.size()));
    }
    return vec;
  };

  auto getRequiredInt = [&](std::string const& key) {
    if (!cfg.contains(key)) {
      throw std::runtime_error("[CAGeometryHostESProducer - getRequiredInt ERROR] Missing required config key '" + key + "'");
    }

    int value;
    try {
      value = cfg.at(key).get<int>();
    } catch (std::exception const& e) {
      throw std::runtime_error("[CAGeometryHostESProducer - getRequiredInt ERROR] Invalid config key '" + key + "': " + e.what());
    }

    return value;
  };

  nLayers_ = getRequiredInt("nLayers");
  nPairs_ = getRequiredInt("nPairs");
  nModules_ = getRequiredInt("nModules");

  thetaCuts_ = getRequiredVector.template operator()<float>("thetaCuts", nLayers_);
  dcaCuts_ = getRequiredVector.template operator()<float>("dcaCuts", nLayers_);
  layerStarts_ = getRequiredVector.template operator()<unsigned int>("layerStarts", nLayers_ + 1);
  phiCuts_ = getRequiredVector.template operator()<short int>("phiCuts", nPairs_);
  minZ_ = getRequiredVector.template operator()<float>("minZ", nPairs_);
  maxZ_ = getRequiredVector.template operator()<float>("maxZ", nPairs_);
  maxR_ = getRequiredVector.template operator()<float>("maxR", nPairs_);
  pairGraph_ = getRequiredVector.template operator()<uint8_t>("pairGraph", nPairs_ * 2);
  startingPairs_ = getRequiredVector.template operator()<uint8_t>("startingPairs", nPairs_);

  for (auto value : startingPairs_) {
    if (value > 1) {
      throw std::runtime_error("[CAGeometryHostESProducer ERROR] startingPairs entries must be 0 or 1");
    }
  }

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
            // << "  - Data directory: " << data_ << "\n"
            << "  - Layers: " << nLayers_
            << ", Pairs: " << nPairs_
            << ", Starting pairs: " << vecToString(startingPairs_) << std::endl;
#endif
}

    void produce(edm::EventSetup& eventSetup);

  private:

    std::filesystem::path defaultPath_ = std::string("data/CAGeometryHostModules") + std::string(TrackerTraits::nameModifier) + ".bin";

    std::filesystem::path data_;
    int nLayers_, nPairs_, nModules_; //int(s) just because the PortableCollection wants so
    
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
    // std::cout << "[GPU_DEBUG] Producing CAGeometryHost for " << TrackerTraits::nameModifier << "\n"
    //           << "  - Reading modules from: " << data_ << std::endl;
    std::cout << "[GPU_DEBUG] Producing CAGeometryHost for " << TrackerTraits::nameModifier << std::endl;
    #endif

    // Construct full host geometry (layers + graph + modules)
    auto caGeometryHost = std::unique_ptr<reco::CAGeometryHost>(
      new reco::CAGeometryHost{{{nLayers_ + 1, nPairs_, nModules_}}, cms::alpakatools::host()});

    auto layers = caGeometryHost->template view<::reco::CALayersSoA>();
    auto graph = caGeometryHost->template view<::reco::CAGraphSoA>();
    // auto modules = caGeometryHost->template view<::reco::CAModulesSoA>();

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
      graph.startingPair(i) = startingPairs_[i] != 0;
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

    // auto in = edm::utils::openInputFile(data_);
    // if (!in.is_open()) {
    //   throw std::runtime_error("CAGeometryHostESProducer: cannot open " + data_.string());
    // }

    // // read number of modules (and check consistency)
    // int nModulesInFile = 0;
    // in.read(reinterpret_cast<char*>(&nModulesInFile), sizeof(uint16_t));


    // if (nModulesInFile < nModules_) {
    //   std::ostringstream msg;
    //   msg << "[CAGeometryHostESProducer ERROR] Module count mismatch when reading file:\n"
    //       << "  File: " << data_ << "\n"
    //       << "  Expected (TrackerTraits::numberOfModules) = " << nModules_ << "\n"
    //       << "  Found in file = " << nModulesInFile << " (too few!)\n";
    // #ifdef GPU_DEBUG
    //   std::cerr << msg.str();
    // #endif
    //   throw std::runtime_error(msg.str());
    // }

    // if (nModulesInFile > nModules_) {
    //   std::ostringstream msg;
    //   msg << "[CAGeometryHostESProducer WARNING] File contains more modules than expected.\n"
    //       << "  File: " << data_ << "\n"
    //       << "  Expected = " << nModules_ << ", Found = " << nModulesInFile << "\n"
    //       << "  Will load only the first " << nModules_ << " modules.\n";
    // #ifdef GPU_DEBUG
    //   std::cerr << msg.str();
    // #else
    //   // Print at least once even without GPU_DEBUG
    //   std::cerr << msg.str();
    // #endif
    //   // Continue loading, but only up to nModules_
    //   nModulesInFile = nModules_;
    // }

    // #ifdef GPU_DEBUG
    // std::cout << "[GPU_DEBUG] Loading " << nModulesInFile << " module frames from file..." << std::endl;
    // #endif

    // for (int i = 0; i < nModulesInFile; ++i) {
    //   Frame frame;
    //   in.read(reinterpret_cast<char*>(&frame), sizeof(frame));
    //   modules.detFrame(i) = frame;
    // }

    // // If file had extra entries, skip to end
    // if (in.peek() != EOF) {
    //   in.ignore(std::numeric_limits<std::streamsize>::max());
    // }

    // in.close();

    // #ifdef GPU_DEBUG
    // std::cout << "[GPU_DEBUG] Finished reading " << nModulesInFile
    //           << " module frames from " << data_.filename() << std::endl;
    // #endif

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

  class CAGeometryHostESProducerPhase2 : public CAGeometryHostESProducer<pixelTopology::Phase2> {
  public:
    using CAGeometryHostESProducer<pixelTopology::Phase2>::CAGeometryHostESProducer;
  };

  class CAGeometryHostESProducerGenericUpgrade : public CAGeometryHostESProducer<pixelTopology::GenericUpgrade> {
  public:
    using CAGeometryHostESProducer<pixelTopology::GenericUpgrade>::CAGeometryHostESProducer;
  };

  class CAGeometryHostESProducerColliderMLPhase1 : public CAGeometryHostESProducer<pixelTopology::ColliderMLPhase1> {
  public:
    using CAGeometryHostESProducer<pixelTopology::ColliderMLPhase1>::CAGeometryHostESProducer;
  };

  class CAGeometryHostESProducerColliderMLPhase2 : public CAGeometryHostESProducer<pixelTopology::ColliderMLPhase2> {
  public:
    using CAGeometryHostESProducer<pixelTopology::ColliderMLPhase2>::CAGeometryHostESProducer;
  };

  class CAGeometryHostESProducerColliderMLPixelPlusShortStripsPhase2 : public CAGeometryHostESProducer<pixelTopology::ColliderMLPixelPlusShortStripsPhase2> {
  public:
    using CAGeometryHostESProducer<pixelTopology::ColliderMLPixelPlusShortStripsPhase2>::CAGeometryHostESProducer;
  };

  class CAGeometryHostESProducerColliderMLAllTrackerPhase2 : public CAGeometryHostESProducer<pixelTopology::ColliderMLAllTrackerPhase2> {
  public:
    using CAGeometryHostESProducer<pixelTopology::ColliderMLAllTrackerPhase2>::CAGeometryHostESProducer;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

  // ---------- Explicit instantiation and registration ----------
  DEFINE_FWK_ALPAKA_EVENTSETUP_MODULE(CAGeometryHostESProducerPhase1);
  DEFINE_FWK_ALPAKA_EVENTSETUP_MODULE(CAGeometryHostESProducerPhase2);
  DEFINE_FWK_ALPAKA_EVENTSETUP_MODULE(CAGeometryHostESProducerGenericUpgrade);
  DEFINE_FWK_ALPAKA_EVENTSETUP_MODULE(CAGeometryHostESProducerColliderMLPhase1);
  DEFINE_FWK_ALPAKA_EVENTSETUP_MODULE(CAGeometryHostESProducerColliderMLPhase2);
  DEFINE_FWK_ALPAKA_EVENTSETUP_MODULE(CAGeometryHostESProducerColliderMLPixelPlusShortStripsPhase2);
  DEFINE_FWK_ALPAKA_EVENTSETUP_MODULE(CAGeometryHostESProducerColliderMLAllTrackerPhase2);