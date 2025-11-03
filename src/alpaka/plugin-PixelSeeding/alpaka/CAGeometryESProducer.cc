#include <memory>
#include <string>
#include <vector>
#include <algorithm>

#include "AlpakaCore/config.h"
#include "Framework/ESProducer.h"
#include "Framework/EventSetup.h"
#include "Framework/ESPluginFactory.h"

#include "AlpakaDataFormats/CAGeometryHost.h"
#include "AlpakaDataFormats/CAGeometrySoA.h"
#include "AlpakaDataFormats/alpaka/CAGeometrySoACollection.h"
#include "AlpakaDataFormats/SOARotation.h"
#include "Geometry/SimplePixelTopology.h"

#define GPU_DEBUG
namespace ALPAKA_ACCELERATOR_NAMESPACE {

  template <typename TrackerTraits>
  class CAGeometryHostESProducer : public edm::ESProducer {
  public:
    explicit CAGeometryHostESProducer(std::string const& datadir) : data_(datadir) {
      // Fill default parameters from TrackerTraits
      constexpr uint32_t nLayers = TrackerTraits::numberOfLayers;
      constexpr uint32_t nPairs = TrackerTraits::nPairsForQuadruplets;

      auto checkSize = [](auto const& vec, size_t expected, std::string const& name) {
      if (vec.size() != expected) {
        std::cerr << "[CAGeometryHostESProducer ERROR] Size mismatch in " << name << ":\n"
                  << "  expected = " << expected << ", actual = " << vec.size() << std::endl;
        std::abort();
      }
    };

#ifdef GPU_DEBUG
      std::cout << "[GPU_DEBUG] Constructing CAGeometryHostESProducer for "
                << TrackerTraits::nameModifier << "\n"
                << "  - Data directory: " << data_ << "\n"
                << "  - Layers: " << nLayers
                << ", Pairs: " << nPairs
                << ", Default start pairs: {0,1,2}" << std::endl;
#endif

      thetaCuts_.assign(TrackerTraits::thetaCuts,
                        TrackerTraits::thetaCuts + nLayers);
      dcaCuts_.assign(TrackerTraits::dcaCuts,
                      TrackerTraits::dcaCuts + nLayers);
      phiCuts_.assign(TrackerTraits::phicuts,
                      TrackerTraits::phicuts + nPairs);
      minZ_.assign(TrackerTraits::minz,
                   TrackerTraits::minz + nPairs);
      maxZ_.assign(TrackerTraits::maxz,
                   TrackerTraits::maxz + nPairs);
      maxR_.assign(TrackerTraits::maxr,
                   TrackerTraits::maxr + nPairs);
      pairGraph_.assign(TrackerTraits::layerPairs,
                        TrackerTraits::layerPairs + (nPairs * 2));
      startingPairs_ = {0u, 1u, 2u}; 

      checkSize(thetaCuts_, nLayers, "thetaCuts_");
      checkSize(dcaCuts_, nLayers, "dcaCuts_");
      checkSize(phiCuts_, nPairs, "phiCuts_");
      checkSize(minZ_, nPairs, "minZ_");
      checkSize(maxZ_, nPairs, "maxZ_");
      checkSize(maxR_, nPairs, "maxR_");
      checkSize(pairGraph_, nPairs * 2, "pairGraph_");

    }

    void produce(edm::EventSetup& eventSetup);

  private:

    std::filesystem::path data_;

    // Geometry parameter data members
    std::vector<float> thetaCuts_;
    std::vector<float> dcaCuts_;
    std::vector<int> phiCuts_;
    std::vector<double> minZ_;
    std::vector<double> maxZ_;
    std::vector<double> maxR_;
    std::vector<unsigned int> pairGraph_;
    std::vector<unsigned int> startingPairs_;
  };

  // ---------- Produce ----------
  template <typename TrackerTraits>
  void CAGeometryHostESProducer<TrackerTraits>::produce(edm::EventSetup& eventSetup) {

    auto filename =
        "CAGeometryHostModules" + std::string(TrackerTraits::nameModifier) + ".bin";
    auto filepath = data_ / filename;

    #ifdef GPU_DEBUG
    std::cout << "[GPU_DEBUG] Producing CAGeometryHost for " << TrackerTraits::nameModifier << "\n"
              << "  - Reading modules from: " << filepath << std::endl;
    #endif

    constexpr uint32_t nLayers  = TrackerTraits::numberOfLayers;
    constexpr uint32_t nPairs   = TrackerTraits::nPairsForQuadruplets;
    constexpr uint32_t nModules = TrackerTraits::numberOfModules;
        
    // Construct full host geometry (layers + graph + modules)
    auto caGeometryHost = std::unique_ptr<reco::CAGeometryHost>(
      new reco::CAGeometryHost{{{nLayers + 1, nPairs, nModules}}, cms::alpakatools::host()});

    auto layers = caGeometryHost->template view<reco::CALayersSoA>();
    auto graph = caGeometryHost->template view<reco::CAGraphSoA>();
    auto modules = caGeometryHost->template view<reco::CAModulesSoA>();

    // Fill layer-level data
    for (uint32_t i = 0; i < thetaCuts_.size(); ++i) {
      layers.layerStarts(i) = i;  // placeholder
      layers.caThetaCut(i) = static_cast<float>(thetaCuts_[i]);
      layers.caDCACut(i) = static_cast<float>(dcaCuts_[i]);
    }

#ifdef GPU_DEBUG
    std::cout << "[GPU_DEBUG] Filled " << nLayers << " layers with theta/DCACuts." << std::endl;
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
    std::cout << "[GPU_DEBUG] Filled " << nPairs << " CA graph pairs." << std::endl;
#endif

    // --- Read CAModulesLayout from binary file ---

    using Rotation = SOARotation<float>;
    using Frame = SOAFrame<float>;

    std::ifstream in(filepath, std::ios::binary);
    if (!in.is_open()) {
      throw std::runtime_error("CAGeometryHostESProducer: cannot open " + filepath.string());
    }

    // read number of modules (and check consistency)
    int nModulesInFile = 0;
    in.read(reinterpret_cast<char*>(&nModulesInFile), sizeof(int));

    const int nModulesExpected = static_cast<int>(TrackerTraits::numberOfModules);

    if (nModulesInFile < nModulesExpected) {
      std::ostringstream msg;
      msg << "[CAGeometryHostESProducer ERROR] Module count mismatch when reading file:\n"
          << "  File: " << filepath << "\n"
          << "  Expected (TrackerTraits::numberOfModules) = " << nModulesExpected << "\n"
          << "  Found in file = " << nModulesInFile << " (too few!)\n";
    #ifdef GPU_DEBUG
      std::cerr << msg.str();
    #endif
      throw std::runtime_error(msg.str());
    }

    if (nModulesInFile > nModulesExpected) {
      std::ostringstream msg;
      msg << "[CAGeometryHostESProducer WARNING] File contains more modules than expected.\n"
          << "  File: " << filepath << "\n"
          << "  Expected = " << nModulesExpected << ", Found = " << nModulesInFile << "\n"
          << "  Will load only the first " << nModulesExpected << " modules.\n";
    #ifdef GPU_DEBUG
      std::cerr << msg.str();
    #else
      // Print at least once even without GPU_DEBUG
      std::cerr << msg.str();
    #endif
      // Continue loading, but only up to nModulesExpected
      nModulesInFile = nModulesExpected;
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
              << " module frames from " << filepath.filename() << std::endl;
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


}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

  // ---------- Explicit instantiation and registration ----------
  DEFINE_FWK_ALPAKA_EVENTSETUP_MODULE(CAGeometryHostESProducerPhase1);