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
#include "AlpakaDataFormats/SOARotation.h"
#include "Geometry/SimplePixelTopology.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  template <typename TrackerTraits>
  class CAGeometryHostESProducer : public edm::ESProducer {
  public:
    explicit CAGeometryHostESProducer(std::string const& datadir) : data_(datadir) {
      // Fill default parameters from TrackerTraits
      constexpr uint32_t nLayers = TrackerTraits::numberOfLayers;
      constexpr uint32_t nPairs = TrackerTraits::nPairsForQuadruplets;

      // Safety runtime assertions
      assert(thetaCuts_.size() == nLayers);
      assert(dcaCuts_.size() == nLayers);
      assert(phiCuts_.size() == nPairs);
      assert(pairGraph_.size() == nPairs * 2);

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
      startingPairs_ = {0u, 1u, 2u};  // CMSSW default

      // --- consistency checks ---
      assert(thetaCuts_.size() == nLayers && "thetaCuts_ size mismatch");
      assert(dcaCuts_.size() == nLayers && "dcaCuts_ size mismatch");
      assert(phiCuts_.size() == nPairs && "phiCuts_ size mismatch");
      assert(minZ_.size()  == nPairs && "minZ_ size mismatch");
      assert(maxZ_.size()  == nPairs && "maxZ_ size mismatch");
      assert(maxR_.size()  == nPairs && "maxR_ size mismatch");
      assert(pairGraph_.size() == nPairs * 2 && "pairGraph_ size mismatch");

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

    // --- Read CAModulesLayout from binary file ---

    using Rotation = SOARotation<float>;
    using Frame = SOAFrame<float>;

    auto filename =
        "CAGeometryHostModules" + std::string(TrackerTraits::nameModifier) + ".bin";
    auto filepath = data_ / filename;

    std::ifstream in(filepath, std::ios::binary);
    if (!in.is_open()) {
      throw std::runtime_error("CAGeometryHostESProducer: cannot open " + filepath.string());
    }

    // read number of modules (and check consistency)
    int nModulesInFile = 0;
    in.read(reinterpret_cast<char*>(&nModulesInFile), sizeof(int));
    if (nModulesInFile != static_cast<int>(TrackerTraits::numberOfModules)) {
      throw std::runtime_error("CAGeometryHostESProducer: mismatch between TrackerTraits::numberOfModules and file");
    }

    // read all module frames
    for (int i = 0; i < nModulesInFile; ++i) {
      Frame frame;
      in.read(reinterpret_cast<char*>(&frame), sizeof(frame));
      modules.detFrame(i) = frame;
    }

    in.close();

    eventSetup.put(std::move(caGeometryHost));
  }

    using CAGeometryHostESProducerPhase1 = CAGeometryHostESProducer<pixelTopology::Phase1>;

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

  // ---------- Explicit instantiation and registration ----------
  DEFINE_FWK_ALPAKA_EVENTSETUP_MODULE(CAGeometryHostESProducerPhase1);