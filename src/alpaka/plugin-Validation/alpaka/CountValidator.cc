#include <atomic>
#include <cmath>
#include <iostream>
#include <mutex>
#include <sstream>
#include <stdexcept>

#include "AlpakaCore/Product.h"
#include "AlpakaCore/ScopedContext.h"
#include "AlpakaCore/config.h"

#include "AlpakaDataFormats/TracksHost.h"
#include "AlpakaDataFormats/ZVertexHost.h"

#include "AlpakaDataFormats/alpaka/SiPixelDigisSoACollection.h"
#include "AlpakaDataFormats/alpaka/SiPixelClustersSoACollection.h"


#include "DataFormats/DigiClusterCount.h"
#include "DataFormats/TrackCount.h"
#include "DataFormats/VertexCount.h"

#include "Framework/EDProducer.h"
#include "Framework/Event.h"
#include "Framework/EventSetup.h"
#include "Framework/PluginFactory.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  class CountValidator : public edm::EDProducer {
  public:
    explicit CountValidator(edm::ProductRegistry& reg);

  private:
    void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override;
    void endJob() override;

    // reference counts
    edm::EDGetTokenT<DigiClusterCount> digiClusterCountToken_;
    edm::EDGetTokenT<TrackCount> trackCountToken_;
    edm::EDGetTokenT<VertexCount> vertexCountToken_;

    // device data
    edm::EDGetTokenT<cms::alpakatools::Product<Queue, SiPixelDigisSoACollection>> digiToken_;
    edm::EDGetTokenT<cms::alpakatools::Product<Queue, SiPixelClustersSoACollection>> clusterToken_;

    // host data
    edm::EDGetTokenT<::reco::TracksHost> trackToken_;
    edm::EDGetTokenT<ZVertexHost> vertexToken_;

    // static counters
    static std::atomic<int> allEvents;
    static std::atomic<int> goodEvents;
    static std::atomic<int> sumVertexDifference;
    static std::mutex sumTrackDifferenceMutex;
    static float sumTrackDifference;
  };

  // static init
  std::atomic<int> CountValidator::allEvents{0};
  std::atomic<int> CountValidator::goodEvents{0};
  std::atomic<int> CountValidator::sumVertexDifference{0};
  std::mutex CountValidator::sumTrackDifferenceMutex;
  float CountValidator::sumTrackDifference = 0.f;

  CountValidator::CountValidator(edm::ProductRegistry& reg)
      : digiClusterCountToken_(reg.consumes<DigiClusterCount>()),
        trackCountToken_(reg.consumes<TrackCount>()),
        vertexCountToken_(reg.consumes<VertexCount>()),
        digiToken_(reg.consumes<cms::alpakatools::Product<Queue, SiPixelDigisSoACollection>>()),
        clusterToken_(reg.consumes<cms::alpakatools::Product<Queue, SiPixelClustersSoACollection>>()),
        trackToken_(reg.consumes<::reco::TracksHost>()),
        vertexToken_(reg.consumes<ZVertexHost>()) {}

  void CountValidator::produce(edm::Event& iEvent, const edm::EventSetup&) {
    constexpr float trackTolerance = 0.012f;
    constexpr int vertexTolerance = 1;
    std::stringstream ss;
    bool ok = true;

    ss << "Event " << iEvent.eventID() << " ";

    // ---- digis + clusters (device)
    auto const& pdigis = iEvent.get(digiToken_);
    cms::alpakatools::ScopedContextProduce<Queue> ctx{pdigis};
    {
      auto const& count = iEvent.get(digiClusterCountToken_);
      auto const& digis = ctx.get(iEvent, digiToken_);
      auto const& clusters = ctx.get(iEvent, clusterToken_);

      if (digis.nModules() != count.nModules()) {
        ss << "\n N(modules) is " << digis.nModules() << " expected " << count.nModules();
        ok = false;
      }
      if (digis.nDigis() != count.nDigis()) {
        ss << "\n N(digis) is " << digis.nDigis() << " expected " << count.nDigis();
        ok = false;
      }
      if (clusters.nClusters() != count.nClusters()) {
        ss << "\n N(clusters) is " << clusters.nClusters() << " expected " << count.nClusters();
        ok = false;
      }
    }

    // ---- tracks (host)
    {
      auto const& count = iEvent.get(trackCountToken_);
      auto const& tracks = iEvent.get(trackToken_);
      auto const view = tracks.view();
      int nTracks = 0; //view.nTracks();
      for (int i = 0; i < static_cast<int>(view.metadata().size()); ++i) {
        if (reco::nHits(view,i) > 0)
          ++nTracks;
      }
      assert(nTracks == view.nTracks());

      float rel = std::abs(float(nTracks - int(count.nTracks())) / count.nTracks());
      if (static_cast<unsigned int>(nTracks) != count.nTracks()) {
        std::lock_guard<std::mutex> guard(sumTrackDifferenceMutex);
        sumTrackDifference += rel;
      }
      if (rel >= trackTolerance) {
        ss << "\n N(tracks) is " << nTracks << " expected " << count.nTracks()
           << ", relative difference " << rel << " outside tolerance " << trackTolerance;
        ok = false;
      }
    }

    // ---- vertices (host)
    {
      auto const& count = iEvent.get(vertexCountToken_);
      auto const& vertices = iEvent.get(vertexToken_);
      int diff = std::abs(int(vertices.view().nvFinal()) - int(count.nVertices()));
      if (diff != 0)
        sumVertexDifference += diff;
      if (diff > vertexTolerance) {
        ss << "\n N(vertices) is " << vertices.view().nvFinal() << " expected " << count.nVertices()
           << ", difference " << diff << " outside tolerance " << vertexTolerance;
        ok = false;
      }
    }

    ++allEvents;
    if (ok) {
      ++goodEvents;
    } else {
      std::cout << ss.str() << std::endl;
    }
  }

  void CountValidator::endJob() {
    if (allEvents == goodEvents) {
      std::cout << "CountValidator: all " << allEvents << " events passed validation\n";
      if (sumTrackDifference != 0.f) {
        std::cout << " Average relative track difference "
                  << sumTrackDifference / allEvents.load() << " (all within tolerance)\n";
      }
      if (sumVertexDifference != 0) {
        std::cout << " Average absolute vertex difference "
                  << float(sumVertexDifference.load()) / allEvents.load()
                  << " (all within tolerance)\n";
      }
    } else {
      std::cout << "CountValidator: " << (allEvents - goodEvents)
                << " events failed validation (see details above)\n";
      throw std::runtime_error("CountValidator failed");
    }
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

DEFINE_FWK_ALPAKA_MODULE(CountValidator);

