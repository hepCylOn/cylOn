#include <atomic>
#include <cmath>
#include <iostream>
#include <mutex>
#include <sstream>
#include <stdexcept>

#include "AlpakaDataFormats/PixelTrackHost.h"
#include "AlpakaDataFormats/ZVertexHost.h"
#include "DataFormats/DigiClusterCount.h"
#include "DataFormats/TrackCount.h"
#include "DataFormats/VertexCount.h"
#include "Framework/EDProducer.h"
#include "Framework/Event.h"
#include "Framework/EventSetup.h"
#include "Framework/PluginFactory.h"

class CountValidatorFromHits : public edm::EDProducer {
public:
  explicit CountValidatorFromHits(edm::ProductRegistry& reg);

private:
  void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override;
  void endJob() override;

  edm::EDGetTokenT<TrackCount> trackCountToken_;
  edm::EDGetTokenT<VertexCount> vertexCountToken_;

  edm::EDGetTokenT<PixelTrackHost> trackToken_;
  edm::EDGetTokenT<ZVertexHost> vertexToken_;

  int allEvents;
  int goodEvents;
  int sumVertexDifference;

  float sumTrackDifference;
};

CountValidatorFromHits::CountValidatorFromHits(edm::ProductRegistry& reg)
      : trackCountToken_(reg.consumes<TrackCount>()),
      vertexCountToken_(reg.consumes<VertexCount>()),
      trackToken_(reg.consumes<PixelTrackHost>()),
      vertexToken_(reg.consumes<ZVertexHost>()) {
        allEvents = 0;
        goodEvents = 0;
        sumVertexDifference = 0;
        sumTrackDifference = 0;
      }

void CountValidatorFromHits::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  constexpr float trackTolerance = 0.012f;  // in 200 runs of 1k events all events are withing this tolerance
  constexpr int vertexTolerance = 1;
  std::stringstream ss;
  bool ok = true;

  ss << "Event " << iEvent.eventID() << " ";

  {
    auto const& count = iEvent.get(trackCountToken_);
    auto const& tracks = iEvent.get(trackToken_);
    int nTracks = 0;
    for (int i = 0; i < tracks->stride(); ++i) {
      if (tracks->nHits(i) > 0) {
        ++nTracks;
      }
    }

    auto rel = std::abs(float(nTracks - int(count.nTracks())) / count.nTracks());
    if (static_cast<unsigned int>(nTracks) != count.nTracks()) {
      sumTrackDifference += rel;
    }
    if (rel >= trackTolerance) {
      ss << "\n N(tracks) is " << nTracks << " expected " << count.nTracks() << ", relative difference " << rel
         << " is outside tolerance " << trackTolerance;
      ok = false;
    }
  }

  {
    auto const& count = iEvent.get(vertexCountToken_);
    auto const& vertices = iEvent.get(vertexToken_);
    auto diff = std::abs(int(vertices->nvFinal) - int(count.nVertices()));
    if (diff != 0) {
      sumVertexDifference += diff;
    }
    if (diff > vertexTolerance) {
      ss << "\n N(vertices) is " << vertices->nvFinal << " expected " << count.nVertices() << ", difference " << diff
         << " is outside tolerance " << vertexTolerance;
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

void CountValidatorFromHits::endJob() {
  if (allEvents == goodEvents) {
    std::cout << "CountValidatorFromHits: all " << allEvents << " events passed validation\n";
    if (sumTrackDifference != 0.f) {
      std::cout << " Average relative track difference " << sumTrackDifference / allEvents
                << " (all within tolerance)\n";
    }
    if (sumVertexDifference != 0) {
      std::cout << " Average absolute vertex difference " << float(sumVertexDifference) / allEvents
                << " (all within tolerance)\n";
    }
  } else {
    std::cout << "CountValidatorFromHits: " << (allEvents - goodEvents) << " events failed validation (see details above)\n";
    throw std::runtime_error("CountValidatorFromHits failed");
  }
}


DEFINE_FWK_MODULE(CountValidatorFromHits);
