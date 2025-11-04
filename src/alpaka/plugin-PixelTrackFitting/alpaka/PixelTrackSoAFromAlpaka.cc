#include <utility>
#include <alpaka/alpaka.hpp>

#include "AlpakaCore/ScopedContext.h"
#include "AlpakaCore/config.h"
#include "AlpakaCore/memory.h"
#include "AlpakaDataFormats/alpaka/TracksSoACollection.h"
#include "AlpakaDataFormats/TracksHost.h"
#include "AlpakaDataFormats/TracksSoA.h"
#include "Framework/EDProducer.h"
#include "Framework/Event.h"
#include "Framework/EventSetup.h"
#include "Framework/PluginFactory.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  // Producer to copy pixel tracks from GPU (device) memory to CPU (host)
  class PixelTrackSoAFromAlpaka : public edm::EDProducer {
  public:
    using TkSoADevice = reco::TracksSoACollection;  // device-side Tracks container
    using TkSoAHost   = ::reco::TracksHost;         // host-side Tracks container

    explicit PixelTrackSoAFromAlpaka(edm::ProductRegistry& reg);
    ~PixelTrackSoAFromAlpaka() override = default;

  private:
    void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override;

    // Tokens for reading GPU product and writing CPU product
    edm::EDGetTokenT<cms::alpakatools::Product<Queue, TkSoADevice>> tokenDevice_;
    edm::EDPutTokenT<TkSoAHost> tokenHost_;
  };

  PixelTrackSoAFromAlpaka::PixelTrackSoAFromAlpaka(edm::ProductRegistry& reg)
      : tokenDevice_(reg.consumes<cms::alpakatools::Product<Queue, TkSoADevice>>()),
        tokenHost_(reg.produces<TkSoAHost>()) {}

  void PixelTrackSoAFromAlpaka::produce(edm::Event& iEvent, edm::EventSetup const& iSetup) {

    // --- 1. Retrieve device-side data product
    cms::alpakatools::Product<Queue, TkSoADevice> const& inputDataWrapped = iEvent.get(tokenDevice_);

    // --- 2. Create a ScopedContext to ensure correct synchronization with GPU stream
    cms::alpakatools::ScopedContextProduce<Queue> ctx{inputDataWrapped};

    // --- 3. Get the device data reference
    auto const& inputData = ctx.get(inputDataWrapped);

#ifdef GPU_DEBUG
    std::cout << "[PixelTrackSoAFromAlpaka::GPU_DEBUG] Starting copy from device to host..." << std::endl;
#endif

    // --- 4. Extract individual SoA views (tracks + hits)
    auto tracks     = inputData.view<::reco::TrackSoA>();
    auto tracksHits = inputData.view<::reco::TrackHitSoA>();

#ifdef GPU_DEBUG
    std::cout << "  Number of tracks SoA entries: " << tracks.metadata().size() << std::endl;
    std::cout << "  Number of track hits SoA entries: " << tracksHits.metadata().size() << std::endl;
#endif

    auto hostSoA = std::unique_ptr<TkSoAHost>(
        new TkSoAHost{{{tracks.metadata().size(), tracksHits.metadata().size()}}, cms::alpakatools::host()});

#ifdef GPU_DEBUG
    std::cout << "  Host buffer allocated, starting async memcpy..." << std::endl;
#endif

    // --- 6. Copy device buffer -> host buffer
    alpaka::memcpy(ctx.stream(), hostSoA->buffer(), inputData.buffer());
    alpaka::wait(ctx.stream());

#ifdef GPU_DEBUG
    std::cout << "  Copy completed. Moving product into event." << std::endl;
#endif

    // --- 7. Move host data product into the event
    iEvent.emplace(tokenHost_, std::move(*hostSoA));

#ifdef GPU_DEBUG
    std::cout << "[PixelTrackSoAFromAlpaka::GPU_DEBUG] Finished successfully.\n" << std::endl;
#endif
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

DEFINE_FWK_ALPAKA_MODULE(PixelTrackSoAFromAlpaka);
