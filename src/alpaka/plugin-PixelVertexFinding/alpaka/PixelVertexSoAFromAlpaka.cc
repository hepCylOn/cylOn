#include <utility>
#include <alpaka/alpaka.hpp>

#include "AlpakaCore/ScopedContext.h"
#include "AlpakaCore/config.h"
#include "AlpakaCore/memory.h"
#include "AlpakaDataFormats/alpaka/ZVertexSoACollection.h"
#include "AlpakaDataFormats/ZVertexHost.h"
#include "Framework/EDProducer.h"
#include "Framework/Event.h"
#include "Framework/EventSetup.h"
#include "Framework/PluginFactory.h"

#define GPU_DEBUG

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  class PixelVertexSoAFromAlpaka : public edm::EDProducer {
  public:
    using VertexDevice = ZVertexSoACollection;  // Device-side vertex data
    using VertexHost   = ZVertexHost;           // Host-side vertex data

    explicit PixelVertexSoAFromAlpaka(edm::ProductRegistry& reg, edm::Config const& cfg);
    ~PixelVertexSoAFromAlpaka() override = default;

  private:
    void produce(edm::Event& iEvent, edm::EventSetup const& iSetup) override;

    edm::EDGetTokenT<cms::alpakatools::Product<Queue, VertexDevice>> tokenDevice_;
    edm::EDPutTokenT<VertexHost> tokenHost_;
  };

  PixelVertexSoAFromAlpaka::PixelVertexSoAFromAlpaka(edm::ProductRegistry& reg, edm::Config const& cfg)
      : tokenDevice_(reg.consumes<cms::alpakatools::Product<Queue, VertexDevice>>()),
        tokenHost_(reg.produces<VertexHost>()) {}

  void PixelVertexSoAFromAlpaka::produce(edm::Event& iEvent, edm::EventSetup const& iSetup) {

    // --- 1. Retrieve vertex data from GPU
    cms::alpakatools::Product<Queue, VertexDevice> const& inputDataWrapped = iEvent.get(tokenDevice_);

    // --- 2. Scoped context for proper stream sync
    cms::alpakatools::ScopedContextProduce<Queue> ctx{inputDataWrapped};
    auto const& inputData = ctx.get(inputDataWrapped);

#ifdef GPU_DEBUG
    std::cout << "[PixelVertexSoAFromAlpaka::GPU_DEBUG] Starting vertex copy from device to host..." << std::endl;
#endif

    // --- 4. Extract individual SoA views (tracks + hits)
    auto vertex     = inputData.view<::reco::ZVertexSoA>();
    auto vertexTr   = inputData.view<::reco::ZVertexTracksSoA>();

    auto hostSoA = std::unique_ptr<VertexHost>(
        new VertexHost{{{vertex.metadata().size(), vertexTr.metadata().size()}}, cms::alpakatools::host()});

#ifdef GPU_DEBUG
    std::cout << "  Host vertex buffer allocated, performing async memcpy..." << std::endl;
#endif

    // --- 4. Perform copy
    alpaka::memcpy(ctx.stream(), hostSoA->buffer(), inputData.buffer());
    alpaka::wait(ctx.stream());

#ifdef GPU_DEBUG
    std::cout << "  Vertex copy complete. Moving product into event." << std::endl;
#endif

    // --- 5. Move to event
    iEvent.emplace(tokenHost_, std::move(*hostSoA));

#ifdef GPU_DEBUG
    std::cout << "[PixelVertexSoAFromAlpaka::GPU_DEBUG] Finished successfully.\n" << std::endl;
#endif
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

DEFINE_FWK_ALPAKA_MODULE(PixelVertexSoAFromAlpaka);
