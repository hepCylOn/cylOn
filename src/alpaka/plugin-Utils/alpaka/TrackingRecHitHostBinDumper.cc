#include <fstream>
#include <filesystem>
#include <iostream>
#include <memory>
#include <vector>
#include <cstring>

#include "AlpakaCore/Product.h"
#include "AlpakaCore/ScopedContext.h"
#include "AlpakaCore/config.h"

#include "AlpakaDataFormats/TrackingRecHitsHost.h"
#include "AlpakaDataFormats/alpaka/TrackingRecHitsSoACollection.h"
#include "AlpakaDataFormats/TrackingRecHitsDevice.h"

#include "Framework/EDProducer.h"
#include "Framework/Event.h"
#include "Framework/PluginFactory.h"
#include "Framework/EventSetup.h"

// #define GPU_DEBUG

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  class TrackingRecHitHostBinDumper : public edm::EDProducer {
  public:
    explicit TrackingRecHitHostBinDumper(edm::ProductRegistry& reg, edm::Config const& cfg);
    void produce(edm::Event& iEvent, edm::EventSetup const& iSetup) override;
    using HitsOnDevice = reco::TrackingRecHitsSoACollection;
    using HitsOnHost = ::reco::TrackingRecHitHost;
    void endJob() override {

        if(!headerWritten_)
            return;

        std::fstream out(outfile_, std::ios::in | std::ios::out | std::ios::binary);
        if (!out)
            return;

        // skip magic + version + endianness
        out.seekp(12, std::ios::beg);
        out.write(reinterpret_cast<const char*>(&nEventsWritten_), sizeof(nEventsWritten_));
#ifdef GPU_DEBUG
        std::cout << "[BinDumper] Updated event count in header: "
                    << nEventsWritten_ << std::endl;
#endif
    }

  private:
    edm::EDGetTokenT<cms::alpakatools::Product<Queue, HitsOnDevice>> tokenHit_;
    std::string outfile_;
    /// TODO: define these somewhere common. It's a mess like this.
    static constexpr uint32_t kVersion = 1;
    static constexpr uint32_t kEndianness = 0x01020304;
    static constexpr char kMagic[4] = {'T', 'R', 'H', '1'};
    bool wroteHeader_ = false;
    static uint32_t nEventsWritten_;
    static bool headerWritten_;
  };

  TrackingRecHitHostBinDumper::TrackingRecHitHostBinDumper(edm::ProductRegistry& reg, edm::Config const& cfg)
      :
      tokenHit_{reg.consumes<cms::alpakatools::Product<Queue, HitsOnDevice>>()}
      {
    outfile_ = "data/hits_from_framework.bin";
    if (std::filesystem::exists(outfile_)) {
    throw std::runtime_error(
        "[TrackingRecHitHostBinDumper] Output file already exists: " + outfile_ +
        "\nI won't remove it or append anything to it. Please delete it before rerunning.");
}
#ifdef GPU_DEBUG
    std::cout << "TrackingRecHitHostBinDumper initialized, output=" << outfile_ << std::endl;
#endif
  }

  void TrackingRecHitHostBinDumper::produce(edm::Event& iEvent, edm::EventSetup const& iSetup) {

    auto const& phits = iEvent.get(tokenHit_);
    cms::alpakatools::ScopedContextProduce<Queue> ctx{phits};
    auto const& dHits = ctx.get(phits);

#if defined ALPAKA_ACC_CPU_B_SEQ_T_SEQ_ENABLED or defined ALPAKA_ACC_CPU_B_TBB_T_SEQ_ENABLED
      ::reco::TrackingRecHitHost const& hits = dHits;
#else
      ::reco::TrackingRecHitHost hits =
          cms::alpakatools::CopyToHost<::reco::TrackingRecHitDevice<Device>>::copyAsync(ctx.stream(), dHits);
      alpaka::wait(ctx.stream());
#endif

    auto view = hits.view<::reco::TrackingRecHitSoA>();
    auto modView = hits.view<::reco::HitModuleSoA>();

    uint32_t nHits = view.metadata().size();
    uint32_t nModules = hits.nModules();

#ifdef GPU_DEBUG
    std::cout << "[BinDumper] Dumping event with " << nHits
              << " hits, " << nModules << " modules\n";
#endif

    std::ofstream out(outfile_, std::ios::binary | std::ios::app);
    if (!out)
      throw std::runtime_error("Cannot open binary output file: " + outfile_);

    if (!wroteHeader_) {
      out.write(kMagic, 4);
      out.write(reinterpret_cast<const char*>(&kVersion), sizeof(kVersion));
      out.write(reinterpret_cast<const char*>(&kEndianness), sizeof(kEndianness));

      // placeholder for number of events (we can't know yet)
      uint32_t dummyEvents = 0;
      out.write(reinterpret_cast<const char*>(&dummyEvents), sizeof(dummyEvents));

      wroteHeader_ = true;

      headerWritten_ |= true; 

#ifdef GPU_DEBUG
      std::cout << "[BinDumper] Wrote header (TRH1 v1)\n";
#endif
    }

    // Write event block
    out.write(reinterpret_cast<const char*>(&nHits), sizeof(nHits));
    out.write(reinterpret_cast<const char*>(&nModules), sizeof(nModules));

    // Write moduleStart (size nModules+1)
    out.write(reinterpret_cast<const char*>(modView.moduleStart().data()),
              (nModules + 1) * sizeof(uint32_t));

    // Write each column (column-major)
    auto write_column = [&](auto const& span) {
      using ElemT = typename std::remove_reference_t<decltype(span)>::value_type;
      out.write(reinterpret_cast<const char*>(span.data()), nHits * sizeof(ElemT));
    };

    write_column(view.xLocal());
    write_column(view.yLocal());
    write_column(view.xerrLocal());
    write_column(view.yerrLocal());
    write_column(view.xGlobal());
    write_column(view.yGlobal());
    write_column(view.zGlobal());
    write_column(view.rGlobal());
    write_column(view.iphi());
    write_column(view.chargeAndStatus());
    write_column(view.clusterSizeX());
    write_column(view.clusterSizeY());
    write_column(view.detectorIndex());

#ifdef GPU_DEBUG
    std::cout << "[BinDumper] Event " << nEventsWritten_
              << " written with " << nHits << " hits\n";
#endif

    ++nEventsWritten_;
    out.close();
  }

  bool TrackingRecHitHostBinDumper::headerWritten_ = false;
  uint32_t TrackingRecHitHostBinDumper::nEventsWritten_ = 0;

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

DEFINE_FWK_ALPAKA_MODULE(TrackingRecHitHostBinDumper);
