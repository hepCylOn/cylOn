#include <fstream>
#include <map>
#include <optional>
#include <string>
#include <utility>

#include <alpaka/alpaka.hpp>

#include <Eigen/Core>

#include "AlpakaCore/Product.h"
#include "AlpakaCore/ScopedContext.h"
#include "AlpakaCore/config.h"
#include "AlpakaCore/memory.h"
#include "AlpakaDataFormats/alpaka/SiPixelDigisSoACollection.h"
#include "AlpakaDataFormats/alpaka/SiPixelClustersSoACollection.h"
#include "AlpakaDataFormats/alpaka/TrackingRecHitsSoACollection.h"
#include "AlpakaDataFormats/alpaka/TracksSoACollection.h"
#include "AlpakaDataFormats/alpaka/ZVertexSoACollection.h"
#include "AlpakaDataFormats/gpuClusteringConstants.h"
#include "AlpakaDataFormats/TrackDefinitions.h"
#include "Framework/EDProducer.h"
#include "Framework/Event.h"
#include "Framework/EventSetup.h"
#include "Framework/PluginFactory.h"

#include "../SimpleAtomicHisto.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  class HistoValidator : public edm::EDProducerExternalWork {
  public:
    explicit HistoValidator(edm::ProductRegistry& reg, edm::Config const& cfg);

  private:
    void acquire(const edm::Event& iEvent,
                 const edm::EventSetup& iSetup,
                 edm::WaitingTaskWithArenaHolder waitingTaskHolder) override;
    void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override;
    void endJob() override;

    edm::EDGetTokenT<cms::alpakatools::Product<Queue, SiPixelDigisSoACollection>> digiToken_;
    edm::EDGetTokenT<cms::alpakatools::Product<Queue, SiPixelClustersSoACollection>> clusterToken_;
    edm::EDGetTokenT<cms::alpakatools::Product<Queue, reco::TrackingRecHitsSoACollection>> hitToken_;
    edm::EDGetTokenT<cms::alpakatools::Product<Queue, reco::TracksSoACollection>> trackToken_;
    edm::EDGetTokenT<cms::alpakatools::Product<Queue, ZVertexSoACollection>> vertexToken_;

    uint32_t nDigis_;
    uint32_t nModules_;
    uint32_t nClusters_;
    uint32_t nHits_;

    static std::map<std::string, SimpleAtomicHisto> histos;
  };

  std::map<std::string, SimpleAtomicHisto> HistoValidator::histos = {
      {"digi_n", SimpleAtomicHisto(100, 0, 1e5)},
      {"digi_adc", SimpleAtomicHisto(250, 0, 5e4)},
      {"module_n", SimpleAtomicHisto(100, 1500, 2000)},
      {"cluster_n", SimpleAtomicHisto(200, 5000, 25000)},
      {"cluster_per_module_n", SimpleAtomicHisto(110, 0, 110)},
      {"hit_n", SimpleAtomicHisto(200, 5000, 25000)},
      {"hit_lx", SimpleAtomicHisto(200, -1, 1)},
      {"hit_ly", SimpleAtomicHisto(800, -4, 4)},
      {"hit_lex", SimpleAtomicHisto(100, 0, 5e-5)},
      {"hit_ley", SimpleAtomicHisto(100, 0, 1e-4)},
      {"hit_gx", SimpleAtomicHisto(200, -20, 20)},
      {"hit_gy", SimpleAtomicHisto(200, -20, 20)},
      {"hit_gz", SimpleAtomicHisto(600, -60, 60)},
      {"hit_gr", SimpleAtomicHisto(200, 0, 20)},
      {"hit_charge", SimpleAtomicHisto(400, 0, 4e6)},
      {"hit_sizex", SimpleAtomicHisto(800, 0, 800)},
      {"hit_sizey", SimpleAtomicHisto(800, 0, 800)},
      {"track_n", SimpleAtomicHisto(150, 0, 15000)},
      {"track_nhits", SimpleAtomicHisto(3, 3, 6)},
      {"track_chi2", SimpleAtomicHisto(100, 0, 40)},
      {"track_pt", SimpleAtomicHisto(400, 0, 400)},
      {"track_eta", SimpleAtomicHisto(100, -3, 3)},
      {"track_phi", SimpleAtomicHisto(100, -3.15, 3.15)},
      {"track_tip", SimpleAtomicHisto(100, -1, 1)},
      {"track_tip_zoom", SimpleAtomicHisto(100, -0.05, 0.05)},
      {"track_zip", SimpleAtomicHisto(100, -15, 15)},
      {"track_zip_zoom", SimpleAtomicHisto(100, -0.1, 0.1)},
      {"track_quality", SimpleAtomicHisto(6, 0, 6)},
      {"vertex_n", SimpleAtomicHisto(60, 0, 60)},
      {"vertex_z", SimpleAtomicHisto(100, -15, 15)},
      {"vertex_chi2", SimpleAtomicHisto(100, 0, 40)},
      {"vertex_ndof", SimpleAtomicHisto(170, 0, 170)},
      {"vertex_pt2", SimpleAtomicHisto(100, 0, 4000)}};

  HistoValidator::HistoValidator(edm::ProductRegistry& reg, edm::Config const& cfg)
      : digiToken_{reg.consumes<cms::alpakatools::Product<Queue, SiPixelDigisSoACollection>>()},
        clusterToken_{reg.consumes<cms::alpakatools::Product<Queue, SiPixelClustersSoACollection>>()},
        hitToken_{reg.consumes<cms::alpakatools::Product<Queue, reco::TrackingRecHitsSoACollection>>()},
        trackToken_{reg.consumes<cms::alpakatools::Product<Queue, reco::TracksSoACollection>>()},
        vertexToken_{reg.consumes<cms::alpakatools::Product<Queue, ZVertexSoACollection>>()} {}

  void HistoValidator::acquire(const edm::Event& iEvent,
                               const edm::EventSetup& iSetup,
                               edm::WaitingTaskWithArenaHolder waitingTaskHolder) {
    auto const& pdigis = iEvent.get(digiToken_);
    cms::alpakatools::ScopedContextAcquire ctx{pdigis, std::move(waitingTaskHolder)};
    auto const& digis = ctx.get(pdigis);
    auto const& clusters = ctx.get(iEvent, clusterToken_);
    auto const& hits = ctx.get(iEvent, hitToken_);

    nDigis_ = digis.nDigis();
    nModules_ = digis.nModules();

    nClusters_ = clusters.nClusters();

    nHits_ = hits.nHits();

  }

  void HistoValidator::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    auto const& pdigis = iEvent.get(digiToken_);
    cms::alpakatools::ScopedContextProduce<Queue> ctx{pdigis};
    auto const& digis = ctx.get(iEvent, digiToken_);
#if defined ALPAKA_ACC_CPU_B_SEQ_T_SEQ_ENABLED or defined ALPAKA_ACC_CPU_B_TBB_T_SEQ_ENABLED
      SiPixelDigisHost const& digis_host = digis;
#else
      SiPixelDigisHost const& digis_host =
          cms::alpakatools::CopyToHost<SiPixelDigisDevice<Device>>::copyAsync(ctx.stream(), digis);
      alpaka::wait(ctx.stream());
#endif
    auto const& clusters = ctx.get(iEvent, clusterToken_);
#if defined ALPAKA_ACC_CPU_B_SEQ_T_SEQ_ENABLED or defined ALPAKA_ACC_CPU_B_TBB_T_SEQ_ENABLED
      SiPixelClustersHost const& clusters_host = clusters;
#else
      SiPixelClustersHost const& clusters_host =
          cms::alpakatools::CopyToHost<SiPixelClustersDevice<Device>>::copyAsync(ctx.stream(), clusters);
      alpaka::wait(ctx.stream());
#endif
    auto const& hits = ctx.get(iEvent, hitToken_);
#if defined ALPAKA_ACC_CPU_B_SEQ_T_SEQ_ENABLED or defined ALPAKA_ACC_CPU_B_TBB_T_SEQ_ENABLED
      ::reco::TrackingRecHitHost const& hits_host = hits;
#else
      ::reco::TrackingRecHitHost const& hits_host =
          cms::alpakatools::CopyToHost<::reco::TrackingRecHitDevice<Device>>::copyAsync(ctx.stream(), hits);
      alpaka::wait(ctx.stream());
#endif
    histos["module_n"].fill(nModules_);
    histos["digi_n"].fill(nDigis_);
    for (uint32_t i = 0; i < nDigis_; ++i) {
      histos["digi_adc"].fill(digis_host.view()[i].adc());
    }

    histos["cluster_n"].fill(nClusters_);
    for (uint32_t i = 0; i < nModules_; ++i) {
      histos["cluster_per_module_n"].fill(clusters_host.view()[i].clusInModule());
    }

    histos["hit_n"].fill(nHits_);
    for (uint32_t i = 0; i < nHits_; ++i) {
      histos["hit_lx"].fill(hits_host.view()[i].xLocal());
      histos["hit_ly"].fill(hits_host.view()[i].yLocal());
      histos["hit_lex"].fill(hits_host.view()[i].xerrLocal());
      histos["hit_ley"].fill(hits_host.view()[i].yerrLocal());
      histos["hit_gx"].fill(hits_host.view()[i].xGlobal());
      histos["hit_gy"].fill(hits_host.view()[i].yGlobal());
      histos["hit_gz"].fill(hits_host.view()[i].zGlobal());
      histos["hit_gr"].fill(hits_host.view()[i].rGlobal());
      histos["hit_sizex"].fill(hits_host.view()[i].clusterSizeX());
      histos["hit_sizey"].fill(hits_host.view()[i].clusterSizeY());
      histos["hit_charge"].fill(hits_host.view()[i].chargeAndStatus().charge);
    }

    {
      auto const& tracks = ctx.get(iEvent, trackToken_);
      int nTracks = 0;
#if defined ALPAKA_ACC_CPU_B_SEQ_T_SEQ_ENABLED or defined ALPAKA_ACC_CPU_B_TBB_T_SEQ_ENABLED
      ::reco::TracksHost const& tracks_host = tracks;
#else
      ::reco::TracksHost const& tracks_host =
          cms::alpakatools::CopyToHost<::reco::TracksDevice<Device>>::copyAsync(ctx.stream(), tracks);
      alpaka::wait(ctx.stream());
#endif
      for (int i = 0; i < tracks_host.view().metadata().size(); ++i) {
        if (nHits(tracks_host.view(),i) > 0 and tracks_host.view()[i].quality() >= ::pixelTrack::qualityByName("loose")) {
          ++nTracks;
          histos["track_nhits"].fill(nHits(tracks_host.view(),i));
          histos["track_chi2"].fill(tracks_host.view()[i].chi2());
          histos["track_pt"].fill(tracks_host.view()[i].pt());
          histos["track_eta"].fill(tracks_host.view()[i].eta());
          histos["track_phi"].fill(tracks_host.view()[i].state()[0]);
          histos["track_tip"].fill(tracks_host.view()[i].state()[1]);
          histos["track_tip_zoom"].fill(tracks_host.view()[i].state()[1]);
          histos["track_zip"].fill(tracks_host.view()[i].state()[4]);
          histos["track_zip_zoom"].fill(tracks_host.view()[i].state()[4]);
          histos["track_quality"].fill(static_cast<uint16_t>(tracks_host.view()[i].quality()));
        }
      }
      histos["track_n"].fill(nTracks);
    }

    {
      auto const& vertices = ctx.get(iEvent, vertexToken_);
#if defined ALPAKA_ACC_CPU_B_SEQ_T_SEQ_ENABLED or defined ALPAKA_ACC_CPU_B_TBB_T_SEQ_ENABLED
      ZVertexHost const& vertices_host = vertices;
#else
      ZVertexHost const& vertices_host =
          cms::alpakatools::CopyToHost<ZVertexDevice<Device>>::copyAsync(ctx.stream(), vertices);
      alpaka::wait(ctx.stream());
#endif
      histos["vertex_n"].fill(vertices_host.view().nvFinal());
      for (uint32_t i = 0; i < vertices_host.view().nvFinal(); ++i) {
        histos["vertex_z"].fill(vertices_host.view()[i].zv());
        histos["vertex_chi2"].fill(vertices_host.view()[i].chi2());
        histos["vertex_ndof"].fill(vertices_host.view<::reco::ZVertexTracksSoA>()[i].ndof());
        histos["vertex_pt2"].fill(vertices_host.view()[i].ptv2());
      }
    }
  }

  void HistoValidator::endJob() {
#if defined ALPAKA_ACC_CPU_B_SEQ_T_SEQ_SYNC_BACKEND
    std::ofstream out("histograms_alpaka_serial.txt");
#elif defined ALPAKA_ACC_CPU_B_TBB_T_SEQ_ASYNC_BACKEND
    std::ofstream out("histograms_alpaka_tbb.txt");
#elif defined ALPAKA_ACC_GPU_CUDA_ASYNC_BACKEND
    std::ofstream out("histograms_alpaka_cuda.txt");
#elif defined ALPAKA_ACC_GPU_HIP_ASYNC_BACKEND
    std::ofstream out("histograms_alpaka_hip.txt");
#else
#error "Support for a new Alpaka backend must be added here"
#endif
    for (auto const& elem : histos) {
      out << elem.first << " " << elem.second << "\n";
    }
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

DEFINE_FWK_ALPAKA_MODULE(HistoValidator);
