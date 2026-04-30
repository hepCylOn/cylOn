#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "AlpakaCore/Product.h"
#include "AlpakaCore/ScopedContext.h"
#include "AlpakaCore/config.h"
#include "AlpakaDataFormats/alpaka/SiPixelClustersSoACollection.h"
#include "AlpakaDataFormats/alpaka/SiPixelDigiErrorsSoACollection.h"
#include "AlpakaDataFormats/alpaka/SiPixelDigisSoACollection.h"
#include "CondFormats/alpaka/SiPixelGainCalibrationForHLTSoACollection.h"
#include "CondFormats/alpaka/SiPixelMappingSoACollection.h"
#include "CondFormats/alpaka/SiPixelMappingUtilities.h"
#include "CondFormats/SiPixelFedIds.h"
#include "DataFormats/FEDNumbering.h"
#include "DataFormats/FEDRawData.h"
#include "DataFormats/FEDRawDataCollection.h"
#include "DataFormats/PixelErrors.h"
#include "Framework/EDProducer.h"
#include "Framework/Event.h"
#include "Framework/EventSetup.h"
#include "Framework/PluginFactory.h"
#include "SiPixelRawToDigi/ErrorChecker.h"

#include "SiPixelRawToClusterKernel.h"

// #define GPU_DEBUG
// #define INPUT_DEBUG

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  template <typename TrackerTraits>
  class SiPixelRawToCluster : public edm::EDProducerExternalWork {
  public:
    explicit SiPixelRawToCluster(edm::ProductRegistry& reg, edm::Config const& cfg);
    ~SiPixelRawToCluster() override = default;

    // static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
    using Algo = pixelDetails::SiPixelRawToClusterKernel<TrackerTraits>;

  private:
    void acquire(const edm::Event& iEvent,
                 const edm::EventSetup& iSetup,
                 edm::WaitingTaskWithArenaHolder waitingTaskHolder) override;
    void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override;
    
    cms::alpakatools::ContextState<Queue> ctxState_;

    edm::EDGetTokenT<FEDRawDataCollection> rawGetToken_;
    edm::EDPutTokenT<cms::alpakatools::Product<Queue,SiPixelDigisSoACollection>> digiPutToken_;
    edm::EDPutTokenT<cms::alpakatools::Product<Queue,SiPixelDigiErrorsSoACollection>> digiErrorPutToken_;
    edm::EDPutTokenT<cms::alpakatools::Product<Queue,SiPixelClustersSoACollection>> clusterPutToken_;


    std::vector<unsigned int> fedIds_;

    Algo Algo_;
    PixelFormatterErrors errors_;

    const bool includeErrors_;
    const bool useQuality_;
    const bool verbose_;
    uint32_t nDigis_;
    const SiPixelClusterThresholds clusterThresholds_;
  };

  template <typename TrackerTraits>
  SiPixelRawToCluster<TrackerTraits>::SiPixelRawToCluster(edm::ProductRegistry& reg, edm::Config const& cfg)
      : rawGetToken_(reg.consumes<FEDRawDataCollection>()),
        digiPutToken_(reg.produces<cms::alpakatools::Product<Queue, SiPixelDigisSoACollection>>()),
        clusterPutToken_(reg.produces<cms::alpakatools::Product<Queue, SiPixelClustersSoACollection>>()),
        includeErrors_(cfg.value("includeErrors", true)),
        useQuality_(cfg.value("useQuality", true)),
        verbose_(cfg.value("verbose", false)),
        clusterThresholds_(
            cfg.value("clusterThreshold_layer1", 4000),
            cfg.value("clusterThreshold_otherLayers", 4000),
            static_cast<float>(cfg.value("VCaltoElectronGain", 47.0)),
            static_cast<float>(cfg.value("VCaltoElectronGain_L1", 50.0)),
            static_cast<float>(cfg.value("VCaltoElectronOffset", -60.0)),
            static_cast<float>(cfg.value("VCaltoElectronOffset_L1", -670.0))) {

    if (includeErrors_) {
      digiErrorPutToken_ = reg.produces<cms::alpakatools::Product<Queue, SiPixelDigiErrorsSoACollection>>();
    }

#ifdef INPUT_DEBUG
  std::cout << "[SiPixelRawToCluster] Configuration:" << std::endl;
  std::cout << "  includeErrors = " << includeErrors_ << std::endl;
  std::cout << "  useQuality    = " << useQuality_ << std::endl;
  std::cout << "  verbose       = " << verbose_ << std::endl;
  // std::cout << "  thresholds L1 = " << clusterThresholds_.clusterThreshold_layer1
  //           << ", other = " << clusterThresholds_.clusterThreshold_otherLayers << std::endl;
#endif
}

  // template <typename TrackerTraits>
  // void SiPixelRawToCluster<TrackerTraits>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //   edm::ParameterSetDescription desc;
  //   desc.add<bool>("IncludeErrors", true);
  //   desc.add<bool>("UseQualityInfo", false);
  //   desc.add<bool>("verbose", false)->setComment("verbose FED / ROC errors output");
  //   // Note: this parameter is obsolete: it is ignored and will have no effect.
  //   // It is kept to avoid breaking older configurations, and will not be printed in the generated cfi.py file.
  //   desc.addOptionalNode(edm::ParameterDescription<uint32_t>("MaxFEDWords", 0, true), false)
  //       ->setComment("This parameter is obsolete and will be ignored.");
  //   desc.add<int32_t>("clusterThreshold_layer1", pixelClustering::clusterThresholdLayerOne);
  //   desc.add<int32_t>("clusterThreshold_otherLayers", pixelClustering::clusterThresholdOtherLayers);
  //   desc.add<double>("VCaltoElectronGain", 47.f);
  //   desc.add<double>("VCaltoElectronGain_L1", 50.f);
  //   desc.add<double>("VCaltoElectronOffset", -60.f);
  //   desc.add<double>("VCaltoElectronOffset_L1", -670.f);

  //   desc.add<edm::InputTag>("InputLabel", edm::InputTag("rawDataCollector"));
  //   {
  //     edm::ParameterSetDescription psd0;
  //     psd0.addOptional<std::vector<edm::InputTag>>("inputs");
  //     psd0.addOptional<std::vector<double>>("deltaPhi");
  //     psd0.addOptional<std::vector<double>>("maxZ");
  //     psd0.addOptional<edm::InputTag>("beamSpot");
  //     desc.add<edm::ParameterSetDescription>("Regions", psd0)
  //         ->setComment("## Empty Regions PSet means complete unpacking");
  //   }
  //   desc.add<std::string>("CablingMapLabel", "")->setComment("CablingMap label");  //Tav
  //   descriptions.addWithDefaultLabel(desc);
  // }

  template <typename TrackerTraits>
  void SiPixelRawToCluster<TrackerTraits>::acquire(const edm::Event& iEvent,
                 const edm::EventSetup& iSetup,
                 edm::WaitingTaskWithArenaHolder waitingTaskHolder) {

#ifdef GPU_DEBUG
    std::cout << "[SiPixelRawToCluster::GPU_DEBUG] acquire started " << std::endl;
#endif
    cms::alpakatools::ScopedContextAcquire<Queue> ctx{iEvent.streamID(), std::move(waitingTaskHolder), ctxState_};

    auto const& hMap = iSetup.get<SiPixelMappingHost>();
    auto const& hGains = iSetup.get<SiPixelGainCalibrationForHLTHost>();

/// TODO: make a helper for this, have the automatic mechamism for copy (later).
#if defined ALPAKA_ACC_CPU_B_SEQ_T_SEQ_ENABLED or defined ALPAKA_ACC_CPU_B_TBB_T_SEQ_ENABLED
    SiPixelMappingSoACollection const& dMap = hMap;
    SiPixelGainCalibrationForHLTSoACollection const& dGains = hGains;
#else
    SiPixelMappingSoACollection dMap = cms::alpakatools::CopyToDevice<SiPixelMappingHost>::copyAsync(ctx.stream(), hMap);
    SiPixelGainCalibrationForHLTSoACollection dGains = cms::alpakatools::CopyToDevice<SiPixelGainCalibrationForHLTHost>::copyAsync(ctx.stream(), hGains);
    alpaka::wait(ctx.stream());
#endif


    // if (hMap.hasQuality() != useQuality_) {
    //   throw std::runtime_error("UseQuality of the module (" + std::to_string(useQuality_) +
    //                            ") differs the one from SiPixelFedCablingMapGPUWrapper. Please fix your configuration.");
    // }

    // // initialize cabling map or update if necessary
    // if (recordWatcher_.check(iSetup)) {
    //   // cabling map, which maps online address (fed->link->ROC->local pixel) to offline (DetId->global pixel)
    //   cablingMap_ = &iSetup.getData(cablingMapToken_);
    //   fedIds_ = cablingMap_->fedIds();
    //   cabling_ = cablingMap_->cablingTree();
    //   LogDebug("map version:") << cablingMap_->version();
    // }

    fedIds_ = iSetup.get<SiPixelFedIds>().fedIds();
    const unsigned char* modulesToUnpack = dMap->modToUnpDefault().data();

    const auto& buffers = iEvent.get(rawGetToken_);

    errors_.clear();

    // GPU specific: Data extraction for RawToDigi GPU
    unsigned int wordCounter = 0;
    unsigned int fedCounter = 0;
    bool errorsInEvent = false;
    std::vector<unsigned int> index(fedIds_.size(), 0);
    std::vector<uint32_t const*> start(fedIds_.size(), nullptr);
    std::vector<ptrdiff_t> words(fedIds_.size(), 0);
    // In CPU algorithm this loop is part of PixelDataFormatter::interpretRawData()
    ErrorChecker errorcheck;
    for (uint32_t i = 0; i < fedIds_.size(); ++i) {
      
      const int fedId = fedIds_[i];

      if (fedId == 40)
        continue;  // skip pilot blade data

      // for GPU
      // first 150 index stores the fedId and next 150 will store the
      // start index of word in that fed
      assert(fedId >= FEDNumbering::MINSiPixeluTCAFEDID);
      fedCounter++;

      // get event data for this fed
      const FEDRawData& rawData = buffers.FEDData(fedId);

      // GPU specific
      int nWords = rawData.size() / sizeof(uint64_t);
      if (nWords == 0) {
        continue;
      }
      // check CRC bit
      const uint64_t* trailer = reinterpret_cast<const uint64_t*>(rawData.data()) + (nWords - 1);
      if (not errorcheck.checkCRC(errorsInEvent, fedId, trailer, errors_)) {
        continue;
      }
      // check headers
      const uint64_t* header = reinterpret_cast<const uint64_t*>(rawData.data());
      header--;
      bool moreHeaders = true;
      while (moreHeaders) {
        header++;
        bool headerStatus = errorcheck.checkHeader(errorsInEvent, fedId, header, errors_);
        moreHeaders = headerStatus;
      }

      // check trailers
      bool moreTrailers = true;
      trailer++;
      while (moreTrailers) {
        trailer--;
        bool trailerStatus = errorcheck.checkTrailer(errorsInEvent, fedId, nWords, trailer, errors_);
        moreTrailers = trailerStatus;
      }

      const cms_uint32_t* bw = (const cms_uint32_t*)(header + 1);
      const cms_uint32_t* ew = (const cms_uint32_t*)(trailer);

      assert(0 == (ew - bw) % 2);
      index[i] = wordCounter;
      start[i] = bw;
      words[i] = (ew - bw);
      wordCounter += (ew - bw);

    }  // end of for loop
    nDigis_ = wordCounter;
    if (nDigis_ == 0)
      return;

#ifdef GPU_DEBUG
    std::cout << "[SiPixelRawToCluster::GPU_DEBUG] Finished FED unpacking. "
              << "FEDs processed: " << fedCounter
              << ", total digis: " << nDigis_ << "\n";
    if (fedCounter > 0) {
      std::cout << "[SiPixelRawToCluster::GPU_DEBUG] Example FED IDs: ";
      for (size_t j = 0; j < std::min<size_t>(fedIds_.size(), 5); ++j)
        std::cout << fedIds_[j] << " ";
      std::cout << "\n";
    }
#endif

    // copy the FED data to a single cpu buffer
    pixelDetails::WordFedAppender wordFedAppender(ctx.stream(), nDigis_);
    for (uint32_t i = 0; i < fedIds_.size(); ++i) {
      wordFedAppender.initializeWordFed(fedIds_[i], index[i], start[i], words[i]);
    }
    Algo_.makePhase1ClustersAsync(ctx.stream(),
                                  clusterThresholds_,
                                  dMap.const_view(),
                                  modulesToUnpack,
                                  dGains.const_view(),
                                  wordFedAppender,
                                  wordCounter,
                                  fedCounter,
                                  useQuality_,
                                  includeErrors_,
                                  verbose_);
  }

  template <typename TrackerTraits>
  void SiPixelRawToCluster<TrackerTraits>::produce(edm::Event& iEvent, edm::EventSetup const& iSetup) {
    cms::alpakatools::ScopedContextProduce ctx{ctxState_};
    // if (nDigis_ == 0) {
    //   // Cannot use the default constructor here, as it would not allocate memory.
    //   // In the case of no digis, clusters_d are not being instantiated, but are
    //   // still used downstream to initialize TrackingRecHitSoADevice. If there
    //   // are no valid pointers to clusters' Collection columns, instantiation
    //   // of TrackingRecHits fail. Example: workflow 11604.0

    //   ctx.emplace(iEvent,digiPutToken_, 0);
    //   ctx.emplace(iEvent,clusterPutToken_, pixelTopology::Phase1::numberOfModules);
    //   if (includeErrors_) {
    //     ctx.emplace(iEvent,digiErrorPutToken_, 0);
    //   }
    //   return;
    // }

    ctx.emplace(iEvent,digiPutToken_, std::move(Algo_.getDigis()));
    ctx.emplace(iEvent,clusterPutToken_, std::move(Algo_.getClusters()));
    if (includeErrors_) {
      ctx.emplace(iEvent,digiErrorPutToken_, std::move(Algo_.getErrors()));
    }
  }

  // using SiPixelRawToClusterPhase1 = SiPixelRawToCluster<pixelTopology::Phase1>;
  // using SiPixelRawToClusterHIonPhase1 = SiPixelRawToCluster<pixelTopology::HIonPhase1>;


  /// FIXME: These are needed to make these plugins visible when building the plugins.txt list
  /// see: src/alpaka/Makefile:204. This is a workaround but it works for the moment.
  class SiPixelRawToClusterPhase1 : public SiPixelRawToCluster<pixelTopology::Phase1> {
  public:
    using SiPixelRawToCluster<pixelTopology::Phase1>::SiPixelRawToCluster;
  };

  class SiPixelRawToClusterHIonPhase1 : public SiPixelRawToCluster<pixelTopology::HIonPhase1> {
  public:
    using SiPixelRawToCluster<pixelTopology::HIonPhase1>::SiPixelRawToCluster;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

// define as framework plugin
DEFINE_FWK_ALPAKA_MODULE(SiPixelRawToClusterPhase1);
DEFINE_FWK_ALPAKA_MODULE(SiPixelRawToClusterHIonPhase1);
