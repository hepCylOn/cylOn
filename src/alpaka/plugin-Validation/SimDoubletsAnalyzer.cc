// -*- C++ -*-
//
// Package:    Validation/TrackingMCTruth
// Class:      SimDoubletsAnalyzer
//

// user include files
#include "SimDoubletsAnalyzer.h"

#include "DataFormats/SimDoublets.h"
#include "DataFormats/ParticleSimpleSoA.h"
#include "DataFormats/TrackingRecHitSimpleSoA.h"
#include "Framework/EDProducer.h"
#include "Framework/Event.h"
#include "Framework/EventSetup.h"
#include "Framework/PluginFactory.h"

#include <cstddef>

namespace simdoublets {
  template <int DEGREE>
constexpr float approx_atan2s_P(float x);

// degree =  3   => absolute accuracy is  53
template <>
constexpr float approx_atan2s_P<3>(float x) {
  auto z = x * x;
  return x * ((-10142.439453125f) + z * 2002.0908203125f);
}
// degree =  5   => absolute accuracy is  7
template <>
constexpr float approx_atan2s_P<5>(float x) {
  auto z = x * x;
  return x * ((-10381.9609375f) + z * ((3011.1513671875f) + z * (-827.538330078125f)));
}
// degree =  7   => absolute accuracy is  2
template <>
constexpr float approx_atan2s_P<7>(float x) {
  auto z = x * x;
  return x * ((-10422.177734375f) + z * (3349.97412109375f + z * ((-1525.589599609375f) + z * 406.64190673828125f)));
}
// degree =  9   => absolute accuracy is 1
template <>
constexpr float approx_atan2s_P<9>(float x) {
  auto z = x * x;
  return x * ((-10428.984375f) + z * (3445.20654296875f + z * ((-1879.137939453125f) +
                                                               z * (888.22314453125f + z * (-217.42669677734375f)))));
}

template <int DEGREE>
constexpr short unsafe_atan2s_impl(float y, float x) {
  constexpr int maxshort = (int)(std::numeric_limits<short>::max()) + 1;
  constexpr short pi4 = short(maxshort / 4);
  constexpr short pi34 = short(3 * maxshort / 4);

  auto r = (std::abs(x) - std::abs(y)) / (std::abs(x) + std::abs(y));
  if (x < 0)
    r = -r;

  auto angle = (x >= 0) ? pi4 : pi34;
  angle += short(approx_atan2s_P<DEGREE>(r));

  return (y < 0) ? -angle : angle;
}

template <int DEGREE>
constexpr short unsafe_atan2s(float y, float x) {
  return unsafe_atan2s_impl<DEGREE>(y, x);
}

double deltaPhi(double phi1, double phi2) {
  double o2pi = 1. / (2. * M_PI);
  if (std::abs((phi1 - phi2)) <= double(M_PI))
    return (phi1 - phi2);
  double n = std::round((phi1 - phi2) * o2pi);
  return (phi1 - phi2) - n * double(2. * M_PI);
}

template<typename T>
void linSpace (const unsigned n, const double a, const double b, std::vector<T> &bins)
{
  T step = (b - a) / ((T) n);

  bins.clear ();
  for (T i = a; i < b + 0.5 * step; i += step)
    bins.push_back (i);
}

template<typename T>
void logSpace (const unsigned n, const double a, const double b, std::vector<T> &bins)
{
  T step = (b - a) / ((T) n);

  bins.clear ();
  for (T i = a; i < b + 0.5 * step; i += step)
    bins.push_back (pow (10.0, i));
}

template<typename T>
void writeHisto(std::vector<uint16_t>& vec, std::ofstream& file, std::vector<T>& bins)
{
  for (int i = 0; i < int(vec.size()); ++i){
    file << vec[i] << ",";
  }
  file << std::endl;
  for (int i = 0; i < int(bins.size()) - 1; ++i) file << (bins[i+1] + bins[i])/2.0 << ",";
  file << std::endl;
}
  // helper function that takes the layerPairId and returns two strings with the
  // inner and outer layer id
  // std::pair<std::string, std::string> getInnerOuterLayerNames(int const layerPairId) {
  //   // make a string from the Id (int)
  //   std::string index = std::to_string(layerPairId);
  //   // determine inner and outer layer name
  //   std::string innerLayerName;
  //   std::string outerLayerName;
  //   if (index.size() < 3) {
  //     innerLayerName = "0";
  //     outerLayerName = index;
  //   } else if (index.size() == 3) {
  //     innerLayerName = index.substr(0, 1);
  //     outerLayerName = index.substr(1, 3);
  //   } else {
  //     innerLayerName = index.substr(0, 2);
  //     outerLayerName = index.substr(2, 4);
  //   }
  //   if (outerLayerName[0] == '0') {
  //     outerLayerName = outerLayerName.substr(1, 2);
  //   }

  //   return {innerLayerName, outerLayerName};
  // }

  // // make bins logarithmic
  // void BinLogX(TH1* h) {
  //   TAxis* axis = h->GetXaxis();
  //   int bins = axis->GetNbins();

  //   float from = axis->GetXmin();
  //   float to = axis->GetXmax();
  //   float width = (to - from) / bins;
  //   std::vector<float> new_bins(bins + 1, 0);

  //   for (int i = 0; i <= bins; i++) {
  //     new_bins[i] = TMath::Power(10, from + i * width);
  //   }
  //   axis->Set(bins, new_bins.data());
  // }

  // // function to produce histogram with log scale on x (taken from MultiTrackValidator)
  // template <typename... Args>
  // dqm::reco::MonitorElement* make1DLogX(dqm::reco::DQMStore::IBooker& ibook, Args&&... args) {
  //   auto h = std::make_unique<TH1F>(std::forward<Args>(args)...);
  //   BinLogX(h.get());
  //   const auto& name = h->GetName();
  //   return ibook.book1D(name, h.release());
  // }

  // // function to produce profile with log scale on x (taken from MultiTrackValidator)
  // template <typename... Args>
  // dqm::reco::MonitorElement* makeProfileLogX(dqm::reco::DQMStore::IBooker& ibook, Args&&... args) {
  //   auto h = std::make_unique<TProfile>(std::forward<Args>(args)...);
  //   BinLogX(h.get());
  //   const auto& name = h->GetName();
  //   return ibook.bookProfile(name, h.release());
  // }

  // function that checks if two vector share a common element
  // template <typename T>
  // bool haveCommonElement(std::vector<T> const& v1, std::vector<T> const& v2) {
  //   return std::find_first_of(v1.begin(), v1.end(), v2.begin(), v2.end()) != v1.end();
  // }

}  // namespace simdoublets

// -------------------------------------------------------------------------------------------------------------
// constructors and destructor
// -------------------------------------------------------------------------------------------------------------

SimDoubletsAnalyzer::SimDoubletsAnalyzer(edm::ProductRegistry& reg)
    : simDoublets_getToken_(reg.consumes<SimDoubletsCollection>()),
      particles_getToken_(reg.consumes<ParticleSimpleSoA>()),
      hits_getToken_(reg.consumes<TrackingRecHitSimpleSoA>()) {//,
      // cellMinz_(iConfig.getParameter<std::vector<double>>("cellMinz")),
      // cellMaxz_(iConfig.getParameter<std::vector<double>>("cellMaxz")),
      // cellPhiCuts_(iConfig.getParameter<std::vector<int>>("cellPhiCuts")),
      // cellMaxr_(iConfig.getParameter<std::vector<double>>("cellMaxr")),
      // // cellMinYSizeB1_(iConfig.getParameter<int>("cellMinYSizeB1")),
      // // cellMinYSizeB2_(iConfig.getParameter<int>("cellMinYSizeB2")),
      // // cellMaxDYSize12_(iConfig.getParameter<int>("cellMaxDYSize12")),
      // // cellMaxDYSize_(iConfig.getParameter<int>("cellMaxDYSize")),
      // // cellMaxDYPred_(iConfig.getParameter<int>("cellMaxDYPred")),
      // cellZ0Cut_(iConfig.getParameter<double>("cellZ0Cut")),
      // cellPtCut_(iConfig.getParameter<double>("cellPtCut")) {
  // These could be taken from the geometry
  // get layer pairs from configuration
  std::vector<int> layerPairs{

    0,  1,  0,  4,  0,  11,  // BPIX1 (3)
    1,  2,  1,  4,  1,  11,  // BPIX2 (6)
    2,  3,  2,  4,  2,  11,  // BPIX3 (9)

    4,  5,  5,  6,  6,  7,  7,  8,  8,  9,  9,  10,  // POS (15)
    11, 12, 12, 13, 13, 14, 14, 15, 15, 16, 16, 17,  // NEG (21)

    0,  2,  0,  5,  0,  12, 0,  6,  0,  13,  // BPIX1 Jump (26)
    1,  3,  1,  5,  1,  12, 1,  6,  1,  13,  // BPIX2 Jump (31)

    4,  6,  5,  7,  6,  8,  7,  9,  8,  10,  // POS Jump (36)
    11, 13, 12, 14, 13, 15, 14, 16, 15, 17,  // NEG Jump (41)

  };

  cellMinz_ = {

    //0,  1     0,  4     0,  11
      2.0*-34.0,    1.0*4.0,      2.0*-48.0, 
    //1,  2     1,  4     1,  11      
      2.0*-44.0,    1.0*6.0,      2.0*-72.0, 
    //2,  3     2,  4     2,  11      
      2.0*-54.0,    1.0*11.0,     2.0*-96.0,
      
    //4,  5     5,  6     6,  7    7,  8    8,  9    9,  10      
      1.0*23.0,     1.0*30.0,     1.0*39.0,    1.0*50.0,    1.0*65.0,    1.0*82.0, 
    //11, 12    12, 13    13, 14   14, 15   15, 16   16, 17      
      2.0*-84.0,    2.0*-105.0,   2.0*-132.0,  2.0*-165.0,  2.0*-210.0,  2.0*-327.0, 
      
    //0,  2     0,  5   0,  12    0,  6    0,  13      
      -17.0,    7.0,    -24.0,    11.0,    -24.0,
    //1,  3     1,  5   1,  12    1,  6    1,  13       
      -17.0,    9.0,    -24.0,    13.0,    -24.0, 
      
    //4,  6     5,  7    6,  8    7,  9    8,  10      
      23.0,     30.0,    39.0,    50.0,    65.0, 
    //11, 13    12, 14   13, 15   14, 16   15, 17      
      -84.0,    -105.0,  -132.0,  -165.0,  -210.0
    
  };
  cellMaxz_ = {

    
    //0,  1     0,  4     0,  11
      2.0*34.0,     2.0*48.0,     1.0*-4.0,
    //1,  2     1,  4     1,  11
      2.0*44.0,     2.0*72.0,     1.0*-6.0,
    //2,  3     2,  4     2,  11
      2.0*54.0,     2.0*96.0,     1.0*-11.0,
      
    //4,  5    5,  6    6,  7    7,  8    8,  9    9,  10
      2.0*84.0,    2.0*105.0,   2.0*132.0,   2.0*165.0,   2.0*210.0,   2.0*327.0, 
    //11, 12   12, 13   13, 14   14, 15   15, 16   16, 17
      1.0*-23.0,   1.0*-30.0,   1.0*-39.0,   1.0*-50.0,   1.0*-65.0,   1.0*-82.0,
      
    //0,  2    0,  5    0,  12   0,  6    0,  13 
      17.0,    24.0,    -7.0,    24.0,    -11.0,
    //1,  3    1,  5    1,  12   1,  6    1,  13
      17.0,    24.0,    -9.0,    24.0,    -13.0,
      
    //4,  6    5,  7    6,  8    7,  9    8,  10
      84.0,    105.0,   132.0,   165.0,   210.0,
    //11, 13   12, 14   13, 15   14, 16   15, 17
      -23.0,   -30.0,   -39.0,   -50.0,   -65.0

  };
  constexpr int16_t phi0p05 = int16_t(1.1*522);
  constexpr int16_t phi0p06 = int16_t(1.1*626);
  constexpr int16_t phi0p07 = int16_t(1.1*730);
  cellPhiCuts_ = {

    //0,  1    0,  4    0,  11
      phi0p05, phi0p05, phi0p05, 
    //1,  2    1,  4    1,  11
      phi0p06, phi0p07, phi0p07,
    //2,  3    2,  4    2,  11
      phi0p06, phi0p07, phi0p07, 
      
    //4,  5    5,  6    6,  7    7,  8    8,  9    9,  10
      phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, 
    //11, 12   12, 13   13, 14   14, 15   15, 16   16, 17
      phi0p05, phi0p05, phi0p05, phi0p05, phi0p05, phi0p05,
      
    //0,  2    0,  5    0,  12   0,  6    0,  13
      phi0p05, phi0p07, phi0p07, phi0p07, phi0p07,
    //1,  3    1,  5    1,  12   1,  6    1,  13 
      phi0p05, phi0p07, phi0p07, phi0p07, phi0p07, 
      
    //4,  6    5,  7    6,  8    7,  9    8,  10
      phi0p07, phi0p07, phi0p07, phi0p07, phi0p07, 
    //11, 13   12, 14   13, 15   14, 16   15, 17
      phi0p07, phi0p07, phi0p07, phi0p07, phi0p07
    
  };
  cellMaxr_ = {
    
  //0,  1   0,  4   0,  11
    6.0,    6.0,    6.0,
  //1,  2   1,  4   1,  11
    7.0,    8.0,    8.0,
  //2,  3   2,  4   2,  11
    7.0,    8.0,    8.0,

  //4,  5   5,  6   6,  7   7,  8   8,  9   9,  10
    6.0,    6.0,    6.0,    6.0,    6.0,    6.0,
  //11, 12  12, 13  13, 14  14, 15  15, 16  16, 17
    6.0,    6.0,    6.0,    6.0,    6.0,    6.0,

  //0,  2    0,  5   0,  12  0,  6   0,  13 
    10.0,    5.0,    5.0,    5.0,    5.0,
  //1,  3    1,  5   1,  12  1,  6   1,  13
    12.0,    8.0,    8.0,    8.0,    8.0,

  //4,  6   5,  7   6,  8   7,  9   8,  10
    9.0,    9.0,    9.0,    8.0,    8.0,
  //11, 13  12, 14  13, 15  14, 16  15, 17
    9.0,    9.0,    9.0,    8.0,    8.0

  };
  cellZ0Cut_ = 1.7*12.0;
  cellPtCut_ = 0.5;

  // number of configured layer pairs
  size_t numLayerPairs = layerPairs.size() / 2;

  // fill the map of layer pairs
  for (size_t i{0}; i < numLayerPairs; i++) {
    int layerPairId = 100 * layerPairs[2 * i] + layerPairs[2 * i + 1];
    layerPairId2Index_.insert({layerPairId, i});
  }

  totalDoublets = 0;
  totalPassedDoublets = 0;

  nBins = 100;

  histoInnerZ.resize(numLayerPairs);
  histoDR.resize(numLayerPairs);
  histoZ0.resize(numLayerPairs);
  histoPT.resize(numLayerPairs);
  histoIPhi.resize(numLayerPairs);

  for (size_t i{0}; i < numLayerPairs; i++) {
    histoInnerZ[i].resize(nBins);
    histoDR[i].resize(nBins);
    histoZ0[i].resize(nBins);
    histoPT[i].resize(nBins);
    histoIPhi[i].resize(nBins);
  }

  // resize all histogram vectors, so that we can fill them according to the
  // layerPairIndex saved in the map that we just created
  // hVector_dr_.resize(numLayerPairs);
  // hVector_dphi_.resize(numLayerPairs);
  // hVector_idphi_.resize(numLayerPairs);
  // hVector_innerZ_.resize(numLayerPairs);
  // hVector_Ysize_.resize(numLayerPairs);
  // hVector_DYsize_.resize(numLayerPairs);
  // hVector_DYPred_.resize(numLayerPairs);
  // hVector_pass_dr_.resize(numLayerPairs);
  // hVector_pass_idphi_.resize(numLayerPairs);
  // hVector_pass_innerZ_.resize(numLayerPairs);
}


SimDoubletsAnalyzer::~SimDoubletsAnalyzer() {
  std::cout << totalDoublets << " -- " << totalPassedDoublets << std::endl;
  file.open("/data/user/borzari/cmssw/pixeltrack-standalone/outputSimDoublets.txt");
  if (file.is_open()) {

    for(int i = 0; i < int(histoInnerZ.size()); ++i) {
      file << "# Inner Z " << i << std::endl;
      simdoublets::writeHisto<double>(histoInnerZ[i], file, binsInnerZ);
      file << "# DR " << i << std::endl;
      simdoublets::writeHisto<double>(histoDR[i], file, binsDR);
      file << "# Z0 " << i << std::endl;
      simdoublets::writeHisto<double>(histoZ0[i], file, binsZ0);
      file << "# PT " << i << std::endl;
      simdoublets::writeHisto<double>(histoPT[i], file, binsPT);
      file << "# IPhi " << i << std::endl;
      simdoublets::writeHisto<double>(histoIPhi[i], file, binsIPhi);
    }

  }
}

// -------------------------------------------------------------------------------------------------------------
// member functions
// -------------------------------------------------------------------------------------------------------------


void SimDoubletsAnalyzer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // get simDoublets
  SimDoubletsCollection const& simDoubletsCollection = iEvent.get(simDoublets_getToken_);
  ParticleSimpleSoA const& particles = iEvent.get(particles_getToken_);
  TrackingRecHitSimpleSoA const& hits = iEvent.get(hits_getToken_);

  // create vectors for inner and outer RecHits of SimDoublets passing all cuts
  std::vector<uint32_t> innerRecHitsPassing;
  std::vector<uint32_t> outerRecHitsPassing;

  // initialize a bunch of variables that we will use in the coming for loops
  [[maybe_unused]] double true_pT, true_eta, inner_r, inner_z, inner_phi, outer_r, outer_z, outer_phi, dz, dr, dphi, z0, curvature, pT;
  [[maybe_unused]] int numSimDoublets, pass_numSimDoublets, inner_iphi, outer_iphi, idphi, layerPairId, layerPairIdIndex, passZWindow, passZ0Cutoff, passPTCut, passIPhiCut;//,
      // innerClusterSizeY;
  // bool doubletGetsCut, subjectToYsizeB1, subjectToYsizeB2, subjectToDYsize, subjectToDYsize12, subjectToDYPred;
  bool doubletGetsCut;

  passZWindow = 0;
  passZ0Cutoff = 0;
  passPTCut = 0;
  passIPhiCut = 0;

  // loop over SimDoublets (= loop over TrackingParticles)
  for (auto const& simDoublets : simDoubletsCollection) {
    // get true pT of the TrackingParticle
    true_pT = particles.pt(simDoublets.simParticle());
    true_eta = particles.eta(simDoublets.simParticle());

    // create the true RecHit doublets of the TrackingParticle
    auto doublets = simDoublets.getSimDoublets();

    // number of SimDoublets of the Tracking Particle
    numSimDoublets = doublets.size();
    // number of SimDoublets of the Tracking Particle passing all cuts
    pass_numSimDoublets = 0;

    // fill histograms for number of SimDoublets
    // h_numSimDoubletsPerTrackingParticle_->Fill(numSimDoublets);
    // h_numLayersPerTrackingParticle_->Fill(simDoublets.numLayers());

    // fill histograms for number of TrackingParticles
    // h_numTPVsPt_->Fill(true_pT);
    // h_numTPVsEta_->Fill(true_eta);

    // clear passing inner and outer RecHits
    innerRecHitsPassing.clear();
    outerRecHitsPassing.clear();

    // loop over those doublets
    for (auto const& doublet : doublets) {
      totalDoublets++;
      // RecHit properties
      inner_r = std::sqrt((doublet.innerGlobalPos(hits)[0]*doublet.innerGlobalPos(hits)[0])+(doublet.innerGlobalPos(hits)[1]*doublet.innerGlobalPos(hits)[1]));
      inner_z = doublet.innerGlobalPos(hits)[2];
      inner_phi = std::atan2(doublet.innerGlobalPos(hits)[1],doublet.innerGlobalPos(hits)[0]);
      inner_iphi = simdoublets::unsafe_atan2s<7>(doublet.innerGlobalPos(hits)[1], doublet.innerGlobalPos(hits)[0]);
      outer_r = std::sqrt((doublet.outerGlobalPos(hits)[0]*doublet.outerGlobalPos(hits)[0])+(doublet.outerGlobalPos(hits)[1]*doublet.outerGlobalPos(hits)[1]));
      outer_z = doublet.outerGlobalPos(hits)[2];
      outer_phi = std::atan2(doublet.outerGlobalPos(hits)[1],doublet.outerGlobalPos(hits)[0]);
      outer_iphi = simdoublets::unsafe_atan2s<7>(doublet.outerGlobalPos(hits)[1], doublet.outerGlobalPos(hits)[0]);

      dz = outer_z - inner_z;
      dr = outer_r - inner_r;
      dphi = simdoublets::deltaPhi(inner_phi, outer_phi);
      idphi = std::min(std::abs(int16_t(outer_iphi - inner_iphi)), std::abs(int16_t(inner_iphi - outer_iphi)));

      if (dr < 0.0) {
        std::cout << doublet.innerLayerId() << " -- " << doublet.outerLayerId() << std::endl;
      }

      // ----------------------------------------------------------
      // general plots (general folder)
      // ----------------------------------------------------------

      // outer layer vs inner layer of SimDoublets
      // h_layerPairs_->Fill(doublet.innerLayerId(), doublet.outerLayerId());

      // number of skipped layers by SimDoublets
      // h_numSkippedLayers_->Fill(doublet.numSkippedLayers());

      // ----------------------------------------------------------
      // layer pair independent cuts (global folder)
      // ----------------------------------------------------------

      // longitudinal impact parameter with respect to the beamspot
      z0 = std::abs(inner_r * outer_z - inner_z * outer_r) / dr;
      // h_z0_->Fill(z0);

      // radius of the circle defined by the two RecHits and the beamspot
      curvature = 1.f / 2.f * std::sqrt((dr / dphi) * (dr / dphi) + (inner_r * outer_r));
      // h_curvatureR_->Fill(curvature);

      // pT that this curvature radius corresponds to
      // pT = curvature / 87.78f;
      pT = curvature / 128.29f;
      // h_pTFromR_->Fill(pT);

      // ----------------------------------------------------------
      // layer pair dependent cuts (sub-folders for layer pairs)
      // ----------------------------------------------------------

      // first, get layer pair Id and exclude layer pairs that are not considered
      layerPairId = doublet.layerPairId();
      if (layerPairId2Index_.find(layerPairId) == layerPairId2Index_.end()) {
        continue;
      }

      // get the position of the layer pair in the vectors of histograms
      layerPairIdIndex = layerPairId2Index_.at(layerPairId);

      // dr = (outer_r - inner_r) histogram
      // hVector_dr_[layerPairIdIndex]->Fill(dr);

      // dphi histogram
      // hVector_dphi_[layerPairIdIndex]->Fill(dphi);
      // hVector_idphi_[layerPairIdIndex]->Fill(idphi);

      // z of the inner RecHit histogram
      // hVector_innerZ_[layerPairIdIndex]->Fill(inner_z);

      // ----------------------------------------------------------
      // cluster size cuts (global + sub-folders for layer pairs)
      // ----------------------------------------------------------

      // cluster size in local y histogram
      // innerClusterSizeY = doublet.innerRecHit()->cluster()->sizeY();
      // hVector_Ysize_[layerPairIdIndex]->Fill(innerClusterSizeY);

      // create bool that indicates if the doublet gets cut
      doubletGetsCut = false;
      // create bools that trace if doublet is subject to any clsuter size cut
      // subjectToYsizeB1 = false;
      // subjectToYsizeB2 = false;
      // subjectToDYsize = false;
      // subjectToDYsize12 = false;
      // subjectToDYPred = false;

      vecInnerZ.push_back(inner_z);
      vecDR.push_back(dr);
      vecZ0.push_back(z0);
      vecPT.push_back(pT);
      vecIPhi.push_back(idphi);
      vecLayerPairId.push_back(layerPairIdIndex);

      // apply all cuts that do not depend on the cluster size
      // z window cut
      if (inner_z < cellMinz_[layerPairIdIndex] || inner_z > cellMaxz_[layerPairIdIndex]) {
        doubletGetsCut = true;
        passZWindow++;
      }
      // else {passZWindow++;}
      // z0cutoff
      if (dr > cellMaxr_[layerPairIdIndex] || dr < 0 || z0 > cellZ0Cut_) {
        doubletGetsCut = true;
        passZ0Cutoff++;
      }
      // else {passZ0Cutoff++;}
      // ptcut
      if (pT < cellPtCut_) {
        doubletGetsCut = true;
        passPTCut++;
      }
      // else {passPTCut++;}
      // iphicut
      if (idphi > cellPhiCuts_[layerPairIdIndex]) {
        doubletGetsCut = true;
        passIPhiCut++;
      }
      // else {passIPhiCut++;}

      // determine the moduleId
      // const GeomDetUnit* geomDetUnit = doublet.innerRecHit()->det();
      // const uint32_t moduleId = geomDetUnit->index();

      // define bools needed to decide on cutting parameters
      // const bool innerInB1 = (doublet.innerLayerId() == 0);
      // const bool innerInB2 = (doublet.innerLayerId() == 1);
      // // const bool isOuterLadder = (0 == (moduleId / 8) % 2);  // check if this even makes sense in Phase-2
      // const bool innerInBarrel = (doublet.innerLayerId() < 4);
      // const bool outerInBarrel = (doublet.outerLayerId() < 4);

      // // histograms for clusterCut
      // // cluster size in local y
      // if (!outerInBarrel) {
      //   if (innerInB1 && isOuterLadder) {
      //     subjectToYsizeB1 = true;
      //     h_YsizeB1_->Fill(innerClusterSizeY);
      //     // apply the cut
      //     if (innerClusterSizeY < cellMinYSizeB1_) {
      //       doubletGetsCut = true;
      //     }
      //   }
      //   if (innerInB2) {
      //     subjectToYsizeB2 = true;
      //     h_YsizeB2_->Fill(innerClusterSizeY);
      //     // apply the cut
      //     if (innerClusterSizeY < cellMinYSizeB2_) {
      //       doubletGetsCut = true;
      //     }
      //   }
      // }

      // // histograms for zSizeCut
      // int DYsize{0}, DYPred{0};
      // if (innerInBarrel) {
      //   if (outerInBarrel) {  // onlyBarrel
      //     DYsize = std::abs(innerClusterSizeY - doublet.outerRecHit()->cluster()->sizeY());
      //     if (innerInB1 && isOuterLadder) {
      //       subjectToDYsize12 = true;
      //       // hVector_DYsize_[layerPairIdIndex]->Fill(DYsize);
      //       h_DYsize12_->Fill(DYsize);
      //       // apply the cut
      //       if (DYsize > cellMaxDYSize12_) {
      //         doubletGetsCut = true;
      //       }
      //     } else if (!innerInB1) {
      //       subjectToDYsize = true;
      //       // hVector_DYsize_[layerPairIdIndex]->Fill(DYsize);
      //       h_DYsize_->Fill(DYsize);
      //       // apply the cut
      //       if (DYsize > cellMaxDYSize_) {
      //         doubletGetsCut = true;
      //       }
      //     }
      //   } else {  // not onlyBarrel
      //     subjectToDYPred = true;
      //     DYPred = std::abs(innerClusterSizeY - int(std::abs(dz / dr) * pixelTopology::Phase2::dzdrFact + 0.5f));
      //     // hVector_DYPred_[layerPairIdIndex]->Fill(DYPred);
      //     h_DYPred_->Fill(DYPred);
      //     // apply the cut
      //     if (DYPred > cellMaxDYPred_) {
      //       doubletGetsCut = true;
      //     }
      //   }
      // }

      // ----------------------------------------------------------
      // all kinds of plots for doublets passing all cuts
      // ----------------------------------------------------------

      // fill the number histograms
      // histogram of all valid doublets
      // h_numVsPt_->Fill(true_pT);
      // h_numVsEta_->Fill(true_eta);

      // if the doublet passes all cuts
      if (!doubletGetsCut) {
        // increment number of SimDoublets passing all cuts
        pass_numSimDoublets++;
        totalPassedDoublets++;

        // fill histogram of doublets that pass all cuts
        // h_pass_layerPairs_->Fill(doublet.innerLayerId(), doublet.outerLayerId());
        // h_pass_numVsPt_->Fill(true_pT);
        // h_pass_numVsEta_->Fill(true_eta);

        // also put the inner/outer RecHit in the respective vector
        innerRecHitsPassing.push_back(doublet.innerRecHit());
        outerRecHitsPassing.push_back(doublet.outerRecHit());

        // fill pass_ histograms
        // h_pass_z0_->Fill(z0);
        // h_pass_pTFromR_->Fill(pT);
        // hVector_pass_dr_[layerPairIdIndex]->Fill(dr);
        // hVector_pass_idphi_[layerPairIdIndex]->Fill(idphi);
        // hVector_pass_innerZ_[layerPairIdIndex]->Fill(inner_z);
        // if (subjectToDYPred) {
        //   h_pass_DYPred_->Fill(DYPred);
        // }
        // if (subjectToDYsize) {
        //   h_pass_DYsize_->Fill(DYsize);
        // }
        // if (subjectToDYsize12) {
        //   h_pass_DYsize12_->Fill(DYsize);
        // }
        // if (subjectToYsizeB1) {
        //   h_pass_YsizeB1_->Fill(innerClusterSizeY);
        // }
        // if (subjectToYsizeB2) {
        //   h_pass_YsizeB2_->Fill(innerClusterSizeY);
        // }
      }

      // std::cout << true_pT << " -- " << true_eta << " -- " << inner_r << " -- " << inner_z << " -- " << inner_phi << " -- " << outer_r << " -- " << outer_z << " -- " << outer_phi << " -- " << dz << " -- " << dr << " -- " << dphi << " -- " << z0 << " -- " << curvature << " -- " << pT << " -- " << inner_iphi << " -- " << outer_iphi << " -- " << idphi << " -- " << layerPairId << " -- " << !doubletGetsCut << std::endl;

    }  // end loop over those doublets

    // // Now check if the TrackingParticle is reconstructable by at least two conencted SimDoublets surviving the cuts
    // if (simdoublets::haveCommonElement<SiPixelRecHitRef>(innerRecHitsPassing, outerRecHitsPassing)) {
    //   h_pass_numTPVsPt_->Fill(true_pT);
    //   h_pass_numTPVsEta_->Fill(true_eta);
    // }

    // // Fill the efficiency profile per Tracking Particle only if the TP has at least one SimDoublet
    // if (numSimDoublets > 0) {
    //   h_effSimDoubletsPerTPVsEta_->Fill(true_eta, pass_numSimDoublets / numSimDoublets);
    //   h_effSimDoubletsPerTPVsPt_->Fill(true_pT, pass_numSimDoublets / numSimDoublets);
    // }

    // std::cout << numSimDoublets << " -- " << pass_numSimDoublets << std::endl;

  }  // end loop over SimDoublets (= loop over TrackingParticles)
  std::cout << passZWindow << " -- " << passZ0Cutoff << " -- " << passPTCut << " -- " << passIPhiCut << std::endl;

  simdoublets::linSpace<double> (nBins, *(std::min_element(vecInnerZ.begin(), vecInnerZ.end())), *(std::max_element(vecInnerZ.begin(), vecInnerZ.end())), binsInnerZ);
  simdoublets::linSpace<double> (nBins, *(std::min_element(vecDR.begin(), vecDR.end())), *(std::max_element(vecDR.begin(), vecDR.end())), binsDR);
  simdoublets::linSpace<double> (nBins, *(std::min_element(vecZ0.begin(), vecZ0.end())), *(std::max_element(vecZ0.begin(), vecZ0.end())), binsZ0);
  simdoublets::logSpace<double> (nBins, -0.5, 2.0, binsPT);
  simdoublets::linSpace<double> (nBins, *(std::min_element(vecIPhi.begin(), vecIPhi.end())), *(std::max_element(vecIPhi.begin(), vecIPhi.end())), binsIPhi);

  for(int i = 0; i < nBins; ++i){
    for(int j = 0; j < int(vecInnerZ.size()); ++j){
      int k = vecLayerPairId[j];
      if(vecInnerZ[j] > binsInnerZ[i] && vecInnerZ[j] < binsInnerZ[i+1]) histoInnerZ[k][i]++;
      if(vecDR[j] > binsDR[i] && vecDR[j] < binsDR[i+1]) histoDR[k][i]++;
      if(vecZ0[j] > binsZ0[i] && vecZ0[j] < binsZ0[i+1]) histoZ0[k][i]++;
      if(vecPT[j] > binsPT[i] && vecPT[j] < binsPT[i+1]) histoPT[k][i]++;
      if(vecIPhi[j] > binsIPhi[i] && vecIPhi[j] < binsIPhi[i+1]) histoIPhi[k][i]++;
    }
  }

}

// // booking the histograms
// 
// void SimDoubletsAnalyzer::bookHistograms(DQMStore::IBooker& ibook,
//                                                         edm::Run const& run,
//                                                         edm::EventSetup const& iSetup) {
//   // set some common parameters
//   int pTNBins = 50;
//   double pTmin = log10(0.01);
//   double pTmax = log10(1000);
//   int etaNBins = 80;
//   double etamin = -4.;
//   double etamax = 4.;

//   // ----------------------------------------------------------
//   // booking general histograms (general folder)
//   // ----------------------------------------------------------

//   ibook.setCurrentFolder(folder_ + "/general");

//   // overview histograms and profiles
//   h_effSimDoubletsPerTPVsPt_ =
//       simdoublets::makeProfileLogX(ibook,
//                                    "efficiencyPerTP_vs_pT",
//                                    "SimDoublets efficiency per TP vs p_{T}; TP transverse momentum p_{T} [GeV]; "
//                                    "Average fraction of SimDoublets per TP passing all cuts",
//                                    pTNBins,
//                                    pTmin,
//                                    pTmax,
//                                    0,
//                                    1,
//                                    " ");
//   h_effSimDoubletsPerTPVsEta_ = ibook.bookProfile("efficiencyPerTP_vs_eta",
//                                                   "SimDoublets efficiency per TP vs #eta; TP transverse momentum #eta; "
//                                                   "Average fraction of SimDoublets per TP passing all cuts",
//                                                   etaNBins,
//                                                   etamin,
//                                                   etamax,
//                                                   0,
//                                                   1,
//                                                   " ");
//   h_layerPairs_ = ibook.book2D("layerPairs",
//                                "Layer pairs in SimDoublets; Inner layer ID; Outer layer ID",
//                                TrackerTraits::numberOfLayers,
//                                -0.5,
//                                -0.5 + TrackerTraits::numberOfLayers,
//                                TrackerTraits::numberOfLayers,
//                                -0.5,
//                                -0.5 + TrackerTraits::numberOfLayers);
//   h_pass_layerPairs_ = ibook.book2D("pass_layerPairs",
//                                     "Layer pairs in SimDoublets passing all cuts; Inner layer ID; Outer layer ID",
//                                     TrackerTraits::numberOfLayers,
//                                     -0.5,
//                                     -0.5 + TrackerTraits::numberOfLayers,
//                                     TrackerTraits::numberOfLayers,
//                                     -0.5,
//                                     -0.5 + TrackerTraits::numberOfLayers);
//   h_numSkippedLayers_ = ibook.book1D(
//       "numSkippedLayers", "Number of skipped layers; Number of skipped layers; Number of SimDoublets", 16, -1.5, 14.5);
//   h_numSimDoubletsPerTrackingParticle_ =
//       ibook.book1D("numSimDoubletsPerTrackingParticle",
//                    "Number of SimDoublets per Tracking Particle; Number of SimDoublets; Number of Tracking Particles",
//                    31,
//                    -0.5,
//                    30.5);
//   h_numLayersPerTrackingParticle_ =
//       ibook.book1D("numLayersPerTrackingParticle",
//                    "Number of layers hit by Tracking Particle; Number of layers; Number of Tracking Particles",
//                    29,
//                    -0.5,
//                    28.5);
//   h_numTPVsPt_ = simdoublets::make1DLogX(
//       ibook,
//       "numTPVsPt",
//       "Total number of TrackingParticles; True transverse momentum p_{T} [GeV]; Total number of TrackingParticles",
//       pTNBins,
//       pTmin,
//       pTmax);
//   h_pass_numTPVsPt_ = simdoublets::make1DLogX(ibook,
//                                               "pass_numTPVsPt",
//                                               "Reconstructable TrackingParticles (two or more connected SimDoublets "
//                                               "pass cuts); True transverse momentum p_{T} [GeV]; "
//                                               "Number of reconstructable TrackingParticles",
//                                               pTNBins,
//                                               pTmin,
//                                               pTmax);
//   h_numTPVsEta_ =
//       ibook.book1D("numTPVsEta",
//                    "Total number of TrackingParticles; True pseudorapidity #eta; Total number of TrackingParticles",
//                    etaNBins,
//                    etamin,
//                    etamax);
//   h_pass_numTPVsEta_ = ibook.book1D("pass_numTPVsEta",
//                                     "Reconstructable TrackingParticles (two or more connected SimDoublets "
//                                     "pass cuts); True pseudorapidity #eta; Number of reconstructable TrackingParticles",
//                                     etaNBins,
//                                     etamin,
//                                     etamax);
//   h_numVsPt_ = simdoublets::make1DLogX(
//       ibook,
//       "numVsPt",
//       "Total number of SimDoublets; True transverse momentum p_{T} [GeV]; Total number of SimDoublets",
//       pTNBins,
//       pTmin,
//       pTmax);
//   h_pass_numVsPt_ = simdoublets::make1DLogX(ibook,
//                                             "pass_numVsPt",
//                                             "Number of passing SimDoublets; True transverse momentum p_{T} [GeV]; "
//                                             "Number of SimDoublets passing all cuts",
//                                             pTNBins,
//                                             pTmin,
//                                             pTmax);
//   h_numVsEta_ = ibook.book1D("numVsEta",
//                              "Total number of SimDoublets; True pseudorapidity #eta; Total number of SimDoublets",
//                              etaNBins,
//                              etamin,
//                              etamax);
//   h_pass_numVsEta_ =
//       ibook.book1D("pass_numVsEta",
//                    "Number of SimDoublets; True pseudorapidity #eta; Number of SimDoublets passing all cuts",
//                    etaNBins,
//                    etamin,
//                    etamax);

//   // -------------------------------------------------------------
//   // booking layer pair independent cut histograms (global folder)
//   // -------------------------------------------------------------

//   ibook.setCurrentFolder(folder_ + "/cutParameters/global");

//   // histogram for z0cutoff  (z0Cut)
//   h_z0_ = ibook.book1D("z0", "z_{0}; Longitudinal impact parameter z_{0} [cm]; Number of SimDoublets", 51, -1, 50);
//   h_pass_z0_ = ibook.book1D(
//       "pass_z0",
//       "z_{0} of SimDoublets passing all cuts; Longitudinal impact parameter z_{0} [cm]; Number of SimDoublets",
//       51,
//       -1,
//       50);

//   // histograms for ptcut  (ptCut)
//   h_curvatureR_ = ibook.book1D(
//       "curvatureR", "Curvature from SimDoublet+beamspot; Curvature radius [cm] ; Number of SimDoublets", 100, 0, 1000);
//   h_pTFromR_ = simdoublets::make1DLogX(
//       ibook,
//       "pTFromR",
//       "Transverse momentum from curvature; Transverse momentum p_{T} [GeV]; Number of SimDoublets",
//       pTNBins,
//       pTmin,
//       pTmax);
//   h_pass_pTFromR_ = simdoublets::make1DLogX(ibook,
//                                             "pass_pTFromR",
//                                             "Transverse momentum from curvature of SimDoublets passing all cuts; "
//                                             "Transverse momentum p_{T} [GeV]; Number of SimDoublets",
//                                             pTNBins,
//                                             pTmin,
//                                             pTmax);

//   // histograms for clusterCut  (minYsizeB1 and minYsizeB2)
//   h_YsizeB1_ = ibook.book1D(
//       "YsizeB1",
//       "Cluster size along z (inner from B1); Size along z of inner cluster [num of pixels]; Number of SimDoublets",
//       51,
//       -1,
//       50);
//   h_YsizeB2_ = ibook.book1D(
//       "YsizeB2",
//       "Cluster size along z (inner not from B1); Size along z of inner cluster [num of pixels]; Number of SimDoublets",
//       51,
//       -1,
//       50);
//   h_pass_YsizeB1_ = ibook.book1D("pass_YsizeB1",
//                                  "Cluster size along z of SimDoublets passing all cuts (inner from B1); Size along z "
//                                  "of inner cluster [num of pixels]; Number of SimDoublets",
//                                  51,
//                                  -1,
//                                  50);
//   h_pass_YsizeB2_ = ibook.book1D("pass_YsizeB2",
//                                  "Cluster size along z of SimDoublets passing all cuts (inner not from B1); Size along "
//                                  "z of inner cluster [num of pixels]; Number of SimDoublets",
//                                  51,
//                                  -1,
//                                  50);

//   // histograms for zSizeCut  (maxDYsize12, maxDYsize and maxDYPred)
//   h_DYsize12_ =
//       ibook.book1D("DYsize12",
//                    "Difference in cluster size along z (inner from B1); Absolute difference in cluster size along z of "
//                    "the two RecHits [num of pixels]; Number of SimDoublets",
//                    31,
//                    -1,
//                    30);
//   h_DYsize_ = ibook.book1D("DYsize",
//                            "Difference in cluster size along z; Absolute difference in cluster size along z of the two "
//                            "RecHits [num of pixels]; Number of SimDoublets",
//                            31,
//                            -1,
//                            30);
//   h_DYPred_ = ibook.book1D("DYPred",
//                            "Difference between actual and predicted cluster size along z of inner cluster; Absolute "
//                            "difference [num of pixels]; Number of SimDoublets",
//                            201,
//                            -1,
//                            200);
//   h_pass_DYsize12_ = ibook.book1D("pass_DYsize12",
//                                   "Difference in cluster size along z of SimDoublets passing all cuts (inner from B1); "
//                                   "Absolute difference in cluster size along z of "
//                                   "the two RecHits [num of pixels]; Number of SimDoublets",
//                                   31,
//                                   -1,
//                                   30);
//   h_pass_DYsize_ =
//       ibook.book1D("pass_DYsize",
//                    "Difference in cluster size along z of SimDoublets passing all cuts; Absolute difference in "
//                    "cluster size along z of the two RecHits [num of pixels]; Number of SimDoublets",
//                    31,
//                    -1,
//                    30);
//   h_pass_DYPred_ =
//       ibook.book1D("pass_DYPred",
//                    "Difference between actual and predicted cluster size along z of inner cluster of SimDoublets "
//                    "passing all cuts; Absolute difference [num of pixels]; Number of SimDoublets",
//                    201,
//                    -1,
//                    200);

//   // -----------------------------------------------------------------------
//   // booking layer pair dependent histograms (sub-folders for layer pairs)
//   // -----------------------------------------------------------------------

//   // loop through valid layer pairs and add for each one booked hist per vector
//   for (auto id = layerPairId2Index_.begin(); id != layerPairId2Index_.end(); ++id) {
//     // get the position of the layer pair in the histogram vectors
//     int layerPairIdIndex = id->second;

//     // get layer names from the layer pair Id
//     auto layerNames = simdoublets::getInnerOuterLayerNames(id->first);
//     std::string innerLayerName = layerNames.first;
//     std::string outerLayerName = layerNames.second;

//     // name the sub-folder for the layer pair "lp_${innerLayerId}_${outerLayerId}"
//     std::string subFolderName = "/cutParameters/lp_" + innerLayerName + "_" + outerLayerName;

//     // layer mentioning in histogram titles
//     std::string layerTitle = "(layers (" + innerLayerName + "," + outerLayerName + "))";

//     // set folder to the sub-folder for the layer pair
//     ibook.setCurrentFolder(folder_ + subFolderName);

//     // histogram for z0cutoff  (maxr)
//     // hVector_dr_.at(layerPairIdIndex) = ibook.book1D(
//         "dr",
//         "dr of RecHit pair " + layerTitle + "; dr between outer and inner RecHit [cm]; Number of SimDoublets",
//         31,
//         -1,
//         30);
//     // hVector_pass_dr_.at(layerPairIdIndex) = ibook.book1D(
//         "pass_dr",
//         "dr of RecHit pair " + layerTitle +
//             " for SimDoublets passing all cuts; dr between outer and inner RecHit [cm]; Number of SimDoublets",
//         31,
//         -1,
//         30);

//     // histograms for iphicut  (phiCuts)
//     // hVector_dphi_.at(layerPairIdIndex) = ibook.book1D(
//         "dphi",
//         "dphi of RecHit pair " + layerTitle + "; d#phi between outer and inner RecHit [rad]; Number of SimDoublets",
//         50,
//         -M_PI,
//         M_PI);
//     // hVector_idphi_.at(layerPairIdIndex) =
//         ibook.book1D("idphi",
//                      "idphi of RecHit pair " + layerTitle +
//                          "; Absolute int d#phi between outer and inner RecHit; Number of SimDoublets",
//                      50,
//                      0,
//                      1000);
//     // hVector_pass_idphi_.at(layerPairIdIndex) = ibook.book1D("pass_idphi",
//                                                             "idphi of RecHit pair " + layerTitle +
//                                                                 " for SimDoublets passing all cuts; Absolute int d#phi "
//                                                                 "between outer and inner RecHit; Number of SimDoublets",
//                                                             50,
//                                                             0,
//                                                             1000);

//     // histogram for z window  (minz and maxz)
//     // hVector_innerZ_.at(layerPairIdIndex) =
//         ibook.book1D("innerZ",
//                      "z of the inner RecHit " + layerTitle + "; z of inner RecHit [cm]; Number of SimDoublets",
//                      100,
//                      -300,
//                      300);
//     // hVector_pass_innerZ_.at(layerPairIdIndex) =
//         ibook.book1D("pass_innerZ",
//                      "z of the inner RecHit " + layerTitle +
//                          " for SimDoublets passing all cuts; z of inner RecHit [cm]; Number of SimDoublets",
//                      100,
//                      -300,
//                      300);

//     // histograms for cluster size and size differences
//     // hVector_DYsize_.at(layerPairIdIndex) =
//         ibook.book1D("DYsize",
//                      "Difference in cluster size along z between outer and inner RecHit " + layerTitle +
//                          "; Absolute difference in cluster size along z of the two "
//                          "RecHits [num of pixels]; Number of SimDoublets",
//                      51,
//                      -1,
//                      50);
//     // hVector_DYPred_.at(layerPairIdIndex) =
//         ibook.book1D("DYPred",
//                      "Difference between actual and predicted cluster size along z of inner cluster " + layerTitle +
//                          "; Absolute difference [num of pixels]; Number of SimDoublets",
//                      51,
//                      -1,
//                      50);
//     // hVector_Ysize_.at(layerPairIdIndex) = ibook.book1D(
//         "Ysize",
//         "Cluster size along z " + layerTitle + "; Size along z of inner cluster [num of pixels]; Number of SimDoublets",
//         51,
//         -1,
//         50);
//   }
// }

// define this as a plug-in
DEFINE_FWK_MODULE(SimDoubletsAnalyzer);
