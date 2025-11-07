#include <iostream>
#include <cmath>
#include <utility>
#include <string>
#include <map>
#include <unordered_map>
#include <set>
#include <memory>

#include "AlpakaDataFormats/TracksHost.h"
#include "AlpakaDataFormats/TrackingRecHitsHost.h"
#include "AlpakaDataFormats/ZVertexHost.h"
#include "AlpakaDataFormats/SimpleMapHost.h"
#include "AlpakaDataFormats/ParticleHost.h"

#include "Framework/EDProducer.h"
#include "Framework/Event.h"
#include "Framework/EventSetup.h"
#include "Framework/PluginFactory.h"
#include "Framework/ConfigRegistry.h"

#include <TFile.h>
#include <TTree.h>

class SimpleTrackValidation : public edm::EDProducer {
public:
  explicit SimpleTrackValidation(edm::ProductRegistry& reg, edm::Config const& cfg);
  ~SimpleTrackValidation() override = default;

private:
  void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override;
  void endJob() override;

  std::pair<int, uint32_t> getMostRepeatingPart(const std::vector<uint32_t>& vec) const;
  double deltaPhi(double phi1, double phi2) const;
  void logSpace(unsigned n, double a, double b, std::vector<double>& bins) const;
  void linSpace(unsigned n, double a, double b, std::vector<double>& bins) const;
  double stdDev(const std::vector<double>& vec) const;

  bool isGoodParticle(const sim::ParticleSoAConstView& parts, uint32_t idx) const;

  // mapping helpers
  uint32_t nTracks(const reco::TrackSoAConstView& tracks) const;
  int trackNHits(const reco::TrackSoAConstView& tracks, int i) const;
  std::pair<uint32_t, uint32_t> trackHitRange(const reco::TrackSoAConstView& tracks, int i) const;

  // === EDM tokens ===
  edm::EDGetTokenT<sim::ParticleHost> tSimpleParticles_;
  edm::EDGetTokenT<reco::TracksHost> tokenTracks_;
  edm::EDGetTokenT<reco::TrackingRecHitHost> tokenHits_;
  edm::EDGetTokenT<ZVertexHost> tokenVertex_; //TODO: put me under reco::
  edm::EDGetTokenT<utils::SimpleMapHost> tokenHitMap_;

  // === configuration (from JSON) ===
  const double minPt_;
  const double maxEta_;
  const unsigned int nBins_;
  const unsigned int minHits_;
  const float purity_;
  const std::string outputFileName_;
  std::set<int> allowedPdgIds_;

  // === global counters ===
  int totalParticles_;
  int trueTracks_;
  int fakeTracks_;
  int duplicateTracks_;
  // === binning ===
  
  std::vector<double> binsPt_;
  std::vector<double> binsEtaPhi_;

  // efficiencies: num/den per bin
  std::vector<double> effPtDen_;
  std::vector<double> effPtNum_;
  std::vector<double> effEtaDen_;
  std::vector<double> effEtaNum_;
  std::vector<double> effPhiDen_;
  std::vector<double> effPhiNum_;

  // fake rates: num/den per bin
  std::vector<double> fakePtDen_;
  std::vector<double> fakePtNum_;
  std::vector<double> fakeEtaDen_;
  std::vector<double> fakeEtaNum_;
  std::vector<double> fakePhiDen_;
  std::vector<double> fakePhiNum_;

  // dup rates: num/den per bin
  std::vector<double> dupPtDen_;
  std::vector<double> dupPtNum_;
  std::vector<double> dupEtaDen_;
  std::vector<double> dupEtaNum_;
  std::vector<double> dupPhiDen_;
  std::vector<double> dupPhiNum_;

  // resolutions: bin index -> vector of residuals
  std::map<int, std::vector<double>> resPt_;
  std::map<int, std::vector<double>> resEta_;
  std::map<int, std::vector<double>> resPhi_;
  std::map<int, std::vector<double>> resD0_;
  std::map<int, std::vector<double>> resDZ_;

//   TODO: use RNTuples 
//   using RNTupleModel = ROOT::Experimental::RNTupleModel;
//   using RNTupleWriter = ROOT::Experimental::RNTupleWriter;

//   std::unique_ptr<RNTupleWriter> ntuple_;
//   std::shared_ptr<int> fieldKind_;      // 0=eff, 1=fake, 2=res
//   std::shared_ptr<int> fieldCoord_;     // 0=pt, 1=eta, 2=phi, 3=d0, 4=dz
//   std::shared_ptr<float> fieldBinCenter_;
//   std::shared_ptr<float> fieldValue_;
//   std::shared_ptr<float> fieldNum_;
//   std::shared_ptr<float> fieldDen_;
    std::unique_ptr<TFile> outFile_;
    TTree* ttree_;

    int kind_;      // 0=eff, 1=fake, 2=res
    int coord_;     // 0=pt, 1=eta, 2=phi, 3=d0, 4=dz
    float binCenter_;
    float value_;
    float num_;
    float den_;
};


namespace {
  enum class kType : int { Efficiency = 0, FakeRate = 1, Duplicates = 2, Resolution = 3 };
  enum class kCoord : int { Pt = 0, Eta = 1, Phi = 2, D0 = 3, DZ = 4 };
}

SimpleTrackValidation::SimpleTrackValidation(edm::ProductRegistry& reg, edm::Config const& cfg)
    : tSimpleParticles_(reg.consumes<sim::ParticleHost>()),
      tokenTracks_(reg.consumes<reco::TracksHost>()),
      tokenHits_(reg.consumes<reco::TrackingRecHitHost>()),
      tokenVertex_(reg.consumes<ZVertexHost>()),
      tokenHitMap_(reg.consumes<utils::SimpleMapHost>()),
      minPt_(cfg.value("minPt", -1.0)),          
      maxEta_(cfg.value("maxEta", 999.)),
      nBins_(cfg.value("nBins", 40)),
      minHits_(cfg.value("minHits",4)),
      purity_(cfg.value("trackPurity",0.75)),
      outputFileName_(cfg.value("outputFile", std::string("pixelTrackValidation_ntuple.root")))
     {

    totalParticles_ = 0;
    trueTracks_ = 0;
    fakeTracks_ = 0;
    duplicateTracks_ = 0;

    {
        if (cfg.contains("pdgIds") && cfg["pdgIds"].is_array()) {
        for (auto const& v : cfg["pdgIds"]) {
            allowedPdgIds_.insert(v.get<int>());
        }
        }
    }


    logSpace(nBins_, -0.1, 2.0, binsPt_);
    linSpace(nBins_, -4.0, 4.0, binsEtaPhi_);

    effPtDen_.assign(nBins_, 0.0);
    effPtNum_.assign(nBins_, 0.0);
    effEtaDen_.assign(nBins_, 0.0);
    effEtaNum_.assign(nBins_, 0.0);
    effPhiDen_.assign(nBins_, 0.0);
    effPhiNum_.assign(nBins_, 0.0);

    fakePtDen_.assign(nBins_, 0.0);
    fakePtNum_.assign(nBins_, 0.0);
    fakeEtaDen_.assign(nBins_, 0.0);
    fakeEtaNum_.assign(nBins_, 0.0);
    fakePhiDen_.assign(nBins_, 0.0);
    fakePhiNum_.assign(nBins_, 0.0);

    dupPtDen_.assign(nBins_, 0.0);
    dupPtNum_.assign(nBins_, 0.0);
    dupEtaDen_.assign(nBins_, 0.0);
    dupEtaNum_.assign(nBins_, 0.0);
    dupPhiDen_.assign(nBins_, 0.0);
    dupPhiNum_.assign(nBins_, 0.0);

    // // --- RNTuple creation ---
    // auto model = RNTupleModel::Create();
    // fieldKind_      = model->MakeField<int>("kind", 0);
    // fieldCoord_     = model->MakeField<int>("coord", 0);
    // fieldBinCenter_ = model->MakeField<float>("binCenter", 0.f);
    // fieldValue_     = model->MakeField<float>("value", 0.f);
    // fieldNum_       = model->MakeField<float>("num", 0.f);
    // fieldDen_       = model->MakeField<float>("den", 0.f);

    // ntuple_ = RNTupleWriter::Recreate(std::move(model),
    //                                     "PixelTrackValidation",
    //                                     outputFileName_);

    // open file and create tree
    outFile_ = std::make_unique<TFile>(outputFileName_.c_str(), "RECREATE");
    ttree_ = new TTree("Validation", "PixelTrack Validation");

    ttree_->Branch("kind", &kind_, "kind/I");
    ttree_->Branch("coord", &coord_, "coord/I");
    ttree_->Branch("binCenter", &binCenter_, "binCenter/F");
    ttree_->Branch("value", &value_, "value/F");
    ttree_->Branch("num", &num_, "num/F");
    ttree_->Branch("den", &den_, "den/F");

}

bool SimpleTrackValidation::isGoodParticle(const sim::ParticleSoAConstView& parts, uint32_t idx) const {
  const float pt  = parts.pt(idx);
  const float eta = parts.eta(idx);
  const int pdg   = parts.pdgID(idx);

  if (pt < minPt_) return false;
  if (std::abs(eta) > maxEta_) return false;

  if (!allowedPdgIds_.empty()) {
    if (allowedPdgIds_.find(pdg) == allowedPdgIds_.end())
      return false;
  }

  return true;
}

// =======================================================================================
// Track / hit indexing helpers for the new SoA
// =======================================================================================

uint32_t SimpleTrackValidation::nTracks(const reco::TrackSoAConstView& tracks) const {
  // metadata().size() is how other code gets the number of entries
  return tracks.metadata().size();
}

int SimpleTrackValidation::trackNHits(const reco::TrackSoAConstView& tracks, int i) const {
  auto [start, end] = trackHitRange(tracks, i);
  return static_cast<int>(end - start);
}

std::pair<uint32_t, uint32_t> SimpleTrackValidation::trackHitRange(
    const reco::TrackSoAConstView& tracks, int i) const {
  auto start = (i == 0) ? 0u : tracks[i - 1].hitOffsets();
  auto end   = tracks[i].hitOffsets();
  return {start, end};
}

// =======================================================================================
// produce()
// =======================================================================================

void SimpleTrackValidation::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  auto const& tracksHost   = iEvent.get(tokenTracks_);
  auto const& hitsHost     = iEvent.get(tokenHits_);
  auto const& vertices     = iEvent.get(tokenVertex_); // currently unused
  auto const& simpleParts  = iEvent.get(tSimpleParticles_);
  auto const& hitPartMap      = iEvent.get(tokenHitMap_); 

  auto tracksView = tracksHost.view<reco::TrackSoA>();
  auto hitsView   = tracksHost.view<reco::TrackHitSoA>(); 
  auto hitMapView = hitPartMap.view();

  auto nParticles = simpleParts.view().metadata().size();
  // total truth particles (after selection)
  for (int j = 0; j < nParticles; ++j) {
    if (!isGoodParticle(simpleParts.view(), j)) continue;
    ++totalParticles_;
  }

  // Denominator for efficiencies: loop over all selected particles
  for (int ib = 0; ib < static_cast<int>(binsPt_.size()) - 1; ++ib) {
    const double ptMin   = binsPt_[ib];
    const double ptMax   = binsPt_[ib + 1];
    const double etaMin  = binsEtaPhi_[ib];
    const double etaMax  = binsEtaPhi_[ib + 1];
    const double phiMin  = binsEtaPhi_[ib];
    const double phiMax  = binsEtaPhi_[ib + 1];

    for (int j = 0; j < nParticles; ++j) {
      if (!isGoodParticle(simpleParts.view(), j)) continue;

      const float pt  = simpleParts.view().pt(j);
      const float eta = simpleParts.view().eta(j);
      const float phi = simpleParts.view().phi(j);

      if (pt  > ptMin  && pt  < ptMax)   effPtDen_[ib]  += 1.0;
      if (eta > etaMin && eta < etaMax)  effEtaDen_[ib] += 1.0;
      if (phi > phiMin && phi < phiMax)  effPhiDen_[ib] += 1.0;
    }
  }

  uint32_t lastMatchedPartInd = std::numeric_limits<uint32_t>::max();

  std::map<uint32_t,uint32_t> simToRecoMap;
  // Loop over tracks
  const auto nTrk = nTracks(tracksView);
  for (int i = 0; i < static_cast<int>(nTrk); ++i) {

    const int nHitsTrk = trackNHits(tracksView, i);
    if (nHitsTrk < int(minHits_))
      continue;

    // quality cut: BAD / DUP are rejected
    const auto q = tracksView[i].quality();
    if (q == pixelTrack::Quality::bad || q == pixelTrack::Quality::dup)
      continue;

    const float trkPt  = tracksView[i].pt();
    const float trkEta = tracksView[i].eta();
    const float trkPhi = reco::phi(tracksView, i);   // helper from TracksSoA
    const float trkTip = reco::tip(tracksView, i);
    const float trkZip = reco::zip(tracksView, i);

    // fake-rate denominators
    for (int ib = 0; ib < static_cast<int>(binsPt_.size()) - 1; ++ib) {
      const double ptMin   = binsPt_[ib];
      const double ptMax   = binsPt_[ib + 1];
      const double etaMin  = binsEtaPhi_[ib];
      const double etaMax  = binsEtaPhi_[ib + 1];
      const double phiMin  = binsEtaPhi_[ib];
      const double phiMax  = binsEtaPhi_[ib + 1];

      if (trkPt  > ptMin  && trkPt  < ptMax)   fakePtDen_[ib]  += 1.0;
      if (trkEta > etaMin && trkEta < etaMax)  fakeEtaDen_[ib] += 1.0;
      if (trkPhi > phiMin && trkPhi < phiMax)  fakePhiDen_[ib] += 1.0;
    }

    std::vector<uint32_t> partIndices;
    auto [hitStart, hitEnd] = trackHitRange(tracksView, i);

    for (uint32_t ih = hitStart; ih < hitEnd; ++ih) {
      uint32_t partIndex = hitMapView[ih].id(); 
      partIndices.push_back(partIndex);
    }

    if (partIndices.empty())
      continue;

    auto [occurrences, bestPartInd] = getMostRepeatingPart(partIndices);

    bool isDuplicate = false;
    if (simToRecoMap.find(bestPartInd) != simToRecoMap.end())
    {
      ++duplicateTracks_;
      isDuplicate = true;
    }
    /// TODO: implement a duplicate check
    // // avoid counting multiple tracks mapped to the same particle in a row
    // if (bestPartInd == lastMatchedPartInd)
    //   continue;
    // lastMatchedPartInd = bestPartInd;

    constexpr auto cutForTriplets = 0.6667;

    const double purity = static_cast<double>(occurrences) / static_cast<double>(nHitsTrk);
    bool good = (purity >= purity_) or bestPartInd > 0; // (nHits < 4 and purity >= cutForTriplets) and 
    // 2/3 hits for triplets + bestPartInd < 0 for noise
    if (good) { 
      

      ++trueTracks_;

      if (!isGoodParticle(simpleParts.view(), bestPartInd))
        continue;

      const float partPt  = simpleParts.view().pt(bestPartInd);
      const float partEta = simpleParts.view().eta(bestPartInd);
      const float partPhi = simpleParts.view().phi(bestPartInd);
      const float partPx  = simpleParts.view().px(bestPartInd);
      const float partPy  = simpleParts.view().py(bestPartInd);
      const float partPz  = simpleParts.view().pz(bestPartInd);
      const float partVx  = simpleParts.view().vx(bestPartInd);
      const float partVy  = simpleParts.view().vy(bestPartInd);

      // d0 and dz from particle
      const double partD0 = (-partVx * partPy + partVy * partPx) / partPt;
      const double partDZ = (-(partVx * partPx + partVy * partPy) / partPt) * (partPz / partPt);

      for (int ib = 0; ib < static_cast<int>(binsPt_.size()) - 1; ++ib) {
        const double ptMin   = binsPt_[ib];
        const double ptMax   = binsPt_[ib + 1];
        const double etaMin  = binsEtaPhi_[ib];
        const double etaMax  = binsEtaPhi_[ib + 1];
        const double phiMin  = binsEtaPhi_[ib];
        const double phiMax  = binsEtaPhi_[ib + 1];

        if (partPt > ptMin && partPt < ptMax)
          effPtNum_[ib] += 1.0;

        if (partEta > etaMin && partEta < etaMax) {
          effEtaNum_[ib] += 1.0;

          // resolutions in bins of eta (as before)
          resPt_[ib].push_back(trkPt  - partPt);
          resEta_[ib].push_back(trkEta - partEta);
          resPhi_[ib].push_back(deltaPhi(trkPhi, partPhi));
          resD0_[ib].push_back(trkTip - partD0);
          resDZ_[ib].push_back(trkZip - partDZ);
        }

        if (partPhi > phiMin && partPhi < phiMax)
          effPhiNum_[ib] += 1.0;
      }
    } else {
      // fake track (not pure enough)
      for (int ib = 0; ib < static_cast<int>(binsPt_.size()) - 1; ++ib) {
        const double ptMin   = binsPt_[ib];
        const double ptMax   = binsPt_[ib + 1];
        const double etaMin  = binsEtaPhi_[ib];
        const double etaMax  = binsEtaPhi_[ib + 1];
        const double phiMin  = binsEtaPhi_[ib];
        const double phiMax  = binsEtaPhi_[ib + 1];

        if (trkPt  > ptMin  && trkPt  < ptMax)   dupPtNum_[ib]  += 1.0;
        if (trkEta > etaMin && trkEta < etaMax)  dupEtaNum_[ib] += 1.0;
        if (trkPhi > phiMin && trkPhi < phiMax)  dupPhiNum_[ib] += 1.0;
      }
      ++fakeTracks_;
    }

    if (isDuplicate)
    {
      for (int ib = 0; ib < static_cast<int>(binsPt_.size()) - 1; ++ib) {
        const double ptMin   = binsPt_[ib];
        const double ptMax   = binsPt_[ib + 1];
        const double etaMin  = binsEtaPhi_[ib];
        const double etaMax  = binsEtaPhi_[ib + 1];
        const double phiMin  = binsEtaPhi_[ib];
        const double phiMax  = binsEtaPhi_[ib + 1];

        if (trkPt  > ptMin  && trkPt  < ptMax)   fakePtNum_[ib]  += 1.0;
        if (trkEta > etaMin && trkEta < etaMax)  fakeEtaNum_[ib] += 1.0;
        if (trkPhi > phiMin && trkPhi < phiMax)  fakePhiNum_[ib] += 1.0;
      }
    }
  }
}

void SimpleTrackValidation::endJob() {
  // Fill efficiencies & fake rates
  auto fillEffOrFake = [&](const std::vector<double>& num,
                           const std::vector<double>& den,
                           const std::vector<double>& bins,
                           kType kind,
                           kCoord coord) {
    for (int i = 0; i < static_cast<int>(num.size()); ++i) {
      const double n = num[i];
      const double d = den[i];
      const double bc = 0.5 * (bins[i] + bins[i + 1]);
      const double v  = (d > 0.0) ? (n / d) : 0.0;

      kind_      = static_cast<int>(kind);
      coord_     = static_cast<int>(coord);
      binCenter_ = static_cast<float>(bc);
      value_     = static_cast<float>(v);
      num_       = static_cast<float>(n);
      den_       = static_cast<float>(d);

      ttree_->Fill();
    }
  };

  // efficiencies
  fillEffOrFake(effPtNum_,  effPtDen_,  binsPt_,     kType::Efficiency, kCoord::Pt);
  fillEffOrFake(effEtaNum_, effEtaDen_, binsEtaPhi_, kType::Efficiency, kCoord::Eta);
  fillEffOrFake(effPhiNum_, effPhiDen_, binsEtaPhi_, kType::Efficiency, kCoord::Phi);

  // fake rates
  fillEffOrFake(fakePtNum_,  fakePtDen_,  binsPt_,     kType::FakeRate, kCoord::Pt);
  fillEffOrFake(fakeEtaNum_, fakeEtaDen_, binsEtaPhi_, kType::FakeRate, kCoord::Eta);
  fillEffOrFake(fakePhiNum_, fakePhiDen_, binsEtaPhi_, kType::FakeRate, kCoord::Phi);

  // dup rates
  fillEffOrFake(dupPtNum_,  fakePtDen_,  binsPt_,     kType::Duplicates, kCoord::Pt);
  fillEffOrFake(dupEtaNum_, fakeEtaDen_, binsEtaPhi_, kType::Duplicates, kCoord::Eta);
  fillEffOrFake(dupPhiNum_, fakePhiDen_, binsEtaPhi_, kType::Duplicates, kCoord::Phi);

  // resolutions
  auto fillResolution = [&](std::map<int, std::vector<double>>& res,
                            const std::vector<double>& bins,
                            kCoord coord) {
    for (int ib = 0; ib < static_cast<int>(bins.size()) - 1; ++ib) {
      const auto it = res.find(ib);
      const double sd = (it != res.end()) ? stdDev(it->second) : 0.0;
      const double bc = 0.5 * (bins[ib] + bins[ib + 1]);

      kind_      = static_cast<int>(kType::Resolution);
      coord_     = static_cast<int>(coord);
      binCenter_ = static_cast<float>(bc);
      value_     = static_cast<float>(sd);
      num_       = 0.f;
      den_       = 0.f;

      ttree_->Fill();
    }
  };

  fillResolution(resPt_,  binsPt_,     kCoord::Pt);
  fillResolution(resEta_, binsEtaPhi_, kCoord::Eta);
  fillResolution(resPhi_, binsEtaPhi_, kCoord::Phi);
  fillResolution(resD0_,  binsEtaPhi_, kCoord::D0);
  fillResolution(resDZ_,  binsEtaPhi_, kCoord::DZ);

  std::cout << "=====================================\n";
  std::cout << "Matched   tracks: " << trueTracks_ << "\n";
  std::cout << "Fake      tracks: " << fakeTracks_ << "\n";
  std::cout << "Duplicate tracks: " << duplicateTracks_ << "\n";
  std::cout << "Total selected particles: " << totalParticles_ << "\n";
  std::cout << "Validation written to: " << outputFileName_ << "\n";

  outFile_->cd();
  ttree_->Write();
  outFile_->Close();  
}

// =======================================================================================
// Small helpers (unchanged logic)
// =======================================================================================

std::pair<int, uint32_t> SimpleTrackValidation::getMostRepeatingPart(
    const std::vector<uint32_t>& vec) const {

  std::unordered_map<uint32_t, int> freq;
  uint32_t mostFrequent = vec.front();
  int maxCount = 0;

  for (auto x : vec) {
    auto& c = freq[x];
    ++c;
    if (c > maxCount) {
      maxCount = c;
      mostFrequent = x;
    }
  }

  return std::make_pair(maxCount, mostFrequent);
}

double SimpleTrackValidation::deltaPhi(double phi1, double phi2) const {
  double o2pi = 1. / (2. * M_PI);
  if (std::abs(phi1 - phi2) <= double(M_PI))
    return (phi1 - phi2);
  double n = std::round((phi1 - phi2) * o2pi);
  return (phi1 - phi2) - n * double(2. * M_PI);
}

void SimpleTrackValidation::logSpace(unsigned n, double a, double b, std::vector<double>& bins) const {
  const double step = (b - a) / static_cast<double>(n);
  bins.clear();
  for (double x = a; x < b + 0.5 * step; x += step)
    bins.push_back(std::pow(10.0, x));
}

void SimpleTrackValidation::linSpace(unsigned n, double a, double b, std::vector<double>& bins) const {
  const double step = (b - a) / static_cast<double>(n);
  bins.clear();
  for (double x = a; x < b + 0.5 * step; x += step)
    bins.push_back(x);
}

double SimpleTrackValidation::stdDev(const std::vector<double>& vec) const {
  const int n = vec.size();
  if (n == 0)
    return 0.0;

  double sum = 0.0;
  for (double v : vec)
    sum += v;
  const double avg = sum / n;

  double sumSq = 0.0;
  for (double v : vec) {
    const double d = v - avg;
    sumSq += d * d;
  }

  return std::sqrt(sumSq / n);
}

DEFINE_FWK_MODULE(SimpleTrackValidation);
