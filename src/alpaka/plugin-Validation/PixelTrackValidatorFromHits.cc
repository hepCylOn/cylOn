#include <iostream>
#include <cmath>
#include <utility>
#include <fstream>
#include <string>
#include <map>

#include "DataFormats/ParticleSimpleSoA.h"
#include "AlpakaDataFormats/PixelTrackHost.h"
#include "Framework/EDProducer.h"
#include "Framework/Event.h"
#include "Framework/EventSetup.h"
#include "Framework/PluginFactory.h"

class ParticleFromSimple : public edm::EDProducer {
public:
  explicit ParticleFromSimple(edm::ProductRegistry& reg);
  ~ParticleFromSimple() override = default;

private:
  void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override;
  void endJob() override;
  std::pair<int,uint32_t> getMostRepeatingPart(std::vector<uint32_t> vec);
  double deltaPhi(double phi1, double phi2);
  void logSpace (const unsigned, const double, const double, std::vector<double> &) const;
  void linSpace (const unsigned, const double, const double, std::vector<double> &) const;
  void writeEffAndFake(std::vector<double>&, std::vector<double>&, std::ofstream &, std::vector<double>&);
  void writeResolution(std::map<int,std::vector<double>>&, std::ofstream&, std::vector<double>&, std::vector<double>&);
  double stdDev(const std::vector<double>&) const;
    
  edm::EDGetTokenT<ParticleSimpleSoA> tSimpleParticles_;
  edm::EDGetTokenT<PixelTrackHost> tokenTrack_;

  std::ofstream file;

  int totalParticles;
  int passedTracks;

  int nBins;
  std::vector<double> binsPt;
  std::vector<double> binsEtaPhi;

  std::vector<double> effPtDen;
  std::vector<double> effPtNum;
  std::vector<double> effEtaDen;
  std::vector<double> effEtaNum;
  std::vector<double> effPhiDen;
  std::vector<double> effPhiNum;

  std::vector<double> fakePtDen;
  std::vector<double> fakePtNum;
  std::vector<double> fakeEtaDen;
  std::vector<double> fakeEtaNum;
  std::vector<double> fakePhiDen;
  std::vector<double> fakePhiNum;

  std::map<int,std::vector<double>> resPt;
  std::map<int,std::vector<double>> resEta;
  std::map<int,std::vector<double>> resPhi;
  std::map<int,std::vector<double>> resD0;
  std::map<int,std::vector<double>> resDZ;

};

ParticleFromSimple::ParticleFromSimple(edm::ProductRegistry& reg)
    : tSimpleParticles_(reg.consumes<ParticleSimpleSoA>()),
      tokenTrack_(reg.consumes<PixelTrackHost>()) {

        totalParticles = 0;
        passedTracks = 0;

        nBins = 40;
        
        logSpace (nBins, -0.1, 2.0, binsPt);
        linSpace (nBins, -4.0, 4.0, binsEtaPhi);

        effPtDen.resize(nBins);
        effPtNum.resize(nBins);
        effEtaDen.resize(nBins);
        effEtaNum.resize(nBins);
        effPhiDen.resize(nBins);
        effPhiNum.resize(nBins);

        fakePtDen.resize(nBins);
        fakePtNum.resize(nBins);
        fakeEtaDen.resize(nBins);
        fakeEtaNum.resize(nBins);
        fakePhiDen.resize(nBins);
        fakePhiNum.resize(nBins);

      }

void ParticleFromSimple::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  auto const& tracks = iEvent.get(tokenTrack_);
  auto const& simpleParticles = iEvent.get(tSimpleParticles_);

  totalParticles = totalParticles + simpleParticles.nParticles();

  for(int i = 0; i < int(binsPt.size()); ++i){
    for(uint32_t j = 0; j < simpleParticles.nParticles(); ++j){
      if (simpleParticles.pt(j) > binsPt[i] and simpleParticles.pt(j) < binsPt[i+1]) effPtDen[i] = effPtDen[i] + 1.0;
      if (simpleParticles.eta(j) > binsEtaPhi[i] and simpleParticles.eta(j) < binsEtaPhi[i+1]) effEtaDen[i] = effEtaDen[i] + 1.0;
      if (simpleParticles.phi(j) > binsEtaPhi[i] and simpleParticles.phi(j) < binsEtaPhi[i+1]) effPhiDen[i] = effPhiDen[i] + 1.0;
    }
  }

  int trackHitIncr = 0;
  for(uint32_t i = 0; i < tracks->m_nTracks; ++i){
    if(tracks->nHits(i) < 4) continue;
    if(tracks->quality(i) == trackQuality::bad or tracks->quality(i) == trackQuality::dup) continue;

    for(int j = 0; j < int(binsPt.size()); ++j){
      if (tracks->pt(i) > binsPt[j] and tracks->pt(i) < binsPt[j+1]) fakePtDen[j] = fakePtDen[j] + 1.0;
      if (tracks->eta(i) > binsEtaPhi[j] and tracks->eta(i) < binsEtaPhi[j+1]) fakeEtaDen[j] = fakeEtaDen[j] + 1.0;
      if (tracks->phi(i) > binsEtaPhi[j] and tracks->phi(i) < binsEtaPhi[j+1]) fakePhiDen[j] = fakePhiDen[j] + 1.0;
    }

    std::vector<uint32_t> partIndices;
    auto start = (i == 0) ? 0 : tracks->partIndices.off[i];
    auto end = tracks->partIndices.off[i + 1];
    for (auto iHit = start; iHit < end; ++iHit){
      partIndices.push_back(tracks->partIndices.bins[trackHitIncr]);
      ++trackHitIncr;
    }

    if (partIndices.size() <= 0) continue;
    std::pair<int,uint32_t> repeatingPart = getMostRepeatingPart(partIndices);
    if (repeatingPart.first/tracks->nHits(i) > 0.75){
      passedTracks = passedTracks + 1;
      auto it = std::find(simpleParticles.partIndVector().begin(), simpleParticles.partIndVector().end(), repeatingPart.second);
      uint32_t particle = 0;
      if (it != simpleParticles.partIndVector().end()) {
        size_t index = std::distance(simpleParticles.partIndVector().begin(), it);
        particle = index;
      }
      for(int j = 0; j < int(binsPt.size()); ++j){
        if (simpleParticles.pt(particle) > binsPt[j] and simpleParticles.pt(particle) < binsPt[j+1]) effPtNum[j] = effPtNum[j] + 1.0;
        if (simpleParticles.eta(particle) > binsEtaPhi[j] and simpleParticles.eta(particle) < binsEtaPhi[j+1]) {

          effEtaNum[j] = effEtaNum[j] + 1.0;

          // Calculating resolutions
          resPt[j].push_back(tracks->pt(i) - simpleParticles.pt(particle));
          resEta[j].push_back(tracks->eta(i) - simpleParticles.eta(particle));
          resPhi[j].push_back(tracks->phi(i) - simpleParticles.phi(particle));
          // (-vx*py + vy*px)/pt for D0
          // (-(vx*px + vy*py)/pt) * (pz/pt) for DZ
          double partD0 = (-simpleParticles.vx(particle)*simpleParticles.py(particle) + simpleParticles.vy(particle)*simpleParticles.px(particle))/simpleParticles.pt(particle);
          double partDZ = (-(simpleParticles.vx(particle)*simpleParticles.px(particle) + simpleParticles.vy(particle)*simpleParticles.py(particle))/simpleParticles.pt(particle)) * (simpleParticles.pz(particle)/simpleParticles.pt(particle));
          resD0[j].push_back(tracks->tip(i) - partD0);
          resDZ[j].push_back(tracks->zip(i) - partDZ);

        }
        if (simpleParticles.phi(particle) > binsEtaPhi[j] and simpleParticles.phi(particle) < binsEtaPhi[j+1]) effPhiNum[j] = effPhiNum[j] + 1.0;

      }
    }
    else{
      for(int j = 0; j < int(binsPt.size()); ++j){
        if (tracks->pt(i) > binsPt[j] and tracks->pt(i) < binsPt[j+1]) fakePtNum[j] = fakePtNum[j] + 1.0;
        if (tracks->eta(i) > binsEtaPhi[j] and tracks->eta(i) < binsEtaPhi[j+1]) fakeEtaNum[j] = fakeEtaNum[j] + 1.0;
        if (tracks->phi(i) > binsEtaPhi[j] and tracks->phi(i) < binsEtaPhi[j+1]) fakePhiNum[j] = fakePhiNum[j] + 1.0;
      }
    }
  }
}

void ParticleFromSimple::endJob() {
  file.open("/data/user/borzari/cmssw/pixeltrack-standalone/output.txt");
  if (file.is_open()) {

    writeEffAndFake(effPtNum, effPtDen, file, binsPt);
    writeEffAndFake(effEtaNum, effEtaDen, file, binsEtaPhi);
    writeEffAndFake(effPhiNum, effPhiDen, file, binsEtaPhi);
    writeEffAndFake(fakePtNum, fakePtDen, file, binsPt);
    writeEffAndFake(fakeEtaNum, fakeEtaDen, file, binsEtaPhi);
    writeEffAndFake(fakePhiNum, fakePhiDen, file, binsEtaPhi);

    writeResolution(resPt, file, binsPt, binsEtaPhi);
    writeResolution(resEta, file, binsEtaPhi, binsEtaPhi);
    writeResolution(resPhi, file, binsEtaPhi, binsEtaPhi);
    writeResolution(resD0, file, binsEtaPhi, binsEtaPhi);
    writeResolution(resDZ, file, binsEtaPhi, binsEtaPhi);

  }
  std::cout << "=====================================" << std::endl;
  std::cout << "The number of matched tracks is: " << passedTracks << std::endl;
  std::cout << "The total number of particles is: " << totalParticles << std::endl;
  std::cout << "Validation ended!!" << std::endl;
  file.close();  // fecha o arquivo
}

std::pair<int,uint32_t> ParticleFromSimple::getMostRepeatingPart( std::vector<uint32_t> vec) {

  std::unordered_map<int, int> freq;
  int mostFrequent = vec[0];
  int maxCount = 0;

  for (int x : vec) {
      freq[x]++;
      if (freq[x] > maxCount) {
          maxCount = freq[x];
          mostFrequent = x;
      }
  }

  return std::make_pair(maxCount,mostFrequent);
}

double ParticleFromSimple::deltaPhi(double phi1, double phi2) {
  double o2pi = 1. / (2. * M_PI);
  if (std::abs((phi1 - phi2)) <= double(M_PI))
    return (phi1 - phi2);
  double n = std::round((phi1 - phi2) * o2pi);
  return (phi1 - phi2) - n * double(2. * M_PI);
}

void ParticleFromSimple::logSpace (const unsigned n, const double a, const double b, std::vector<double> &bins) const
{
  double step = (b - a) / ((double) n);

  bins.clear ();
  for (double i = a; i < b + 0.5 * step; i += step)
    bins.push_back (pow (10.0, i));
}

void ParticleFromSimple::linSpace (const unsigned n, const double a, const double b, std::vector<double> &bins) const
{
  double step = (b - a) / ((double) n);

  bins.clear ();
  for (double i = a; i < b + 0.5 * step; i += step)
    bins.push_back (i);
}

void ParticleFromSimple::writeEffAndFake(std::vector<double>& num, std::vector<double>& den, std::ofstream& file, std::vector<double>& bins)
{
  for (int i = 0; i < int(num.size()); ++i){
    if (den[i] == 0) file << den[i] << ",";
    else{
      double x = num[i]/den[i];
      file << x << ",";
    }
  }
  file << std::endl;
  for (int i = 0; i < int(bins.size()) - 1; ++i) file << (bins[i+1] + bins[i])/2.0 << ",";
  file << std::endl;
}

void ParticleFromSimple::writeResolution(std::map<int,std::vector<double>>& res, std::ofstream& file, std::vector<double>& binsVar, std::vector<double>& binsEta)
{
  for(int j = 0; j < int(binsVar.size()) - 1; ++j){
        double stdDevRes = stdDev(res[j]);
        file << stdDevRes << ",";
      }
      file << std::endl;
      for (int i = 0; i < int(binsEta.size()) - 1; ++i) file << (binsEta[i+1] + binsEta[i])/2.0 << ",";
      file << std::endl;
}

double ParticleFromSimple::stdDev(const std::vector<double>& vec) const
{
    int n = vec.size();
    if (n == 0) return 0.0;

    double sum = 0.0;
    for (double v : vec)
        sum += v;
    double avg = sum / n;

    double sumSquare = 0.0;
    for (double v : vec)
        sumSquare += (v - avg) * (v - avg);

    double var = sumSquare / n;

    return std::sqrt(var);
}



DEFINE_FWK_MODULE(ParticleFromSimple);
