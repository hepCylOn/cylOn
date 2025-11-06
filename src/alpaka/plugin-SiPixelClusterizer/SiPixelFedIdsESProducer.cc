#include <filesystem>
#include <fstream>
#include <ios>
#include <memory>
#include <utility>
#include <vector>

#include "CondFormats/SiPixelFedIds.h"
#include "Framework/ESPluginFactory.h"
#include "Framework/ESProducer.h"
#include "Framework/EventSetup.h"
#include "Framework/ConfigRegistry.h"

class SiPixelFedIdsESProducer : public edm::ESProducer {
public:
  explicit SiPixelFedIdsESProducer(edm::Config const& cfg) : data_(static_cast<std::string>(cfg.value("data", defaultPath_))) {}
  void produce(edm::EventSetup& eventSetup);

private:
  std::filesystem::path data_;
  std::filesystem::path defaultPath_ = "data/fedIds.bin";
};

void SiPixelFedIdsESProducer::produce(edm::EventSetup& eventSetup) {
  std::ifstream in(data_, std::ios::binary);
  in.exceptions(std::ifstream::badbit | std::ifstream::failbit | std::ifstream::eofbit);
  unsigned int nfeds;
  in.read(reinterpret_cast<char*>(&nfeds), sizeof(unsigned));
  std::vector<unsigned int> fedIds(nfeds);
  in.read(reinterpret_cast<char*>(fedIds.data()), sizeof(unsigned int) * nfeds);
  eventSetup.put(std::make_unique<SiPixelFedIds>(std::move(fedIds)));
}

DEFINE_FWK_EVENTSETUP_MODULE(SiPixelFedIdsESProducer);
