#include <filesystem>
#include <fstream>
#include <ios>
#include <iostream>  // for std::cout, std::endl
#include <memory>
#include <utility>

#include "DataFormats/BeamSpotPOD.h"
#include "Framework/ESPluginFactory.h"
#include "Framework/ESProducer.h"
#include "Framework/EventSetup.h"

#define GPU_DEBUG

class BeamSpotESProducer : public edm::ESProducer {
public:
  explicit BeamSpotESProducer(std::filesystem::path const& datadir) : data_(datadir) {
#ifdef GPU_DEBUG
    std::cout << "[GPU_DEBUG] BeamSpotESProducer constructed with data path: "
              << data_ << std::endl;
#endif
  }

  void produce(edm::EventSetup& eventSetup);

private:
  std::filesystem::path data_;
};

void BeamSpotESProducer::produce(edm::EventSetup& eventSetup) {
#ifdef GPU_DEBUG
  std::cout << "[GPU_DEBUG] BeamSpotESProducer::produce() called" << std::endl;
#endif

  auto bs = std::make_unique<BeamSpotPOD>();

  const auto filePath = data_ / "beamspot.bin";
#ifdef GPU_DEBUG
  std::cout << "[GPU_DEBUG] Attempting to open file: " << filePath << std::endl;
#endif

  try {
    std::ifstream in(filePath, std::ios::binary);
    in.exceptions(std::ifstream::badbit | std::ifstream::failbit | std::ifstream::eofbit);

    in.read(reinterpret_cast<char*>(bs.get()), sizeof(BeamSpotPOD));

#ifdef GPU_DEBUG
    std::cout << "[GPU_DEBUG] Successfully read BeamSpotPOD (" << sizeof(BeamSpotPOD)
              << " bytes) from " << filePath << std::endl;
    std::cout << "[GPU_DEBUG] BeamSpot values: "
              << "x=" << bs->x << "  y=" << bs->y << "  z=" << bs->z
              << "  sigmaZ=" << bs->sigmaZ << std::endl;
#endif

  } catch (std::exception const& e) {
#ifdef GPU_DEBUG
    std::cout << "[GPU_DEBUG] ERROR reading BeamSpot file: " << e.what() << std::endl;
#endif
    throw;  // rethrow to preserve framework error handling
  }

  eventSetup.put(std::move(bs));

#ifdef GPU_DEBUG
  std::cout << "[GPU_DEBUG] BeamSpotESProducer::produce() completed successfully"
            << std::endl;
#endif
}

DEFINE_FWK_EVENTSETUP_MODULE(BeamSpotESProducer);
