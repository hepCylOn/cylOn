#include <iostream>
#include <filesystem>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>

#include "DataFormats/BeamSpotPODInp.h"
#include "AlpakaDataFormats/BeamSpotPOD.h"
#include "Framework/ESPluginFactory.h"
#include "Framework/ESProducer.h"
#include "Framework/EventSetup.h"
#include "Framework/ConfigRegistry.h"

// #define GPU_DEBUG

namespace {
  void readBeamSpot(std::filesystem::path const& path, BeamSpotPOD& beamSpot) {
    auto const extension = path.extension();
    if (extension != ".bin" && extension != ".txt") {
      throw std::runtime_error("Unsupported beam spot file extension: " + extension.string());
    }

    std::ifstream input(path, extension == ".bin" ? std::ios::binary : std::ios::in);
    if (!input) {
      throw std::runtime_error("Cannot open beam spot file: " + path.string());
    }

    BeamSpotPODInp values{};
    if (extension == ".bin") {
      input.read(reinterpret_cast<char*>(&values), sizeof(values));
      if (input.gcount() != static_cast<std::streamsize>(sizeof(values))) {
        throw std::runtime_error("Beam spot binary must contain exactly 11 floats: " + path.string());
      }
      char extra;
      if (input.read(&extra, 1)) {
        throw std::runtime_error("Beam spot binary has trailing data: " + path.string());
      }
    } else {
      if (!(input >> values.x >> values.y >> values.z >> values.sigmaZ >> values.beamWidthX >> values.beamWidthY >>
            values.dxdz >> values.dydz >> values.emittanceX >> values.emittanceY >> values.betaStar)) {
        throw std::runtime_error("Beam spot text must contain 11 whitespace-separated floats: " + path.string());
      }
      std::string extra;
      if (input >> extra) {
        throw std::runtime_error("Beam spot text has extra data after the 11 expected floats: " + path.string());
      }
    }

    beamSpot = BeamSpotPOD{values.x,
                           values.y,
                           values.z,
                           values.sigmaZ,
                           values.beamWidthX,
                           values.beamWidthY,
                           values.dxdz,
                           values.dydz,
                           values.emittanceX,
                           values.emittanceY,
                           values.betaStar};
  }
}  // namespace

class BeamSpotESProducer : public edm::ESProducer {
public:
  explicit BeamSpotESProducer(edm::Config const& config)
      : beamSpot_{0.f, 0.f, 5e-05f, 4.f, 0.0015f, 0.0015f, 0.f, 0.f, 0.f, 0.f, 0.f} {
    auto const data = config.value("data", std::string{});
    if (!data.empty()) {
      readBeamSpot(data, beamSpot_);
#ifdef GPU_DEBUG
      std::cout << "[GPU_DEBUG] BeamSpotESProducer loaded beam spot from " << data << std::endl;
#endif
    } else {
#ifdef GPU_DEBUG
      std::cout << "[GPU_DEBUG] BeamSpotESProducer using default beam spot values" << std::endl;
#endif
    }
  }

  void produce(edm::EventSetup& eventSetup);

private:
  BeamSpotPOD beamSpot_;
};

void BeamSpotESProducer::produce(edm::EventSetup& eventSetup) {
#ifdef GPU_DEBUG
  std::cout << "[GPU_DEBUG] BeamSpotESProducer::produce() called" << std::endl;
#endif

  auto bs = std::make_unique<BeamSpotPOD>(beamSpot_);

#ifdef GPU_DEBUG
    std::cout << "[GPU_DEBUG] BeamSpot values: "
              << "x=" << bs->x << "  y=" << bs->y << "  z=" << bs->z
              << "  sigmaZ=" << bs->sigmaZ << std::endl;
    std::cout << "[GPU_DEBUG] BeamSpot values: "
              << "beamWidthX=" << bs->beamWidthX << "  beamWidthY=" << bs->beamWidthY
              << " dxdz=" << bs->dxdz << "  dydz=" << bs->dydz
              << " emittanceX=" << bs->emittanceX << "  emittanceY=" << bs->emittanceY
              << "  betaStar=" << bs->betaStar << std::endl;
#endif

  eventSetup.put(std::move(bs));

#ifdef GPU_DEBUG
  std::cout << "[GPU_DEBUG] BeamSpotESProducer::produce() completed successfully"
            << std::endl;
#endif
}

DEFINE_FWK_EVENTSETUP_MODULE(BeamSpotESProducer);
