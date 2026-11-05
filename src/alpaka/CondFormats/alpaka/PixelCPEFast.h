#ifndef CondFormats_alpaka_PixelCPEFast_h
#define CondFormats_alpaka_PixelCPEFast_h

#include <utility>
#include <iostream>
#include <iomanip>

#include "AlpakaCore/ESProduct.h"
#include "AlpakaCore/config.h"
#include "AlpakaCore/memory.h"
#include "CondFormats/pixelCPEforDevice.h"

#define GPU_DEBUG
// #define DUMPDETS

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  template <typename TrackerTraits>
  class PixelCPEFast {
    using ParamsOnDevice = pixelCPEforDevice::ParamsOnDevice;

  public:
    PixelCPEFast(std::string const &path)
        : m_commonParamsGPU(cms::alpakatools::make_host_buffer<pixelCPEforDevice::CommonParams, Platform>()) {
#ifdef GPU_DEBUG
      std::cout << "[GPU_DEBUG] Loading PixelCPEFast<" << TrackerTraits::nameModifier
                << "> from file: " << path << std::endl;
#endif

      std::ifstream in(path, std::ios::binary);
      in.exceptions(std::ifstream::badbit | std::ifstream::failbit | std::ifstream::eofbit);

      // --- CommonParams
      in.read(reinterpret_cast<char *>(m_commonParamsGPU.data()), sizeof(pixelCPEforDevice::CommonParams));
#ifdef GPU_DEBUG
      auto const *cp = m_commonParamsGPU.data();
      std::cout << "[GPU_DEBUG] CommonParams loaded:"
                << " theThicknessB=" << cp->theThicknessB
                << " theThicknessE=" << cp->theThicknessE
                << " maxModuleStride=" << cp->maxModuleStride
                << " nLadders=" << static_cast<int>(cp->numberOfLaddersInBarrel)
                << std::endl;
#endif

      // --- DetParams
      unsigned int ndetParams = 0;
      in.read(reinterpret_cast<char *>(&ndetParams), sizeof(unsigned int));
      m_detParamsGPU.resize(ndetParams);
      in.read(reinterpret_cast<char *>(m_detParamsGPU.data()),
              ndetParams * sizeof(pixelCPEforDevice::DetParams));

#ifdef GPU_DEBUG
      std::cout << "[GPU_DEBUG] Read " << ndetParams << " DetParams entries" << std::endl;
      unsigned int nPrint = std::min<unsigned int>(ndetParams, 3);
      for (unsigned int i = 0; i < nPrint; ++i) {
        const auto &d = m_detParamsGPU[i];
        std::cout << std::fixed << std::setprecision(4);
        std::cout << "  [" << i << "]"
                  << " rawId=" << d.rawId
                  << " layer=" << d.layer
                  << " isBarrel=" << d.isBarrel
                  << " isPosZ=" << d.isPosZ
                  << " shiftX=" << d.shiftX
                  << " shiftY=" << d.shiftY
                  << " chargeWidthX=" << d.chargeWidthX
                  << " chargeWidthY=" << d.chargeWidthY
                  << " pitchX=" << d.thePitchX
                  << " pitchY=" << d.thePitchY
                  << " x0=" << d.x0
                  << " y0=" << d.y0
                  << " z0=" << d.z0
                  << std::endl;
      }
#endif

#ifdef GPU_DEBUG
      std::cout << "[GPU_DEBUG] (no AverageGeometry or LayerGeometry blocks in this simplified format)" << std::endl;
#endif

#ifdef DUMPDETS
      std::ofstream out(std::string(TrackerTraits::cpeModules) + ".bin", std::ios::binary);
      if (!out) throw std::runtime_error("Cannot open file for writing");
      std::cout << "=== Writing " << ndetParams << " modules" << std::endl;
      for (auto const &v : m_detParamsGPU) {
        auto f = v.frame;
        out.write(reinterpret_cast<const char *>(&f), sizeof(f));
      }
      out.close();
#endif
    }

    ~PixelCPEFast() = default;

    const ParamsOnDevice *getGPUProductAsync(Queue &queue) const {
      const auto &data = gpuData_.dataForDeviceAsync(queue, [this](Queue &queue) {
        unsigned int ndetParams = m_detParamsGPU.size();
        GPUData gpuData(queue, ndetParams);

#ifdef GPU_DEBUG
        std::cout << "[GPU_DEBUG] Copying PixelCPEFast<" << TrackerTraits::nameModifier
                  << "> data to GPU..." << std::endl;
#endif
        alpaka::memcpy(queue, gpuData.d_commonParams, m_commonParamsGPU);
        gpuData.h_paramsOnGPU->m_commonParams = gpuData.d_commonParams.data();

        auto detParams_h = cms::alpakatools::make_host_view(m_detParamsGPU.data(), ndetParams);
        alpaka::memcpy(queue, gpuData.d_detParams, detParams_h);
        gpuData.h_paramsOnGPU->m_detParams = gpuData.d_detParams.data();

        alpaka::memcpy(queue, gpuData.d_paramsOnGPU, gpuData.h_paramsOnGPU);
#ifdef GPU_DEBUG
        std::cout << "[GPU_DEBUG] GPU buffers populated for " << ndetParams << " modules" << std::endl;
#endif
        return gpuData;
      });
      return data.d_paramsOnGPU.data();
    }

  private:
    std::vector<pixelCPEforDevice::DetParams> m_detParamsGPU;
    cms::alpakatools::host_buffer<pixelCPEforDevice::CommonParams> m_commonParamsGPU;

    struct GPUData {
      GPUData() = delete;
      GPUData(Queue &queue, unsigned int ndetParams)
          : h_paramsOnGPU{cms::alpakatools::make_host_buffer<ParamsOnDevice>(queue)},
            d_paramsOnGPU{cms::alpakatools::make_device_buffer<ParamsOnDevice>(queue)},
            d_commonParams{cms::alpakatools::make_device_buffer<pixelCPEforDevice::CommonParams>(queue)},
            d_detParams{cms::alpakatools::make_device_buffer<pixelCPEforDevice::DetParams[]>(queue, ndetParams)} {};
      ~GPUData() = default;

      cms::alpakatools::host_buffer<ParamsOnDevice> h_paramsOnGPU;
      cms::alpakatools::device_buffer<Device, ParamsOnDevice> d_paramsOnGPU;
      cms::alpakatools::device_buffer<Device, pixelCPEforDevice::CommonParams> d_commonParams;
      cms::alpakatools::device_buffer<Device, pixelCPEforDevice::DetParams[]> d_detParams;
    };

    cms::alpakatools::ESProduct<Queue, GPUData> gpuData_;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif  // CondFormats_alpaka_PixelCPEFast_h
