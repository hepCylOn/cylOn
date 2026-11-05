#ifndef CondFormats_alpaka_CAGeometry_h
#define CondFormats_alpaka_CAGeometry_h

#include <utility>

#include "AlpakaCore/ESProduct.h"
#include "AlpakaCore/config.h"
#include "AlpakaCore/memory.h"
#include "CondFormats/CAGeometrySoA.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  class CAGeometry {
  public:
    CAGeometry(std::string const &path)
    {
  
      // std::cout << "Building CA Geometry" << std::endl;
      std::ifstream in(path, std::ios::binary);
      if (!in) {
          throw std::runtime_error("Failed to open file: " + path);
      }

      // Reading sizes as uint16_t, to match CMSPhase1Geometry.bin
      uint16_t nModules, nLayers, nPairs;
      in.read(reinterpret_cast<char*>(&nModules), sizeof(nModules));
      in.read(reinterpret_cast<char*>(&nLayers),  sizeof(nLayers));
      in.read(reinterpret_cast<char*>(&nPairs),   sizeof(nPairs));

      // Resize vectors
      m_caModules.resize(nModules);
      m_caLayers.resize(nLayers + 1);  // layer vector has nLayers+1 elements
      m_caPairs.resize(nPairs);

      // to be used on host
      m_sizes.nModules = nModules;
      m_sizes.nLayers = nLayers;
      m_sizes.nPairs = nPairs;
      // Read modules
      in.read(reinterpret_cast<char*>(m_caModules.data()), nModules * sizeof(caGeometry::CAModule));

      // Read layers (nLayers + 1)
      in.read(reinterpret_cast<char*>(m_caLayers.data()), (nLayers + 1) * sizeof(caGeometry::CALayer));

      // Read pairs
      in.read(reinterpret_cast<char*>(m_caPairs.data()), nPairs * sizeof(caGeometry::CAPair));

      if (!in) {
          throw std::runtime_error("Error reading binary file: " + path);
      }
    }

    ~CAGeometry() = default;

    caGeometry::CASizes const& sizes() const { return m_sizes;}
    // The return value can only be used safely in kernels launched on
    // the same cudaStream, or after cudaStreamSynchronize.
    const caGeometry::CAGeometrySoA *getGPUProductAsync(Queue &queue) const {
      const auto &data = gpuData_.dataForDeviceAsync(queue, [this](Queue &queue) {

        unsigned int nModules = m_caModules.size();
        unsigned int nLayers  = m_caLayers.size();
        unsigned int nPairs   = m_caPairs.size();

        GPUData gpuData(queue, nModules, nLayers, nPairs);

        auto modules_h = cms::alpakatools::make_host_view(m_caModules.data(), nModules);
        alpaka::memcpy(queue, gpuData.d_caModules, modules_h);
        gpuData.h_geometry->m_modules = gpuData.d_caModules.data();

        auto layers_h = cms::alpakatools::make_host_view(m_caLayers.data(), nLayers);
        alpaka::memcpy(queue, gpuData.d_caLayers, layers_h);
        gpuData.h_geometry->m_layers = gpuData.d_caLayers.data();

        auto pairs_h = cms::alpakatools::make_host_view(m_caPairs.data(), nPairs);
        alpaka::memcpy(queue, gpuData.d_caPairs, pairs_h);
        gpuData.h_geometry->m_pairs = gpuData.d_caPairs.data();
        
        gpuData.h_geometry->m_nModules = nModules;
        gpuData.h_geometry->m_nLayers = nLayers;
        gpuData.h_geometry->m_nPairs = nPairs;
        
        alpaka::memcpy(queue, gpuData.d_geometry, gpuData.h_geometry);
        
        return gpuData;
      });
      return data.d_geometry.data();
    }

  private:
    // allocate it with posix malloc to be compatible with cpu wf
    std::vector<caGeometry::CAModule> m_caModules;
    std::vector<caGeometry::CALayer> m_caLayers;
    std::vector<caGeometry::CAPair> m_caPairs;
    caGeometry::CASizes m_sizes;

    struct GPUData {
      // not needed if not used on CPU...
    public:
      GPUData() = delete;
      GPUData(Queue &queue, unsigned int nModules, unsigned int mLayers , unsigned int nPairs)
          : h_geometry{cms::alpakatools::make_host_buffer<caGeometry::CAGeometrySoA>(queue)},
            d_geometry{cms::alpakatools::make_device_buffer<caGeometry::CAGeometrySoA>(queue)},
            d_caModules{cms::alpakatools::make_device_buffer<caGeometry::CAModule[]>(queue, nModules)},
            d_caLayers{cms::alpakatools::make_device_buffer<caGeometry::CALayer[]>(queue, mLayers)},
            d_caPairs{cms::alpakatools::make_device_buffer<caGeometry::CAPair[]>(queue, nPairs)} {};
      ~GPUData() = default;

    public:
      cms::alpakatools::host_buffer<caGeometry::CAGeometrySoA> h_geometry;
      cms::alpakatools::device_buffer<Device, caGeometry::CAGeometrySoA> d_geometry;  // copy of the above on the Device
      
      cms::alpakatools::device_buffer<Device, caGeometry::CAModule[]> d_caModules;
      cms::alpakatools::device_buffer<Device, caGeometry::CALayer[]> d_caLayers;
      cms::alpakatools::device_buffer<Device, caGeometry::CAPair[]> d_caPairs;
    };

    cms::alpakatools::ESProduct<Queue, GPUData> gpuData_;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif  // CondFormats_alpaka_CAGeometry_h
