#ifndef AlpakaDataFormats_alpaka_TrackingRecHit2DAlpaka_h
#define AlpakaDataFormats_alpaka_TrackingRecHit2DAlpaka_h

#include <memory>
#include <optional>

#include <alpaka/alpaka.hpp>

#include "AlpakaCore/config.h"
#include "AlpakaCore/memory.h"
#include "AlpakaDataFormats/TrackingRecHit2DSoAView.h"
#include "DataFormats/TrackingRecHitSimpleSoA.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  class TrackingRecHit2DAlpaka {
  public:

    TrackingRecHit2DAlpaka() = delete;  // alpaka buffers are not default-constructible
    
    explicit TrackingRecHit2DAlpaka(const TrackingRecHitSimpleSoA& simpleHits, Queue& queue):
      m_nHits(simpleHits.nHits())
    {
      auto nModules = simpleHits.moduleStartVec().size();
      m_store32 = cms::alpakatools::make_device_buffer<uint32_t[]>(queue, m_nHits * static_cast<uint32_t>(Fields32::size_) + nModules + 1);
      m_store16 = cms::alpakatools::make_device_buffer<uint16_t[]>(queue, m_nHits * static_cast<uint32_t>(Fields16::size_));

      m_view = cms::alpakatools::make_device_buffer<TrackingRecHit2DSoAView>(queue);
      m_view_h = cms::alpakatools::make_host_buffer<TrackingRecHit2DSoAView>(queue);

      field32ToDeviceAsync(Fields32::xl, simpleHits.xlVector().data(), m_nHits, queue);
      field32ToDeviceAsync(Fields32::yl, simpleHits.ylVector().data(), m_nHits, queue);
      field32ToDeviceAsync(Fields32::xerr, simpleHits.xerrVector().data(), m_nHits, queue);
      field32ToDeviceAsync(Fields32::yerr, simpleHits.yerrVector().data(), m_nHits, queue);
      field32ToDeviceAsync(Fields32::xg, simpleHits.xgVector().data(), m_nHits, queue);
      field32ToDeviceAsync(Fields32::yg, simpleHits.ygVector().data(), m_nHits, queue);
      field32ToDeviceAsync(Fields32::zg, simpleHits.zgVector().data(), m_nHits, queue);
      field32ToDeviceAsync(Fields32::rg, simpleHits.rgVector().data(), m_nHits, queue);
      field16ToDeviceAsync(Fields16::iphi, simpleHits.iphiVector().data(), m_nHits, queue);
      field32ToDeviceAsync(Fields32::charge, simpleHits.chargeVector().data(), m_nHits, queue);
      field16ToDeviceAsync(Fields16::xsize, simpleHits.xsizeVector().data(), m_nHits, queue);
      field16ToDeviceAsync(Fields16::ysize, simpleHits.ysizeVector().data(), m_nHits, queue);
      field16ToDeviceAsync(Fields16::detInd, simpleHits.detIndVector().data(), m_nHits, queue);
      field32ToDeviceAsync(Fields32::partInd, simpleHits.partIndVector().data(), m_nHits, queue);
      
      m_nHits = (*m_view_h)->m_nHits = uint32_t(m_nHits);
      m_hitsModuleStart =  (*m_view_h)->m_hitsModuleStart = endOf32();

      // pointers to data owned by this TrackingRecHit2DAlpaka object:
      (*m_view_h)->m_xl     = const_cast<float*>(simpleHits.xlVector().data());
      (*m_view_h)->m_yl     = const_cast<float*>(simpleHits.ylVector().data());
      (*m_view_h)->m_xerr   = const_cast<float*>(simpleHits.xerrVector().data());
      (*m_view_h)->m_yerr   = const_cast<float*>(simpleHits.yerrVector().data());
      (*m_view_h)->m_xg     = const_cast<float*>(simpleHits.xgVector().data());
      (*m_view_h)->m_yg     = const_cast<float*>(simpleHits.ygVector().data());
      (*m_view_h)->m_zg     = const_cast<float*>(simpleHits.zgVector().data());
      (*m_view_h)->m_rg     = const_cast<float*>(simpleHits.rgVector().data());
      (*m_view_h)->m_iphi   = const_cast<short int*>(simpleHits.iphiVector().data());
      (*m_view_h)->m_charge = const_cast<int*>(simpleHits.chargeVector().data());
      (*m_view_h)->m_xsize  = const_cast<short int*>(simpleHits.xsizeVector().data());
      (*m_view_h)->m_ysize  = const_cast<short int*>(simpleHits.ysizeVector().data());
      (*m_view_h)->m_detInd = const_cast<short int*>(simpleHits.detIndVector().data());
      (*m_view_h)->m_partInd = const_cast<uint32_t*>(simpleHits.partIndVector().data());
      (*m_view_h)->m_hitsModuleStart = const_cast<uint32_t*>(simpleHits.moduleStartVec().data());

      // copy the SoA view to the device
      alpaka::memcpy(queue, *m_view, *m_view_h);

      alpaka::wait(queue);

    };

    explicit TrackingRecHit2DAlpaka(
      size_t nHits, size_t nModules,
      const float* xl, const float* yl, const float* xerr, const float* yerr,
      const float* xg, const float* yg, const float* zg, const float* rg, const int16_t* iphi,
      const int32_t* charge, const int16_t* xsize, const int16_t* ysize, const int16_t* detInd,
      const uint32_t* modStart, Queue& queue)
      : m_nHits(nHits)
  {
      m_store32 = cms::alpakatools::make_device_buffer<uint32_t[]>(queue, m_nHits * static_cast<uint32_t>(Fields32::size_) + nModules + 1);
      m_store16 = cms::alpakatools::make_device_buffer<uint16_t[]>(queue, m_nHits * static_cast<uint32_t>(Fields16::size_));

      m_view = cms::alpakatools::make_device_buffer<TrackingRecHit2DSoAView>(queue);
      m_view_h = cms::alpakatools::make_host_buffer<TrackingRecHit2DSoAView>(queue);
      
      field32ToDeviceAsync(Fields32::yl, xl, m_nHits, queue);
      field32ToDeviceAsync(Fields32::yl, yl, m_nHits, queue);
      field32ToDeviceAsync(Fields32::xerr, xerr, m_nHits, queue);
      field32ToDeviceAsync(Fields32::yerr, yerr, m_nHits, queue);
      field32ToDeviceAsync(Fields32::xg, xg, m_nHits, queue);
      field32ToDeviceAsync(Fields32::yg, yg, m_nHits, queue);
      field32ToDeviceAsync(Fields32::zg, zg, m_nHits, queue);
      field32ToDeviceAsync(Fields32::rg, rg, m_nHits, queue);
      field32ToDeviceAsync(Fields32::charge, charge, m_nHits, queue);
      field16ToDeviceAsync(Fields16::xsize, xsize, m_nHits, queue);
      field16ToDeviceAsync(Fields16::ysize, ysize, m_nHits, queue);
      field16ToDeviceAsync(Fields16::detInd, detInd, m_nHits, queue);

      m_nHits = (*m_view_h)->m_nHits = uint32_t(m_nHits);
      m_hitsModuleStart =  (*m_view_h)->m_hitsModuleStart = endOf32();

      alpaka::memcpy(queue, *m_view_h, *m_view);
      alpaka::wait(queue);

  }

    explicit TrackingRecHit2DAlpaka(uint32_t nHits,
                                    const uint32_t* hitsModuleStart,
                                    Queue& queue)
        : m_nHits(nHits),
          // non-owning device pointers
          m_hitsModuleStart(reinterpret_cast<const uint32_t*>(hitsModuleStart)),
          // TODO replace with Eric's SoA
          // 32-bit SoA data members packed in a single buffer
          m_store32{cms::alpakatools::make_device_buffer<uint32_t[]>(
              queue, nHits * static_cast<uint32_t>(Fields32::size_))},
          // 16-bit SoA data members packed in a single buffer
          m_store16{
              cms::alpakatools::make_device_buffer<uint16_t[]>(queue, nHits * static_cast<uint32_t>(Fields16::size_))},
          // SoA view
          m_view{cms::alpakatools::make_device_buffer<TrackingRecHit2DSoAView>(queue)},
          m_view_h{cms::alpakatools::make_host_buffer<TrackingRecHit2DSoAView>(queue)} {
      // the hits are actually accessed in order only in building
      // if ordering is relevant they may have to be stored phi-ordered by layer or so
      // this will break 1to1 correspondence with cluster and module locality
      // so unless proven VERY inefficient we keep it ordered as generated

      // copy all the pointers
      (*m_view_h)->m_nHits = nHits;

      // pointer to data already owned in the event by SiPixelClusterAlpaka object:
      (*m_view_h)->m_hitsModuleStart = reinterpret_cast<const uint32_t*>(hitsModuleStart);

      // pointers to data owned by this TrackingRecHit2DAlpaka object:
      (*m_view_h)->m_xl = xl();
      (*m_view_h)->m_yl = yl();
      (*m_view_h)->m_xerr = xerr();
      (*m_view_h)->m_yerr = yerr();
      (*m_view_h)->m_xg = xg();
      (*m_view_h)->m_yg = yg();
      (*m_view_h)->m_zg = zg();
      (*m_view_h)->m_rg = rg();
      (*m_view_h)->m_iphi = iphi();
      (*m_view_h)->m_charge = charge();
      (*m_view_h)->m_xsize = xsize();
      (*m_view_h)->m_ysize = ysize();
      (*m_view_h)->m_detInd = detInd();

      // copy the SoA view to the device
      alpaka::memcpy(queue, *m_view, *m_view_h);

    }

    ~TrackingRecHit2DAlpaka() = default;

    TrackingRecHit2DAlpaka(const TrackingRecHit2DAlpaka&) = delete;
    TrackingRecHit2DAlpaka& operator=(const TrackingRecHit2DAlpaka&) = delete;

    TrackingRecHit2DAlpaka(TrackingRecHit2DAlpaka&&) = default;
    TrackingRecHit2DAlpaka& operator=(TrackingRecHit2DAlpaka&&) = default;

    TrackingRecHit2DSoAView* view() { return (*m_view).data(); }
    TrackingRecHit2DSoAView const* view() const { return (*m_view).data(); }

    TrackingRecHit2DSoAView* view_h() { return (*m_view_h).data(); }
    TrackingRecHit2DSoAView const* view_h() const { return (*m_view_h).data(); }

    uint32_t nHits() const { return m_nHits; }

    uint32_t const* hitsModuleStart() const { return m_hitsModuleStart; }

    auto xlToHostAsync(Queue& queue) const {
      auto device_view = cms::alpakatools::make_device_view(alpaka::getDev(queue), xl(), nHits());
      auto host_buffer = cms::alpakatools::make_host_buffer<float[]>(queue, nHits());
      alpaka::memcpy(queue, host_buffer, device_view);
      return host_buffer;
    }
    auto ylToHostAsync(Queue& queue) const {
      auto device_view = cms::alpakatools::make_device_view(alpaka::getDev(queue), yl(), nHits());
      auto host_buffer = cms::alpakatools::make_host_buffer<float[]>(queue, nHits());
      alpaka::memcpy(queue, host_buffer, device_view);
      return host_buffer;
    }
    auto xerrToHostAsync(Queue& queue) const {
      auto device_view = cms::alpakatools::make_device_view(alpaka::getDev(queue), xerr(), nHits());
      auto host_buffer = cms::alpakatools::make_host_buffer<float[]>(queue, nHits());
      alpaka::memcpy(queue, host_buffer, device_view);
      return host_buffer;
    }
    auto yerrToHostAsync(Queue& queue) const {
      auto device_view = cms::alpakatools::make_device_view(alpaka::getDev(queue), yerr(), nHits());
      auto host_buffer = cms::alpakatools::make_host_buffer<float[]>(queue, nHits());
      alpaka::memcpy(queue, host_buffer, device_view);
      return host_buffer;
    }
    auto xgToHostAsync(Queue& queue) const {
      auto device_view = cms::alpakatools::make_device_view(alpaka::getDev(queue), xg(), nHits());
      auto host_buffer = cms::alpakatools::make_host_buffer<float[]>(queue, nHits());
      alpaka::memcpy(queue, host_buffer, device_view);
      return host_buffer;
    }
    auto ygToHostAsync(Queue& queue) const {
      auto device_view = cms::alpakatools::make_device_view(alpaka::getDev(queue), yg(), nHits());
      auto host_buffer = cms::alpakatools::make_host_buffer<float[]>(queue, nHits());
      alpaka::memcpy(queue, host_buffer, device_view);
      return host_buffer;
    }
    auto zgToHostAsync(Queue& queue) const {
      auto device_view = cms::alpakatools::make_device_view(alpaka::getDev(queue), zg(), nHits());
      auto host_buffer = cms::alpakatools::make_host_buffer<float[]>(queue, nHits());
      alpaka::memcpy(queue, host_buffer, device_view);
      return host_buffer;
    }
    auto rgToHostAsync(Queue& queue) const {
      auto device_view = cms::alpakatools::make_device_view(alpaka::getDev(queue), rg(), nHits());
      auto host_buffer = cms::alpakatools::make_host_buffer<float[]>(queue, nHits());
      alpaka::memcpy(queue, host_buffer, device_view);
      return host_buffer;
    }
    auto chargeToHostAsync(Queue& queue) const {
      auto device_view = cms::alpakatools::make_device_view(alpaka::getDev(queue), charge(), nHits());
      auto host_buffer = cms::alpakatools::make_host_buffer<int32_t[]>(queue, nHits());
      alpaka::memcpy(queue, host_buffer, device_view);
      return host_buffer;
    }
    auto xsizeToHostAsync(Queue& queue) const {
      auto device_view = cms::alpakatools::make_device_view(alpaka::getDev(queue), xsize(), nHits());
      auto host_buffer = cms::alpakatools::make_host_buffer<int16_t[]>(queue, nHits());
      alpaka::memcpy(queue, host_buffer, device_view);
      return host_buffer;
    }
    auto ysizeToHostAsync(Queue& queue) const {
      auto device_view = cms::alpakatools::make_device_view(alpaka::getDev(queue), ysize(), nHits());
      auto host_buffer = cms::alpakatools::make_host_buffer<int16_t[]>(queue, nHits());
      alpaka::memcpy(queue, host_buffer, device_view);
      return host_buffer;
    }

#ifdef TODO
    // only the local coord and detector index
    cms::alpakatools::host_buffer<uint16_t[]> detIndexToHostAsync(Queue& queue) const;
    cms::alpakatools::host_buffer<uint32_t[]> hitsModuleStartToHostAsync(Queue& queue) const;
#endif

    // non-const accessors
    float* xl() { return getField32<float>(Fields32::xl); }
    float* yl() { return getField32<float>(Fields32::yl); }
    float* xerr() { return getField32<float>(Fields32::xerr); }
    float* yerr() { return getField32<float>(Fields32::yerr); }

    float* xg() { return getField32<float>(Fields32::xg); }
    float* yg() { return getField32<float>(Fields32::yg); }
    float* zg() { return getField32<float>(Fields32::zg); }
    float* rg() { return getField32<float>(Fields32::rg); }

    int32_t* charge() { return getField32<int32_t>(Fields32::charge); }
    int16_t* xsize() { return getField16<int16_t>(Fields16::xsize); }
    int16_t* ysize() { return getField16<int16_t>(Fields16::ysize); }

    int16_t* iphi() { return getField16<int16_t>(Fields16::iphi); }
    int16_t* detInd() { return getField16<int16_t>(Fields16::detInd); }

    uint32_t* partInd() { return getField32<uint32_t>(Fields32::partInd); }

    // const accessors
    float const* xl() const { return getField32<float>(Fields32::xl); }
    float const* yl() const { return getField32<float>(Fields32::yl); }
    float const* xerr() const { return getField32<float>(Fields32::xerr); }
    float const* yerr() const { return getField32<float>(Fields32::yerr); }

    float const* xg() const { return getField32<float>(Fields32::xg); }
    float const* yg() const { return getField32<float>(Fields32::yg); }
    float const* zg() const { return getField32<float>(Fields32::zg); }
    float const* rg() const { return getField32<float>(Fields32::rg); }

    int32_t const* charge() const { return getField32<int32_t>(Fields32::charge); }
    uint32_t const* endOf32() const { return getField32<uint32_t>(Fields32::size_); }

    int16_t const* xsize() const { return getField16<int16_t>(Fields16::xsize); }
    int16_t const* ysize() const { return getField16<int16_t>(Fields16::ysize); }

    int16_t const* iphi() const { return getField16<int16_t>(Fields16::iphi); }
    int16_t const* detInd() const { return getField16<int16_t>(Fields16::detInd); }

    uint32_t const* partInd() const { return getField32<uint32_t>(Fields32::partInd); }

    // explicitly const accessors
    int16_t const* c_iphi() const { return getField16<int16_t>(Fields16::iphi); }

  private:
    uint32_t m_nHits;

    // non-owning device pointers
    // m_hitsModuleStart data is already owned by SiPixelClusterAlpaka, let's not abuse of shared_ptr!!
    uint32_t const* m_hitsModuleStart;  // needed for legacy, this is on GPU!

    // TODO replace with Eric's SoA
    static_assert(sizeof(uint32_t) == sizeof(float));  // just stating the obvious
    // 32-bit SoA data members (float, int32_t or uint32_t) packed in a single buffer
    enum class Fields32 : uint32_t { xl, yl, xerr, yerr, xg, yg, zg, rg, charge, partInd, size_ };
    // 16-bit SoA data members (int16_t or uint16_t) packed in a single buffer
    enum class Fields16 : uint32_t { iphi, detInd, xsize, ysize, size_ };

    template <typename T>
    T* getField32(Fields32 field) {
      return reinterpret_cast<T*>((*m_store32).data() + static_cast<uint32_t>(field) * m_nHits);
    }

    template <typename T>
    T const* getField32(Fields32 field) const {
      return reinterpret_cast<T const*>((*m_store32).data() + static_cast<uint32_t>(field) * m_nHits);
    }

    template <typename T>
    T* getField16(Fields16 field) {
      return reinterpret_cast<T*>((*m_store16).data() + static_cast<uint32_t>(field) * m_nHits);
    }

    template <typename T>
    T const* getField16(Fields16 field) const {
      return reinterpret_cast<T const*>((*m_store16).data() + static_cast<uint32_t>(field) * m_nHits);
    }

    // 32-bit SoA data members packed in a single buffer
    std::optional<cms::alpakatools::device_buffer<Device, uint32_t[]>> m_store32;
    // 16-bit SoA data members packed in a single buffer
    std::optional<cms::alpakatools::device_buffer<Device, uint16_t[]>> m_store16;

    // This is a SoA view which itself gathers non-owning pointers to the data owned above (in TrackingRecHit2DAlpaka instance).
    // This is used to access and modify data on GPU in a SoA format (TrackingRecHit2DSoAView),
    // while the data itself is owned here in the TrackingRecHit2DAlpaka instance.
    std::optional<cms::alpakatools::device_buffer<Device, TrackingRecHit2DSoAView>> m_view;

    // Keep a host copy of the device view alive during the asynchronous copy
    std::optional<cms::alpakatools::host_buffer<TrackingRecHit2DSoAView>> m_view_h;

  public:

    template<typename T>
    void field32ToDeviceAsync(Fields32 field, const T* inputs, unsigned int n, Queue& queue) {
      auto device_view = cms::alpakatools::make_device_view(alpaka::getDev(queue), getField32<T>(field), n);
      auto host_view = cms::alpakatools::make_host_view(inputs, n);
      alpaka::memcpy(queue, device_view, host_view);
    }

    template<typename T>
    void field16ToDeviceAsync(Fields16 field, const T* inputs, unsigned int n, Queue& queue) {
      auto device_view = cms::alpakatools::make_device_view(alpaka::getDev(queue), getField16<T>(field), n);
      auto host_view = cms::alpakatools::make_host_view(inputs, n);
      alpaka::memcpy(queue, device_view, host_view);
    }


  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif  // AlpakaDataFormats_alpaka_TrackingRecHit2DAlpaka_h
