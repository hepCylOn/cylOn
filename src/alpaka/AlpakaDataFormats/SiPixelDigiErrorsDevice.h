#ifndef AlpakaDataFormats_SiPixelDigiErrorsDevice_h
#define AlpakaDataFormats_SiPixelDigiErrorsDevice_h

#include <cstdint>

#include <alpaka/alpaka.hpp>

#include "Framework/Uninitialized.h"
#include "Portable/PortableDeviceCollection.h"
#include "AlpakaDataFormats/SiPixelDigiErrorsSoA.h"
#include "DataFormats/PixelErrors.h"
#include "AlpakaCore/SimpleVector.h"
#include "AlpakaCore/config.h"

template <typename TDev>
class SiPixelDigiErrorsDevice : public PortableDeviceCollection<SiPixelDigiErrorsSoA, TDev> {
public:
  SiPixelDigiErrorsDevice(edm::Uninitialized)
      : PortableDeviceCollection<SiPixelDigiErrorsSoA, TDev>{edm::kUninitialized} {}

  template <typename TQueue>
  explicit SiPixelDigiErrorsDevice(size_t maxFedWords, TQueue queue)
      : PortableDeviceCollection<SiPixelDigiErrorsSoA, TDev>(maxFedWords, queue), maxFedWords_(maxFedWords) {}

  // Constructor which specifies the SoA size
  explicit SiPixelDigiErrorsDevice(size_t maxFedWords, TDev const& device)
      : PortableDeviceCollection<SiPixelDigiErrorsSoA, TDev>(maxFedWords, device) {}

  auto& error_data() const { return (*this->view().pixelErrors()); }
  auto maxFedWords() const { return maxFedWords_; }

private:
  int maxFedWords_;
};

#endif  // AlpakaDataFormats_SiPixelDigiErrorsDevice_h