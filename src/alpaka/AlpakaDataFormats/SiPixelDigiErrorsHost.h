
#ifndef AlpakaDataFormats_SiPixelDigiErrorsHost_h
#define AlpakaDataFormats_SiPixelDigiErrorsHost_h

#include <utility>

#include <alpaka/alpaka.hpp>

#include "Framework/Uninitialized.h"
#include "Portable/PortableHostCollection.h"
#include "AlpakaDataFormats/SiPixelDigiErrorsSoA.h"
#include "DataFormats/PixelErrors.h"
#include "AlpakaCore/SimpleVector.h"
#include "AlpakaCore/memory.h"

class SiPixelDigiErrorsHost : public PortableHostCollection<SiPixelDigiErrorsSoA> {
public:
  SiPixelDigiErrorsHost(edm::Uninitialized) : PortableHostCollection<SiPixelDigiErrorsSoA>{edm::kUninitialized} {}

  template <typename TQueue>
  explicit SiPixelDigiErrorsHost(int maxFedWords, TQueue queue)
      : PortableHostCollection<SiPixelDigiErrorsSoA>(maxFedWords, queue), maxFedWords_(maxFedWords) {}

  int maxFedWords() const { return maxFedWords_; }

  auto& error_data() { return (*view().pixelErrors().data()); }
  auto const& error_data() const { return (*view().pixelErrors().data()); }

private:
  int maxFedWords_ = 0;
};

#endif  // AlpakaDataFormats_SiPixelDigiErrorsHost_h