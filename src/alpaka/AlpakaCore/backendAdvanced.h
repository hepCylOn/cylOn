#ifndef AlpakaCore_backend_h
#define AlpakaCore_backend_h

#include <string_view>


/// TODO: align to cmssw (needs FWCore/Utilities/interface/Exception.h) 
namespace cms::alpakatools {
    
  // Enumeration whose value EDModules can put in the event
  enum class Backend : unsigned short { SerialSync = 0, TbbAsync = 1, CudaAsync = 2, ROCmAsync = 3, size };

  Backend toBackend(std::string_view name);
  std::string_view toString(Backend backend);

  inline std::string const& backendName(Backend backend) {
    static const std::string names[] = {"serial_sync", "tbb_async", "cuda_async", "rocm_async"};
    return names[static_cast<int>(backend)];
  }

  template <typename T>
  inline T& operator<<(T& out, Backend backend) {
    out << backendName(backend);
    return out;
  }

}

#endif  // AlpakaCore_backend_h
