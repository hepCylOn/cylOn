#ifndef AlpakaCore_backend_h
#define AlpakaCore_backend_h

#include <string_view>

namespace cms::alpakatools {
    
  // Enumeration whose value EDModules can put in the event
  enum class Backend : unsigned short { SerialSync = 0, TbbAsync = 1, CudaAsync = 2, ROCmAsync = 3, size };

  Backend toBackend(std::string_view name);
  std::string_view toString(Backend backend);

}

#endif  // AlpakaCore_backend_h
