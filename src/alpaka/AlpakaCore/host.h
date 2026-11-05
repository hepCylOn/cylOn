#ifndef AlpakaCore_host_h
#define AlpakaCore_host_h

#include <alpaka/alpaka.hpp>
#include "AlpakaCore/common.h"

namespace cms::alpakatools {

  // alpaka host platform and device

  // return the alpaka host platform
  inline alpaka::PlatformCpu const& host_platform() {
    static const auto platform = alpaka_common::PlatformHost{};
    return platform;
  }

  // return the alpaka host device
  inline alpaka::DevCpu const& host() {
    static const auto host = alpaka::getDevByIdx(host_platform(), 0u);
    // assert on the host index ?
    return host;
  }


}  // namespace cms::alpakatools

#endif  // AlpakaCore_host_h
