#ifndef AlpakaCore_concepts_h
#define AlpakaCore_concepts_h

#include <type_traits>

#include <alpaka/alpaka.hpp>

namespace cms::alpakatools {
  template <typename T>
  concept NonCPUQueue = alpaka::isQueue<T> and not std::is_same_v<alpaka::Dev<T>, alpaka::DevCpu>;
}  // namespace cms::alpakatools

#endif
