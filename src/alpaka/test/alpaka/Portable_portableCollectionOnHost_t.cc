#include <catch2/catch_test_macros.hpp>

#include "Portable/PortableCollection.h"
#include "Portable/PortableHostCollection.h"
#include "SoATemplate/SoACommon.h"
#include "SoATemplate/SoALayout.h"

namespace {
  GENERATE_SOA_LAYOUT(TestLayout, SOA_COLUMN(double, x), SOA_COLUMN(int32_t, id))

  using TestSoA = TestLayout<>;

  constexpr auto s_tag = "[PortableCollection]";
}  // namespace

// This test is currently mostly about the code compiling
TEST_CASE("Use of PortableCollection<T, TDev> on host code", s_tag) {
  auto const size = 10;
  PortableCollection<TestSoA, alpaka::DevCpu> coll(size, cms::alpakatools::host());

  SECTION("Tests") { REQUIRE(coll->metadata().size() == size); }

  static_assert(std::is_same_v<PortableCollection<TestSoA, alpaka::DevCpu>, PortableHostCollection<TestSoA>>);
}
