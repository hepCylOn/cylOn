#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <random>

#include "AlpakaCore/alpaka/devices.h"
#include "AlpakaCore/host.h"
#include "AlpakaCore/memory.h"
#include "AlpakaCore/initialise.h"
#include "AlpakaCore/workdivision.h"

#include "AlpakaDataFormats/GeometrySoA.h"
#include "AlpakaDataFormats/GeometryHost.h"
#include "AlpakaDataFormats/alpaka/GeometryAlpaka.h"

using namespace cms::alpakatools;
using namespace ALPAKA_ACCELERATOR_NAMESPACE;

int main() {
  initialise();
  const Device device = devices<Platform>().at(0);
  Queue queue(device);
  return 0;
}
