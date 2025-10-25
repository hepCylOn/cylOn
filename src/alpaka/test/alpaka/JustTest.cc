#include <Eigen/Core>

#include "SoATemplate/SoACommon.h"
#include "SoATemplate/SoALayout.h"
#include "AlpakaCore/concepts.h"
#include "Portable/PortableHostObject.h"
#include "Portable/PortableObject.h"
#include "Portable/PortableDeviceObject.h"

#include "Portable/PortableCollection.h"
#include "Portable/PortableDeviceCollection.h"
#include "Portable/PortableHostCollection.h"
#include "Portable/PortableCollectionCommon.h"

#include "AlpakaCore/CopyToHost.h"
#include "AlpakaCore/CopyToDevice.h"
 
#include "AlpakaDataFormats/VertexSoA.h"
#include "AlpakaDataFormats/VertexHost.h"
#include "AlpakaDataFormats/VertexDevice.h"
#include "AlpakaDataFormats/alpaka/VertexSoACollection.h"

#include "Geometry/SimplePixelTopology.h"

#include "AlpakaDataFormats/TracksSoA.h"
#include "AlpakaDataFormats/TracksHost.h"
#include "AlpakaDataFormats/TracksDevice.h"
#include "AlpakaDataFormats/alpaka/TracksSoACollection.h"

//Dummy test just to compile

int main() {
  return 0;
}