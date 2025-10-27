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

#include "AlpakaDataFormats/CAGeometrySoA.h"
#include "AlpakaDataFormats/CAGeometryHost.h"
#include "AlpakaDataFormats/CAGeometryDevice.h"
#include "AlpakaDataFormats/alpaka/CAGeometrySoACollection.h"

#include "AlpakaDataFormats/CAPairSoA.h"
#include "AlpakaDataFormats/CAPairHost.h"
#include "AlpakaDataFormats/CAPairDevice.h"
#include "AlpakaDataFormats/alpaka/CAPairSoACollection.h"

#include "AlpakaDataFormats/TrackingRecHitsSoA.h"
#include "AlpakaDataFormats/TrackingRecHitsHost.h"
#include "AlpakaDataFormats/TrackingRecHitsDevice.h"
#include "AlpakaDataFormats/alpaka/TrackingRecHitsSoACollection.h"

#include "AlpakaDataFormats/SiPixelClustersSoA.h"
#include "AlpakaDataFormats/SiPixelClustersHost.h"
#include "AlpakaDataFormats/SiPixelClustersDevice.h"
#include "AlpakaDataFormats/alpaka/SiPixelClustersSoACollection.h"

#include "plugin-PixelSeeding/CAStructures.h"

#include "CondFormats/SiPixelMappingSoA.h"
#include "CondFormats/alpaka/SiPixelMappingSoACollection.h"
#include "CondFormats/SiPixelMappingHost.h"
#include "CondFormats/SiPixelMappingDevice.h"
#include "CondFormats/alpaka/SiPixelMappingUtilities.h"
#include "CondFormats/SiPixelROCsStatusAndMappingConstants.h"

#include "CondFormats/SiPixelGainCalibrationForHLTSoA.h"
#include "CondFormats/alpaka/SiPixelGainCalibrationForHLTSoACollection.h"
#include "CondFormats/SiPixelGainCalibrationForHLTHost.h"
#include "CondFormats/SiPixelGainCalibrationForHLTDevice.h"

#include "plugin-PixelSeeding/alpaka/CACell.h"
//Dummy test just to compile

int main() {
  return 0;
}