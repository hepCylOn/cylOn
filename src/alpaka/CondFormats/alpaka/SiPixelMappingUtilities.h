// NO REGIONALITY for the moment

// #ifndef CondFormats_interface_alpaka_SiPixelMappingUtilities_h
// #define CondFormats_interface_alpaka_SiPixelMappingUtilities_h

// #include <set>
// #include <vector>

// #include <alpaka/alpaka.hpp>

// #include "CondFormats/CablingPathToDetUnit.h"
// #include "CondFormats/PixelROC.h"
// #include "CondFormats/SiPixelFedCablingMap.h"
// #include "CondFormats/SiPixelFedCablingTree.h"
// #include "CondFormats/SiPixelMappingSoA.h"
// #include "CondFormats/SiPixelROCsStatusAndMappingConstants.h"
// #include "AlpakaCore/config.h"
// #include "AlpakaCore/memory.h"

// namespace ALPAKA_ACCELERATOR_NAMESPACE {

//   struct SiPixelMappingUtilities {
//     ALPAKA_FN_HOST_ACC ALPAKA_FN_ACC ALPAKA_FN_INLINE static bool hasQuality(const SiPixelMappingSoAConstView& view) {
//       return view.hasQuality();
//     }

//     ALPAKA_FN_HOST_ACC ALPAKA_FN_ACC ALPAKA_FN_INLINE static cms::alpakatools::device_buffer<Device, unsigned char[]>
//     getModToUnpRegionalAsync(std::set<unsigned int> const& modules,
//                              const SiPixelFedCablingTree* cabling,
//                              std::vector<unsigned int> const& fedIds,
//                              Queue& queue) {
//       auto modToUnpDevice = cms::alpakatools::make_device_buffer<unsigned char[]>(queue, pixelgpudetails::MAX_SIZE);
//       auto modToUnpHost = cms::alpakatools::make_host_buffer<unsigned char[]>(queue, pixelgpudetails::MAX_SIZE);

//       unsigned int startFed = fedIds.front();
//       unsigned int endFed = fedIds.back() - 1;

//       sipixelobjects::CablingPathToDetUnit path;
//       int index = 1;

//       for (unsigned int fed = startFed; fed <= endFed; fed++) {
//         for (unsigned int link = 1; link <= pixelgpudetails::MAX_LINK; link++) {
//           for (unsigned int roc = 1; roc <= pixelgpudetails::MAX_ROC; roc++) {
//             path = {fed, link, roc};
//             const sipixelobjects::PixelROC* pixelRoc = cabling->findItem(path);
//             if (pixelRoc != nullptr) {
//               modToUnpHost[index] = (not modules.empty()) and (modules.find(pixelRoc->rawId()) == modules.end());
//             } else {  // store some dummy number
//               modToUnpHost[index] = true;
//             }
//             index++;
//           }
//         }
//       }

//       alpaka::memcpy(queue, modToUnpDevice, modToUnpHost);

//       return modToUnpDevice;
//     }
//   };
// }  // namespace ALPAKA_ACCELERATOR_NAMESPACE
// #endif  //CondFormats_interface_alpaka_SiPixelMappingUtilities_h
