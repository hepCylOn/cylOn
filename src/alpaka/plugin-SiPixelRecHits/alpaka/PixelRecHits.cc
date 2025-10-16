#ifdef GPU_DEBUG
#include <iostream>
#endif

#include "AlpakaCore/AllocatorPolicy.h"
#include "AlpakaCore/config.h"
#include "CondFormats/pixelCPEforGPU.h"
#include "DataFormats/TrackingRecHitSimpleSoA.h"
#include "PixelRecHits.h"
#include "gpuPixelRecHits.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  namespace {
    struct setHitsLayerStart {
      template <typename TAcc>
      ALPAKA_FN_ACC void operator()(const TAcc& acc,
                                    uint32_t const* __restrict__ hitsModuleStart,
                                    pixelCPEforGPU::ParamsOnGPU const* cpeParams,
                                    uint32_t* hitsLayerStart) const {
        ALPAKA_ASSERT_ACC(0 == hitsModuleStart[0]);

        cms::alpakatools::for_each_element_in_grid(acc, 11, [&](uint32_t i) {
          hitsLayerStart[i] = hitsModuleStart[cpeParams->layerGeometry().layerStart[i]];
#ifdef GPU_DEBUG
          printf("LayerStart %d %d: %d\n", i, cpeParams->layerGeometry().layerStart[i], hitsLayerStart[i]);
#endif
        });
      }
    };
  }  // namespace

  namespace pixelgpudetails {

    TrackingRecHit2DAlpaka PixelRecHitGPUKernel::makeHitsAsync(SiPixelDigisAlpaka const& digis_d,
                                                               SiPixelClustersAlpaka const& clusters_d,
                                                               BeamSpotAlpaka const& bs_d,
                                                               pixelCPEforGPU::ParamsOnGPU const* cpeParams,
                                                               Queue& queue) const {


       {
        // std::cout << "Testing hits from simple hits"<< std::endl;
        size_t nh = 3;
        float xl[3] = {1.0f, 2.0f, 3.0f};
        float yl[3] = {4.0f, 5.0f, 6.0f};
        float xerr[3] = {0.1f, 0.2f, 0.3f};
        float yerr[3] = {0.4f, 0.5f, 0.6f};
        float xg[] = {7.0f, 8.0f, 9.0f};
        float yg[] = {10.0f, 11.0f, 12.0f};
        float zg[] = {13.0f, 14.0f, 15.0f};
        float rg[] = {16.0f, 17.0f, 18.0f};
        int16_t iphi[] = {1, 2, 3};
        int32_t charge[] = {10, 20, 30};
        int16_t xsize[] = {4, 5, 6};
        int16_t ysize[] = {7, 8, 9};
        int16_t detInd[] = {100, 200, 300};
        uint32_t partInd[] = {100, 200, 300};
        uint32_t modStart[] = {0,1,3};

        TrackingRecHitSimpleSoA soa(
            nh, xl, yl, xerr, yerr, xg, yg, zg, rg,
            iphi, charge, xsize, ysize, detInd, partInd, modStart);

        TrackingRecHit2DAlpaka hits_test(soa, queue);
      }
      auto nHits = clusters_d.nClusters();
      TrackingRecHit2DAlpaka hits_d(nHits, clusters_d.clusModuleStart(), queue);

      const int threadsPerBlockOrElementsPerThread = 128;
      const int blocks = digis_d.nModules();  // active modules (with digis)
      const auto getHitsWorkDiv = cms::alpakatools::make_workdiv<Acc1D>(blocks, threadsPerBlockOrElementsPerThread);

#ifdef GPU_DEBUG
      std::cout << "launching getHits kernel for " << blocks << " blocks" << std::endl;
#endif
      if (blocks) {  // protect from empty events
        alpaka::enqueue(queue,
                        alpaka::createTaskKernel<Acc1D>(getHitsWorkDiv,
                                                        gpuPixelRecHits::getHits(),
                                                        cpeParams,
                                                        bs_d.data(),
                                                        digis_d.view(),
                                                        digis_d.nDigis(),
                                                        clusters_d.view(),
                                                        hits_d.view()));
      }

#ifdef GPU_DEBUG
      alpaka::wait(queue);
#endif

#if defined ALPAKA_ACC_CPU_B_TBB_T_SEQ_ASYNC_BACKEND && defined ALPAKA_DISABLE_CACHING_ALLOCATOR
      // FIXME this is required to keep the host buffer inside hits_d alive; it could be removed once the host buffers are also stream-ordered
      alpaka::wait(queue);
#endif

      return hits_d;
    }

  }  // namespace pixelgpudetails

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
