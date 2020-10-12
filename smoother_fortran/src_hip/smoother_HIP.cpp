
#include "precision.h"
#include <hip/hip_runtime.h>

__global__ void applySmoother_gpu( real *weights_dev, real *f_dev, real *smoothF_dev, int nX, int nY, int nW )
{
  size_t i = hipThreadIdx_x + hipBlockIdx_x*hipBlockDim_x + nW;
  size_t j = hipThreadIdx_y + hipBlockIdx_y*hipBlockDim_y + nW;
  int iel, ism;

  if( i >= nW && i < nX-nW && j >= nW && j< nY-nW){
    real smLocal = 0.0;
    for( int jj=-nW; jj <= nW; jj++ ){
      for( int ii=-nW; ii <= nW; ii++ ){
        iel = (i+ii)+(j+jj)*nX;
        ism = (ii+nW) + (jj+nW)*(2*nW+1);
        smLocal += f[iel]*weights_dev[ism];
      }
    }
    iel = i+j*nX;
    smoothF[iel] = smLocal;
  }
}

extern "C"
{
  void applySmoother_HIP(real **weights_dev, real **f_dev, real **smoothF_dev, int nX, int nY, int nW)
  {
    int threadsPerBlock = 16;
    int gridDimX = (nX-2*nW)/threadsPerBlock + 1;
    int gridDimY = (nY-2*nW)/threadsPerBlock + 1;

    hipLaunchKernelGGL((applySmoother_gpu), dim3(gridDimX,gridDimY,1), dim3(threadsPerBlock,threadsPerBlock,1), 0, 0, *weights_dev, *f_dev, *smoothF_dev, nX, nY, nW);
  }
}

