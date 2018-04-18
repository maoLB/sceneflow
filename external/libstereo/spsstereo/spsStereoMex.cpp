#include "mex.h"
#include <iostream>
#include <string>
// #include "Image.h"
#include "SPSStereo.h"
#include "defParameter.h"

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  
  if (nrhs!=2 && nrhs!=3 && nrhs!=4 && nrhs!=5)
    mexErrMsgTxt("2-5 inputs required: I1 (left image), I2 (right image), [data,smoothness,sgm] parameters");
  if (nlhs!=2) 
    mexErrMsgTxt("2 outputs required: D1 (left disparities), D2 (right disparities)");  
  if (!mxIsUint8(prhs[0]))// || mxGetNumberOfDimensions(prhs[0])!=3)
    mexErrMsgTxt("Input I1 (left image) must be a uint8 image.");
  if (!mxIsUint8(prhs[1]))// || mxGetNumberOfDimensions(prhs[1])!=3)
    mexErrMsgTxt("Input I2 (right image) must be a uint8 image.");
  
  // get input pointers
  unsigned char*   I1 = (unsigned char*)mxGetPr(prhs[0]);
  unsigned char*   I2 = (unsigned char*)mxGetPr(prhs[1]);
  const int32_t *dims = mxGetDimensions(prhs[0]);
  
  // width/height
  int w = dims[1];
  int h = dims[0];

  // create input images
  //rev::Image<unsigned char> leftImage(w,h,1);
  //rev::Image<unsigned char> rightImage(w,h,1);
  png::image<png::rgb_pixel> leftImage(w,h);
	png::image<png::rgb_pixel> rightImage(w,h);
  
  // copy input
  for (int u=0; u<w; u++) {
    for (int v=0; v<h; v++) {
      leftImage.set_pixel(u,v,png::rgb_pixel(I1[u*h+v],I1[u*h+v+w*h],I1[u*h+v+2*(w*h)]));
      rightImage.set_pixel(u,v,png::rgb_pixel(I2[u*h+v],I2[u*h+v+w*h],I2[u*h+v+2*(w*h)]));
    }
  }
  //leftImage.write("test_l.png");
  //rightImage.write("test_r.png");
  
  int disp_num = 256;

  // run sgm
//   SgmStereo sgmStereo;
//   sgmStereo.setDisparityTotal(disp_num);
//   sgmStereo.setOutputDisparityFactor(disp_num);
  
  // run sps
  SPSStereo sps;
  
//   if (nrhs>=3) {
//     double* p = (double*)mxGetPr(prhs[2]);
//     sgmStereo.setDataCostParameters((int)p[0],(int)p[1],(double)p[2],(int)p[3]);
//   }
//   
//   if (nrhs>=4) {
//     double* p = (double*)mxGetPr(prhs[3]);
//     sgmStereo.setSmoothnessCostParameters((int)p[0],(int)p[1]);
//   }
//   
//   if (nrhs>=5) {
//     double* p = (double*)mxGetPr(prhs[4]);
//     sgmStereo.setSgmParameters((bool)p[0],(int)p[1]);
//   }
  
//   rev::Image<unsigned short> leftDisparityImage;
//   rev::Image<unsigned short> rightDisparityImage;
//   sgmStereo.computeLeftRight(leftImage, rightImage, leftDisparityImage, rightDisparityImage);
  
  png::image<png::gray_pixel_16> segmentImage;
	png::image<png::gray_pixel_16> disparityImage;
	std::vector< std::vector<double> > disparityPlaneParameters;
	std::vector< std::vector<int> > boundaryLabels;
	sps.compute(superpixelTotal, leftImage, rightImage, segmentImage, disparityImage, disparityPlaneParameters, boundaryLabels);
  
  //disparityImage.write("test.png");
  
  // create outputs
  plhs[0]      = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS,mxREAL);
  double*   D1 = (double*)mxGetPr(plhs[0]);
  plhs[1]      = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS,mxREAL);
  double*   L  = (double*)mxGetPr(plhs[1]);


 // copy to output
  for (int u=0; u<w; u++) {
    for (int v=0; v<h; v++) {
      D1[u*h+v] = (double)disparityImage.get_pixel(u,v)/(double)disp_num;
      L [u*h+v] = (double)segmentImage.get_pixel(u,v);
      if (D1[u*h+v]<0.001) D1[u*h+v] = -1;
//       if (D2[u*h+v]<0.001) D2[u*h+v] = -1;
    }
  }
}
 
