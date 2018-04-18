#include <mex.h>
#include <vector>
#include <math.h>   
#include <cmath>     
#include <assert.h>
#include <iostream>

#ifdef MM_POPCNT
  #include <nmmintrin.h>
#endif

#include <stdint.h>

#include "../../roadplane/matrix.h"
#include <bitset>


using namespace std;

// extracts calibration matrix form mxArray respecting col-major ordering
Matrix getCalibrationMatrix(const double* K_ptr) {
  
  Matrix K = Matrix::eye(3);

  // focal lengths
  K.val[0][0] = K_ptr[0];
  K.val[1][1] = K_ptr[4];
  
  // principal point
  K.val[0][2] = K_ptr[6];
  K.val[1][2] = K_ptr[7];

  return K;
}

// Prints all object particles
void printObjectParticles(int nObj, int nObjParticles, double* objParticles) {

  int nObjParameters = 6;
  
  // all object particles
  for( int iOP=0; iOP<nObjParticles; iOP++ ) { 

    // all objects
    for( int iO=0; iO<nObj; iO++ ) {

      double rx = objParticles[iO + 0*nObj + iOP*nObj*nObjParameters];
      double ry = objParticles[iO + 1*nObj + iOP*nObj*nObjParameters];
      double rz = objParticles[iO + 2*nObj + iOP*nObj*nObjParameters];
      double tx = objParticles[iO + 3*nObj + iOP*nObj*nObjParameters];
      double ty = objParticles[iO + 4*nObj + iOP*nObj*nObjParameters];
      double tz = objParticles[iO + 5*nObj + iOP*nObj*nObjParameters];

      cout << "\nObject " << iO << ", particle " << iOP << endl;
      cout << "rx: " << rx << endl;
      cout << "ry: " << ry << endl;
      cout << "rz: " << rz << endl;
      cout << "tx: " << tx << endl;
      cout << "ty: " << ty << endl;
      cout << "tz: " << tz << endl;

    }
  }
}  

int abcFromAlphaBetaGamma(double* alpha,double* abc, double cu, double cv) {
  
  abc[0] = alpha[0];
  abc[1] = alpha[1];
  abc[2] = alpha[2] - alpha[0] * cu - alpha[1] * cv;
  
  return 0;
}

// computes 3D scaled normal from disparity plane
int dispPlaneToNormal(double* dispPlane, const Matrix& K, const double& base_m, Matrix& n) {

// Parameters
//    dispPlane   coefficients of the disparity plane in (a,b,c)-representation
//    K           3x3 camera matrix containing principal distance and point
//    base_m      (positive) length of the stereo baseline

    // camera parameters
    double f  = K.val[0][0];
    double cu = K.val[0][2];
    double cv = K.val[1][2];
    
    // parameters of the disparity plane  
    // disp = a*u + b*v + c
    double a = dispPlane[0];
    double b = dispPlane[1];
    double c = dispPlane[2];
        
    // transform disparity plane to scaled normal
    n.val[0][0] = -a/base_m;
    n.val[0][1] = -b/base_m;
    n.val[0][2] = -(c + cu*a + cv*b) / (f*base_m);
    
    return 0;
}

inline double disparity (double u, double v, double* abc) {
  return abc[0]*u + abc[1]*v + abc[2];
}


double energyCensus(int nImCo, 
        const Matrix& u1, const Matrix& v1,
        const Matrix& u2, const Matrix& v2,
        const double& oCens, const double& wCens, const double& tCens, 
        const int& w, const int& h, const int& strideCensus,
        const uint32_t* Cl_1, const uint32_t* Cr_1) 
{

    double E_Cens = 0.0;
    for (int i=0; i<nImCo; i+=strideCensus) {

      // penalize 'invisible' pixels
      if (round(u2.val[0][i]-1.0) < 0.0 || round(u2.val[0][i]-1) >= w ||
          round(v2.val[0][i]-1.0) < 0.0 || round(v2.val[0][i]-1) >= h) {
        E_Cens += oCens;
      }
      else {

        // index in left and right census image
        int idxL = (round(u1.val[0][i]-1)*h + round(v1.val[0][i])-1);
        int idxR = (round(u2.val[0][i]-1)*h + round(v2.val[0][i])-1);

        // just to make sure...
        if (idxL < 0 || idxL > w*h)
          mexErrMsgTxt("out-of-bounds access to Cl_1");
        if (idxR < 0 || idxR > w*h){
          cout << "u2: " << round(u2.val[0][i]-1) << ", v2: " << round(v2.val[0][i])-1 << endl;
          mexErrMsgTxt("out-of-bounds access to Cl_2");
        }

        // get local binary patterns
        unsigned int cl = Cl_1[ idxL ];
        unsigned int cr = Cr_1[ idxR ];

        // compute Hamming distance
        uint32_t xor_c = cl^cr;

        #ifdef MM_POPCNT
          int hammingDistance = _mm_popcnt_u32(xor_c);
        #else
          int hammingDistance = __builtin_popcount(xor_c);
        #endif

        // update energy
        E_Cens += min(double(hammingDistance)/24.0, tCens); // truncated penalty
      }
    }
    return E_Cens;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  
  // Check number of inputs   
  if (nrhs!=12)
    mexErrMsgTxt("10 inputs required: S, K, m, ShapeParticles, ObjectParticles, parameters, D1, Fl, FC1, Cl, Cr, verbosityLevel ");
  if (nlhs<1)
    mexErrMsgTxt("At least 1 output required (E_data)");
  
  // Check data types
  if (!mxIsStruct(prhs[0]))
    mexErrMsgTxt("S must be a struct containing the superpixels");
  if (!mxIsDouble(prhs[1]))
    mexErrMsgTxt("K must be a double precision camera matrix");
  if (!mxIsDouble(prhs[2]))
    mexErrMsgTxt("m must be a double matrix");
  if (!mxIsStruct(prhs[5]))
    mexErrMsgTxt("parameters must be a struct");
  if (!mxIsDouble(prhs[6]))
    mexErrMsgTxt("D must be a double matrix");
  
  if (!mxIsCell(prhs[9]))
    mexErrMsgTxt("Cl must be a cell containing two census transformed images");
  if (!mxIsUint32( mxGetCell(prhs[9],0) ))
    mexErrMsgTxt("census transformed images must be of type uint32");
  if (!mxIsUint32( mxGetCell(prhs[9],1) ))
    mexErrMsgTxt("census transformed images must be of type uint32");
  if (!mxIsCell(prhs[10]))
    mexErrMsgTxt("Cr must be a cell containing two census transformed images");
  
  // Check dimensions
  if (mxGetNumberOfDimensions(prhs[7])!=3)
    mexErrMsgTxt("Fl must be a three-dimensional matrix");
  if (mxGetNumberOfDimensions(prhs[8])!=3)
    mexErrMsgTxt("FC1 must be a three-dimensional matrix");
  
  
  // ====== get input =======
  
  // 0) struct array of superpixels
  int nS = mxGetNumberOfElements(prhs[0]); // number of superpixels 
  const mxArray* S = prhs[0];
   
  // 1) calibration matrix
  Matrix K = getCalibrationMatrix( mxGetPr(prhs[1]) );
  Matrix K_inv( Matrix::inv(K) );
  
  // 2) translation vector
  double *m_ptr = (double*)mxGetPr(prhs[2]);
  Matrix m(3,1,m_ptr);

  // stereo base in m
  double base_m = -1.0*m.val[0][0] / K.val[0][0];
  
  // 3) shape particles
  const mwSize *shpDims = mxGetDimensions(prhs[3]);
  
  double *shpParticles = (double*)mxGetPr(prhs[3]);
  int dimShpParticles  = mxGetNumberOfDimensions(prhs[3]);
  
  int nShp           = shpDims[0];
  int nShpParameters = shpDims[1];
  int nShpParticles;
  
  if (dimShpParticles > 2)
    nShpParticles  =  shpDims[2];
  else
    nShpParticles  = 1;
  
  // 4) object particles
  const mwSize *objDims = mxGetDimensions(prhs[4]);
  
  double *objParticles = (double*)mxGetPr(prhs[4]);
  int dimObjParticles  = mxGetNumberOfDimensions(prhs[4]);
  
  int nObj           = objDims[0];
  int nObjParameters = objDims[1];
  int nObjParticles;
  
  if (dimObjParticles > 2)
    nObjParticles  = objDims[2];
  else
    nObjParticles  = 1;
      
  // 5) parameters struct
  const mxArray* param_ = prhs[5];
  
  // 6) disparity map at t0
  int w = mxGetN(prhs[6]);
  int h = mxGetM(prhs[6]);
  double *D1 = (double*)mxGetPr(prhs[6]);
  
  // 7) optical flow left
  const mwSize *f1Dims = mxGetDimensions(prhs[7]);
  double *F1 = (double*)mxGetPr(prhs[7]);
  
  // 8) optical flow left t0 to right t1 (cross term 1)
  double *FC1 = (double*)mxGetPr(prhs[8]);
  
  // 9) cell containing left census images
  const mxArray* Cl_1_ = mxGetCell(prhs[9],0);
  unsigned int* Cl_1   = (unsigned int*)mxGetPr(Cl_1_);
  
  const mxArray* Cl_2_ = mxGetCell(prhs[9],1);
  unsigned int* Cl_2   = (unsigned int*)mxGetPr(Cl_2_);
  
  // 10) cell containing right census images
  const mxArray* Cr_1_ = mxGetCell(prhs[10],0);
  unsigned int* Cr_1   = (unsigned int*)mxGetPr(Cr_1_);
  
  const mxArray* Cr_2_ = mxGetCell(prhs[10],1);
  unsigned int* Cr_2   = (unsigned int*)mxGetPr(Cr_2_);
  
  // 11) verbosity level
  double verbose = mxGetScalar(prhs[11]);
  
  if (verbose>0) {
    cout << endll << "Verbosity level: " << verbose << endll;
    cout << "\n Number of superpixels: " << nS << endll;
    
    cout << "  n shp:            " << nShp           << endl;
    cout << "  n shp parameters: " << nShpParameters << endl;
    cout << "  n shp particles:  " << nShpParticles  << endl;
    
    cout << "  n obj:            " << nObj           << endl;
    cout << "  n obj parameters: " << nObjParameters << endl;
    cout << "  n obj particles:  " << nObjParticles  << endl;

    cout << "\nSize of disp map (wxh): " << w << " x " << h << endll;    
  }
  
  // Query parameter struct
  if(mxGetField(param_,0,"sigmaSGM") == NULL)
    mexErrMsgTxt("no field 'sigmaSGM'");
  
  const mxArray* S1_SGM_  = mxGetField(param_,0,"S1_SGM");
  const mxArray* S1_Cens_ = mxGetField(param_,0,"S1_Census");
  const mxArray* F1_Fl_   = mxGetField(param_,0,"F1_spFlow");
  const mxArray* F1_Cens_ = mxGetField(param_,0,"F1_Census");
  const mxArray* C1_Fl_   = mxGetField(param_,0,"C1_FlowDisp");
  const mxArray* C1_Cens_ = mxGetField(param_,0,"C1_Census");
  
  bool S1_SGM  = mxGetScalar(S1_SGM_)  > 0.0;
  bool S1_Cens = mxGetScalar(S1_Cens_) > 0.0;
  bool F1_Fl   = mxGetScalar(F1_Fl_)   > 0.0;
  bool F1_Cens = mxGetScalar(F1_Cens_) > 0.0;
  bool C1_Fl   = mxGetScalar(C1_Fl_)   > 0.0;
  bool C1_Cens = mxGetScalar(C1_Cens_) > 0.0;
  
  const mxArray* wDisp_ = mxGetField(param_,0,"weightSGM");
  const mxArray* sDisp_ = mxGetField(param_,0,"sigmaSGM");
  const mxArray* wFlow_ = mxGetField(param_,0,"weightSpFlow");
  const mxArray* sFlow_ = mxGetField(param_,0,"sigmaSpFlow");
  const mxArray* tCens_ = mxGetField(param_,0,"threshCensus");
  const mxArray* wCens_ = mxGetField(param_,0,"weightCensus");
  const mxArray* strCs_ = mxGetField(param_,0,"strideCensus");
  const mxArray* oCens_ = mxGetField(param_,0,"occCensus");
  
  double wDisp = mxGetScalar(wDisp_);
  double sDisp = mxGetScalar(sDisp_);
  double wFlow = mxGetScalar(wFlow_);
  double sFlow = mxGetScalar(sFlow_);
  double tCens = mxGetScalar(tCens_);
  double wCens = mxGetScalar(wCens_);
  double oCens = mxGetScalar(oCens_);
  
  int strideCensus = (int)mxGetScalar(strCs_);
  
  if (verbose>0) {
  	cout << endl;
    cout <<   " wDisp:   " << wDisp << endl;
    cout <<   " sDisp:   " << sDisp << endl;
    cout <<   " wFlow:   " << wFlow << endl;
    cout <<   " sFlow:   " << sFlow << endl;
    cout <<   " tauCens: " << tCens << endl;
    cout <<   " wCens:   " << wCens << endl;
    cout <<   " oCens:   " << oCens << endl;
  }
  
  // create output matrix
  plhs[0]   = mxCreateDoubleMatrix(nS,nShpParticles*nObjParticles,mxREAL);
  double* E = mxGetPr(plhs[0]);
     
  // create variables for the loop 
  int iO;             // index of the assigned object
  int nImCo;          // number of image coordinates
  
  double rx_deg, rx_rad;
  double ry_deg, ry_rad;
  double rz_deg, rz_rad;
  double d, d1, d2;
  double cu, cv;
  double f;
  double fu,fu1,fu2,fu3;
  double fv,fv1,fv2,fv3;
  double dfu,dfv;
  
  double abc[3];      // disparity plane
  double alpha[3];    // disparity plane
  
  Matrix n(1,3);      // scaled 3D normal
  
  Matrix R(3,3);      // rotation matrix
  Matrix t(3,1);      // translation vector
  
  Matrix H2(3,3);     // the homographies
  Matrix H3(3,3);
  Matrix H4(3,3);
    
  // all superpixels
  for (int iS=0; iS<nS; iS++) {
    
    // query image coordinates
    mxArray* u_ = mxGetField(S,iS,"u");
    mxArray* v_ = mxGetField(S,iS,"v");
    
    // query centroid 
    mxArray* cu_ = mxGetField(S,iS,"cu");
    mxArray* cv_ = mxGetField(S,iS,"cv");
    
    cu = mxGetScalar(cu_);
    cv = mxGetScalar(cv_);
    
    // number of image coordinates
    nImCo = mxGetNumberOfElements(u_);
       
    Matrix u1(1,nImCo,(double*)mxGetPr(u_));
    Matrix v1(1,nImCo,(double*)mxGetPr(v_));
    
    // matrix of homogeneous image coordinates
    Matrix uv_h = Matrix::ones(3,nImCo);
    uv_h.setMat(u1,0,0);
    uv_h.setMat(v1,1,0);
    
    // init matrices for projected image coordinates
    Matrix uv2_h = Matrix::ones(3,nImCo);
    Matrix uv3_h = Matrix::ones(3,nImCo);
    Matrix uv4_h = Matrix::ones(3,nImCo);
    
    if (verbose > 2) {
      cout << "\n======== Superpixel " << iS << " ========\n";
      cout << "\nuv_h\n" << uv_h << endl;
    }
    
    // for all shape particles
    for (int iSP=0; iSP<nShpParticles; iSP++) {
      
      // compute scaled normal from disp plane of the shape particle
      alpha[0] = shpParticles[iS + 0*nS + iSP*nS*nShpParameters];
      alpha[1] = shpParticles[iS + 1*nS + iSP*nS*nShpParameters];
      alpha[2] = shpParticles[iS + 2*nS + iSP*nS*nShpParameters];
        
      abcFromAlphaBetaGamma(alpha,abc,cu,cv);
      dispPlaneToNormal(abc, K, base_m, n);
      
      // compute homography from left t0 to right t0
      // (it is independent of the motion parameters)
      H2 = (K - m*n) * K_inv;

      // project image coordinates
      uv2_h = H2*uv_h;
      
      Matrix u2 = uv2_h.getMat(0,0,0,-1) / uv2_h.getMat(2,0,2,-1);
      Matrix v2 = uv2_h.getMat(1,0,1,-1) / uv2_h.getMat(2,0,2,-1);
      
      // compute data term from disp map at t0
      double E_S1_Disp = 0.0;
      
      if (S1_SGM) {
        for (int i=0; i<nImCo; i++) {

          d = D1[ int(round(u1.val[0][i]-1)*h + round(v1.val[0][i])-1) ];
          if (d<0.0) continue; // ignore invalid disparities

          d1 = disparity(u1.val[0][i],v1.val[0][i],abc);

          E_S1_Disp += min( abs(d-d1), sDisp ); // truncated penalty
        }
      }
      
      // compute census data term at t0
      double E_S1_Cens = 0.0;
      
      if (S1_Cens) {
        E_S1_Cens = energyCensus(nImCo, u1, v1, u2, v2, oCens, wCens, tCens, w, h, strideCensus, Cl_1, Cr_1);
      }
      
      // for all object particles
      for (int iOP=0; iOP<nObjParticles; iOP++) {
             
        // get assigned object from the particle
        iO = shpParticles[iS + 3*nS + iSP*nS*nShpParameters] - 1; 
        
        // get rotation and translation from object particle
        rx_deg = objParticles[iO + 0*nObj + iOP*nObj*nObjParameters];
        ry_deg = objParticles[iO + 1*nObj + iOP*nObj*nObjParameters];
        rz_deg = objParticles[iO + 2*nObj + iOP*nObj*nObjParameters];
        
        rx_rad = rx_deg / 180.0 * M_PI;
        ry_rad = ry_deg / 180.0 * M_PI;
        rz_rad = rz_deg / 180.0 * M_PI;

        R = Matrix::rotMatX(rx_rad) * Matrix::rotMatY(ry_rad) * Matrix::rotMatZ(rz_rad);
        
        t.val[0][0] = objParticles[iO + 3*nObj + iOP*nObj*nObjParameters];
        t.val[1][0] = objParticles[iO + 4*nObj + iOP*nObj*nObjParameters];
        t.val[2][0] = objParticles[iO + 5*nObj + iOP*nObj*nObjParameters];
        
        // compute homographies across time for specified terms
        if (F1_Fl || F1_Cens) 
          H3 = K*(R - t*n) * K_inv; // left t0 to left  t1
        
        if (C1_Fl || C1_Cens) 
          H4 = (K*R - (K*t+m)*n) * K_inv; // left t0 to right t1

        // project image coordinates
        Matrix u3,v3;
        Matrix u4,v4;
        
        if (F1_Fl || F1_Cens) {
          uv3_h = H3*uv_h;
          u3 = uv3_h.getMat(0,0,0,-1) / uv3_h.getMat(2,0,2,-1);
          v3 = uv3_h.getMat(1,0,1,-1) / uv3_h.getMat(2,0,2,-1);
        }
          
        if (C1_Fl || C1_Cens) {
          uv4_h = H4*uv_h;
          u4 = uv4_h.getMat(0,0,0,-1) / uv4_h.getMat(2,0,2,-1);
          v4 = uv4_h.getMat(1,0,1,-1) / uv4_h.getMat(2,0,2,-1);
        }
        
        // compute data term from left optical flow
        double E_F1_Flow = 0.0;

        if (F1_Fl) {
          for (int i=0; i<nImCo; i++) {

            // check whether there is a valid observation (third channel of flow map)
            f = F1[ int(round(u1.val[0][i]-1)*h + round(v1.val[0][i])-1) + 2*w*h ];
            if (f<1.0) continue; // ignore invalid observations

            // get flow observations
            fu = F1[ int(round(u1.val[0][i]-1)*h + round(v1.val[0][i])-1) + 0*w*h ];
            fv = F1[ int(round(u1.val[0][i]-1)*h + round(v1.val[0][i])-1) + 1*w*h ];

            // compute flow from particles
            fu1 = u3.val[0][i] - u1.val[0][i];
            fv1 = v3.val[0][i] - v1.val[0][i];

            // differences of u and v components
            dfu = fu1 - fu;
            dfv = fv1 - fv;

            // sum robust penalty
            E_F1_Flow += min( hypot(dfu,dfv), sFlow ); 
          }
        }        
        
        // compute census data for left flow
        double E_F1_Cens = 0.0;
        if (F1_Cens) {
          E_F1_Cens = energyCensus(nImCo, u1, v1, u3, v3, oCens, wCens, tCens, w, h, strideCensus, Cl_1, Cl_2);
        }
        
        // compute data term from concatenated disp + flow (cross term 1)
        double E_C1_Flow = 0.0;

        if (C1_Fl) {
          for (int i=0; i<nImCo; i++) {

            // check whether there is a valid observation (third channel of flow map)
            f = FC1[ int(round(u1.val[0][i]-1)*h + round(v1.val[0][i])-1) + 2*w*h ];
            if (f<1.0) continue; // ignore invalid observations

            // get flow observations
            fu = FC1[ int(round(u1.val[0][i]-1)*h + round(v1.val[0][i])-1) + 0*w*h ];
            fv = FC1[ int(round(u1.val[0][i]-1)*h + round(v1.val[0][i])-1) + 1*w*h ];

            // compute flow from particles
            fu1 = u4.val[0][i] - u1.val[0][i];
            fv1 = v4.val[0][i] - v1.val[0][i];

            // differences of u and v components
            dfu = fu1 - fu;
            dfv = fv1 - fv;

            // sum robust penalty
            E_C1_Flow += min( hypot(dfu,dfv), sFlow ); 
          }
        }
        
        // compute census cross term
        double E_C1_Cens = 0.0;
        if (C1_Cens) {
          E_C1_Cens = energyCensus(nImCo, u1, v1, u4, v4, oCens, wCens, tCens, w, h, strideCensus, Cl_1, Cr_2);
        }
        
        // set energy
        double E_S1 = wDisp * E_S1_Disp + wCens * strideCensus * E_S1_Cens;
        double E_F1 = wFlow * E_F1_Flow + wCens * strideCensus * E_F1_Cens;
        double E_C1 = wFlow * E_C1_Flow + wCens * strideCensus * E_C1_Cens;
        E[ (iSP*nObjParticles+iOP)*nS+iS ] = E_S1 + E_F1 + E_C1;
        
      } // all object particles
    } // all shape particles
  } // all superpixels

}
