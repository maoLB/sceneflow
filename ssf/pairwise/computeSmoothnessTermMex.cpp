#include <mex.h>
#include <vector>
#include <iostream>
#include <algorithm> 
#include <cmath>     

#include "../../roadplane/matrix.h"

using namespace std;

#define endll endl << endl 

// extracts calibration matrix form mxArray
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

// extracts boundary pixels to vector<int> (0-based)
inline vector<int> boundaryPixels (const mxArray *SP_bnd, int32_t idx) {
  mxArray *SP_bnd_idx = mxGetCell(SP_bnd,idx);
  int32_t n = mxGetM(SP_bnd_idx);
  double* vals = (double*)mxGetPr(SP_bnd_idx);
  vector<int> bp; bp.resize(n);
  for (int i=0; i<n; i++)
    bp[i] = vals[i]-1;
  return bp;
}

inline double disparity (double u, double v, double* abc) {
  return abc[0]*u + abc[1]*v + abc[2];
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


// Prints all object particles
void printShpParticles(int nS, int nShpParameters, int nShpParticles, double* shpParticles) {
  
  // all shape particles
  for (int iSP=0; iSP<nShpParticles; iSP++) { 

    // all superpixels
    for (int iS=0; iS<nS; iS++) {
      
      cout << "\nSuperpixel " << iS << ", particle " << iSP << endl;
      
      // all parameters
      for (int iShpParam = 0; iShpParam<nShpParameters; iShpParam++) {
        cout << shpParticles[iS + iShpParam*nS + iSP*nS*nShpParameters] << " ";
      }

    }
  }
  cout << endll;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  
  // check number of arguments
  if (nrhs!=10)
    mexErrMsgTxt("10 inputs required: S, SP_adj, SP_bnd, shapeParticles, nStates, img_height, K, m, parameters, verbosityLevel");
  if (nlhs!=4)
    mexErrMsgTxt("4 outputs required: E_disp, E_normal, E_Potts, S_pairs");
  
  // check data types
  if (!mxIsStruct(prhs[0]))
    mexErrMsgTxt("S must be a struct containing the superpixels");
  if (!mxIsDouble(prhs[1]))
    mexErrMsgTxt("SP_adj must be a double matrix");
  if (!mxIsCell(prhs[2]))
    mexErrMsgTxt("SP_bnd must be a cell array");
  if (!mxIsDouble(prhs[3]))
    mexErrMsgTxt("shapeParticles must be a double matrix");
  if (!mxIsDouble(prhs[6]))
    mexErrMsgTxt("K must be a double precision camera matrix");
  if (!mxIsStruct(prhs[8]))
    mexErrMsgTxt("parameters must be a struct");
  
  // get number of superpixels
  int nS = mxGetNumberOfElements(prhs[0]); 
  
  // get pointers to input
  
  // 0) superpixels
  const mxArray* S = prhs[0];
  
  // 1) adjacency matrix
  double* SP_adj = mxGetPr(prhs[1]);
  
  // 2) superpixel boundaries
  const mxArray* SP_bnd = prhs[2];
  
  // 3) shape particles
  double* shpParticles  = mxGetPr(prhs[3]);
  
  const mwSize *shpDims = mxGetDimensions(prhs[3]);
  int dimShpParticles   = mxGetNumberOfDimensions(prhs[3]);
  int nShp              = shpDims[0];
  int nShpParameters    = shpDims[1];
  int nShpParticles;
  
  if (dimShpParticles > 2)
    nShpParticles  =  shpDims[2];
  else
    nShpParticles  = 1;
  
  // 4) number of shape particles
  int nStates           = (int)mxGetScalar(prhs[4]);
  if (nStates!=nShpParticles)
    mexErrMsgTxt("Inconsistent number of shape particles");
  
  // 5) image height (for indexing)
  int img_height        = (int)mxGetScalar(prhs[5]);
  
  // 6) calibration matrix
  Matrix K = getCalibrationMatrix( mxGetPr(prhs[6]) );
  Matrix K_inv( Matrix::inv(K) );
  
  // 7) translation vector
  double *m_ptr = (double*)mxGetPr(prhs[7]);
  Matrix m(3,1,m_ptr);
  
  // stereo base in m
  double base_m = -1.0*m.val[0][0] / K.val[0][0];
  
  // 8) parameters struct
  const mxArray* param_ = prhs[8];
  
  // 9) verbosity level
  int verbose           = (int)mxGetScalar(prhs[9]);
  
  if (verbose>0) {
    cout << endll << "Verbosity level: " << verbose << endll;
    cout << "\n Number of superpixels: " << nS << endll;
    
    cout << "  n shp:            " << nShp           << endl;
    cout << "  n shp parameters: " << nShpParameters << endl;
    cout << "  n shp particles:  " << nShpParticles  << endl;
    
    cout << "\n  image height:     " << img_height << endll;    
    
    cout << "\nK\n" << K << endll << "inv(K)\n" << K_inv << endll;
  }
  
  // query parameter struct
  const mxArray* Disp_    = mxGetField(param_,0,"Disp");
  const mxArray* Norm_    = mxGetField(param_,0,"Normal");
  const mxArray* Potts_   = mxGetField(param_,0,"Potts");
  
  bool Disp   = mxGetScalar(Disp_)    > 0.1;
  bool Norm   = mxGetScalar(Norm_)    > 0.1;
  bool Potts  = mxGetScalar(Potts_)   > 0.1;
  
  const mxArray* sDisp_   = mxGetField(param_,0,"sigmaDisp");
  const mxArray* sNorm_   = mxGetField(param_,0,"sigmaNormal");
  const mxArray* lPotts_  = mxGetField(param_,0,"lambdaPotts");
  
  double sDisp  = mxGetScalar(sDisp_);
  double sNorm  = mxGetScalar(sNorm_);
  double lPotts = mxGetScalar(lPotts_);
  
  if (verbose>0) {
    cout << "\nParameters:" << endl; 
    cout << " Disp:      " << Disp << endl;
    cout << " Norm:      " << Norm << endl;
    cout << " Potts:     " << Potts << endll;
    cout << " sigmaDisp: " << sDisp << endl;
    cout << " sigmaNorm: " << sNorm << endl;
    cout << " lPotts:    " << lPotts << endl;
  }
  
  // create output cell arrays (pairwise energies)
  int E_pair_dims[2] = {nS,nS};
  plhs[0] = mxCreateCellArray(2,E_pair_dims); // disparity differences
  plhs[1] = mxCreateCellArray(2,E_pair_dims); // 3D normals
  plhs[2] = mxCreateCellArray(2,E_pair_dims); // Potts
  
  // collect adjacent superpixel pairs
  vector<int> S1_idx,S2_idx;
  
  // loop variables
  double cu1, cv1, cu2, cv2;
  double alpha1[3], alpha2[3], abc1[3], abc2[3];
  Matrix n1(1,3), n2(1,3);
  
  // for all superpixels do
  for (int iS=0; iS<nS; iS++) {
    
    // query centroid 
    mxArray* cu1_ = mxGetField(S,iS,"cu");
    mxArray* cv1_ = mxGetField(S,iS,"cv");
    
    cu1 = mxGetScalar(cu1_);
    cv1 = mxGetScalar(cv1_);
    
    // get boundary pixel indices for superpixel i
    vector<int> bp_i = boundaryPixels(SP_bnd,iS);
    
    // for all other superpixels do
    for (int jS=iS+1; jS<nS; jS++) {
      
      // query centroid 
      mxArray* cu2_ = mxGetField(S,jS,"cu");
      mxArray* cv2_ = mxGetField(S,jS,"cv");
    
      cu2 = mxGetScalar(cu2_);
      cv2 = mxGetScalar(cv2_);
      
      // if superpixels are adjacent
      if (SP_adj[jS*nS+iS]>0.5) {

        // store superpixel pair
        S1_idx.push_back(iS+1);
        S2_idx.push_back(jS+1);
        
        // get boundary pixel indices for superpixel j
        vector<int> bp_j = boundaryPixels(SP_bnd,jS);
        
        // compute indices of intersecting boundary pixels
        // note: this assumes that the inputs bp_i and bp_j are sorted!!
        vector<int> bp_int(bp_i.size()+bp_j.size());
        vector<int>::iterator it;
        it = set_intersection (bp_i.begin(),bp_i.end(),bp_j.begin(),bp_j.end(),bp_int.begin());
        bp_int.resize(it-bp_int.begin());
        int bp_num = bp_int.size();
       
        // create output for disp differences
        int E_dims[2] = {nShpParticles,nShpParticles};
        mxArray* E_disp_ = mxCreateNumericArray(2,E_dims,mxDOUBLE_CLASS,mxREAL);
        double*  E_disp  = (double*)mxGetPr(E_disp_);

        // create output for normals
        mxArray* E_norm_ = mxCreateNumericArray(2,E_dims,mxDOUBLE_CLASS,mxREAL);
        double*  E_norm  = (double*)mxGetPr(E_norm_);
        
        // create output for Potts
        mxArray* E_pott_ = mxCreateNumericArray(2,E_dims,mxDOUBLE_CLASS,mxREAL);
        double*  E_pott  = (double*)mxGetPr(E_pott_);
        
        // all shapeParticles in the first superpixel
        for (int iSP=0; iSP<nShpParticles; iSP++) {
          
          // extract plane parameters for first superpixel
          alpha1[0] = shpParticles[iS + 0*nS + iSP*nS*nShpParameters];
          alpha1[1] = shpParticles[iS + 1*nS + iSP*nS*nShpParameters];
          alpha1[2] = shpParticles[iS + 2*nS + iSP*nS*nShpParameters];
          
          abcFromAlphaBetaGamma(alpha1,abc1,cu1,cv1);
          
          dispPlaneToNormal(abc1, K, base_m, n1);
           
          // for all planes in second superpixel do
          for (int jSP=0; jSP<nShpParticles; jSP++) {
           
            alpha2[0] = shpParticles[jS + 0*nS + jSP*nS*nShpParameters];
            alpha2[1] = shpParticles[jS + 1*nS + jSP*nS*nShpParameters];
            alpha2[2] = shpParticles[jS + 2*nS + jSP*nS*nShpParameters];
          
            abcFromAlphaBetaGamma(alpha2,abc2,cu2,cv2);
            dispPlaneToNormal(abc2, K, base_m, n2);

            // compute accumulated truncated l1 disparity error
            double e   = 0.0;
            double sse = 0.0;

            // for all boundary pixels (pixels of intersection) do
            for (int k=0; k<bp_num; k++) {

              // extract u and v image coordinates
              // note: 1-based index is calculated as planes are represented
              //       using MATLAB 1-based indices
              double u = bp_int[k]/img_height+1;
              double v = bp_int[k]%img_height+1;

              // compute disparity values given both plane hypotheses
              double d1 = disparity(u,v,abc1);
              double d2 = disparity(u,v,abc2);

              // accumulate disparity error
              double d = abs(d1-d2);
              
              sse += d*d; 
              
              e += min(d,sDisp); // truncated L1 
            }

            if (Disp) {
              E_disp[jSP*nShpParticles+iSP] = e;
            }          
            
            // Compute cos similarity
            double cos_sim = (n1.val[0][0]*n2.val[0][0] + n1.val[0][1]*n2.val[0][1] + n1.val[0][2]*n2.val[0][2]) / 
                             (sqrt(n1.val[0][0]*n1.val[0][0] + n1.val[0][1]*n1.val[0][1] + n1.val[0][2]*n1.val[0][2]) * 
                              sqrt(n2.val[0][0]*n2.val[0][0] + n2.val[0][1]*n2.val[0][1] + n2.val[0][2]*n2.val[0][2]));
            
            if (Norm) {
              // Compute penalty for difference in normals
              double dNorm   = 0.5*(1.0-cos_sim);
              E_norm[jSP*nShpParticles+iSP] = min(dNorm,sNorm);
            }
            
            if (Potts) {
              // Potts model penalizing differences in k
              double k1 = shpParticles[iS + 3*nS + iSP*nS*nShpParameters];
              double k2 = shpParticles[jS + 3*nS + jSP*nS*nShpParameters];
              E_pott[jSP*nShpParticles+iSP] = k1==k2 ? 0 : exp( (-1.0*lPotts) * sse/(bp_num+1.0)) * abs(cos_sim);
            }
          }
        }

        // set energy values in cell corresponding to superpixel (iS,jS)
        int subs[2] = {iS,jS};
        int idx_1 = mxCalcSingleSubscript(plhs[0],2,subs);
        mxSetCell(plhs[0],idx_1,E_disp_);
        
        int idx_2 = mxCalcSingleSubscript(plhs[1],2,subs);
        mxSetCell(plhs[1],idx_2,E_norm_);
        
        int idx_3 = mxCalcSingleSubscript(plhs[2],2,subs);
        mxSetCell(plhs[2],idx_3,E_pott_);
      }
    }
  }
          
  // create output for superpixel pairs
  plhs[3]         = mxCreateDoubleMatrix(S1_idx.size(),2,mxREAL);
  double* S_pairs = mxGetPr(plhs[3]);
  
  for (int i=0; i<S1_idx.size(); i++) {
    S_pairs[i] = S1_idx[i];
    S_pairs[i+S1_idx.size()] = S2_idx[i];
  }
  
}
