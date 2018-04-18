#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <time.h>

#define SWAP(a,b) {temp=a;a=b;b=temp;}

using namespace std;

struct disp {
  float u,v,d;
  disp(float u,float v,float d) : u(u),v(v),d(d) {}
};

inline uint32_t getAddressOffsetImage (const int32_t& u,const int32_t& v,const int32_t& width) {
  return v*width+u;
}

float **allocateMatrix(int32_t nrow,int32_t ncol) {
  float **m;
  m    = (float**)malloc(nrow*sizeof(float*));
  m[0] = (float*)calloc(nrow*ncol,sizeof(float));
  for(int32_t i=1; i<nrow; i++) m[i]=m[i-1]+ncol;
  return m;
}

void freeMatrix(float **m) {
  free(m[0]);
  free(m);
}

void zeroMatrix(float** m, int32_t nrow,int32_t ncol) {
  for (int32_t i=0; i<nrow; i++)
    for (int32_t j=0; j<ncol; j++)
      m[i][j] = 0;
}

void printMatrix(float** m, int32_t nrow,int32_t ncol) {
  for (int32_t i=0; i<nrow; i++) {
    for (int32_t j=0; j<ncol; j++)
      cout << m[i][j] << " ";
    cout << endl;
  }
}

// Adopted by "Numerical Recipies in C"
// Linear equation solution by Gauss-Jordan elimination, equation (2.1.1) above. a[1..n][1..n]
// is the input matrix. b[1..n][1..m] is input containing the m right-hand side vectors. On
// output, a is replaced by its matrix inverse, and b is replaced by the corresponding set
// of solution vectors.
bool gaussJordanElimination(float **a, int n, float **b, int m, float eps=1e-8) {
  
  // index vectors for bookkeeping on the pivoting
  int32_t indxc[n];
  int32_t indxr[n];
  int32_t ipiv[n];
  
  // loop variables
  int32_t i,icol,irow,j,k,l,ll;
  float big,dum,pivinv,temp; 
  
  // initialize pivots to zero
  for (j=0;j<n;j++) ipiv[j]=0;
  
  // main loop over the columns to be reduced
  for (i=0;i<n;i++) {
    
    big=0.0;
    
    // search for a pivot element
    for (j=0;j<n;j++) 
      if (ipiv[j]!=1)
        for (k=0;k<n;k++)
          if (ipiv[k]==0)
            if (fabs(a[j][k])>=big) {
              big=fabs(a[j][k]);
              irow=j;
              icol=k;
            }
    ++(ipiv[icol]);

    // We now have the pivot element, so we interchange rows, if needed, to put the pivot
    // element on the diagonal. The columns are not physically interchanged, only relabeled:
    // indxc[i], the column of the ith pivot element, is the ith column that is reduced, while
    // indxr[i] is the row in which that pivot element was originally located. If indxr[i] !=
    // indxc[i] there is an implied column interchange. With this form of bookkeeping, the
    // solution b’s will end up in the correct order, and the inverse matrix will be scrambled
    // by columns.
    if (irow != icol) {
      for (l=0;l<n;l++) SWAP(a[irow][l],a[icol][l])
      for (l=0;l<m;l++) SWAP(b[irow][l],b[icol][l])
    }
        
    indxr[i]=irow; // We are now ready to divide the pivot row by the
    indxc[i]=icol; // pivot element, located at irow and icol.
        
    // check for singularity
    if (fabs(a[icol][icol]) < eps) {
      return false;
    }
    
    pivinv=1.0/a[icol][icol];
    a[icol][icol]=1.0;
    for (l=0;l<n;l++) a[icol][l] *= pivinv;
    for (l=0;l<m;l++) b[icol][l] *= pivinv;
    
    // Next, we reduce the rows except for the pivot one
    for (ll=0;ll<n;ll++) 
    if (ll!=icol) {
      dum = a[ll][icol];
      a[ll][icol] = 0.0;
      for (l=0;l<n;l++) a[ll][l] -= a[icol][l]*dum;
      for (l=0;l<m;l++) b[ll][l] -= b[icol][l]*dum;
    }
  }
  
  // This is the end of the main loop over columns of the reduction. It only remains to unscramble
  // the solution in view of the column interchanges. We do this by interchanging pairs of
  // columns in the reverse order that the permutation was built up.
  for (l=n-1;l>=0;l--) {
    if (indxr[l]!=indxc[l])
      for (k=0;k<n;k++)
        SWAP(a[k][indxr[l]],a[k][indxc[l]])
  }
  
  // success
  return true;
}

vector<disp> sparseDisparityGrid (float* D,const int *dims,int32_t* roi,int32_t step_size) {
  
  // get image width and height
  const int width  = dims[0];
  const int height = dims[1];
  
  // init list
  vector<disp> d_list;
  
  // loop through disparity image
  for (int32_t u=max(roi[0],0); u<=min(roi[1],width-1); u+=step_size) {
    for (int32_t v=max(roi[2],0); v<=min(roi[3],height-1); v+=step_size) {
      float d = *(D+getAddressOffsetImage(u,v,width));
      if (d>=1)
        d_list.push_back(disp(u,v,d));
    }
  }
  
  // return list
  return d_list;
}

void leastSquarePlane(vector<disp> &d_list,vector<int32_t> &ind,float *plane) {
  
  int32_t n = 3; int32_t m = 1;
  float** A = allocateMatrix(n,n);
  float** b = allocateMatrix(n,m);
  
  // find parameters
  for (vector<int32_t>::iterator it=ind.begin(); it!=ind.end(); it++) {
    float u = d_list[*it].u;
    float v = d_list[*it].v;
    float d = d_list[*it].d;
    A[0][0] += u*u;
    A[0][1] += u*v;
    A[0][2] += u;
    A[1][1] += v*v;
    A[1][2] += v;
    A[2][2] += 1;
    b[0][0] += u*d;
    b[1][0] += v*d;
    b[2][0] += d;
  }
  A[1][0] = A[0][1];
  A[2][0] = A[0][2];
  A[2][1] = A[1][2];

  if (gaussJordanElimination(A,3,b,1)) {
    plane[0] = b[0][0];
    plane[1] = b[1][0];
    plane[2] = b[2][0];
  } else {
    plane[0] = 0;
    plane[1] = 0;
    plane[2] = 0;
  }
  
  freeMatrix(A);
  freeMatrix(b);
}

void drawRandomPlaneSample (vector<disp> &d_list,float *plane) {
  
  int32_t num_data = d_list.size();
  vector<int32_t> ind;
  
  // draw 3 measurements
  int32_t k=0;
  while (ind.size()<3 && k<1000) {
    
    // draw random measurement
    int32_t curr_ind = rand()%num_data;
    
    // first observation
    if (ind.size()==0) {
      
      // simply add
      ind.push_back(curr_ind);
      
    // second observation
    } else if(ind.size()==1) {
      
      // check distance to first point
      float diff_u = d_list[curr_ind].u-d_list[ind[0]].u;
      float diff_v = d_list[curr_ind].v-d_list[ind[0]].v;
      if (sqrt(diff_u*diff_u+diff_v*diff_v)>50)
        ind.push_back(curr_ind);
      
    // third observation
    } else {
      
      // check distance to line between first and second point
      float vu   = d_list[ind[1]].u-d_list[ind[0]].u;
      float vv   = d_list[ind[1]].v-d_list[ind[0]].v;
      float norm = sqrt(vu*vu+vv*vv);
      float nu   = +vv/norm;
      float nv   = -vu/norm;
      float ru   = d_list[curr_ind].u-d_list[ind[0]].u;
      float rv   = d_list[curr_ind].v-d_list[ind[0]].v;
      if (fabs(nu*ru+nv*rv)>50)
        ind.push_back(curr_ind);
    }
    
    k++;
  }
  
  // return zero plane on error
  if (ind.size()==0) {
    plane[0] = 0;
    plane[1] = 0;
    plane[2] = 0;
    return;
  }
  
  // find least squares solution
  leastSquarePlane(d_list,ind,plane);
}

void computeRoadPlaneEstimate (float* D,const int *dims,float* plane,int32_t* roi,int32_t num_samples,float d_threshold) {
  
  // get list with disparities
  vector<disp> d_list = sparseDisparityGrid(D,dims,roi,5);
  
  // loop variables
  vector<int32_t> curr_inlier;
  vector<int32_t> best_inlier;
  
  for (int32_t i=0; i<num_samples; i++) {
    
    // draw random samples and compute plane
    drawRandomPlaneSample(d_list,plane);
    
    // find inlier
    curr_inlier.clear();
    for (int32_t i=0; i<d_list.size(); i++)
      if (fabs(plane[0]*d_list[i].u+plane[1]*d_list[i].v+plane[2]-d_list[i].d)<d_threshold)
        curr_inlier.push_back(i);

    // is this a better solution? (=more inlier)
    if (curr_inlier.size()>best_inlier.size()) {
      best_inlier = curr_inlier;
    }
  }
  
  // reoptimize plane with inliers only
  if (curr_inlier.size()>3)
    leastSquarePlane(d_list,best_inlier,plane);
}

void mexFunction (int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]) {
	   
    // check for proper number of arguments
    if(nrhs!=4) 
        mexErrMsgTxt("Four input required (D,roi,num_samples,d_threshold).");
    if(nlhs!=1) 
        mexErrMsgTxt("One output required (plane).");
    
    // check for proper argument types and sizes
    if(!mxIsSingle(prhs[0]) || mxGetNumberOfDimensions(prhs[0])!=2)
      mexErrMsgTxt("Input D must be a float disparity image.");
    if(!mxIsInt32(prhs[1]) || mxGetM(prhs[1])!=1 || mxGetN(prhs[1])!=4)
      mexErrMsgTxt("Input roi must be a int32 1x4 matrix.");
    if(!mxIsInt32(prhs[2]) || mxGetN(prhs[2])*mxGetM(prhs[2])!=1)
      mexErrMsgTxt("Input num_samples must be a int32 scalar.");
    if(!mxIsSingle(prhs[3]) || mxGetN(prhs[3])*mxGetM(prhs[3])!=1)
      mexErrMsgTxt("Input d_threshold must be a float scalar.");
    
    // get pointers
    float*   D           =     (float*)mxGetPr(prhs[0]);
    int32_t* roi         =   (int32_t*)mxGetPr(prhs[1]);
    int32_t  num_samples = *((int32_t*)mxGetPr(prhs[2]));
    float    d_threshold =   *((float*)mxGetPr(prhs[3]));
    const int *dims      =  mxGetDimensions(prhs[0]);

    // create output
    const int plane_dims[] = {1,3};
    plhs[0]      = mxCreateNumericArray(2,plane_dims,mxSINGLE_CLASS,mxREAL);
    float* plane = (float*)mxGetPr(plhs[0]);
    
    srand(0);

    // do computation
    computeRoadPlaneEstimate(D,dims,plane,roi,num_samples,d_threshold);
}
