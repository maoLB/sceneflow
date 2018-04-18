#include <cstdlib>
#include <stdio.h>
#include <limits>
#include <time.h>
#include <iostream>

#include "TRW_S-v1.3/MRFEnergy.h"

#include "mex.h"

#define MATLAB_ASSERT(expr,msg) if (!(expr)) { mexErrMsgTxt(msg);}

#if !defined(MX_API_VER) || MX_API_VER < 0x07030000
typedef int mwSize;
typedef int mwIndex;
#endif

using namespace std;

MRFEnergy<TypeGeneral>::Options getOptions (const mxArray *oInPtr) {
    
    //prepare default options
    MRFEnergy<TypeGeneral>::Options options;
    options.m_eps = 1e-3;
    options.m_iterMax = 500;
    options.m_printIter = 5;
    options.m_printMinIter = 10;
    
    //get options structure
    if(oInPtr != NULL){
        MATLAB_ASSERT(mxIsStruct(oInPtr), "Expected structure array for options");
        MATLAB_ASSERT(mxGetNumberOfElements(oInPtr) == 1, "Wrong size of options structure: expected 1");
        mxArray *curField = NULL;
        if((curField = mxGetField(oInPtr, 0, "method")) != NULL){
            MATLAB_ASSERT(mxGetClassID(curField) == mxCHAR_CLASS, "Wrong structure type for options: expected STRING for field <<method>>");
            mwSize buflen = mxGetN(curField)*sizeof(mxChar)+1;
            char *buf = (char*)mxMalloc(buflen);
            if(!mxGetString(curField, buf, buflen)){
                if(!strcmp(buf, "trw-s")) options.m_method = 0;
                if(!strcmp(buf, "bp")) options.m_method = 1;
            }
            mxFree(buf);
        }
        if((curField = mxGetField(oInPtr, 0, "maxIter")) != NULL){
            MATLAB_ASSERT(mxGetClassID(curField) == mxDOUBLE_CLASS, "Wrong structure type for options: expected DOUBLE for field <<maxIter>>");
            MATLAB_ASSERT(mxGetNumberOfElements(curField) == 1, "Wrong structure type for options: expected 1 number for field <<maxIter>>");
            options.m_iterMax = (int)(*(double*)mxGetData(curField));
            MATLAB_ASSERT(options.m_iterMax >= 1, "Wrong value for options.maxIter: expected value is >= 1");
        }
        if((curField = mxGetField(oInPtr, 0, "verbosity")) != NULL){
            MATLAB_ASSERT(mxGetClassID(curField) == mxDOUBLE_CLASS, "Wrong structure type for options: expected DOUBLE for field <<verbosity>>");
            MATLAB_ASSERT(mxGetNumberOfElements(curField) == 1, "Wrong structure type for options: expected 1 number for field <<verbosity>>");
            options.m_verbosityLevel = (int)(*(double*)mxGetData(curField));
            MATLAB_ASSERT(options.m_verbosityLevel == 0 || options.m_verbosityLevel == 1 || options.m_verbosityLevel == 2, "Wrong value for options.verbosity: expected value is 0, 1, or 2");
        }
        if((curField = mxGetField(oInPtr, 0, "funcEps")) != NULL){
            MATLAB_ASSERT(mxGetClassID(curField) == mxDOUBLE_CLASS, "Wrong structure type for options: expected DOUBLE for field <<funcEps>>");
            MATLAB_ASSERT(mxGetNumberOfElements(curField) == 1, "Wrong structure type for options: expected 1 number for field <<funcEps>>");
            options.m_eps = *(double*)mxGetData(curField);
        }
        if((curField = mxGetField(oInPtr, 0, "printMinIter")) != NULL){
            MATLAB_ASSERT(mxGetClassID(curField) == mxDOUBLE_CLASS, "Wrong structure type for options: expected DOUBLE for field <<printMinIter>>");
            MATLAB_ASSERT(mxGetNumberOfElements(curField) == 1, "Wrong structure type for options: expected 1 number for field <<printMinIter>>");
            options.m_printMinIter = (int)(*(double*)mxGetData(curField));
            MATLAB_ASSERT(options.m_printMinIter >= 0, "Wrong value for options.printMinIter: expected value is >= 0");
        }
        if((curField = mxGetField(oInPtr, 0, "printIter")) != NULL){
            MATLAB_ASSERT(mxGetClassID(curField) == mxDOUBLE_CLASS, "Wrong structure type for options: expected DOUBLE for field <<printIter>>");
            MATLAB_ASSERT(mxGetNumberOfElements(curField) == 1, "Wrong structure type for options: expected 1 number for field <<printIter>>");
            options.m_printIter = (int)(*(double*)mxGetData(curField));
            MATLAB_ASSERT(options.m_printIter >= 1, "Wrong value for options.printIter: expected value is >= 1");
        }
    }
    
    return options;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    MATLAB_ASSERT(nrhs >= 2 , "Not enough input arguments, expected 2 or 3" );
    MATLAB_ASSERT(nrhs <= 3, "Too many input arguments, expected 2 or 3");
    
    //Fix input parameter order:
    const mxArray *vInPtr = (nrhs > 0) ? prhs[0] : NULL; // variables
    const mxArray *fInPtr = (nrhs > 1) ? prhs[1] : NULL; // factors
    const mxArray *oInPtr = (nrhs > 2) ? prhs[2] : NULL; // options
    
    int32_t numVariables = mxGetN(vInPtr);
    int32_t numFactors = mxGetN(fInPtr);
    double* dim_mex = (double*)mxGetPr(vInPtr);    
    
    //Fix output parameter order:
    mxArray **sOutPtr  = (nlhs > 0) ? &plhs[0] : NULL; // solution
    mxArray **eOutPtr  = (nlhs > 1) ? &plhs[1] : NULL; // energy
    mxArray **lbOutPtr = (nlhs > 2) ? &plhs[2] : NULL; // lowerbound
    
    // get options
    MRFEnergy<TypeGeneral>::Options options = getOptions(oInPtr);

    // create MRF object
    MRFEnergy<TypeGeneral>* mrf = new MRFEnergy<TypeGeneral>(TypeGeneral::GlobalSize());
    MRFEnergy<TypeGeneral>::NodeId* nodes = new MRFEnergy<TypeGeneral>::NodeId[numVariables];
    TypeGeneral::REAL energy, lowerBound;
    
    // init nodes to 0
    for (int32_t v_id=0; v_id<numVariables; v_id++)
        nodes[v_id] = 0;
    
    // add unary factors
    for (int32_t f=0; f<numFactors; f++) {
        
        mxArray *f_cell  = mxGetCell(fInPtr,f);
        mxArray *v_field = mxGetField(f_cell,0,"v"); // variables
        mxArray *e_field = mxGetField(f_cell,0,"e"); // dense
        double* v_mex = (double*)mxGetPr(v_field);
        double* e_mex = (double*)mxGetPr(e_field);
        int32_t n = mxGetN(e_field);
        int32_t n_var_fac = mxGetN(v_field);
        
        // unary factor
        if (n_var_fac==1) {          
            int32_t v_id = (int32_t)v_mex[0]-1;
            int32_t numLabels = dim_mex[v_id];
            if (!nodes[v_id])
                nodes[v_id] = mrf->AddNode(TypeGeneral::LocalSize(numLabels),TypeGeneral::NodeData(e_mex));
        }
    }
    
    // add uninitialized unary nodes
    for (int32_t v_id=0; v_id<numVariables; v_id++) {
        if (!nodes[v_id]) {
            int32_t numLabels = dim_mex[v_id];
            double* e_mex = new double[numLabels];
            for (int32_t i=0; i<numLabels; i++)
                e_mex[i] = 0;
            nodes[v_id] = mrf->AddNode(TypeGeneral::LocalSize(numLabels),TypeGeneral::NodeData(e_mex));
            delete [] e_mex;
        }
    }
    
    // add pairwise factors
    for (int32_t f=0; f<numFactors; f++) {
        
        mxArray *f_cell  = mxGetCell(fInPtr,f);
        mxArray *v_field = mxGetField(f_cell,0,"v"); // variables
        mxArray *e_field = mxGetField(f_cell,0,"e"); // dense
        double* v_mex = (double*)mxGetPr(v_field);
        double* e_mex = (double*)mxGetPr(e_field);
        int32_t n = mxGetN(e_field);
        int32_t n_var_fac = mxGetN(v_field);

        // pairwise factor
        if (n_var_fac==2) {
            int32_t v_id1 = (int32_t)v_mex[0]-1;
            int32_t v_id2 = (int32_t)v_mex[1]-1;
            int32_t nl1 = dim_mex[v_id1];
            int32_t nl2 = dim_mex[v_id2];
            TypeGeneral::REAL *P = new TypeGeneral::REAL[nl1*nl2];
            for (int32_t i=0; i<nl1*nl2; i++)
                P[i] = e_mex[i];
            mrf->AddEdge(nodes[v_id1],nodes[v_id2],TypeGeneral::EdgeData(TypeGeneral::GENERAL,P));
            delete [] P;
        }
    }
    
    // verbosity
    if (options.m_verbosityLevel < 2)
        options.m_printMinIter = options.m_iterMax + 2;
    
    // run TRW-S
    clock_t tStart = clock();
    mrf->SetAutomaticOrdering();
    if(options.m_method == 0) {
        mrf->Minimize_TRW_S(options, lowerBound, energy);        
        if(options.m_verbosityLevel >= 1)
            printf("TRW-S finished. Time: %f\n", (clock() - tStart) * 1.0 / CLOCKS_PER_SEC);
    }
    else {   
        mrf->Minimize_BP(options, energy);
        lowerBound = std::numeric_limits<double>::signaling_NaN();
        if(options.m_verbosityLevel >= 1)
            printf("BP finished. Time: %f\n", (clock() - tStart) * 1.0 / CLOCKS_PER_SEC);
    }
    
    // output the best solution
    if(sOutPtr != NULL)	{
        *sOutPtr = mxCreateNumericMatrix(1, numVariables, mxDOUBLE_CLASS, mxREAL);
        double* segment = (double*)mxGetData(*sOutPtr);
        for(int i = 0; i < numVariables; ++i)
            segment[i] = (double)(mrf -> GetSolution(nodes[i])) + 1;
    }
    
    // output the best energy value
    if(eOutPtr != NULL)	{
        *eOutPtr = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
        *(double*)mxGetData(*eOutPtr) = (double)energy;
    }
    
    // output the best lower bound
    if(lbOutPtr != NULL) {
        *lbOutPtr = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
        *(double*)mxGetData(*lbOutPtr) = (double)lowerBound;
    }
    
    // done
    delete [] nodes;
    delete mrf;
}

