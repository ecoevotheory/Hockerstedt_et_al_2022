/**************************************************************************
 * [infections] = inf_function(S,IP,Q,mutmatrix,G,connections,muP,rho,loci,n,M);
 *
 * Calculates the infection rates for the whole metapopulation
 *************************************************************************/

#include <mex.h>

/*************************************
 * FUNCTION PROTOTYPES
 *************************************/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *S, *IP, *Q, *mutmatrix, *G, *connections, *W, *Y, *parameter, muP, rho, *infections, loci;
    int i, j, k, n, M;
    mwSize dims[3];
            
    /* Allocate inputs */
    S = mxGetPr(prhs[0]);
    IP = mxGetPr(prhs[1]);
    Q = mxGetPr(prhs[2]);
    mutmatrix = mxGetPr(prhs[3]);
    G = mxGetPr(prhs[4]);
    connections = mxGetPr(prhs[5]);
    parameter= mxGetPr(prhs[6]);
    muP= *parameter;
    parameter= mxGetPr(prhs[7]);
    rho= *parameter;
    parameter= mxGetPr(prhs[8]);
    loci= *parameter;
    parameter= mxGetPr(prhs[9]);
    n= (int)*parameter;
    parameter= mxGetPr(prhs[10]);
    M= (int)*parameter;
    
    dims[0] = M;
    dims[1] = n;
    dims[2] = n;
    
    /* Create temporary arrays */
    W = (double *)malloc(n*M*sizeof(double));
    Y = (double *)malloc(n*M*sizeof(double));
    
    /* Create output */
    plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
    infections = mxGetPr(plhs[0]);
    
    /* Mutations */
    for(j=0;j<n;j++){
        for(k=0;k<M;k++){
            W[k + j*M] = IP[k + j*M]*(1-muP*loci);
            for(i=0;i<n;i++){
                if(mutmatrix[i + j*n]>0) W[k + j*M] += muP*IP[k + i*M];
            }
        }
    }
    
    /* Dispersal */
    for(j=0;j<n;j++){
        for(k=0;k<M;k++){
            Y[k + j*M] = W[k + j*M]*(1-rho*connections[k]);
            for(i=0;i<M;i++){
                if(G[i + k*M]>0) Y[k + j*M] += rho*W[i + j*M];
            }
        }
    }
    
    /* Infections */
    for(j=0;j<n;j++){
        for(i=0;i<n;i++){
            for(k=0;k<M;k++){
                infections[k + i*M + j*M*n] = Q[i + j*n]*S[k + i*M]*Y[k + j*M];
            }
        }
    }
    
    /* Free memory */
    free(W);
    free(Y);
    
    return;
}

