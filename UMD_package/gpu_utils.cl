#define CST 2.0f

// Calculation of umd_pdist with half (16-bit floating-point) precision 
__kernel void half_pdist(const FINT M, const __global half* X,
                         __global half* res, const __global half* coeff)
{
    // thread identifiers
    const FINT globalRow = get_global_id(0); // row index of res in {0, ..., M-1}
    const FINT globalCol = get_global_id(1); // col index of res in {0, ..., M-1}
    // 2d upper triangular indices mapped to 1d
    const FINT k_index = M *(M-1)/2 - (M-globalRow)*((M-globalRow)-1)/2 +
                         globalCol - globalRow - 1;

    // store coeffs in local memory
    __local FREAL lcoeff[K];
    lcoeff[0] = vload_half(0, coeff);
    lcoeff[1] = vload_half(1, coeff);
    lcoeff[2] = vload_half(2, coeff);

    FREAL acc = 0.0f, dist = 0.0f;
    int shift = 0;
    if ((globalRow < globalCol) && (globalRow < M) && (globalCol < M)) {

        // dimension-wise discrepancy 
        dist = vload_half(globalRow * K, X) - vload_half(globalCol * K, X);
        shift = CST * dist / lcoeff[0];
        dist += - shift * lcoeff[0];
        acc += dist * dist;
        // dimension-wise discrepancy 
        dist = vload_half(globalRow * K + 1, X) - vload_half(globalCol * K + 1, X);
        shift = CST * dist / lcoeff[1];
        dist += - shift * lcoeff[1];
        acc += dist * dist;
        // dimension-wise discrepancy 
        dist = vload_half(globalRow * K + 2, X) - vload_half(globalCol * K + 2, X);
        shift = CST * dist / lcoeff[2];
        dist += - shift * lcoeff[2]; 
        acc += dist * dist;
        // store the result
        vstore_half(sqrt(acc), k_index, res);
    }
}



// Calculation of umd_pdist with float or double precision.
__kernel void pdist(const FINT M, const __global FREAL* X,
                    __global FREAL* res, const __global FREAL* coeff)
{    
    // thread identifiers
    const FINT globalRow = get_global_id(0); // row index of res in {0, ..., M-1}
    const FINT globalCol = get_global_id(1); // col index of res in {0, ..., M-1}
    // 2d upper triangular indices mapped to 1d
    const FINT k_index = M *(M-1)/2 - (M-globalRow)*((M-globalRow)-1)/2 +
                         globalCol - globalRow - 1;
                         
    // store coeffs in local memory
    __local FREAL lcoeff[K];
    lcoeff[0] = coeff[0];
    lcoeff[1] = coeff[1];
    lcoeff[2] = coeff[2];

    FREAL acc = 0.0f, dist = 0.0f;
    int shift = 0;
    if ((globalRow < globalCol) && (globalRow < M) && (globalCol < M)) {

        // dimension-wise discrepancy
        dist = X[globalRow * K] - X[globalCol * K ];
        shift = CST * dist / lcoeff[0];
        dist += - shift * lcoeff[0];
        acc += dist * dist;
        // dimension-wise discrepancy 
        dist = X[globalRow * K + 1] - X[globalCol * K + 1];
        shift = CST * dist / lcoeff[1];
        dist += - shift * lcoeff[1];
        acc += dist * dist;
        // dimension-wise discrepancy 
        dist = X[globalRow * K + 2] - X[globalCol * K + 2];
        shift = CST * dist / lcoeff[2];
        dist += - shift * lcoeff[2]; 
        acc += dist * dist;
        
        // store the result
        res[k_index] = sqrt(acc);
    }
}
