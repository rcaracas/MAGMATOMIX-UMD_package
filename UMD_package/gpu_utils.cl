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

    FREAL acc = 0.0f, dist = 0.0f, k_coeff = 0.0f;
    int shift = 0;
    if ((globalRow < globalCol) && (globalRow < M) && (globalCol < M)) {
        // loop over coordinates space
        for (FINT k = 0; k < K; k++) {
            // dimension-wise discrepancy 
            dist = vload_half(globalRow * K + k, X) - vload_half(globalCol * K + k, X);
            // substract the shift
            k_coeff = vload_half(k, coeff);
            shift = CST * dist / k_coeff;
            dist += - shift * k_coeff;
            // square
            acc += dist * dist;
        }
        // store the result
        vstore_half(sqrt(acc), k_index, res);
    }
}


// Calculation of umd_pdist with half (16-bit floating-point) precision and with
// local memory storage.
__kernel void half_tiled_pdist(const FINT M, const __global half* X,
                               __global half* res, const __global half* coeff)
{    
    // thread identifiers
    const FINT row = get_local_id(0); // local row index in {0, ..., TS-1}
    const FINT col = get_local_id(1); // local col index in {0, ..., TS-1}
    const FINT globalRow = TS * get_group_id(0) + row; // row index of X in {0, ..., M-1}
    const FINT globalCol = TS * get_group_id(1) + col; // col index of X in {0, ..., M-1}
     // local memory to fill tiles of size TS * K
    __local FREAL_IN Asub[TS][K];
    __local FREAL_IN tAsub[TS][K];
    // 2d upper triangular indices mapped to 1dg
    const FINT k_index = M * (M - 1) / 2 - (M - globalRow)*((M - globalRow) - 1) / 2 +
                            globalCol - globalRow - 1;

    FREAL acc = 0.0f, dist = 0.0f, k_coeff = 0.0f;
    int shift = 0;
    
    if ((globalRow < globalCol) && (globalRow < M) && (globalCol < M)) {
    
        for (FINT k = 0; k < K; k++) {
        
            Asub[row][k] = vload_half(globalRow * K + k, X);
            tAsub[col][k] = vload_half(globalCol * K + k, X);
        }

        // loop over coordinates space
        for (FINT k = 0; k < K; k++) {
            // dimension-wise discrepancy 
            dist = Asub[row][k] - tAsub[col][k];
            // substract the shift
            k_coeff = vload_half(k, coeff);
            shift = CST * dist / k_coeff;
            dist += - shift * k_coeff;
            // square
            acc += dist * dist;
        }
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

    FREAL acc = 0.0f, dist = 0.0f, k_coeff = 0.0f;
    int shift = 0;
    if ((globalRow < globalCol) && (globalRow < M) && (globalCol < M)) {
        // loop over coordinates space
        for (FINT k = 0; k < K; k++) {

            // dimension-wise discrepancy 
            dist = X[globalRow * K + k] - X[globalCol * K + k];
            // substract the shift
            k_coeff = coeff[k];
            shift = CST * dist / k_coeff;
            dist += - shift * k_coeff;
            // square
            acc += dist * dist;
        }
        // store the result
        res[k_index] = sqrt(acc);
    }
}


// Calculation of umd_pdist with float or double precision and with
// local memory storage.
__kernel void tiled_pdist(const FINT M, const __global FREAL* X,
                          __global FREAL* res, const __global FREAL* coeff)
{    
    // thread identifiers
    const FINT row = get_local_id(0); // local row index in {0, ..., TS-1}
    const FINT col = get_local_id(1); // local col index in {0, ..., TS-1}
    const FINT globalRow = TS * get_group_id(0) + row; // row index of X in {0, ..., M-1}
    const FINT globalCol = TS * get_group_id(1) + col; // col index of X in {0, ..., M-1}
     // local memory to fill tiles of size TS * K
    __local FREAL_IN Asub[TS][K];
    __local FREAL_IN tAsub[TS][K];
    // 2d upper triangular indices mapped to 1dg
    const FINT k_index = M * (M - 1) / 2 - (M - globalRow)*((M - globalRow) - 1) / 2 +
                            globalCol - globalRow - 1;

    FREAL acc = 0.0f, dist = 0.0f, k_coeff = 0.0f;
    int shift = 0;
    
    if ((globalRow < globalCol) && (globalRow < M) && (globalCol < M)) {
    
        for (FINT k = 0; k < K; k++) {
        
            Asub[row][k] = X[globalRow * K + k];
            tAsub[col][k] = X[globalCol * K + k];
        }

        // loop over coordinates space
        for (int k = 0; k < K; k++) {
            // dimension-wise discrepancy 
            dist = Asub[row][k] - tAsub[col][k];
            // substract the shift
            k_coeff = coeff[k];
            shift = CST * dist / k_coeff;
            dist += - shift * k_coeff;
            // square
            acc += dist * dist;
        }
        // store the result
        res[k_index] = sqrt(acc);
    }
}
