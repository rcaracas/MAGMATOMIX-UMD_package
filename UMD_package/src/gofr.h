/*
 * Author: Laurent Gilquin
 * Date:   September 2, 2020
 */ 


static NPY_INLINE npy_intp
        ij_to_k(const npy_intp M, const npy_intp row, const npy_intp col)
{
        npy_intp pos;
        if (row <= col){  
            pos = M * (M + 1) / 2 - (M - row) * (M - row + 1) / 2 + col - row;
        } else {
            pos = M * (M + 1) / 2 - (M - col) * (M - col + 1) / 2 + row - col;
        }

        return pos;
}



static NPY_INLINE npy_intp
        k_to_ij(const npy_intp M, const npy_intp k_index, npy_intp* indexes)
{

        npy_intp k_bis, p_index;

        // recover the 2d (i,j) index
        k_bis = M * (M - 1) / 2 - k_index - 1;
        p_index = floor((sqrt(1.0 + 8.0 * k_bis) - 1.0) / 2.0);
        indexes[0] = M - 2 - p_index;
        indexes[1] = k_index + M - M * (M - 1) / 2 + p_index * (p_index + 1) / 2;
        
        return 0;
}



static NPY_INLINE double
umd_elem(const double *coeff, const double *u, const double *v,
         const npy_intp ncols)
{
    double sum = 0.0;
    npy_intp it;

    for (it = 0; it < ncols; ++it) {
        double dist = u[it] - v[it];
        npy_intp shift = 2.0 * dist / coeff[it];
        dist += - coeff[it] * shift;
        sum += dist * dist;
    }
    return sqrt(sum);
}


static NPY_INLINE int
umd_distance(const double *X, const double *coeff, double *res,
             const npy_intp nrows, const npy_intp ncols)
{
    npy_intp i, j;

    for (i = 0; i < nrows; ++i) {
    
        const double *u = X + (ncols * i);
        for (j = i + 1; j < nrows; ++j) {
        
            const double *v = X + (ncols * j);
            *res = umd_elem(coeff, u, v, ncols);
            ++res;
        }
    }
    return 0;
}


static NPY_INLINE int
compute_gofr(const double *X, npy_intp *res, const double *coeff, const npy_intp *types,
             const double maxlength, const double discrete, const npy_intp ntypes, 
             const npy_intp Xrows, const npy_intp Xcols, const npy_intp Rcols)
{

    npy_intp i, j, k_index, bin;

    for (i = 0; i < Xrows; ++i) {
    
        const double *u = X + (Xcols * i);
        for (j = i + 1; j < Xrows; ++j) {
        
            const double *v = X + (Xcols * j);
            double dist = umd_elem(coeff, u, v, Xcols);
            if (dist < (0.5 * maxlength)) {
            
                k_index = ij_to_k(ntypes, types[i], types[j]);
                bin = dist / discrete;
                
                if (types[i] == types[j]) {
                
                    // we double since we work only with the upper triangular
                    *(res + Rcols * k_index + bin) += 2;
                } else {
                    *(res + Rcols * k_index + bin) += 1;
                }
            }
        }
    }
        
    return 0;
}


static NPY_INLINE int
compute_gofrM(const double *dist, npy_intp *res, const npy_intp *types, const double maxlength,
              const double discrete, const npy_intp ntypes, const npy_intp Xrows,
              const npy_intp Trows, const npy_intp Rcols)
{

    npy_intp i, j, kk, k_index, bin;

    for (kk = 0; kk < Xrows; ++kk) {
    
        // recover the 2d (i,j) index
        npy_intp indexes[2];
        k_to_ij(Trows, kk, indexes);
        i = indexes[0];
        j = indexes[1];
    
        if (dist[kk] < (0.5 * maxlength)) {
            
            k_index = ij_to_k(ntypes, types[i], types[j]);
            bin = dist[kk] / discrete;
                
            if (types[i] == types[j]) {
                
                // we double since we work only with the upper triangular
                *(res + Rcols * k_index + bin) += 2;
            } else {
                *(res + Rcols * k_index + bin) += 1;
            }
        }
    }
   
    return 0;
}


static NPY_INLINE int
print_gofr(const npy_intp *gofr, const npy_intp *types, const npy_intp ntypes,
           const double maxlength, const double discrete, const double normalization,
           const npy_intp Gcols, FILE *fp)
{

    // allocate memory for intgofr
    npy_intp count = ntypes * ntypes;
    double *intgofr = calloc(count, sizeof (double));
    
    // loop over bins
    npy_intp kk, itype, jtype, k_index;
    double smallr, bigr, volshell, volcell, temp;
    for (kk =  0; kk < Gcols; ++kk) {
    
        smallr = kk * discrete;
        bigr = smallr + discrete;
        volshell = 4.0 * M_PI * (pow(bigr, 3) - pow(smallr, 3)) / 3.0;
        volcell = pow(maxlength, 3);
        // write first entry to file
        fprintf(fp, "%g\t", smallr + discrete / 2.0);
        
        for (itype = 0; itype < ntypes; ++itype) {
        
            for (jtype = 0; jtype < ntypes; ++jtype) {
            
                // recover 1D index
                k_index = ij_to_k(ntypes, itype, jtype);
                // evaluate temp value
                temp = *(gofr + k_index * Gcols + kk) * volcell /
                (volshell * *(types + itype) * *(types + jtype) * normalization);
                // cumulative integral of gofr in spherical coordinates
                double *pos = intgofr + itype * ntypes + jtype;
                *pos += temp * pow(smallr, 2) * discrete * 4.0 * M_PI * *(types + jtype) / volcell;
                // add temp and intogfr to string
                fprintf(fp, "%.15g\t%.15g\t", temp, *pos);
            }
        }
        // add line break
        fputs("\n", fp);
    }

    // de-allocate memory
    free(intgofr);
        
    return 0;
}