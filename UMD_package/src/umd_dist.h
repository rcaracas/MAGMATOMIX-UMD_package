/*
 * Author: Laurent Gilquin
 * Date:   August 26, 2020
 */ 

static NPY_INLINE double
umd_elem(const double *coeff, const double *u, const double *v,
         const npy_intp ncols)
{
    double sum = 0.0;
    npy_intp it;

    for (it = 0; it < ncols; ++it) {
        double dist = u[it] - v[it];
        int shift = 2.0 * dist / coeff[it];
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
        for (j = i + 1; j < nrows; ++j, ++res) {
        
            const double *v = X + (ncols * j);
            *res = umd_elem(coeff, u, v, ncols);
        }
    }
    return 0;
}