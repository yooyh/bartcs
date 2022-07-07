#include <Rcpp.h>

using namespace Rcpp;

#ifdef _OPENMP
    #include <omp.h>
#endif

//' Count number of threads for parallel computation
//'
//' \code{count_omp_thread} counts the number of threads for
//' parallel computation using openMP.
//' If it returns 1, openMP is not viable.
//'
//' @export
// [[Rcpp::export]]
int count_omp_thread() 
{
    int n = 0;

    #ifdef _OPENMP
        #pragma omp parallel reduction(+:n)
    #endif

    n += 1;

    return n;
}
