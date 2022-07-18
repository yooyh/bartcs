#include <Rcpp.h>

using namespace Rcpp;

#ifdef _OPENMP
    #include <omp.h>
#endif

//' Count the number of OpenMP threads for parallel computation
//'
//' \code{count_omp_thread()} counts the number of OpenMP threads for
//' parallel computation.
//' If it returns 1, OpenMP is not viable.
//'
//' @return
//' Number of OpenMP thread(s).
//'
//' @examples
//' count_omp_thread()
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
