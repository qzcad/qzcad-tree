#ifndef CONJUGATEGRADIENT_HPP
#define CONJUGATEGRADIENT_HPP

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <boost/numeric/ublas/operation.hpp>

#include <iostream>


using namespace boost::numeric::ublas;

/** \brief compute residual r = b - Ax
 *
 */
template<class V, class T1, class IA1, class TA1, class E2>
void
residual ( const compressed_matrix<T1, row_major, 0, IA1, TA1> & A,
           const vector_expression<E2> &x,
           const vector_expression<E2> &b,
           V &r, row_major_tag) {

    BOOST_UBLAS_CHECK( x().size() >= A.size2(), external_logic() );
    BOOST_UBLAS_CHECK( b().size() >= A.size1(), external_logic() );
    BOOST_UBLAS_CHECK( r.size() >= A.size1(), external_logic() );

    typedef typename V::size_type size_type;
    typedef typename V::value_type value_type;

    for (size_type i = 0; i < A.size1 (); ++ i) {
        size_type begin = A.index1_data () [i];
        size_type end = A.index1_data () [i + 1];
        value_type t (b() (i));
        for (size_type j = begin; j < end; ++ j)
            t -= A.value_data () [j] * x () (A.index2_data () [j]);
        r (i) = t;
    }
    return;
}

template<class V, class T1, class IA1, class TA1, class E2>
void
sps_prod ( const compressed_matrix<T1, row_major, 0, IA1, TA1> & A,
           const vector_expression<E2> &x,
           V &r, row_major_tag) {

    BOOST_UBLAS_CHECK( x().size() >= A.size2(), external_logic() );
    BOOST_UBLAS_CHECK( r.size() >= A.size1(), external_logic() );

    typedef typename V::size_type size_type;
    typedef typename V::value_type value_type;
#if USE_OMP
#pragma omp parallel for
#endif
    for (size_type i = 0; i < A.size1 (); ++ i) {
        size_type begin = A.index1_data () [i];
        size_type end = A.index1_data () [i + 1];
        value_type t (0);
        for (size_type j = begin; j < end; ++ j)
            t += A.value_data () [j] * x () (A.index2_data () [j]);
        r (i) = t;
    }
    return;
}

template<class T1, class IA1, class TA1, class E2>
int conjugateGradient(const compressed_matrix<T1, row_major, 0, IA1, TA1> & A,
             vector<E2> &x,
             const vector<E2> &b,
             row_major_tag,
             double tol = 0.000001,
             int niter = 30000,
             bool messages = true)
{
    int i;
    double alpha, beta, residn;

    vector<E2> resid (b.size()) ;
    vector<E2> d ;            // search direction
    vector<E2> resid_old;
    vector<E2> temp (b.size()) ;

    for (typename vector<E2>::size_type i = 0; i < x.size(); i++) x[i] = 0.0;

    residual(A, x, b,resid, row_major_tag());

    d = resid;
    if (messages) std::cout << "Невязка: " << norm_2(resid) << std::endl;
    // CG loop
    for(i = 1; i <= niter; i++)
    {
        sps_prod (A, d, temp, row_major_tag());
        alpha = inner_prod(resid,resid)/inner_prod(d,temp);
        x += (d*alpha);
        resid_old = resid;
        resid -= (temp*alpha);
        residn = norm_2(resid);
        if(residn <= tol)
        {
            break;
        }
        if (messages && (i % 100 == 0)) std::cout << "Невязка: " << residn << std::endl;
        beta = inner_prod(resid,resid)/inner_prod(resid_old,resid_old);
        d = resid + d*beta;
    }

    return i;
}

#endif // CONJUGATEGRADIENT_HPP
