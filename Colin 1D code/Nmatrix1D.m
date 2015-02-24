% evaluates shape functions (in physical coordinates) at point xt
 function N = Nmatrix1D(psi)
    N = 1/2 * [1-psi   1+psi];

