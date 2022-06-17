# Optimized 4th Kind Chebyshev Polynomials

## 2D Poisson Results

2D Poisson, $E_x = 2, E_y = 32, N=7, \Omega=[0,1]^2$ with Dirichlet Boundary Conditions.
Two-level multigrid with $N=1$ Galerkin coarse grid solve used to precondition KSP (CG + keeping all projection vectors, so should be O.K. for non-symmetric WASM smoother).
Post smoothing is used, so the Jacobi-based multigrid methods are fully SPD.

### 1st Kind Results, Weighted Additive Schwarz Method (WASM) smoothing:

Smoother| Iterations| Matvecs| Smoothers| Initial Residual| Final Relative Residual  |
--------|-----------|--------|----------|-----------------|--------------------------|
WASM,(7,1), $\omega=0.8413$|54|269|108|6915.9255|7.2393e-11|
Cheby-WASM(1),(7,1)        |54|269|108|6915.9255|7.2393e-11|
Cheby-WASM(2),(7,1)        |31|216|124|6915.9255|6.1281e-11|
Cheby-WASM(3),(7,1)        |23|206|138|6915.9255|6.5132e-11|
Cheby-WASM(4),(7,1)        |20|219|160|6915.9255|2.7355e-11|
Cheby-WASM(5),(7,1)        |18|233|180|6915.9255|2.9309e-11|
Cheby-WASM(6),(7,1)        |16|239|192|6915.9255|3.8236e-11|
Cheby-WASM(7),(7,1)        |15|254|210|6915.9255|3.1146e-11|

### Optimal 4th Kind Results, Weighted Additive Schwarz Method (WASM) smoothing:

Smoother| Iterations| Matvecs| Smoothers| Initial Residual| Final Relative Residual|
--------|-----------|--------|----------|-----------------|--------------------------|
WASM,(7,1), $\omega=0.8413$|54|269|108|6915.9255|7.2393e-11|
Opt. 4th Kind Cheby-WASM(1),(7,1)|59|294|118|6915.9255|6.0384e-11|
Opt. 4th Kind Cheby-WASM(2),(7,1)|33|230|132|6915.9255|8.8919e-11|
Opt. 4th Kind Cheby-WASM(3),(7,1)|23|206|138|6915.9255|6.5069e-11|
Opt. 4th Kind Cheby-WASM(4),(7,1)|18|197|144|6915.9255|6.4798e-11|
Opt. 4th Kind Cheby-WASM(5),(7,1)|15|194|150|6915.9255|2.0337e-11|
Opt. 4th Kind Cheby-WASM(6),(7,1)|13|194|156|6915.9255|1.4384e-11|
Opt. 4th Kind Cheby-WASM(7),(7,1)|11|186|154|6915.9255|2.2537e-11|

### 1st Kind Results, Jacobi smoothing:

Smoother| Iterations| Matvecs| Smoothers| Initial Residual| Final Relative Residual|
--------|-----------|--------|----------|-----------------|--------------------------|
Jacobi,(7,1), $\omega=0.68547$|297|1484|594|6915.9255 | 9.803e-11|
Cheby-Jacobi(1),(7,1)|297|1484|594|6915.9255| 9.803e-11|
Cheby-Jacobi(2),(7,1)|170|1189|680|6915.9255| 9.5788e-11|
Cheby-Jacobi(3),(7,1)|130|1169|780|6915.9255| 8.9169e-11|
Cheby-Jacobi(4),(7,1)|110|1209|880|6915.9255| 8.0867e-11|
Cheby-Jacobi(5),(7,1)|98|1273|980|6915.9255 | 9.4575e-11|
Cheby-Jacobi(6),(7,1)|89|1334|1068|6915.9255| 7.3354e-11|
Cheby-Jacobi(7),(7,1)|82|1393|1148|6915.9255| 7.2088e-11|

### Optimal 4th Kind Results, Jacobi smoothing:

Smoother| Iterations| Matvecs| Smoothers| Initial Residual| Final Relative Residual|
--------|-----------|--------|----------|-----------------|--------------------------|
Jacobi,(7,1), $\omega=0.68547$|297|1484|594|6915.9255     |9.803e-11|
Opt. 4th Kind Cheby-Jacobi(1),(7,1)|328|1639|656|6915.9255|9.7937e-11|
Opt. 4th Kind Cheby-Jacobi(2),(7,1)|182|1273|728|6915.9255|9.8686e-11|
Opt. 4th Kind Cheby-Jacobi(3),(7,1)|129|1160|774|6915.9255|9.3539e-11|
Opt. 4th Kind Cheby-Jacobi(4),(7,1)|100|1099|800|6915.9255|7.8685e-11|
Opt. 4th Kind Cheby-Jacobi(5),(7,1)|83|1078|830|6915.9255 |7.728e-11|
Opt. 4th Kind Cheby-Jacobi(6),(7,1)|72|1079|864|6915.9255 |7.583e-11|
Opt. 4th Kind Cheby-Jacobi(7),(7,1)|63|1070|882|6915.9255 |6.9955e-11|

## nekRS case: 67 pebbles (tet-to-hex mesh, E=122,284, N=7, n=43,385,013 pressure DOFS)

Solver is GMRES(15), A-conjugate residual projection is used with $L=10$.
Tolerance is $10^{-4}$.
$\Delta t = 10^{-4}$, single subcycling step with CFL ~ 2.

Smoother            |  Avg. Iterations, 1st Kind | Avg. Iterations, Opt. 4th Kind | Ratio   |
--------------------|----------------------------|--------------------------------|---------|
Cheb-ASM(1),(7,3,1) |    42.499$\dagger$                 |      45.048$\dagger$                   | 0.943   |
Cheb-ASM(2),(7,3,1) |    29.0455                 |      27.4085                   | 1.06    |
Cheb-ASM(3),(7,3,1) |    24.956                  |      20.6125                   | 1.21    |
Cheb-ASM(4),(7,3,1) |    TBD                     |      TBD                       | TBD     |

Avg. iterations over 2000 steps, 67 pebble case

$\dagger$ only over first 1000 steps

Asymptotically, expect 18% improvement as $k \rightarrow \infty$ (https://arxiv.org/pdf/2202.08830.pdf, bottom of page 15)

Pros:
  - weights can be pre-tabulated, and only depend on k
  - minimal change to existing implementation
  - no additional matvecs, etc. -- can directly compare like-with-like
  - good for relatively large k (k >= 2?)
  - might be able to choose a more aggressive coarsening strategy with a heavier weight smoother(?)
Cons:
  - may be relatively poor for small k

We could probably ammend the issue for small k by choosing the method with the lowest iteration count during the first solve.

Experimental, not production ready nekRS implementation: https://github.com/MalachiTimothyPhillips/nekRS/tree/optimal-chebyshev-polynomials
Code for generating weights (with fix for W), stolen from ArXiv paper: https://github.com/MalachiTimothyPhillips/optimal-fourth-kind-chebyshev-weights

Next steps?

AMG Solver integration:
1. Trilinos/MueLu (https://github.com/trilinos/Trilinos/blob/0bd35b1d0298b79ab6e18619e61fa364c84a6267/packages/ifpack2/src/Ifpack2_Details_Chebyshev_def.hpp#L1236)
2. Hypre: (https://github.com/hypre-space/hypre/blob/ac9d7d0d7b43cd3d0c7f24ec5d64b58fbf900097/src/parcsr_ls/par_cheby.c#L219)
