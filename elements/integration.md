---
title: "Numerical Integration"
---

# Numerical Integration

The numerical integration schemes can be found in Stenger, F. Approximate Calculation of Multiple Integrals (A. H. Stroud). SIAM Rev. 1973, 15 (1), 234–235. https://doi.org/10.1137/1015023. 

The properties of an integration scheme includes the following: 

1. region of integration

2. degree of the scheme

3. number of integration points

## Region of integration

We follow the notations for regions used in Stenger's book. These are the notations for some commonly used regions of integration in finite element analysis. 

| Regions of integration | Notation |
| --------------------------- | -------- |
| Line | `C1` |
| Triangle | `T2` |
| Quadrilateral | `C2` |
| Tetrahedron | `T3` |
| Wedge | `T2+C1` |
| Hexahedron | `C3` |

For the wedge region, the integration is treated in a hybrid form: the 3D integral is decoupled into a 2D integration over the triangular cross-section `T2` and a 1D integration along the line `C1` in the out-of-plane direction.

## Degree of the scheme and number of integration points

The degree of an integration scheme is the highest polynomial degree it integrates exactly. A scheme of degree n yields exact results for all polynomial integrands of total degree ≤ n. 

Achieving a higher degree of accuracy in general requires more integration points, because more information about the integrand must be sampled. However, increasing the number of points also increases the computational cost of the finite element analysis. Therefore, in practice we aim to choose, for a given region of integration, an integration scheme that achieves the required degree of accuracy with as few integration points as possible.

For example, for 2nd-order 3D continuum elements, the external pressure loading is integrated over the element surface. Expressed in parametric coordinates $(\xi, \eta)$:
$$F^{ext} = \iint p \, N^T \left(\frac{\partial \mathbf{x}}{\partial \xi} \times \frac{\partial \mathbf{x}}{\partial \eta}\right) d\xi \, d\eta$$
where $\mathbf{x}(\xi,\eta) = N_I(\xi,\eta)\,\mathbf{x}_I$ is the position interpolation. For 2nd-order elements, $N$ is quadratic (degree 2), and the surface normal term involves products of shape function derivatives (degree 1 each), giving a total integrand degree of 4. Thus, the integration scheme must have degree $\geq 4$. 

Similarly, the mass matrix for higher-order elements requires integrating products of shape functions ($N^T N$), so for 2nd-order elements the integrand is of degree 4, again requiring a scheme of degree $\geq 4$.
