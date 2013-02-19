
Elliptic: elliptic integral special functions
=============================================

This julia package provides the following functions:

Elliptic Integrals
------------------

```
F(phi, m)

    Incomplete elliptic integral of the first kind:

    F(\phi|m) = \int_0^\phi d\theta (1-m\sin^2\theta)^{-1/2}.

K(m)

    Complete elliptic integral of the first kind: K(m) = F(pi/2|m).

E(phi, m)

    Incomplete elliptic integral of the second kind:

    E(\phi|m) = \int_0^\phi d\theta (1 - m\sin^2\theta)^{1/2}.

E(m)

    Complete elliptic integral of the second kind: E(m) = E(pi/2|m).

```

Jacobi Elliptic Functions
-------------------------

All 12 `pq(u|m)` functions where p, q are chosen without replacement from
the letters s, c, d, n. For example,

```jlcon
julia> import Elliptic

julia> Elliptic.sn(0.672, 0.36)
0.6095196917919022
```

Matlab-Compatible Functions
---------------------------

For convenience, the matlab-compatible `ellipj` and `ellipke` functions are
also provided.
