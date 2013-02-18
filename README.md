
Elliptic: elliptic integral special functions
=============================================

This julia package provides the following functions:

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
