Elliptic Special Functions for Julia
====================================

[![Build Status](https://travis-ci.org/nolta/Elliptic.jl.svg?branch=master)](https://travis-ci.org/nolta/Elliptic.jl)

This julia package provides the following:

Elliptic Integrals
------------------

<table>
  <tr>
    <td>K(m)</td>
    <td>Complete elliptic integral of the first kind, K(m) = F(π/2|m)</td>
  </tr>
  <tr>
    <td>F(phi, m)</td>
    <td>Incomplete elliptic integral of the first kind, F(φ|m)</td>
  </tr>
  <tr>
    <td>E(m)</td>
    <td>Complete elliptic integral of the second kind, E(m) = E(π/2|m)</td>
  </tr>
  <tr>
    <td>E(phi, m)</td>
    <td>Incomplete elliptic integral of the second kind, E(φ|m)</td>
  </tr>
  <tr>
    <td>Pi(n, phi, m)</td>
    <td>Incomplete elliptic integral of the third kind, Π(n;φ|m)</td>
  </tr>
</table>

The parameter `m = k^2 = sin(α)^2` where `α` is the modular angle and `k` is the modulus.

```jlcon
julia> import Elliptic

julia> Elliptic.K(0.5)
1.854074677301372
```

Jacobi Elliptic Functions
-------------------------

<table>
  <tr>
    <td>am(u, m)</td>
    <td>Jacobi amplitude, defined by u = F(am(u|m)|m)</td>
  </tr>
  <tr>
    <td>sn(u, m)</td>
    <td>Jacobi elliptic function, sn(u|m) = sin(am(u|m))</td>
  </tr>
  <tr>
    <td>cn(u, m)</td>
    <td>Jacobi elliptic function, cn(u|m) = cos(am(u|m))</td>
  </tr>
  <tr>
    <td>dn(u, m)</td>
    <td>Jacobi elliptic function, dn(u|m) = sqrt(1 - m sn(u|m)^2)</td>
  </tr>

  <tr>
    <td>cd(u, m)</td>
    <td>Jacobi elliptic function, cd(u|m)</td>
  </tr>
  <tr>
    <td>sd(u, m)</td>
    <td>Jacobi elliptic function, sd(u|m)</td>
  </tr>
  <tr>
    <td>nd(u, m)</td>
    <td>Jacobi elliptic function, nd(u|m) = 1/dn(u|m)</td>
  </tr>

  <tr>
    <td>dc(u, m)</td>
    <td>Jacobi elliptic function, dc(u|m)</td>
  </tr>
  <tr>
    <td>nc(u, m)</td>
    <td>Jacobi elliptic function, nc(u|m) = 1/cn(u|m)</td>
  </tr>
  <tr>
    <td>sc(u, m)</td>
    <td>Jacobi elliptic function, sc(u|m)</td>
  </tr>

  <tr>
    <td>ns(u, m)</td>
    <td>Jacobi elliptic function, ns(u|m) = 1/sn(u|m)</td>
  </tr>
  <tr>
    <td>ds(u, m)</td>
    <td>Jacobi elliptic function, ds(u|m)</td>
  </tr>
  <tr>
    <td>cs(u, m)</td>
    <td>Jacobi elliptic function, cs(u|m)</td>
  </tr>
</table>

```jlcon
julia> import Elliptic.Jacobi

julia> Jacobi.sn(2, 9)
-0.15028246569211734
```

Matlab Compatibility
--------------------

<table>
  <tr>
    <td>ellipj(u, m)</td>
    <td>returns (sn(u,m), cn(u,m), dn(u,m))</td>
  </tr>
  <tr>
    <td>ellipke(m)</td>
    <td>returns (K(m), E(m))</td>
  </tr>
</table>

For convenience, the matlab compatible `ellipj` and `ellipke` routines are
also provided. `ellipj(u,m)` is equivalent to `sn(u,m), cn(u,m), dn(u,m)`,
but faster if you want all three. Likewise, `ellipke(m)` is equivalent to
`K(m), E(m)`, but faster if you want both.

```jlcon
julia> import Elliptic

julia> k,e = Elliptic.ellipke(0.5)
(1.854074677301372,1.3506438810476757)

julia> sn,cn,dn = Elliptic.ellipj(0.672, 0.36)
(0.6095196917919022,0.792770928653356,0.9307281387786907)
```

Installation
------------

```jlcon
julia> Pkg.update()

julia> Pkg.add("Elliptic")
```

Definitions
-----------

![F(\phi|m) = \int_0^\phi d\theta (1 - m\sin^2\theta)^{-1/2}](http://mathurl.com/av9eou5.png)

![E(\phi|m) = \int_0^\phi d\theta (1 - m\sin^2\theta)^{1/2}](http://mathurl.com/al2zsok.png)

![\Pi(n;\varphi|m) = \int_0^\varphi d\theta\, (1-n\sin^2\theta)^{-1}(1 - m\sin^2\theta)^{-1/2}](http://mathurl.com/bzsx5tw.png)
