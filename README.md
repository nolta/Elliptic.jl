
Elliptic Special Functions
==========================

This julia package provides the following:

Elliptic Integrals
------------------

<table>
  <tr>
    <td>K(m)</td>
    <td>Complete elliptic integral of the first kind, K(m) = F(π/2,m)</td>
  </tr>
  <tr>
    <td>F(phi, m)</td>
    <td>Incomplete elliptic integral of the first kind</td>
  </tr>
  <tr>
    <td>E(m)</td>
    <td>Complete elliptic integral of the second kind, E(m) = E(π/2,m)</td>
  </tr>
  <tr>
    <td>E(phi, m)</td>
    <td>Incomplete elliptic integral of the second kind</td>
  </tr>
  <tr>
    <td>Pi(n, phi, m)</td>
    <td>Incomplete elliptic integral of the third kind</td>
  </tr>
</table>

The parameter `m = k^2 = sin(α)^2` where `α` is the modular angle.

```jlcon
julia> import Elliptic

julia> julia> Elliptic.K(0.5)
1.854074677301372
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

Matlab Compatibility
--------------------

For convenience, the matlab compatible `ellipj` and `ellipke` routines are
also provided. `ellipj(u,m)` is equivalent to `sn(u,m), cn(u,m), dn(u,m)`,
but faster if you want all three. Likewise, `ellipke(m)` is equivalent to
`K(m), E(M)`, but faster if you want both.

```jlcon
julia> import Elliptic

julia> k,e = Elliptic.ellipke(0.5)
(1.854074677301372,1.3506438810476757)

julia> sn,cn,dn = Elliptic.ellipj(0.672, 0.36)
0.6095196917919022,0.792770928653356,0.9307281387786907)
```

Installation
------------

```jlcon
julia> Pkg.update()

julia> Pkg.add("Elliptic")
```

Details
-------

![F(\phi|m) = \int_0^\phi d\theta (1 - m\sin^2\theta)^{-1/2}](http://mathurl.com/akv49po.png)

![E(\phi|m) = \int_0^\phi d\theta (1 - m\sin^2\theta)^{1/2}](http://mathurl.com/amde52p.png)

![\Pi(n;\varphi|m) = \int_0^\varphi d\theta\, (1-n\sin^2\theta)^{-1}(1 - m\sin^2\theta)^{-1/2}](http://mathurl.com/bzsx5tw.png)
