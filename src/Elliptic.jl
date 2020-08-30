module Elliptic
# elliptic integrals of 1st/2nd/3rd kind
export E, F, K, Pi

# jacobi elliptic functions
export Jacobi

# matlab compatible
export ellipj, ellipke

include("jacobi.jl")
include("slatec.jl")

function E(phi::Float64, m::Float64)
    if isnan(phi) || isnan(m) return NaN end
    if m < 0. || m > 1. throw(DomainError(m, "argument m not in [0,1]")) end
    if abs(phi) > pi/2
        phi2 = phi + pi/2
        return 2*fld(phi2,pi)*E(m) - _E(cos(mod(phi2,pi)), m)
    end
    _E(sin(phi), m)
end
function _E(sinphi::Float64, m::Float64)
    sinphi2 = sinphi^2
    cosphi2 = 1. - sinphi2
    y = 1. - m*sinphi2
    drf,ierr1 = SLATEC.DRF(cosphi2, y, 1.)
    drd,ierr2 = SLATEC.DRD(cosphi2, y, 1.)
    if ierr1 == ierr2 == 0
        return sinphi*(drf - m*sinphi2*drd/3)
    elseif ierr1 == ierr2 == 2
        # 2 - (1+m)*sinphi2 < tol
        return sinphi
    end
    NaN
end
E(phi::Real, m::Real) = E(Float64(phi), Float64(m))

"""
    ellipke(M)

Complete elliptic integrals of first and second kind. Each element of `M`, say `m`, should satisfy `m ≤ 1.0`.
A tuple `(K, E)` is returned. If `M` is an array, then `K` and `E` share the same size with `M`.
"""
function ellipke end

function ellipke(m::Float64)
    if m < 1.
        y = 1. - m
        drf,ierr1 = SLATEC.DRF(0., y, 1.)
        drd,ierr2 = SLATEC.DRD(0., y, 1.)
        @assert ierr1 == 0 && ierr2 == 0
        return drf, drf - m*drd/3
    elseif m == 1.
        return Inf, 1.
    elseif isnan(m)
        return NaN, NaN
    else
        throw(DomainError(m, "argument m not <= 1"))
    end
end
ellipke(x::Real) = ellipke(Float64(x))

function ellipke(M::AbstractArray{T}) where T <: Real
    K = similar(M, Float64)
    E = similar(M, Float64)
    for i in eachindex(M)
        @inbounds K[i], E[i] = ellipke(M[i])
    end
    return K, E
end

E(m::Float64) = ellipke(m)[2]
E(x::Float32) = Float32(E(Float64(x)))
E(x::Real) = E(Float64(x))

# assumes 0 ≤ m ≤ 1
function rawF(sinphi::Float64, m::Float64)
    if abs(sinphi) == 1. && m == 1. return sign(sinphi)*Inf end
    sinphi2 = sinphi^2
    drf,ierr = SLATEC.DRF(1. - sinphi2, 1. - m*sinphi2, 1.)
    @assert ierr == 0
    sinphi*drf
end

function F(phi::Float64, m::Float64)
    if isnan(phi) || isnan(m) return NaN end
    if m < 0. || m > 1. throw(DomainError(m, "argument m not in [0,1]")) end
    if abs(phi) > pi/2
        # Abramowitz & Stegun (17.4.3)
        phi2 = phi + pi/2
        return 2*fld(phi2,pi)*K(m) - rawF(cos(mod(phi2,pi)), m)
    end
    rawF(sin(phi), m)
end
F(phi::Real, m::Real) = F(Float64(phi), Float64(m))

function K(m::Float64)
    if m < 1.
        drf,ierr = SLATEC.DRF(0., 1. - m, 1.)
        @assert ierr == 0
        return drf
    elseif m == 1.
        return Inf
    elseif isnan(m)
        return NaN
    else
        throw(DomainError(m, "argument m not <= 1"))
    end
end
K(x::Float32) = Float32(K(Float64(x)))
K(x::Real) = K(Float64(x))

function Pi(n::Float64, phi::Float64, m::Float64)
    if isnan(n) || isnan(phi) || isnan(m) return NaN end
    if m < 0. || m > 1. throw(DomainError(m, "argument m not in [0,1]")) end
    sinphi = sin(phi)
    sinphi2 = sinphi^2
    cosphi2 = 1. - sinphi2
    y = 1. - m*sinphi2
    drf,ierr1 = SLATEC.DRF(cosphi2, y, 1.)
    drj,ierr2 = SLATEC.DRJ(cosphi2, y, 1., 1. - n*sinphi2)
    if ierr1 == 0 && ierr2 == 0
        return sinphi*(drf + n*sinphi2*drj/3)
    elseif ierr1 == 2 && ierr2 == 2
        # 2 - (1+m)*sinphi2 < tol
        return Inf
    elseif ierr1 == 0 && ierr2 == 2
        # 1 - n*sinphi2 < tol
        return Inf
    end
    NaN
end
Pi(n::Real, phi::Real, m::Real) = Pi(Float64(n), Float64(phi), Float64(m))
Π = Pi

"""
    ellipj(U, M)

Compute the Jacobi elliptic functions `SN`, `CN`, and `DN` evaluated for corresponding 
elements of argument `U` and parameter `M`.
Note that `M` can only take values in range `[0, 1]` (both ends closed).
If both `U` and `M` are arrays, they must be the same size, and each element of the returned tuple `(SN, CN, DN)`
also has the same size. If either `U` or `M` is an array, then the other (a scalar) is broadcasted, and 
each element of the returned tuple `(SN, CN, DN)` has the same size as the array.
"""
function ellipj end

function ellipj(u::Float64, m::Float64, tol::Float64)
    phi = Jacobi.am(u, m, tol)
    s = sin(phi)
    c = cos(phi)
    d = sqrt(1. - m*s^2)
    s, c, d
end
ellipj(u::Float64, m::Float64) = ellipj(u, m, eps(Float64))
ellipj(u::Real, m::Real) = ellipj(Float64(u), Float64(m))

function ellipj(U::AbstractArray{TU}, M::AbstractArray{TM}) where {TU <: Real,TM <: Real}
    @assert size(U) == size(M) "U and M must be the same size"
    SN = similar(U, Float64)
    CN = similar(U, Float64)
    DN = similar(U, Float64)
    for i in eachindex(U)
        @inbounds SN[i], CN[i], DN[i] = ellipj(U[i], M[i])
    end
    return SN, CN, DN
end
ellipj(U::AbstractArray{TU}, m::Real) where TU <: Real = ellipj(U, fill(Float64(m), size(U)))
ellipj(u::Real, M::AbstractArray{TM}) where TM <: Real = ellipj(fill(Float64(u), size(M)), M)
end # module
