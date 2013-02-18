module Elliptic

export E, F, K

include("slatec.jl")

function E(phi::Float64, m::Float64)
    if m < 0. || m > 1. throw(DomainError()) end
    sinphi = sin(phi)
    sinphi2 = sinphi^2
    cosphi2 = 1. - sinphi2
    y = 1. - m*sinphi2
    drf,ierr1 = SLATEC.DRF(cosphi2, y, 1.)
    drd,ierr2 = SLATEC.DRD(cosphi2, y, 1.)
    @assert ierr1 == 0 && ierr2 == 0
    sinphi*(drf - m*sinphi2*drd/3)
end
E(phi::Real, m::Real) = E(float64(phi), float64(m))
@vectorize_2arg Real E

function E(m::Float64)
    if m < 0. || m > 1. throw(DomainError()) end
    if m == 1. return 1. end
    y = 1. - m
    drf,ierr1 = SLATEC.DRF(0., y, 1.)
    drd,ierr2 = SLATEC.DRD(0., y, 1.)
    @assert ierr1 == 0 && ierr2 == 0
    drf - m*drd/3
end
E(x::Float32) = float32(E(float64(x)))
E(x::Real) = E(float64(x))
@vectorize_1arg Real E

function F(phi::Float64, m::Float64)
    if m < 0. || m > 1. throw(DomainError()) end
    sinphi = sin(phi)
    drf,ierr = SLATEC.DRF(cos(phi)^2, 1. - m*sinphi^2, 1.)
    @assert ierr == 0
    sinphi*drf
end
F(phi::Real, m::Real) = F(float64(phi), float64(m))
@vectorize_2arg Real F

function K(m::Float64)
    if m < 0. || m > 1. throw(DomainError()) end
    if m == 1. return Inf end
    drf,ierr = SLATEC.DRF(0., 1. - m, 1.)
    @assert ierr == 0
    drf
end
K(x::Float32) = float32(K(float64(x)))
K(x::Real) = K(float64(x))
@vectorize_1arg Real K

end # module
