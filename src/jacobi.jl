module Jacobi

export am, sn, cn, dn, cd, sd, nd, dc, nc, sc, ns, ds, cs

# Abramowitz & Stegun, section 16.4, p571
const _ambuf = Array(Float64, 10)
function _am(u::Float64, m::Float64, tol::Float64)
    if u == 0. return 0. end

    sqrt_tol = sqrt(tol)
    if m < sqrt_tol
        # A&S 16.13.4
        return u - 0.25*m*(u - 0.5*sin(2.*u))
    end
    m1 = 1. - m
    if m1 < sqrt_tol
        # A&S 16.15.4
        t = tanh(u)
        return asin(t) + 0.25*m1*(t - u*(1. - t^2))*cosh(u)
    end

    a,b,c,n = 1., sqrt(m1), sqrt(m), 0
    while abs(c) > tol
        @assert n < 10
        a,b,c,n = 0.5*(a+b), sqrt(a*b), 0.5*(a-b), n+1
        _ambuf[n] = c/a
    end

    phi = ldexp(a*u, n)
    for i = n:-1:1
        phi = 0.5*(phi + asin(_ambuf[i]*sin(phi)))
    end
    phi
end
_am(u::Float64, m::Float64) = _am(u, m, eps(Float64))
function am(u::Float64, m::Float64, tol::Float64)
    if m < 0. || m > 1. throw(DomainError()) end
    _am(u, m, tol)
end
am(u::Float64, m::Float64) = am(u, m, eps(Float64))
am(u::Real, m::Real) = am(float64(u), float64(m))
@vectorize_2arg Real am

for (f,a,b,c) in ((:sn, :(sin(phi)),                :(sqrtmu1*s), :(sqrt(mu)*sin(phi))),
                  (:cn, :(cos(phi)),                :(cos(phi)),  :(sqrt(1. - mu*sin(phi)^2))),
                  (:dn, :(sqrt(1. - m*sin(phi)^2)), :(1.),        :(cos(phi))))
    @eval begin
        function ($f)(u::Float64, m::Float64)
            # Abramowitz & Stegun, section 16.10, p573
            lt0 = m < 0.
            gt1 = m > 1.
            if !(lt0 || gt1)
                phi = _am(u,m)
                return $a
            elseif lt0
                mu1 = 1./(1. - m)
                mu = -m*mu1
                sqrtmu1 = sqrt(mu1)
                v = u/sqrtmu1
                phi = _am(v,mu)
                s = sin(phi)
                return ($b)/sqrt(1. - mu*s^2)
            elseif gt1
                mu = 1/m
                v = u*sqrt(m)
                phi = _am(v,mu)
                return $c
            end
        end
    end
end

xn = Dict{Symbol, Any}(:s => :(sn(u,m)), :c => :(cn(u,m)), :d => :(dn(u,m)), :n => :(1.))
for (p,num) in xn, (q,den) in xn
    if p == q continue end
    f = symbol(string(p,q))
    if q != :n
        @eval ($f)(u::Float64, m::Float64) = ($num)/($den)
    end
    @eval begin
        ($f)(u::Real, m::Real) = ($f)(float64(u), float64(m))
        @vectorize_2arg Real $f
    end
end

end # module
