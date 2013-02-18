using Test

include("Elliptic.jl")
using Elliptic

@test_approx_eq E(0) pi/2
@test_approx_eq E(1) 1
@test_approx_eq E(0,0.5) 0
@test_approx_eq E(0.5,0) 0.5
@test_approx_eq E(0.5,1) sin(0.5)

@test_approx_eq F(0,0.5) 0
@test_approx_eq F(0.5,0) 0.5
@test_approx_eq F(0.5,1) log(sec(0.5) + tan(0.5))

@test_approx_eq K(0) pi/2
@test_approx_eq K(0.5) (0.25/sqrt(pi))*gamma(0.25)^2
@test           K(1) == Inf

# values from Abramowitz & Stegun, Table 17.1 (p608-609)
@test_approx_eq_eps E(0.1) 1.53075_7637 1e-9
@test_approx_eq_eps K(0.1) 1.61244_13487_20219 1e-15
K20 = 1.62002_58991_24204
m20 = sind(20)^2
@test_approx_eq_eps K(m20) K20 1e-15

# values from Abramowitz & Stegun, Table 17.5 (p613-615)
F2020 = 0.34988_016
for i = -2:2
    @test_approx_eq_eps F(degrees2radians(20+180i), m20) F2020 + 2i*K20 1e-8
end

# values from Abramowitz & Stegun, Table 17.6 (p616-618)
@test_approx_eq_eps E(degrees2radians(20), sind(20)^2) 0.34825_492 1e-8

# values from Abramowitz & Stegun, Table 17.9 (p625-626)
@test_approx_eq_eps Pi(0.2, degrees2radians(30), sind(30)^2) 0.53896 1e-5
