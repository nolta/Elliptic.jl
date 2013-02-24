using Test

include("Elliptic.jl")
using Elliptic
using Elliptic.Jacobi

# values from Abramowitz & Stegun, Table 17.1 (p608-609)
#    m            K(m)             E(m)
table17p1 = [
    0.00  1.57079_63267_94897  1.57079_6327;
    0.01  1.57474_55615_17356  1.56686_1942;
    0.10  1.61244_13487_20219  1.53075_7637;
    0.20  1.65962_35986_10528  1.48903_5058;
    0.50  1.85407_46773_01372  1.35064_3881;
    0.90  2.57809_21133_48173  1.10477_4733;
    0.99  3.69563_73629_89875  1.01599_3546;
]

for i = 1:size(table17p1,1)
    m = table17p1[i,1]
    k = table17p1[i,2]
    e = table17p1[i,3]
    @test_approx_eq_eps K(m) k 1e-15
    @test_approx_eq_eps E(m) e 1e-9
    @test_approx_eq_eps F(pi/2,m) k 1e-15
    @test_approx_eq_eps E(pi/2,m) e 1e-9
    a,b = ellipke(m)
    @test_approx_eq_eps a k 1e-15
    @test_approx_eq_eps b e 1e-9

end

@test_approx_eq K(0) pi/2
@test_approx_eq K(0.5) (0.25/sqrt(pi))*gamma(0.25)^2
@test           K(1) == Inf

@test_approx_eq E(0) pi/2
@test_approx_eq E(1) 1

@test_approx_eq E(0,0.5) 0
@test_approx_eq E(0.5,0) 0.5
@test_approx_eq E(0.5,1) sin(0.5)

@test_approx_eq F(0,0.5) 0
@test_approx_eq F(0.5,0) 0.5
@test_approx_eq F(0.5,1) log(sec(0.5) + tan(0.5))

# values from Abramowitz & Stegun, Table 17.2 (p610-611)
m20 = sind(20)^2
K20 = 1.62002_58991_24204
E20 = 1.52379_92052_59774
@test_approx_eq_eps K(m20) K20 1e-15
@test_approx_eq_eps E(m20) E20 1e-15

# values from Abramowitz & Stegun, Table 17.5 (p613-615)
# values from Abramowitz & Stegun, Table 17.6 (p616-618)
E2020 = 0.34825_492
F2020 = 0.34988_016
for i = -2:2
    @test_approx_eq_eps E(degrees2radians(20+180i), m20) E2020 + 2i*E20 1e-8
    @test_approx_eq_eps F(degrees2radians(20+180i), m20) F2020 + 2i*K20 1e-8
end

# values from Abramowitz & Stegun, Table 17.9 (p625-626)
@test_approx_eq_eps Pi(0.2, degrees2radians(30), sind(30)^2) 0.53896 1e-5

# values from Abramowitz & Stegun, Table 16.1 (p582-583)
sn(0.672,0.36)
u20 = 20*K20/90
θs2020 = 0.35274_9211
θc2020 = 0.96935_0025/sqrt(secd(20))
θd2020 = 1.02789_45992/sqrt(secd(20))
θn2020 = 1.00369_53131
# sn, cn, dn
@test_approx_eq_eps sn(u20,m20) θs2020/θn2020 1e-9
@test_approx_eq_eps cn(u20,m20) θc2020/θn2020 1e-9
@test_approx_eq_eps dn(u20,m20) θd2020/θn2020 1e-9
# cd, sd, nd
@test_approx_eq_eps cd(u20,m20) θc2020/θd2020 1e-9
@test_approx_eq_eps sd(u20,m20) θs2020/θd2020 1e-9
@test_approx_eq_eps nd(u20,m20) θn2020/θd2020 1e-9
# dc, nc, sc
@test_approx_eq_eps dc(u20,m20) θd2020/θc2020 1e-9
@test_approx_eq_eps nc(u20,m20) θn2020/θc2020 1e-9
@test_approx_eq_eps sc(u20,m20) θs2020/θc2020 1e-9
# ns, ds, cs
@test_approx_eq_eps ns(u20,m20) θn2020/θs2020 4e-9
@test_approx_eq_eps ds(u20,m20) θd2020/θs2020 4e-9
@test_approx_eq_eps cs(u20,m20) θc2020/θs2020 4e-9
# ellipj
s,c,d = ellipj(u20,m20)
@test_approx_eq_eps s θs2020/θn2020 1e-9
@test_approx_eq_eps c θc2020/θn2020 1e-9
@test_approx_eq_eps d θd2020/θn2020 1e-9

for u = -1.:0.21:1.
    @test_approx_eq am(u,0) u
    @test_approx_eq sn(u,0) sin(u)
    @test_approx_eq cn(u,0) cos(u)
    @test_approx_eq dn(u,0) 1.
    @test_approx_eq cd(u,0) cos(u)
    @test_approx_eq sd(u,0) sin(u)
    @test_approx_eq nd(u,0) 1.
    @test_approx_eq dc(u,0) sec(u)
    @test_approx_eq nc(u,0) sec(u)
    @test_approx_eq sc(u,0) tan(u)
    @test_approx_eq ns(u,0) csc(u)
    @test_approx_eq ds(u,0) csc(u)
    @test_approx_eq cs(u,0) cot(u)

    @test_approx_eq am(u,1) atan(sinh(u))
    @test_approx_eq sn(u,1) tanh(u)
    @test_approx_eq cn(u,1) sech(u)
    @test_approx_eq dn(u,1) sech(u)
    @test_approx_eq cd(u,1) 1.
    @test_approx_eq sd(u,1) sinh(u)
    @test_approx_eq nd(u,1) cosh(u)
    @test_approx_eq dc(u,1) 1.
    @test_approx_eq nc(u,1) cosh(u)
    @test_approx_eq sc(u,1) sinh(u)
    @test_approx_eq ns(u,1) coth(u)
    @test_approx_eq ds(u,1) csch(u)
    @test_approx_eq cs(u,1) csch(u)
end
