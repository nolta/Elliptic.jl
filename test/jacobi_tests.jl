#
# Tests for Jacobi sub-module
#

@testset "Jacobi" begin
    using Jacobi

    @testset "Abramowitz & Stegun, Table 16.1 (p582-583)" begin
        # from Abramowitz & Stegun, Table 17.2 (p610-611)
        K20 = 1.62002_58991_24204

        u20 = 20*K20/90
        θs2020 = 0.35274_9211
        θc2020 = 0.96935_0025/sqrt(secd(20))
        θd2020 = 1.02789_45992/sqrt(secd(20))
        θn2020 = 1.00369_53131
        # sn, cn, dn, nn
        @test sn(u20,m20) ≈ θs2020/θn2020 atol=1e-9
        @test cn(u20,m20) ≈ θc2020/θn2020 atol=1e-9
        @test dn(u20,m20) ≈ θd2020/θn2020 atol=1e-9
        @test nn(u20, m20) ≈ 1.0

        # cd, sd, nd, dd
        @test cd(u20,m20) ≈ θc2020/θd2020 atol=1e-9
        @test sd(u20,m20) ≈ θs2020/θd2020 atol=1e-9
        @test nd(u20,m20) ≈ θn2020/θd2020 atol=1e-9
        @test dd(u20, m20) ≈ 1.0

        # dc, nc, sc, cc
        @test dc(u20,m20) ≈ θd2020/θc2020 atol=1e-9
        @test nc(u20,m20) ≈ θn2020/θc2020 atol=1e-9
        @test sc(u20,m20) ≈ θs2020/θc2020 atol=1e-9
        @test cc(u20, m20) ≈ 1.0

        # ns, ds, cs, ss
        @test ns(u20,m20) ≈ θn2020/θs2020 atol=4e-9
        @test ds(u20,m20) ≈ θd2020/θs2020 atol=4e-9
        @test cs(u20,m20) ≈ θc2020/θs2020 atol=4e-9
        @test ss(u20, m20) ≈ 1.0

        # ellipj
        s,c,d = ellipj(u20,m20)
        @test s ≈ θs2020/θn2020 atol=1e-9
        @test c ≈ θc2020/θn2020 atol=1e-9
        @test d ≈ θd2020/θn2020 atol=1e-9
    end

    @testset "u = $u" for u in -1.:0.21:1.
        @test am(u,0) ≈ u
        @test sn(u,0) ≈ sin(u)
        @test cn(u,0) ≈ cos(u)
        @test dn(u,0) ≈ 1.
        @test nn(u,0) ≈ 1.
        @test cd(u,0) ≈ cos(u)
        @test sd(u,0) ≈ sin(u)
        @test nd(u,0) ≈ 1.
        @test dd(u,0) ≈ 1.
        @test dc(u,0) ≈ sec(u)
        @test nc(u,0) ≈ sec(u)
        @test sc(u,0) ≈ tan(u)
        @test cc(u,0) ≈ 1.
        @test ns(u,0) ≈ csc(u)
        @test ds(u,0) ≈ csc(u)
        @test cs(u,0) ≈ cot(u)
        @test ss(u,0) ≈ 1.

        @test am(u,1) ≈ atan(sinh(u))
        @test sn(u,1) ≈ tanh(u)
        @test cn(u,1) ≈ sech(u)
        @test dn(u,1) ≈ sech(u)
        @test nn(u,1) ≈ 1.
        @test cd(u,1) ≈ 1.
        @test sd(u,1) ≈ sinh(u)
        @test nd(u,1) ≈ cosh(u)
        @test dd(u,1) ≈ 1.
        @test dc(u,1) ≈ 1.
        @test nc(u,1) ≈ cosh(u)
        @test sc(u,1) ≈ sinh(u)
        @test cc(u,1) ≈ 1.
        @test ns(u,1) ≈ coth(u)
        @test ds(u,1) ≈ csch(u)
        @test cs(u,1) ≈ csch(u)
        @test ss(u,1) ≈ 1.
    end

    @testset "errors" begin
        u = 0.5 # random value for testing
        @test_throws DomainError am(u, -0.1)
        @test_throws DomainError am(u, -eps(0.0))

        @test_throws DomainError am(u, 1 + eps(1.0))
        @test_throws DomainError am(u, 1.1)
    end
end
