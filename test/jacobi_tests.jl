#
# Tests for Jacobi sub-module
#

@testset "Jacobi" begin
    @testset "Abramowitz & Stegun, Table 16.1 (p582-583)" begin
        # from Abramowitz & Stegun, Table 17.2 (p610-611)
        K20 = 1.62002_58991_24204

        u20 = 20*K20/90
        θs2020 = 0.35274_9211
        θc2020 = 0.96935_0025/sqrt(secd(20))
        θd2020 = 1.02789_45992/sqrt(secd(20))
        θn2020 = 1.00369_53131
        # sn, cn, dn
        @test Jacobi.sn(u20,m20) ≈ θs2020/θn2020 atol=1e-9
        @test Jacobi.cn(u20,m20) ≈ θc2020/θn2020 atol=1e-9
        @test Jacobi.dn(u20,m20) ≈ θd2020/θn2020 atol=1e-9
        # Jacobi.cd, sd, nd
        @test Jacobi.cd(u20,m20) ≈ θc2020/θd2020 atol=1e-9
        @test Jacobi.sd(u20,m20) ≈ θs2020/θd2020 atol=1e-9
        @test Jacobi.nd(u20,m20) ≈ θn2020/θd2020 atol=1e-9
        # Jacobi.dc, nc, sc
        @test Jacobi.dc(u20,m20) ≈ θd2020/θc2020 atol=1e-9
        @test Jacobi.nc(u20,m20) ≈ θn2020/θc2020 atol=1e-9
        @test Jacobi.sc(u20,m20) ≈ θs2020/θc2020 atol=1e-9
        # Jacobi.ns, ds, cs
        @test Jacobi.ns(u20,m20) ≈ θn2020/θs2020 atol=4e-9
        @test Jacobi.ds(u20,m20) ≈ θd2020/θs2020 atol=4e-9
        @test Jacobi.cs(u20,m20) ≈ θc2020/θs2020 atol=4e-9
        # ellipj
        s,c,d = ellipj(u20,m20)
        @test s ≈ θs2020/θn2020 atol=1e-9
        @test c ≈ θc2020/θn2020 atol=1e-9
        @test d ≈ θd2020/θn2020 atol=1e-9
    end

    @testset "u = $u" for u in -1.:0.21:1.
        @test Jacobi.am(u,0) ≈ u
        @test Jacobi.sn(u,0) ≈ sin(u)
        @test Jacobi.cn(u,0) ≈ cos(u)
        @test Jacobi.dn(u,0) ≈ 1.
        @test Jacobi.cd(u,0) ≈ cos(u)
        @test Jacobi.sd(u,0) ≈ sin(u)
        @test Jacobi.nd(u,0) ≈ 1.
        @test Jacobi.dc(u,0) ≈ sec(u)
        @test Jacobi.nc(u,0) ≈ sec(u)
        @test Jacobi.sc(u,0) ≈ tan(u)
        @test Jacobi.ns(u,0) ≈ csc(u)
        @test Jacobi.ds(u,0) ≈ csc(u)
        @test Jacobi.cs(u,0) ≈ cot(u)

        @test Jacobi.am(u,1) ≈ atan(sinh(u))
        @test Jacobi.sn(u,1) ≈ tanh(u)
        @test Jacobi.cn(u,1) ≈ sech(u)
        @test Jacobi.dn(u,1) ≈ sech(u)
        @test Jacobi.cd(u,1) ≈ 1.
        @test Jacobi.sd(u,1) ≈ sinh(u)
        @test Jacobi.nd(u,1) ≈ cosh(u)
        @test Jacobi.dc(u,1) ≈ 1.
        @test Jacobi.nc(u,1) ≈ cosh(u)
        @test Jacobi.sc(u,1) ≈ sinh(u)
        @test Jacobi.ns(u,1) ≈ coth(u)
        @test Jacobi.ds(u,1) ≈ csch(u)
        @test Jacobi.cs(u,1) ≈ csch(u)
    end

    @testset "errors" begin
        u = 0.5 # random value for testing
        @test_throws DomainError Jacobi.am(u, -0.1)
        @test_throws DomainError Jacobi.am(u, -eps(0.0))

        @test_throws DomainError Jacobi.am(u, 1 + eps(1.0))
        @test_throws DomainError Jacobi.am(u, 1.1)
    end
end
