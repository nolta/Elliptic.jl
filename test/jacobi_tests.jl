#
# Tests for Jacobi sub-module
#

@testset "Jacobi" begin
    @testset "Abramowitz & Stegun Table 16.1" begin
        dataloc = "data/ab"
        # table 16.1
        t161, _ = readdlm(joinpath(dataloc, "table_16_1.csv"), ',', header=true)
        # table 17.2
        t172, _ = readdlm(joinpath(dataloc, "table_17_2.csv"), ',', header=true)

        K_lut = Dict(zip(t172[:, 1], t172[:, 2]))

        # data vary first by ϵ ∈ 0:5:90 then α ∈ 0:5:85
        αs = 0:5:85
        ϵs = 0:5:90
        θss = reshape(t161[:, 3], length(ϵs), length(αs))
        θns = reshape(t161[:, 4], length(ϵs), length(αs))

        @testset "α = $α" for (i, α) in enumerate(αs)
            K = K_lut[α]
            m = sind(α)^2
            denom = sqrt(secd(α))

            @testset "ϵ = $ϵ" for (j, ϵ) in enumerate(ϵs)
                j₁ = length(ϵs) - j + 1
                ϵ₁ = ϵs[j₁]
                u = ϵ * K / 90

                θs = θss[j, i]
                θn = θns[j, i]
                θc = θss[j₁, i]/denom
                θd = θns[j₁, i]/denom

                @test Jacobi.sn(u, m) ≈ θs / θn atol=2.5e-9
                @test Jacobi.cn(u, m) ≈ θc / θn atol=2.5e-9
                @test Jacobi.dn(u, m) ≈ θd / θn atol=2.5e-9
                @test Jacobi.nn(u, m) == 1.0

                @test Jacobi.sd(u, m) ≈ θs / θd atol=2.5e-9
                @test Jacobi.cd(u, m) ≈ θc / θd atol=2.5e-9
                @test Jacobi.dd(u, m) == 1.0
                @test Jacobi.nd(u, m) ≈ θn / θd atol=2.5e-9

                @test Jacobi.cc(u, m) == 1.0
                if ϵ != 90
                    # very sensitive around u = K,
                    # estimate of K(0) = pi/2 + 4e-9, so cosine causes errors
                    # also, errors build up in ϵ ≥75°, so lower tolerences for that
                    @test Jacobi.sc(u, m) ≈ θs / θc atol=1e-7
                    @test Jacobi.dc(u, m) ≈ θd / θc atol=1e-7
                    @test Jacobi.nc(u, m) ≈ θn / θc atol=1e-7
                end

                @test Jacobi.ss(u, m) == 1.0
                @test Jacobi.cs(u, m) ≈ θc / θs atol=5.5e-8
                @test Jacobi.ds(u, m) ≈ θd / θs atol=5.5e-8
                @test Jacobi.ns(u, m) ≈ θn / θs atol=5.5e-8

                # ellipj
                s, c, d = ellipj(u, m)
                @test s ≈ θs / θn atol=1e-9
                @test c ≈ θc / θn atol=1e-9
                @test d ≈ θd / θn atol=1e-9
            end
        end
    end

    @testset "u = $u" for u in -1.:0.21:1.0
        @test Jacobi.am(u,0) ≈ u
        @test Jacobi.sn(u,0) ≈ sin(u)
        @test Jacobi.cn(u,0) ≈ cos(u)
        @test Jacobi.dn(u,0) ≈ 1.
        @test Jacobi.nn(u,0) ≈ 1.
        @test Jacobi.cd(u,0) ≈ cos(u)
        @test Jacobi.sd(u,0) ≈ sin(u)
        @test Jacobi.nd(u,0) ≈ 1.
        @test Jacobi.dd(u,0) ≈ 1.
        @test Jacobi.dc(u,0) ≈ sec(u)
        @test Jacobi.nc(u,0) ≈ sec(u)
        @test Jacobi.sc(u,0) ≈ tan(u)
        @test Jacobi.cc(u,0) ≈ 1.
        @test Jacobi.ns(u,0) ≈ csc(u)
        @test Jacobi.ds(u,0) ≈ csc(u)
        @test Jacobi.cs(u,0) ≈ cot(u)
        @test Jacobi.ss(u,0) ≈ 1.

        @test Jacobi.am(u,1) ≈ atan(sinh(u))
        @test Jacobi.sn(u,1) ≈ tanh(u)
        @test Jacobi.cn(u,1) ≈ sech(u)
        @test Jacobi.dn(u,1) ≈ sech(u)
        @test Jacobi.nn(u,1) ≈ 1.
        @test Jacobi.cd(u,1) ≈ 1.
        @test Jacobi.sd(u,1) ≈ sinh(u)
        @test Jacobi.nd(u,1) ≈ cosh(u)
        @test Jacobi.dd(u,1) ≈ 1.
        @test Jacobi.dc(u,1) ≈ 1.
        @test Jacobi.nc(u,1) ≈ cosh(u)
        @test Jacobi.sc(u,1) ≈ sinh(u)
        @test Jacobi.cc(u,1) ≈ 1.
        @test Jacobi.ns(u,1) ≈ coth(u)
        @test Jacobi.ds(u,1) ≈ csch(u)
        @test Jacobi.cs(u,1) ≈ csch(u)
        @test Jacobi.ss(u,1) ≈ 1.
    end

    @testset "errors" begin
        u = 0.5 # random value for testing
        @test_throws DomainError Jacobi.am(u, -0.1)
        @test_throws DomainError Jacobi.am(u, -eps(0.0))

        @test_throws DomainError Jacobi.am(u, 1 + eps(1.0))
        @test_throws DomainError Jacobi.am(u, 1.1)
    end
end
