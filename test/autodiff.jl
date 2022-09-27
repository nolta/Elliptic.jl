using Zygote
using ForwardDiff
using SpecialFunctions

@testset "Zygote and ForwardDiff" begin

    num_trials = 1
    ks = rand(num_trials)
    ϕs = rand(num_trials) .* 2π
    ns = rand(num_trials)

    @show ks, ϕs, ns

    # Test several known derivative identities across a wide range of values of ϕ and m
    # in order to verify that derivatives work correctly
    for (k, ϕ, n) in zip(ks, ϕs, ns)
        m = k^2
        dk_dm = 0.5 / k

        # I. Tests for complete integrals, from https://en.wikipedia.org/wiki/Elliptic_integral
        # 1. K'(m) == K'(k) * dk/dm == E(k) / (k * (1 - k^2)) - K(k)/k
        @test Zygote.gradient(K, m)[1] ≈ (E(m) / (k * (1 - k^2)) - K(m) / k) * dk_dm
        @test ForwardDiff.derivative(K, m) ≈ (E(m) / (k * (1 - k^2)) - K(m) / k) * dk_dm

        # 2. E'(m) == E'(k) * dk/dm == (E(k) - K(k))/k
        @test Zygote.gradient(E, m)[1] ≈ (E(m) - K(m))/k * dk_dm
        @test ForwardDiff.derivative(E, m) ≈ (E(m) - K(m))/k * dk_dm

        # II. Tests for incomplete integrals, from https://functions.wolfram.com/EllipticIntegrals/EllipticF/introductions/IncompleteEllipticIntegrals/ShowAll.html
        # 3. ∂ϕ(F(ϕ, m)) == 1 / √(1 - m*sin(ϕ)^2)
        @test Zygote.gradient(ϕ -> F(ϕ, m), ϕ)[1] ≈ 1 / √(1 - m*sin(ϕ)^2)
        @test ForwardDiff.derivative(ϕ -> F(ϕ, m), ϕ) ≈ 1 / √(1 - m*sin(ϕ)^2)

        # 4. ∂m(F(ϕ, m)) == E(ϕ, m) / (2 * m * (1 - m)) - F(ϕ, m) / 2m - sin(2ϕ) / (4 * (1-m) * √(1 - m * sin(ϕ)^2))
        @test Zygote.gradient(m -> F(ϕ, m), m)[1] ≈
            E(ϕ, m) / (2 * m * (1 - m)) -
            F(ϕ, m) / 2 / m -
            sin(2*ϕ) / (4 * (1 - m) * √(1 - m * sin(ϕ)^2))
        @test ForwardDiff.derivative(m -> F(ϕ, m), m) ≈
            E(ϕ, m) / (2 * m * (1 - m)) -
            F(ϕ, m) / 2 / m -
            sin(2*ϕ) / (4 * (1 - m) * √(1 - m * sin(ϕ)^2))

        # 5. ∂ϕ(E(ϕ, m)) == √(1 - m * sin(ϕ)^2)
        @test Zygote.gradient(ϕ -> E(ϕ, m), ϕ)[1] ≈ √(1 - m * sin(ϕ)^2)
        @test ForwardDiff.derivative(ϕ -> E(ϕ, m), ϕ) ≈ √(1 - m * sin(ϕ)^2)

        # 6. ∂m(E(ϕ, m)) == (E(ϕ, m) - F(ϕ, m)) / 2m
        @test Zygote.gradient(m -> E(ϕ, m), m)[1] ≈ (E(ϕ, m) - F(ϕ, m)) / 2m
        @test ForwardDiff.derivative(m -> E(ϕ, m), m) ≈ (E(ϕ, m) - F(ϕ, m)) / 2m

    end
end
