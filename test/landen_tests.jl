#
# Tests for Landen Sequences
#

@testset "Landen" begin
    using Elliptic.Landen

    areseqequal(x, y; args...) = all(isapprox.(x, y; args...))

    @testset "basics" begin
        @test LandenSeq(1, descending=true) isa NonConvergedLandenSeq
        @test LandenSeq(0, descending=false) isa NonConvergedLandenSeq
        @test LandenSeq(-eps(0.0)) isa NonConvergedLandenSeq
        @test LandenSeq(1 + eps()) isa NonConvergedLandenSeq
        @test LandenSeq(NaN) isa NonConvergedLandenSeq
        @test LandenSeq(NaN, descending=false) isa NonConvergedLandenSeq
        @test LandenSeq((0.5, ), descending=false) isa LandenSeq{0, Base.Bottom}
        @test LandenSeq("0.5", descending=false) isa LandenSeq{0, Base.Bottom}

        l0 = LandenSeq(0)
        l1 = LandenSeq(1)


        @test l0.ks == l1.k′s == (0.0, )
        @test l1.ks == l0.k′s == (1.0, )
    end

    @testset "$(f)" for f in ["table1", "table2"]
        # tables go to very high tol, so reduce eps here to match
        ks = first( readdlm(joinpath("data/papers/ow", "$f.csv"), ',', header=true) )[:, 2]
        k′s = Landen._k′.(ks)
        k = ks[1]
        k′ = Landen._k′(k)
        lc = LandenSeq(k, descending=true, ktol=eps())
        lc′ = LandenSeq(k′, descending=false, ktol=eps())

        # test _k′
        @test areseqequal(hypot.(lc.ks, lc.k′s), 1, rtol=1e-9, atol=eps())
        @test areseqequal(hypot.(lc′.ks, lc′.k′s), 1, rtol=1e-9, atol=eps())
        # descending seq is correct
        @test areseqequal(lc.ks, ks, rtol=1e-9, atol=eps())
        @test areseqequal(lc.k′s, k′s, rtol=1e-9, atol=eps())
        # ascending seq works
        @test areseqequal(lc′.k′s, ks, rtol=1e-9, atol=eps())
        @test areseqequal(lc′.ks, k′s, rtol=1e-9, atol=eps())
        # technically these two could be off by 2*rtol ...
        @test areseqequal(lc′.k′s, lc.ks, rtol=1e-9, atol=eps())
        @test areseqequal(lc′.ks, lc.k′s, rtol=1e-9, atol=eps())

        @test areseqequal(landenseq(k, descending=true, ktol=eps()), lc.ks, rtol=1e-15)
        @test areseqequal(landenseq(k′, descending=false, ktol=eps()), lc.k′s, rtol=1e-15)

        N = length(ks)-1
        @test LandenSeq(k, N=N, descending=true, ktol=eps()) == lc
        @test LandenSeq(k, N=N-1, descending=true, ktol=eps()) isa NonConvergedLandenSeq
    end
end
