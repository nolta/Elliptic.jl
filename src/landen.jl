#
# Compute (& store) landen sequence, and other fun stuff
#
module Landen

export LandenSeq, landenseq, NonConvergedLandenSeq

const Maybe{T} = Union{Nothing, T}

#
# Landen Sequence
#

#Returns a tuple of kₙs in descending (ascending) Landen sequence where K′ₙ/Kₙ doubles (halves)
"""
    LandenSeq(k, [k′]; N::Int=10, ktol=√eps(k), descending=(k ≤ 1/√2))

Return the elliptic moduli, `{(kₙ, k′ₙ)}`, in a Landen sequence starting with `k₀ = k`

Compute and store `LandenSeq` containing `{(kₙ, k′ₙ)}` of a descending (ascending) Landen
sequence starting from `k₀ = k` until kₙ converges to 0 (1) via criteria `kₙ ≤ ktol`
(`k′ₙ ≤ ktol` or `kₙ ≥ 1 - ktol`).
Returns a `landen <: LandenSeq{0} == NonConvergedLandenSeq` if `(k ∉ ℝ & k ∉ [0, 1])`,
`(k ∉ ℂ & abs(k) ∉ [0, 1])`, or the sequence did not converge in `N` iterations.
"""
struct LandenSeq{N, T<:Union{AbstractFloat, Complex{<:AbstractFloat}}}
    # basically a SArray{(2, N), T, 2*N}
    ks::NTuple{N, T}
    k′s::NTuple{N, T}

    LandenSeq(U::Type) = new{0, U}((), ())
    function LandenSeq(k::Number, k′::Maybe{Number}=nothing; N=10, ktol=_default_tol(k),
            descending=(k ≤ 1/sqrt(2)))
        k = float(k)
        T = typeof(k)
        # abs for complex moduli
        ka = isreal(k) ? k : abs(k)
        (0 ≤ ka ≤ 1) || return LandenSeq(T)

        # cant compute before checking value
        isnothing(k′) && (k′ = _k′(k))

        # check if already within range
        if (descending && (abs(k) ≤ ktol)) || (!descending && (abs(k′) ≤ ktol))
            return new{1, T}((k, ), (k′, ))
        end

        descending || ((k, k′) = (k′, k))

        ks = zeros(T, N+1)
        k′s = zeros(T, N+1)

        ks[1], k′s[1] = k, k′

        # landen_kernel = descending ? _desc_landen_kernel : _asc_landen_kernel
        # k_conv = descending ? ktol : (1 - ktol)

        n = 0
        for i in 2:(N+1)
            k, k′ = _desc_landen_kernel(k, k′)
            # dirty hack
            # k, k′ = (k, k′) ./ hypot(k, k′)
            ks[i], k′s[i] = k, k′

            if abs(k) ≤ ktol # _converged(k, k′, ktol, descending)
                n = i
                break
            end
        end

        descending || ((ks, k′s) = (k′s, ks))
        return new{n, T}((ks[1:n]..., ), (k′s[1:n]..., ))
    end
end

const NonConvergedLandenSeq = LandenSeq{0}
LandenSeq(k; args...) = LandenSeq(Base.Bottom)

Base.length(::LandenSeq{N}) where N = N
Base.firstindex(::LandenSeq) = 1
Base.lastindex(landen::LandenSeq) = length(landen)
Base.getindex(landen::LandenSeq, i::Int) = (landen.ks[i], landen.k′s[i])
Base.getindex(landen::LandenSeq, I) = [landen[i] for i in I]

ktype(::LandenSeq{N, T}) where {N,T} = T
Base.eltype(landen::LandenSeq)= NTuple{2, ktype(landen)}
Base.iterate(landen::LandenSeq, state=1) = state ≤ length(landen) ?
        (landen[state], state+1) : nothing

# landen can be used on complex moduli, and k can be an Int
# √(eps) because 1 - k² = 1 when k² ≤ eps(k), so m is indistinguishable from 0
@inline  _default_tol(k) = √(eps(float(real(one(k)))))

const _kmid = 1/sqrt(2)

 # TODO improve accuraccy for k ≈ 1
 _k′(k) = sqrt(1-k)*sqrt(1+k)

## TODO eval stablity and error accumulation here
# esp descending


# descending landen iteration k → 0
# not accuarte for k ≈ 1
@inline function _desc_landen_kernel(k, k′=_k′(k))
    if k ≤ eps(typeof(k))^(1/5)
        kk = k^2
        # taylor series expansions around k = 0
        return (1/4*kk*(1 + 1/2*kk*(1 + 5/8*kk*(1 + 7/2*kk))),
                (1 - kk^2/32*(1 - kk*(1 - 57/64*kk))))
    else
        return (k^2/(1+k′)^2, 2*sqrt(k′)/(1+k′))
    end
end

# ascending landen iteration k → 1
#  ignored second argument to match `desc_landen_kernel` call
@inline _asc_landen_kernel(k, k′) = reverse(_desc_landen_kernel(k′, k))

"""
    landenseq(k::Number; N::Int=10, ktol=eps(k), descending=true)

Return a tuple of `{kₙ}` in a Landen sequence starting from `k₀ = k`

Return `()::Tuple{}` for invalid inputs or if the sequence did not converge.
See [`LandenSeq`](@ref) for options or more details.
"""
landenseq(args...; kargs...) = LandenSeq(args...; kargs...).ks

end # module
