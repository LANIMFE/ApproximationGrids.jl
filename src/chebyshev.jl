# Structs
struct ChebyshevGrid{T <: AbstractFloat, P} <: ApproximationGrid
    a::T
    b::T
    n::Int
    p::P
end

# Constructors
function ChebyshevGrid(a::A, b::B, n) where {A, B}
    T = promote_type(float(A), float(B))
    p = FFTW.plan_r2r!(fill(zero(T), n), FFTW.REDFT10)
    return ChebyshevGrid{T, typeof(p)}(a, b, n, p)
end

# Nodes
function nodes(g::ChebyshevGrid{T}) where {T}
    n = g.n
    a = T(0.5) * (g.b - g.a)
    b = T(0.5) * (g.b + g.a)
    return [a * sinpi((2k + one(T) - n) / 2n) + b for k = 0:(n - 1)]
end

# Weights
function weights(g::ChebyshevGrid)
    n = g.n
    a = zeros(n)
    a[1:2:n] = (2 / n) ./ (1 .- ((1:2:n) .- 1).^2)
    return FFTW.r2r(a, FFTW.REDFT01)
end

# Jacobian
jacobian(g::ChebyshevGrid{T}) where {T} = T(0.5) * (g.b - g.a)

# Interpolation
function interpolate(g::ChebyshevGrid, v::AbstractVector, x′)
    @assert (g.a ≤ x′ ≤ g.b)
    x = 2 * (x′ - (g.b + g.a) / 2) / (g.b - g.a)

    n = length(v)

    c = reverse_plan(g.p, v, n)
    c[1] /= 2
    c   ./= n

    if n == 0
        return zero(eltype(c))
    elseif n == 1 # avoid issues with NaN x
        return first(c) * one(x)
    end

    x = 2x
    b₁ = b₂ = zero(eltype(c))
    @inbounds for k = n:-1:2
        b₂, b₁ = b₁, muladd(x, b₁, c[k] - b₂)
    end

    return muladd(x / 2, b₁, c[1] - b₂)
end

reverse_plan(p, v::AbstractVector{Float64}, n) = p * reverse(v)
function reverse_plan(p, v::AbstractVector{S}, n) where {S}
    r = reinterpret(eltype(S), reverse(v))
    m = length(r)
    d = div(m, n)
    w = permutedims(reshape(r, d, n))

    for i in 1:d
        p * view(w, :, i) # `p` is an in-place `FFTW` plan which mutates `w`
        r[i:d:m] = view(w, :, i)
    end

    return reinterpret(S, r)
end
