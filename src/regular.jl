struct RegularGrid{T <: AbstractFloat} <: ApproximationGrid
    a::T
    b::T
    n::Int
end

function RegularGrid(a::A, b::B, n) where {A, B}
    T = promote_type(float(A), float(B))
    return RegularGrid{T}(a, b, n)
end

function nodes(g::RegularGrid{T}) where {T}
    n = g.n
    h = (g.b - g.a) / (n - 1)
    return g.a .+ h * (0:(n - 1))
end

function weights(g::RegularGrid)
    @assert g.n ≥ 4
    m = g.n % 3

    w = fill(oftype(g.a, 9 // 8), g.n)
    w[1] = 3 // 8
    for i = 4:3:(g.n - 3)
        w[i] = 3 // 4
    end
    if m == 1
        w[end    ] =  3 // 8
    elseif m == 0
        w[end - 2] = 17 // 24
        w[end - 1] =  4 // 3
        w[end    ] =  1 // 3
    elseif m == 2
        w[end - 3] =  7 // 6
        w[end - 2] = 11 // 12
        w[end - 1] =  7 // 6
        w[end    ] =  3 // 8
    end

    return w
end

jacobian(g::RegularGrid) = (g.b - g.a) / (g.n - 1)

function interpolate(g::RegularGrid, v::AbstractVector, x′)
    @assert (g.a ≤ x′ ≤ g.b)
    r = x′ .- range(g.a, stop = g.b, length = g.n)

    i = findfirst(x -> x ≤ 0, r)
    if r[i] == 0
        return v[i] / one(r[i]) # in case there are type instabilities
    end
    k = 3 * div(i - 1, 3)
    k = (k + 4 > g.n) ? (g.n - 4) : k

    # Barycentric interpolation weights
    w₁ =  1 / r[k + 1]
    w₂ = -3 / r[k + 2]
    w₃ =  3 / r[k + 3]
    w₄ = -1 / r[k + 4]
    # Interpolation function reference values
    v₁, v₂, v₃, v₄ = v[k + 1], v[k + 2], v[k + 3], v[k + 4]

    return (w₁ * v₁ + w₂ * v₂ + w₃ * v₃ + w₄ * v₄) / (w₁ + w₂ + w₃ + w₄)
end
