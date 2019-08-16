module ApproximationGrids

### Imports
using FFTW


### Exports
export ChebyshevGrid, RegularGrid,
       interpolate, jacobian, nodes, weights


### Implemetation
abstract type ApproximationGrid end

# Objects' string representations
function Base.show(io::IO, g::T) where {T <: ApproximationGrid}
    print(io, "$(nameof(T))(($(g.a),$(g.b)), $(g.n))")
end

include("regular.jl")
include("chebyshev.jl")


end # module
