using Documenter, ApproximationGrids

makedocs(;
    modules=[ApproximationGrids],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://gitlab.com/lanimfe/ApproximationGrids.jl/blob/{commit}{path}#L{line}",
    sitename="ApproximationGrids.jl",
    authors="Pablo Zubieta",
    assets=String[],
)

deploydocs(;
    repo="gitlab.com/lanimfe/ApproximationGrids.jl",
)
