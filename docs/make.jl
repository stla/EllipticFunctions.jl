push!(LOAD_PATH,"../src/")

using Documenter, EllipticFunctions

makedocs(sitename="EllipticFunctions", modules = [EllipticFunctions])

deploydocs(
    repo = "github.com/stla/EllipticFunctions.jl.git",
)