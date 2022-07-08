push!(LOAD_PATH,"../src/")

using Documenter, EllipticFunctions

makedocs(sitename="EllipticFunctions", modules = [EllipticFunctions])