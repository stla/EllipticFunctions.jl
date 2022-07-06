push!(LOAD_PATH,"../src/")

using Documenter, Jacobi

makedocs(sitename="Jacobi", modules = [Jacobi])