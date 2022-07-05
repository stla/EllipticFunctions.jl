push!(LOAD_PATH,"../src/")

using Documenter, Jacobi

makedocs(sitename="My Documentation", modules = [Jacobi])