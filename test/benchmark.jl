using BenchmarkTools, PkgBenchmark, EllipticFunctions

pkgpath = dirname(dirname(pathof(EllipticFunctions)))
# move it out of the repository so that you can check out different branches
script = tempname() * ".jl"
benchpath = joinpath(pkgpath, "benchmark", "benchmarks.jl")
cp(benchpath, script)

# compare to master by default
commit = get(ENV, "BASELINE_COMMIT", "master")

# start
j = judge(pkgpath, commit, f=mean, retune=true, script=script)

println("BASELINE ", commit)
println(j.baseline_results)

println("THIS BRANCH")
println(j.target_results)

println("DIFFGROUP")
println(j.benchmarkgroup)

println("MARKDOWN")
export_markdown(stdout, j, export_invariants=true)
