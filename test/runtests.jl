using Jacobi
using SpecialFunctions
using Test

@testset "Some values of jtheta3" begin
  @test isapprox(
    jtheta3(0, 1im),
    pi^(1/4) / gamma(3/4)
  )
  @test isapprox(
    jtheta3(1+1im, 1im),
    0.86456184935441778 - 0.28488586703507289im
  )
end