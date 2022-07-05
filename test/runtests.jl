using Jacobi
using SpecialFunctions
using Test

@testset "Some values of the jtheta functions" begin
  @test isapprox(
    jtheta3(0, 1im),
    pi^(1/4) / gamma(3/4)
  )
  @test isapprox(
    jtheta1(1+1im, 1im),
    1.1816128551455719 + 0.59589712760417439im
  )
  @test isapprox(
    jtheta2(1+1im, 1im),
    0.74328632006610539 - 0.904159309718008im
  )
  @test isapprox(
    jtheta3(1+1im, 1im),
    0.86456184935441778 - 0.28488586703507289im
  )
  @test isapprox(
    jtheta4(1+1im, 1im),
    1.1351891564632007 + 0.28517396444192509im
  )
end