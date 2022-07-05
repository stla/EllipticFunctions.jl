using Jacobi
using SpecialFunctions
using Test

@testset "Some values of the jtheta functions." begin
  @test isapprox(jtheta3(0, 1im), pi^(1 / 4) / gamma(3 / 4))
  @test isapprox(
    jtheta1(1 + 1im, 1im),
    1.1816128551455719 + 0.59589712760417439im
  )
  @test isapprox(
    jtheta2(1 + 1im, 1im),
    0.74328632006610539 - 0.904159309718008im
  )
  @test isapprox(
    jtheta3(1 + 1im, 1im),
    0.86456184935441778 - 0.28488586703507289im
  )
  @test isapprox(
    jtheta4(1 + 1im, 1im),
    1.1351891564632007 + 0.28517396444192509im
  )
end

@testset "A value of jtheta1dash." begin
  @test isapprox(
    jtheta1dash(1 + 1im, 1im),
    0.81117649363854416 - 0.89452803853474627im
  )
end

@testset "Some values of etaDedekind." begin
  @test isapprox(etaDedekind(1im / 2), gamma(1 / 4) / 2^(7 / 8) / pi^(3 / 4))
  @test isapprox(etaDedekind(2im), gamma(1 / 4) / 2^(11 / 8) / pi^(3 / 4))
end

@testset "Lambda modular relation." begin
  x = 2.0
  @test isapprox(lambda(1im * sqrt(x)) + lambda(1im * sqrt(1 / x)), 1)
end

@testset "A value of kleinj." begin
  @test isapprox(
    kleinj(2im),
    66^3
  )
end

@testset "A value of ellipticK." begin
  @test isapprox(
    ellipticK(0.5),
    8 * pi^(3/2) / gamma(-1/4)^2
  )
end
