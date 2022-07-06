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

@testset "A value of ellipticE." begin
  @test isapprox(
    ellipticE(0.5),
    (2 * gamma(3/4)^4 + pi^2) / (4 * sqrt(pi) * gamma(3/4)^2)
  )
end

@testset "A value of agm." begin
  @test isapprox(
    agm(1, sqrt(2)),
    2 * pi^(3/2) * sqrt(2) / gamma(1/4)^2
  )
end

@testset "kleinjinv works." begin
  j = 0.1 + 0.2im
  @test isapprox(
    kleinj(kleinjinv(j)),
    j
  )
end

@testset "EisensteinE2 development." begin
  q = 0.005 + 0.005im
  @test isapprox(
    EisensteinE2(q),
    1 - 24 * (q/(1-q) + 2*q^2/(1-q^2) + 3*q^3/(1-q^3) + 4*q^4/(1-q^4) + 5*q^5/(1-q^5))
  )
end

@testset "Sum of the e_i is zero." begin
  omega1 = 1.4 - 1im
  omega2 = 1.6 + 0.5im
  omega = (omega1, omega2)
  e1 = wp(omega1; omega = omega)
  e2 = wp(omega2; omega = omega)
  e3 = wp(-omega1-omega2; omega = omega)
  @test isapprox(
    e1 + e2 + e3, 0; atol = 1e-10
  )
end

@testset "A value of wsigma." begin
  omega1 = gamma(1/4)^2 / 4 / sqrt(pi)
  omega2 = 1im * omega1
  omega = (omega1, omega2)
  @test isapprox(
    wsigma(omega1; omega = omega),
    exp(pi/8) * 2^(1/4)
  )
end

@testset "A value of wzeta given omega." begin
  omega1 = gamma(1/4)^2 / 4 / sqrt(pi)
  omega2 = 1im * omega1
  omega = (omega1, omega2)
  @test isapprox(
    wzeta(omega1; omega = omega),
    pi / 4 / omega1
  )
end

@testset "A value of wzeta given g." begin
  g = (5 + 3im, 5 + 3im)
  @test isapprox(
    wzeta(1 + 1im; g = g),
    0.802084165492408 - 0.381791358666872im
  )
end
