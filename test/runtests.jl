using EllipticFunctions
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
  @test real(jtheta1(1-1im, 1.e-13*im)) == -Inf
  @test imag(jtheta1(1-1im, 1.e-13*im)) == Inf
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

@testset "CarlsonRJ homogeneity." begin
  x = 1 + im
  y = -2 + 3im
  z = -3
  p = 4im
  kappa = 2
  @test isapprox(
    CarlsonRJ(x, y, z, p) / kappa / sqrt(kappa),
    CarlsonRJ(kappa*x, kappa*y, kappa*z, kappa*p)
  )
end

@testset "CarlsonRJ(x, y, y, p)." begin
  x = 1 + im
  y = -2 + 3im
  p = 4im
  @test isapprox(
    CarlsonRJ(x, y, y, p),
    3 * (CarlsonRC(x, y) - CarlsonRC(x, p)) / (p - y)
  )
end

@testset "A value of ellipticK." begin
  @test isapprox(
    ellipticK(0.5),
    8 * pi^(3/2) / gamma(-1/4)^2
  )
end

@testset "ellipticK and CarlsonRF." begin
  m = 2 - 3im
  @test isapprox(
    ellipticK(m),
    CarlsonRF(0, 1-m, 1)
  )
end

@testset "Complete ellipticE and CarlsonRG." begin
  m = 2 - 3im
  @test isapprox(
    ellipticE(m),
    2*CarlsonRG(0, 1-m, 1)
  )
end

@testset "Complete ellipticE and CarlsonRD." begin
  m = 2 - 3im
  @test isapprox(
    ellipticE(m),
    (1-m) * (CarlsonRD(0, 1-m, 1) + CarlsonRD(0, 1, 1-m)) / 3
  )
end

@testset "A value of complete ellipticE." begin
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

@testset "Differential equation wp." begin
  g2 = 1.4 - 1im
  g3 = 1.6 + 0.5im
  g = (g2, g3)
  z = 1 + 1im
  p = wp(z; g = g)
  pdash = wp(z; g = g, derivative = 1)
  @test isapprox(
    pdash^2,
    4 * p^3 - g2 * p - g3
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

@testset "Neville theta functions." begin
  z = 2.5
  m = 0.3
  @test isapprox(
    thetaC(z; m = m),
    -0.65900466676738154967
  )
  @test isapprox(
    thetaD(z; m = m),
    0.95182196661267561994
  )
  @test isapprox(
    thetaN(z; m = m),
    1.0526693354651613637
  )
  @test isapprox(
    thetaS(z; m = m),
    0.82086879524530400536
  )
end

@testset "jellip functions." begin
  u = 2 + 2im
  tau = 1im
  @test isapprox(
    jellip("cn", u; tau = tau)^2 + jellip("sn", u; tau = tau)^2,
    1
  )
end

@testset "am function." begin
  phi = 1 + 1im
  m = 2 - 1im
  u = ellipticF(phi, m)
  v = 2 - 2im
  psi = am(v, m)
  @test isapprox(
    am(u, m),
    phi
  )
  @test isapprox(
    ellipticF(psi, m),
    v
  )
end

