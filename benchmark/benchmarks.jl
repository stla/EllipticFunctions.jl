using BenchmarkTools
using EllipticFunctions
using SpecialFunctions

import EllipticFunctions: argtheta3, calctheta3, dologtheta3

if isempty(methods(argtheta3,  (Any,Any,Any)))
  # I had another implementation previously, but they turns out to cause regessions
  f1(z,tau) = argtheta3(z,tau)
else
  f1(z,tau) = argtheta3(z,tau,0)
end
if isempty(methods(dologtheta3, (Any,Any,Any)))
  f2(z,tau) = dologtheta3(z,tau)
else
  f2(z,tau) = dologtheta3(z,tau,0)
end
bench_argtheta3(m) = map(splat(f1), m)
bench_calctheta3(m) = map(splat(calctheta3), m)
bench_dologtheta3(m) = map(splat(f2), m)


const SUITE = @benchmarkset "Elliptic functions" begin
  @benchmarkset "Aid functions of jacobi theta" begin
    # make sure current implementation is not poor for some specific z or tau
    data = ((x+im*y,p+im*q) for x in 0:0.1:4 for y in 0:0.1:4 for p in 0:0.1:4 for q in 0.1:0.1:4)
    @case "argtheta3" bench_argtheta3($data)
    @case "calctheta3" bench_calctheta3($data)
    @case "dologtheta3" bench_dologtheta3($data)
  end
  @benchmarkset "Jacobi thetas" begin
    q1, q2, q3, q4 = qfromtau.((1+1im, 1im, 0.1im, 2+0.3im))
    @case "jtheta1" jtheta1(20, $q1)
    @case "jtheta2" jtheta2(1+1im, $q2)
    @case "jtheta3" jtheta3(2+2im, $q3)
    @case "jtheta4" jtheta4(1+1im, $q4)
    @case "jtheta1dash" jtheta1dash(1 + 1im, $q2)
  end

  @benchmarkset "Log-Jacobi thetas" begin
    z   = 0.1 + 0.4im
    tau = 6.2 + 6.3im
    q   = qfromtau(tau)
    @case "ljtheta1" ljtheta1($z, $q)
    @case "ljtheta2" ljtheta1($z, $q)
    @case "ljtheta3" ljtheta1($z, $q)
    @case "ljtheta4" ljtheta1($z, $q)
  end

  @benchmarkset "etaDedekind" begin
    @case "on 1im / 2" etaDedekind(1im / 2)
    @case "on 2im" etaDedekind(2im)
  end

  @benchmarkset "Lambda" begin
    @case "on 1im * sqrt(x)" lambda(1im * sqrt(2.0))
    @case "on 1im * sqrt(1/x)" lambda(1im * sqrt(1 / 2.0))
  end

  @benchmarkset "Klein J" begin
    @case "kleinj(2im)" kleinj(2im)
    @case "kleinjinv(0.1+0.2im)" kleinjinv(0.1+0.2im)
  end

  @benchmarkset "Carlson functions" begin
    x = 1 + im
    y = -2 + 3im
    z = -3
    p = 4im
    @case "RC" CarlsonRC($y*$y, $y*$y-$x*$x)
    @case "RJ" CarlsonRJ($x, $y, $z, $p)
    @case "RF" CarlsonRF(0, $y, 1)
    @case "RG" CarlsonRG(0, $y, 1)
    @case "RD" CarlsonRD(0, $y, 1)
  end

  @benchmarkset "Various elliptic's" begin
    @case "ellipticK on 0.5" ellipticK(0.5)
    @case "ellipticK on big\"0.5\"" ellipticK(big"0.5")
    @case "ellipticE on 2-3im" ellipticE(2-3im)
    @case "ellipticE on 0.5" ellipticE(0.5)
    @case "ellipticZ(-5+3im,-4-9im)" ellipticZ(-5+3im,-4-9im)
    @case "ellipticPI(1+im, 1, 2-im)" ellipticPI(1+im, 1, 2-im)
    @case "ellipticF(big\"0.15\", big\"0.81\")" ellipticF(big"0.15", big"0.81")
  end

  @benchmarkset "Eisenstein series" begin
    @case "E2 on 0.005+0.005im" EisensteinE2(0.005+0.005im)
  end

  @benchmarkset "Others" begin
    omega1 = gamma(1/3)^3 / 4 / pi
    z0     = omega1 * (1 + 1im/sqrt(3))
    omega2 = 1im * omega1
    omega = (omega1, omega2)

    @case "wp" wp($z0; g=(0, 1))
    @case "agm" agm(1, sqrt(2))
    @case "wsigma" wsigma($omega1, omega=$omega)
    @case "wzeta 1" wzeta($omega1, omega=$omega)
    @case "wzeta 2" wzeta(1 + 1im; g=(5 + 3im, 5 + 3im))
    @case "jellip" jellip("cn", 2+2im, tau=1im)

    x = gamma(1/4)^2 / 4 / sqrt(pi)
    @case "elliptic Invariants" ellipticInvariants($x, im * $x)
  end

  @benchmarkset "Neville theta functions" begin
    z = 2.5
    m = 0.3
    bigz = big(z)
    @case "thetaC" thetaC($z, m=$m)
    @case "thetaD" thetaD($z, m=$m)
    @case "thetaN" thetaN($z, m=$m)
    @case "thetaS" thetaS($z, m=$m)
    @case "thetaC on big" thetaC($bigz, m=$bigz)
  end
end
