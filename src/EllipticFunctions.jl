module EllipticFunctions

import SpecialFunctions

export qfromtau
export taufromq
export ljtheta1
export jtheta1
export ljtheta2
export jtheta2
export ljtheta3
export jtheta3
export ljtheta4
export jtheta4
export jtheta_ab
export jtheta1dash
export etaDedekind
export lambda
export kleinj
export kleinjinv
export CarlsonRF
export CarlsonRC
export CarlsonRD
export CarlsonRG
export CarlsonRJ
export ellipticE
export ellipticF
export ellipticK
export ellipticZ
export ellipticPI
export agm
export EisensteinE2
export EisensteinE4
export EisensteinE6
export halfPeriods
export ellipticInvariants
export wp
export wsigma
export wzeta
export thetaC
export thetaD
export thetaN
export thetaS
export jellip
export am

function xcispi(x)
    return exp(1im * (pi * x))
end

function csqrt(x::Number)
  return sqrt(Complex(x))
end

function areclose(z1::S, z2::S) where {T <: Real, S <: Union{T, Complex{T}}}
  z1 == z2 && (return true)
  eps2 = eps(T)^2
  mod2_z2 = abs2(z2)
  maxmod2 = (mod2_z2 < eps2) ? 1.0 : max(abs2(z1), mod2_z2)
  return abs2(z1 - z2) < 4.0 * eps2 * maxmod2
end

function modulo(a::Real, p::Real)
  i = a > 0 ? floor(a / p) : ceil(a / p)
  return a - i * p
end

function calctheta3(z::Number, tau::Complex)
  out = Complex(one(z) * one(tau))
  n = 0
  while true
    n += 1
    qweight =
      exp(n * 1im * (pi * (n * tau + 2 * z))) +
      exp(n * 1im * (pi * (n * tau - 2 * z)))
    out += qweight
    modulus = abs(out)
    if isnan(modulus)
      error("NaN has occured in the summation.")
    elseif isinf(modulus)
      error("Infinity has occured in the summation.")
    elseif n >= 3 && areclose(out + qweight, out)
      break
    end
  end
  return log(out)
end

function argtheta3(z::Number, tau::Complex, passes::Int64)
  passes += 1
  if passes > 1000
    error("Reached 1000 iterations (argtheta3).")
  end
  z_img = imag(z)
  tau_img = imag(tau)
  h = tau_img / 2.0
  zuse = complex(modulo(real(z), 1), z_img)
  if z_img < -h
    out = argtheta3(-zuse, tau, passes)
  elseif z_img >= h
    quotient = floor(z_img / tau_img + 0.5)
    zmin = zuse - quotient * tau
    out = -2im * quotient * zmin * pi +
          argtheta3(zmin, tau, passes) - 
          1im * tau * quotient * quotient * pi
  else
    out = calctheta3(zuse, tau)
  end
  return out
end

function dologtheta4(z::Number, tau::Complex, passes::Int64) 
  return dologtheta3(z + 0.5, tau, passes + 1)
end

function dologtheta3(z::Number, tau::Complex, passes::Int64)
  passes += 1
  tau_rl = real(tau)
  if tau_rl > 0.6
    tau2_rl = modulo(tau_rl + 1.0, 2) - 1.0
  else 
    tau2_rl = modulo(tau_rl - 1.0, 2) + 1.0
  end
  tau2_img = imag(tau)
  tau2 = complex(tau2_rl, tau2_img)
  if abs(tau2) < 0.98 && tau2_img < 0.98
    tauprime = -1.0 / tau2
    out =
      1im * tauprime * z * z * pi +
      dologtheta3(z * tauprime, tauprime, passes) - 
      log(csqrt(tau2) / csqrt(one(z) * one(tau) * 1im))
  elseif tau2_rl >= 0.6
    out = dologtheta4(z, tau2 - 1.0, passes)
  elseif tau2_rl <= -0.6
    out = dologtheta4(z, tau2 + 1.0, passes)
  else
    out = argtheta3(z, tau2, 0)
  end
  return out
end

function M(z::Number, tau::Complex)
  return 1im * (z + tau / 4) * pi
end

function _ljtheta2_raw(z::Number, tau::Complex)
  return M(z, tau) + dologtheta3(z + 0.5 * tau, tau, 0)
end

function _jtheta2_raw(z::Number, tau::Complex)
  return exp(_ljtheta2_raw(z, tau))
end

function _ljtheta1_raw(z::Number, tau::Complex)
  return _ljtheta2_raw(z - 0.5, tau)
end

function _jtheta1_raw(z::Number, tau::Complex)
  return exp(_ljtheta1_raw(z, tau))
end

function _ljtheta3_raw(z::Number, tau::Complex)
  return dologtheta3(z, tau, 0)
end

function _jtheta3_raw(z::Number, tau::Complex) 
  return exp(_ljtheta3_raw(z, tau))
end

function _ljtheta4_raw(z::Number, tau::Complex)
  return dologtheta4(z, tau, 0)
end

function _jtheta4_raw(z::Number, tau::Complex)
  return exp(_ljtheta4_raw(z, tau))
end

function _jtheta1(z::Number, tau::Complex) 
  return _jtheta1_raw(z/pi, tau);
end

function _jtheta2(z::Number, tau::Complex) 
  return _jtheta2_raw(z/pi, tau);
end

function _jtheta3(z::Number, tau::Complex) 
  return _jtheta3_raw(z/pi, tau);
end

function _jtheta4(z::Number, tau::Complex) 
  return _jtheta4_raw(z/pi, tau);
end

function principal_log_branch(z)
  y = imag(z)
  (y ≤ π) && (-y < π) && return z
  twopi = 2 * one(y) * π
  y = mod(y, twopi)
  y > π && (y -= twopi)
  return complex(real(z), y)
end

function _ljtheta1(z::Number, tau::Complex)
  return principal_log_branch(_ljtheta1_raw(z/pi, tau))
end

function _ljtheta2(z::Number, tau::Complex)
  return principal_log_branch(_ljtheta2_raw(z/pi, tau))
end

function _ljtheta3(z::Number, tau::Complex)
  return principal_log_branch(_ljtheta3_raw(z/pi, tau))
end

function _ljtheta4(z::Number, tau::Complex)
  return principal_log_branch(_ljtheta4_raw(z/pi, tau))
end

function _jtheta_ab(a::Number, b::Number, z::Number, tau::Complex) 
  alpha = a * tau
  beta  = b + z/pi
  C = xcispi(a * (alpha + 2*beta)) 
  return C * _jtheta3_raw(alpha + beta, tau)
end

function _jtheta1dash(z::Number, tau::Complex)
  q = xcispi(tau)
  out = complex(zero(q))
  alt = -one(q)
  q² = q * q
  q²ⁿ = one(q)
  qⁿ⁽ⁿ⁺¹⁾ = one(q)
  for n = 0:3000
    if n > 0
      q²ⁿ *= q²
      qⁿ⁽ⁿ⁺¹⁾ *= q²ⁿ
    end
    alt = -alt
    k = 2 * n + one(q)
    outnew = out + alt * qⁿ⁽ⁿ⁺¹⁾ * k * cos(k * z)
    if areclose(out, outnew)
      return 2 * sqrt(csqrt(q)) * out
    end
    out = outnew
  end
  error("Reached 3000 iterations.")
end

function _etaDedekind(tau::Complex)
  return xcispi(-1 / tau / 12.0) *
    _jtheta3_raw((-1 / tau + 1.0) / 2.0, -3.0 / tau) / sqrt(-1im * tau)
end

function isvector(x)
  return length(size(x)) == 1
end

function _EisensteinE2(tau::Complex)
  j3 = _jtheta3_raw(0, tau)
  lbd = (_jtheta2_raw(0, tau) / j3)^4
  j3sq = j3^2
  return 6.0/pi * ellipticE(lbd) * j3sq - j3sq^2 - _jtheta4_raw(0, tau)^4
end

function _jtheta1dash0(tau::Complex)
  return -2im * _jtheta_ab(1/(6*one(tau)), 1/2, 0, 3*tau)^3
  #return exp(_ljtheta2(0.0, tau) + _ljtheta3(0.0, tau) + _ljtheta4(0.0, tau))
end

function _jtheta1dashdashdash0(tau::Complex)
  return -_jtheta1dash(0, tau) * _EisensteinE2(tau)
end

function _dljtheta1(z::Number, tau::Complex)
  if z == 0
    return _jtheta1dash0(tau) / _jtheta1_raw(0.0, tau)
  end
  return _jtheta1dash(z, tau) / _jtheta1(z, tau)
end

function _E4(tau::Complex)
  return (
    _jtheta2_raw(0, tau)^8 + _jtheta3_raw(0, tau)^8 + _jtheta4_raw(0, tau)^8
  ) / 2
end

function _E6(tau::Complex)
  j2 = _jtheta2_raw(0, tau)
  j3 = _jtheta3_raw(0, tau)
  j4 = _jtheta4_raw(0, tau)
  x3 = j3^4
  x4 = j4^4
  return (x3^3 + x4^3 - 3.0 * j2^8 * (x3 + x4)) / 2.0
end

function _omega1_and_tau(g)
  g2, g3 = g
  if g2 == 0
    omega1 = SpecialFunctions.gamma(1/3)^3 / 4 / pi / g3^(1/6)
    tau    = 0.5 + 1im * sqrt(3)/2
  else
    g2cube = g2*g2*g2
    j      = 1728 * g2cube / (g2cube - 27*g3*g3)
    if isinf(j)
      return (-1im*pi/2/sqrt(3), complex(Inf, Inf))
    end
    tau = kleinjinv(j)
    if g3 == 0
      omega1 = 1im * pi * sqrt(sqrt(1.0 / g2 / 12 * _E4(tau)))
    else
      G6_over_G4 = 2.0 * pi * pi / 21.0 * _E6(tau) / _E4(tau)
      omega1     = csqrt(7.0 * G6_over_G4 * g2 / (12.0 * g3)) 
    end
    #omega1 = 1im * pi * sqrt(sqrt(1.0 / g2 / 12 * _E4(tau)))
  end
  return (omega1, tau)
end

function _g2_from_omega1_and_tau(omega1::Number, tau::Complex)
  j2 = _jtheta2_raw(0, tau)
  j3 = _jtheta3_raw(0, tau)
  return 4/3 * (pi/2/omega1)^4 * (j2^8 - (j2*j3)^4 + j3^8)
end

function _wpFromTau(z, tau::Complex)
  j2 = _jtheta2_raw(0, tau)
  j3 = _jtheta3_raw(0, tau)
  j1 = _jtheta1_raw.(z, tau)
  j4 = _jtheta4_raw.(z, tau)
  return (pi * j2 * j3 * j4 ./ j1)^2 .- (pi^2 * (j2^4 + j3^4) / 3.0)
end

function _wpDerivative(z, omega1::Number, tau::Complex)
  w1 = 2 * omega1 / pi
  z1 = - z / 2 / omega1
  j1 = _jtheta1_raw.(z1, tau)
  j2 = _jtheta2_raw.(z1, tau)
  j3 = _jtheta3_raw.(z1, tau)
  j4 = _jtheta4_raw.(z1, tau)
  f = _jtheta1dash0(tau)^3 /
    (_jtheta2_raw(0, tau) * _jtheta3_raw(0, tau) *
       _jtheta4_raw(0, tau) * j1.^3)
  2/(w1*w1*w1) * j2 .* j3 .* j4 .* f
end

function _thetaS(z, tau::Complex)
  j3sq = _jtheta3_raw(zero(tau), tau)^2
  zprime = z / j3sq / pi
  return j3sq * _jtheta1_raw.(zprime, tau) / _jtheta1dash0(tau)
end

function _thetaC(z, tau::Complex)
  zprime = z / _jtheta3_raw(zero(tau), tau)^2 / pi
  return _jtheta2_raw.(zprime, tau) / _jtheta2_raw(zero(tau), tau)
end

function _thetaN(z, tau::Complex)
  zprime = z / _jtheta3_raw(zero(tau), tau)^2 / pi
  return _jtheta4_raw.(zprime, tau) / _jtheta4_raw(zero(tau), tau)
end

function _thetaD(z, tau::Complex)
  j3 = _jtheta3_raw(zero(tau), tau)
  zprime = z / j3^2 / pi
  return _jtheta3_raw.(zprime, tau) / j3
end

function _tau_from_m(m::Number)
  1im * ellipticK(one(m) - m) / ellipticK(m)
end

function _check_and_get_tau_from_m(tau::Union{Missing,Number}, m::Union{Missing,Number})
  nmissing = ismissing(tau) + ismissing(m)
  @assert nmissing == 1 ArgumentError("You must supply either `tau` or `m`.")
  if !ismissing(tau)
    @assert imag(tau) > 0 ArgumentError("The imaginary part of `tau` must be nonnegative.")
  else
    tau = _tau_from_m(m)
    @assert imag(tau) > 0 ArgumentError("Invalid value of `m`.")
  end
  return tau
end


# exports ####

"""
    qfromtau(tau)

The nome `q` given the `tau` parameter.

# Arguments
- `tau`: complex number with nonnegative imaginary part
""" 
function qfromtau(tau::Complex)
  @assert imag(tau) > 0 ArgumentError("Invalid `tau`.")
  return xcispi(tau)
end

"""
    taufromq(q)

The `tau` parameter given the nome `q`.

# Arguments
- `q`: complex number, the nome; it must not be a negative real number and its modulus must be strictly smaller than 1
""" 
function taufromq(q::Number)
  @assert abs(q) < 1 ArgumentError("Invalid `q`.")
  @assert imag(q) != 0 || real(q) > 0 ArgumentError("Invalid `q`.")
  return -im * (log(q) / pi)
end

"""
    ljtheta1(z, q)

Logarithm of the first Jacobi theta function.

# Arguments
- `z`: real or complex number or array of numbers
- `q`: the nome
"""
function ljtheta1(z, q::Number)
  @assert abs(q) < 1 ArgumentError("Invalid `q`.")
  @assert imag(q) != 0 || real(q) > 0 ArgumentError("Invalid `q`.")
  return _ljtheta1.(z, taufromq(q))
end

"""
    jtheta1(z, q)

First Jacobi theta function.

# Arguments
- `z`: real or complex number or array of numbers
- `q`: the nome
"""
function jtheta1(z, q::Number)
  @assert abs(q) < 1 ArgumentError("Invalid `q`.")
  @assert imag(q) != 0 || real(q) > 0 ArgumentError("Invalid `q`.")
  z,tau = promote(z, taufromq(q))
  return _jtheta1.(z, tau)
end

"""
    ljtheta2(z, q)

Logarithm of the second Jacobi theta function.

# Arguments
- `z`: real or complex number or array of numbers
- `q`: the nome
"""
function ljtheta2(z, q::Number)
  @assert abs(q) < 1 ArgumentError("Invalid `q`.")
  @assert imag(q) != 0 || real(q) > 0 ArgumentError("Invalid `q`.")
  return _ljtheta2.(z, taufromq(q))
end

"""
    jtheta2(z, q)

Second Jacobi theta function.

# Arguments
- `z`: real or complex number or array of numbers
- `q`: the nome
"""
function jtheta2(z, q::Number)
  @assert abs(q) < 1 ArgumentError("Invalid `q`.")
  @assert imag(q) != 0 || real(q) > 0 ArgumentError("Invalid `q`.")
  z,tau = promote(z, taufromq(q))
  return _jtheta2.(z, tau)
end

"""
    ljtheta3(z, q)

Logarithm of the third Jacobi theta function.

# Arguments
- `z`: real or complex number or array of complex numbers
- `q`: the nome
"""
function ljtheta3(z, q::Number)
  @assert abs(q) < 1 ArgumentError("Invalid `q`.")
  @assert imag(q) != 0 || real(q) > 0 ArgumentError("Invalid `q`.")
  return _ljtheta3.(z, taufromq(q))
end

"""
    jtheta3(z, q)

Third Jacobi theta function.

# Arguments
- `z`: real or complex number or array of numbers
- `q`: the nome
"""
function jtheta3(z, q::Number)
  @assert abs(q) < 1 ArgumentError("Invalid `q`.")
  @assert imag(q) != 0 || real(q) > 0 ArgumentError("Invalid `q`.")
  z,tau = promote(z, taufromq(q))
  return _jtheta3.(z, tau)
end

"""
    ljtheta4(z, q)

Logarithm of the fourth Jacobi theta function.

# Arguments
- `z`: real or complex number or array of complex numbers
- `q`: the nome
"""
function ljtheta4(z, q::Number)
  @assert abs(q) < 1 ArgumentError("Invalid `q`.")
  @assert imag(q) != 0 || real(q) > 0 ArgumentError("Invalid `q`.")
  return _ljtheta4.(z, taufromq(q))
end

"""
    jtheta4(z, q)

Fourth Jacobi theta function.

# Arguments
- `z`: real or complex number or array of numbers
- `q`: the nome
"""
function jtheta4(z, q::Number)
  @assert abs(q) < 1 ArgumentError("Invalid `q`.")
  @assert imag(q) != 0 || real(q) > 0 ArgumentError("Invalid `q`.")
  z,tau = promote(z, taufromq(q))
  return _jtheta4.(z, tau)
end

"""
    jtheta_ab(a, b, z, q)

Jacobi theta function with characteristics. This is a family of functions
parameterized by `a` and `b`, which contains the opposite of the first Jacobi 
theta function (`a=b=0.5`), the second Jacobi theta function (`a=0.5,b=0`), 
the third Jacobi theta function (`a=b=0`), and the fourth Jacobi theta 
function (`a=0,b=0.5`).

# Arguments
- `a`: first characteristic, a real or complex number
- `b`: second characteristic, a real or complex number
- `z`: real or complex number or array of numbers
- `q`: the nome
"""
function jtheta_ab(a::Number, b::Number, z, q::Number)
  @assert abs(q) < 1 ArgumentError("Invalid `q`.")
  @assert imag(q) != 0 || real(q) > 0 ArgumentError("Invalid `q`.")
  a,b,z,tau = promote(a, b, z, taufromq(q))
  return _jtheta_ab.(a, b, z, tau)
end

"""
    jtheta1dash(z, q)

Derivative of the first Jacobi theta function.

# Arguments
- `z`: real or complex number or array of numbers
- `q`: the nome
"""
function jtheta1dash(z, q::Number)
  @assert abs(q) < 1 ArgumentError("Invalid `q`.")
  @assert imag(q) != 0 || real(q) > 0 ArgumentError("Invalid `q`.")
  z,tau = promote(z, taufromq(q))
  return _jtheta1dash.(z, tau)
end

"""
    etaDedekind(tau)

Dedekind eta function.

# Arguments
- `tau`: complex number with nonnegative imaginary part
"""
function etaDedekind(tau::Complex)
  @assert imag(tau) > 0 ArgumentError("Invalid `tau`.")
  return _etaDedekind(tau)
end

"""
    lambda(tau)

Lambda modular function.

# Arguments
- `tau`: complex number with nonnegative imaginary part
"""
function lambda(tau::Complex)
  @assert imag(tau) > 0 ArgumentError("Invalid `tau`.")
  return (_jtheta2_raw(0, tau) / _jtheta3_raw(0, tau))^4
end

"""
    kleinj(tau)

Klein j-invariant function.

# Arguments
- `tau`: complex number with nonnegative imaginary part
"""
function kleinj(tau::Complex)
  @assert imag(tau) > 0 ArgumentError("Invalid `tau`.")
  lbd = (_jtheta2_raw(0, tau) / _jtheta3_raw(0, tau))^4
  x = lbd * (1.0 - lbd)
  return 256 * (1/x - 1)^2 * (1 - x)  #256 * (1-x)^3 / x^2
end

"""
    CarlsonRF(x, y, z)

Carlson 'RF' integral.

# Arguments
- `x`,`y`,`z`: real or complex numbers; at most one of them can be zero
"""
function CarlsonRF(x::Number, y::Number, z::Number)
  local A
  (xzero, yzero, zzero) = iszero.((x, y, z))
  (xx, yy, zz, _) = promote(x, y, z, 1.0)
  T = real(typeof(xx))
  xzero + yzero + zzero ≤ 1 || throw(ArgumentError("At most one of `x`, `y`, `z` can be 0."))
  dx = dy = dz = typemax(T)
  epsilon = eps(T)^(4/3)
  while dx > epsilon || dy > epsilon || dz > epsilon
    lambda = csqrt(xx)*csqrt(yy) + csqrt(yy)*csqrt(zz) + csqrt(zz)*csqrt(xx)
    xx = (xx + lambda) / 4.0
    yy = (yy + lambda) / 4.0
    zz = (zz + lambda) / 4.0
    A = (xx + yy + zz) / 3.0
    dx = abs2(1.0 - xx/A)
    dy = abs2(1.0 - yy/A)
    dz = abs2(1.0 - zz/A)
  end
  E2 = sqrt(dx*dy) + sqrt(dy*dz) + sqrt(dz*dx)
  E3 = sqrt(dy*dx*dz)
  return (1 - E2/10 + E3/14 + E2*E2/24 - 3*E2*E3/44 - 5*E2*E2*E2/208 +
    3*E3*E3/104 + E2*E2*E3/16) / sqrt(A)
end

"""
    CarlsonRC(x, y)

Carlson 'RC' integral.

# Arguments
- `x`,`y`: real or complex numbers; `y` cannot be zero
"""
function CarlsonRC(x::Number, y::Number)
  @assert y != 0 ArgumentError("`y` cannot be 0.")
  return CarlsonRF(x, y, y)
end

"""
    CarlsonRD(x, y, z)

Carlson 'RD' integral.

# Arguments
- `x`,`y`,`z`: real or complex numbers; at most one of them can be zero
"""
function CarlsonRD(x::Number, y::Number, z::Number)
  local A
  (xzero, yzero, zzero) = iszero.((x, y, z))
  (xx, yy, zz, _) = promote(x, y, z, 1.0)
  T = real(typeof(xx))
  xzero + yzero + zzero ≤ 1 || throw(ArgumentError("At most one of `x`, `y`, `z` can be 0."))
  dx = typemax(T)
  dy = typemax(T)
  dz = typemax(T)
  epsilon = eps(T)^(4/3)
  s = complex(zero(T), zero(T))
  fac = complex(one(T), zero(T))
  while dx > epsilon || dy > epsilon || dz > epsilon
    lambda = csqrt(xx)*csqrt(yy) + csqrt(yy)*csqrt(zz) + csqrt(zz)*csqrt(xx)
    s = s + fac/(csqrt(zz) * (zz + lambda))
    fac = fac / 4.0
    xx = (xx + lambda) / 4.0
    yy = (yy + lambda) / 4.0
    zz = (zz + lambda) / 4.0
    A = (xx + yy + 3*zz) / 5.0
    dx = abs2(1.0 - xx/A)
    dy = abs2(1.0 - yy/A)
    dz = abs2(1.0 - zz/A)
  end
  dx = sqrt(dx)
  dy = sqrt(dy)
  dz = sqrt(dz)
  E2 = dx * dy + dy * dz + 3 * dz * dz + 2 * dz * dx +
              dx * dz + 2 * dy * dz
  E3 = dz * dz * dz + dx * dz * dz + 3 * dx * dy * dz +
              2 * dy * dz * dz + dy * dz * dz + 2 * dx * dz * dz
  E4 = dy * dz * dz * dz + dx * dz * dz * dz + dx * dy * dz * dz +
              2 * dx * dy * dz * dz
  E5 = dx * dy * dz * dz * dz
  return 3 * s + fac * (1 - 3 * E2/14 + E3/6 + 9 * E2 * E2/88 - 3 * E4/22 -
          9 * E2 * E3/52 + 3 * E5/26 - E2 * E2 * E2/16 +
          3 * E3 * E3/40 + 3 * E2 * E4/20 + 45 * E2 * E2 * E3/272 -
          9 * (E3 * E4 + E2 * E5)/68) / A / sqrt(A)
end

"""
    CarlsonRG(x, y, z)

Carlson 'RG' integral.

# Arguments
- `x`,`y`,`z`: real or complex numbers
"""
function CarlsonRG(x::Number, y::Number, z::Number)
  local A
  (xzero, yzero, zzero) = iszero.((x, y, z))
  nzeros = xzero + yzero + zzero
  (xx, yy, zz, _) = promote(x, y, z, 1.0)
  T = real(typeof(xx))
  if nzeros == 3
    return complex(zero(T), zero(T))
  end
  if nzeros == 2
    return csqrt(x + y + z) / 2
  end
  if zzero
    return CarlsonRG(y, z, x)
  end
  return (z * CarlsonRF(x, y, z) - 
    (x - z) * (y - z) * CarlsonRD(x, y, z) / 3 + 
    csqrt(x) * csqrt(y) / csqrt(z)) / 2
end

"""
    CarlsonRJ(x, y, z, p)

Carlson 'RJ' integral.

# Arguments
- `x`,`y`,`z`,`p`: real or complex numbers; at most one of them can be zero
"""
function CarlsonRJ(x::Number, y::Number, z::Number, p::Number)
  (xzero, yzero, zzero, pzero) = iszero.((x, y, z, p))
  nzeros = xzero + yzero + zzero + pzero
  nzeros ≤ 1 || throw(ArgumentError("At most one of `x`, `y`, `z`, `p` can be 0."))
  (xx, yy, zz, pp, _) = promote(x, y, z, p, 1.0)
  T = real(typeof(xx))
  A0 = (xx + yy + zz + pp + pp) / 5
  A = A0
  delta = (pp - xx) * (pp - yy) * (pp - zz)
  f = 1
  fac = 1
  d = Vector{Complex}(undef, 0)
  e = Vector{Complex}(undef, 0)
  epsilon = 1 * eps(T)^3
  Q = (4 / epsilon)^(1/3) * max(abs2(A-xx), abs2(A-yy), abs2(A-zz), abs2(A-pp))
  xx = Complex(xx)
  yy = Complex(yy)
  zz = Complex(zz)
  pp = Complex(pp)
  while abs2(A) <= Q
    sqrt_x = sqrt(xx)
    sqrt_y = sqrt(yy)
    sqrt_z = sqrt(zz)
    sqrt_p = sqrt(pp)
    dnew = (sqrt_p + sqrt_x) * (sqrt_p + sqrt_y) * (sqrt_p + sqrt_z)
    d = vcat(d, dnew * f)
    e = vcat(e, fac * delta / dnew / dnew)
    f = f * 4
    fac = fac / 64
    lambda = sqrt_x*sqrt_y + sqrt_y*sqrt_z + sqrt_z*sqrt_x
    xx = (xx + lambda) / 4
    yy = (yy + lambda) / 4
    zz = (zz + lambda) / 4
    pp = (pp + lambda) / 4
    A = (A + lambda) / 4
    Q = Q / 16
  end
  M_1_fA = 1 / f / A
  X = (A0-xx) * M_1_fA
  Y = (A0-yy) * M_1_fA
  Z = (A0-zz) * M_1_fA
  P = -(X+Y+Z) / 2
  E2 = X*Y + X*Z + Y*Z - 3*P*P
  E3 = X*Y*Z + 2*E2*P + 4*P*P*P
  E4 = P*(2*X*Y*Z + E2*P + 3*P*P*P)
  E5 = X*Y*Z*P*P
  g = (1 - 3*E2/14 + E3/6 + 9*E2*E2/88 - 3*E4/22 - 9*E2*E3/52 + 3*E5/26) /
    f / A / sqrt(A)
  return length(e) > 1 ? 
    (6 * sum(ifelse.(e .== 0, Complex(one(T)), atan.(sqrt.(e)) ./ sqrt.(e) ) ./ d)) : 
    Complex(zero(T))
end

"""
    ellipticF(phi, m)

Incomplete elliptic integral of the first kind.

# Arguments
- `phi`: real or complex number, the amplitude
- `m`: real or complex number, the squared modulus
"""
function ellipticF(phi::Number, m::Number)
  local k
  if phi == 0 || m == Inf || m == -Inf
    return complex(0.0, 0.0)
  end
  rphi = real(phi)
  if rphi == 0 && imag(phi) == Inf && imag(m) == 0 && real(m) > 0 && real(m) < 1
    return sign(imag(phi)) *
        (ellipticF(pi/2, m) - ellipticF(pi/2, 1/m) / csqrt(m))
  end
  rphiopi = rphi / pi
  if abs(rphiopi) == 1/2 && m == 1
    return complex(NaN, NaN)
  end
  if rphiopi >= -1/2 && rphiopi <= 1/2
    if m == 1 && abs(rphiopi) < 1/2
      return atanh(sin(phi))
    end
    if m == 0
      return phi
    end
    sine = sin(phi)
    if isinf(sine)
      error("`sin(phi)` is not finite.");
    end
    sine2 = sine * sine
    cosine2 = 1.0 - sine2
    oneminusmsine2 = 1.0 - m * sine2
    return sine * CarlsonRF(cosine2, oneminusmsine2, complex(1.0, 0.0))
  end
  if rphiopi > 1/2
    k = ceil(rphiopi-1/2)
    phi = phi - k * pi
  else
    k = -floor(1/2-rphiopi)
    phi = phi - k * pi
  end
  return 2*k*ellipticF(pi/2, m) + ellipticF(phi, m)
end

"""
    ellipticK(m)

Complete elliptic integral of the first kind.

# Arguments
- `m`: real or complex number, the squared modulus
"""
function ellipticK(m::Number)
  return ellipticF(pi * one(m) / 2, m)
end

"""
    ellipticE(phi, m)

Incomplete elliptic integral of the second kind.

# Arguments
- `phi`: real or complex number, the amplitude
- `m`: real or complex number, the squared modulus
"""
function ellipticE(phi::Number, m::Number)
  local k
  if phi == 0
    return zero(phi)
  end
  if real(m) == Inf && imag(m) == 0
    return complex(NaN, NaN)
  end
  rphiopi = real(phi) / pi
  if rphiopi >= -1/2 && rphiopi <= 1/2
    if m == 0
      return phi
    end
    if m == 1
      return sin(phi)
    end
    sine = sin(phi)
    if isinf(sine)
      error("`sin(phi)` is not finite.")
    end
    sine2 = sine*sine
    cosine2 = 1.0 - sine2
    oneminusmsine2 = 1.0 - m * sine2;
    return sine * (CarlsonRF(cosine2, oneminusmsine2, 1.0) -
            m * sine2 * CarlsonRD(cosine2, oneminusmsine2, 1.0) / 3)
  end
  if rphiopi > 1/2
    k = ceil(rphiopi-1/2)
    phi = phi - k * pi
  else
    k = -floor(1/2-rphiopi)
    phi = phi - k * pi
  end
  return 2 * k * ellipticE(pi/2, m) + ellipticE(phi, m)
end

"""
    ellipticE(m)

Complete elliptic integral of the second kind.

# Arguments
- `m`: real or complex number, the squared modulus
"""
function ellipticE(m::Number)
  return ellipticE(pi*one(m)/2, m)
end

"""
    ellipticZ(phi, m)

Jacobi Zeta function.

# Arguments
- `phi`: real or complex number, the amplitude
- `m`: real or complex number, the squared modulus
"""
function ellipticZ(phi::Number, m::Number)
  if isinf(real(m)) && imag(m) == 0
    return NaN
  end
  if m == 1
    rl = real(phi)
    if abs(rl) <= pi/2
      return sin(phi)
    end
    if rl > pi/2
      k = ceil(rl/pi - 0.5)
      return sin(phi - k*pi)
    end
    k = -floor(0.5 - rl/pi)
    return sin(phi - k*pi)
  end
  return ellipticE(phi, m) - ellipticE(m)/ellipticK(m) * ellipticF(phi, m)
end

"""
    ellipticPI(phi, n, m)

Incomplete elliptic integral of first kind.

# Arguments
- `phi`: real or complex number, the amplitude
- `n`: real or complex number, the characteristic
- `m`: real or complex number, the squared modulus
"""
function ellipticPI(phi::Number, n::Number, m::Number)
  if phi == 0 || (isinf(real(m)) && imag(m) == 0) ||
      (isinf(real(n)) && imag(n) == 0)
    return 0.0
  end
  pio2 = pi * one(phi) / 2
  if phi == pio2 && m == 1 && imag(n) == 0 && n != 1
    return real(n) > 1 ? -Inf : Inf
  end
  if phi == pio2 && n == 1
    return NaN
  end
  if phi == pio2 && m == 0 
    return pio2 / csqrt(1-n)
  end
  if phi == pio2 && n == m
    return ellipticE(m) / (1-m)
  end
  if phi == pio2 && n == 0
    return ellipticK(m)
  end
  rphiopi = real(phi) / pi
  if rphiopi >= -1/2 && rphiopi <= 1/2
    sine = sin(phi)
    if isinf(sine)
      error("`sin(phi)` is not finite.")
    end
    sine2 = sine*sine
    cosine2 = 1 - sine2
    oneminusmsine2 = 1 - m*sine2
    return sine * (CarlsonRF(cosine2, oneminusmsine2, 1) +
             n * sine2 * CarlsonRJ(cosine2, oneminusmsine2, 1, 1-n*sine2) / 3)
  end
  if rphiopi > 1/2
    k = ceil(rphiopi - 0.5)
    phi = phi - k*pi
    return 2*k*ellipticPI(pi/2, n, m) + ellipticPI(phi, n, m)
  end
  k = -floor(0.5 - rphiopi)
  phi = phi - k*pi
  return 2*k*ellipticPI(pi/2, n, m) + ellipticPI(phi, n, m)
end

"""
    agm(x, y)

Arithmetic-geometric mean.

# Arguments
- `x`,`y`: real or complex numbers
"""
function agm(x::Number, y::Number)
  if x + y == 0 || x == 0 || y == 0
    return zero(promote_type(typeof(x), typeof(y)))
  end
  return pi * (x + y) / 4 / ellipticK(((x-y)/(x+y))^2)
end

"""
    kleinjinv(j)

Inverse of the Klein j-invariant function.

# Arguments
- `j`: real or complex number
"""
function kleinjinv(j::Number)
  local x
  if isinf(j)
    x = complex(zero(j))
  else
    j2 = j * j
    j3 = j2 * j
    t = (-j3 + 2304 * j2 + 12288 *
          csqrt(3 * (1728 * j2 - j3)) - 884736 * j)^(1/3)
    x = 1/768 * t - (1536 * j - j2) / (768 * t) + (1 - j/768)
  end
  lbd = -(-1 - sqrt(1 - 4*x)) / 2
  return 1im * agm(1, sqrt(1-lbd)) / agm(1, sqrt(lbd))
end

"""
    EisensteinE2(q)

Eisenstein E-series of weight 2.

# Arguments
- `q`: nome, real or complex number; it must not be a negative real number and its modulus must be strictly smaller than 1
"""
function EisensteinE2(q::Number)
  abs(q) < 1 || throw(ArgumentError("Invalid `q`."))
  imag(q) != 0 || real(q) > 0 || throw(ArgumentError("Invalid `q`."))
  tau = -1im * log(q) / pi / 2.0
  return _EisensteinE2(tau)
end

"""
    EisensteinE4(q)

Eisenstein E-series of weight 4.

# Arguments
- `q`: nome, real or complex number; it must not be a negative real number and its modulus must be strictly smaller than 1
"""
function EisensteinE4(q::Number)
  abs(q) < 1 || throw(ArgumentError("Invalid `q`."))
  imag(q) != 0 || real(q) > 0 || throw(ArgumentError("Invalid `q`."))
  tau = -1im * log(q) / pi / 2
  return _E4(tau)
end

"""
    EisensteinE6(q)

Eisenstein E-series of weight 6.

# Arguments
- `q`: nome, real or complex number; it must not be a negative real number and its modulus must be strictly smaller than 1
"""
function EisensteinE6(q::Number)
  abs(q) < 1 || throw(ArgumentError("Invalid `q`."))
  imag(q) != 0 || real(q) > 0 || throw(ArgumentError("Invalid `q`."))
  tau = -1im * log(q) / pi / 2.0
  return _E6(tau)
end

"""
    halfPeriods(g2, g3)

Half-periods ``\\omega_1`` and ``\\omega_2`` from the elliptic invariants.

# Arguments
- `g2`,`g3`: the Weierstrass elliptic invariants, real or complex numbers
"""
function halfPeriods(g2::Number, g3::Number)
  if g2 == 0
    omega1 = SpecialFunctions.gamma(1/3)^3 / 4 / pi / g3^(1/6)
    tau    = 0.5 + 1im * sqrt(3.0)/2.0
  else
    g2cube = g2*g2*g2
    j = 1728 * g2cube / (g2cube - 27*g3*g3)
    if isinf(j)
      return (-1im*pi/2/sqrt(3), complex(Inf, Inf))
    end
    tau = kleinjinv(j)
    if g3 == 0
      omega1 = 1im * pi * sqrt(csqrt(1.0 / g2 / 12 * _E4(tau)))
    else 
      G6_over_G4 = 2.0 * pi * pi / 21.0 * _E6(tau) / _E4(tau)
      omega1     = csqrt(7.0 * G6_over_G4 * g2 / (12.0 * g3))
    end
    #omega1 = 1im * pi * sqrt(csqrt(1.0 / g2 / 12 * _E4(tau)))
  end
  return (omega1, tau*omega1)
end

"""
    ellipticInvariants(omega1, omega2)

Weierstrass elliptic invariants ``g_2`` and ``g_3`` from the half-periods.

# Arguments
- `omega1`,`omega2`: the Weierstrass half periods, real or complex numbers
"""
function ellipticInvariants(omega1::Number, omega2::Number)
  tau = omega2 / omega1
  imag(tau) > 0 || throw(ArgumentError("Invalid pair `(omega1, omega2)`."))
  j2 = _jtheta2_raw(0, tau)
  j3 = _jtheta3_raw(0, tau)
  g2 = 4/3 * (pi/2/omega1)^4 * (j2^8 - (j2*j3)^4 + j3^8)
  g3 = 8/27 * (pi/2/omega1)^6 * (j2^12 - (
      (3/2 * j2^8 * j3^4) + (3/2 * j2^4 * j3^8)
    ) + j3^12)
  return (g2, g3)
end

"""
    wp(z; tau, omega, g, derivative=0)

Weierstrass p-function. One and only one of the parameters `tau`, `omega` or `g` must be given.

# Arguments
- `z`: real or complex number or array of numbers
- `tau`: half-periods ratio, complex number with non negative imaginary part
- `omega`: half-periods, a pair (tuple) of complex numbers
- `g`: elliptic invariants, a pair (tuple) of complex numbers
- `derivative`: order of differentiation, 0, 1, 2 or 3
"""
function wp(z; tau::Union{Missing,Number}=missing, omega::Union{Missing,Tuple{Number,Number}}=missing, g::Union{Missing,Tuple{Number,Number}}=missing, derivative::Int64=0)
  local omega1, weier, weierPrime
  0 ≤ derivative ≤ 3 || throw(ArgumentError("`derivative` must be between 0 and 3."))
  nmissing = ismissing(tau) + ismissing(omega) + ismissing(g)
  nmissing == 2 || throw(ArgumentError("You must supply either `tau`, `omega` or `g`."))
  if !ismissing(tau)
    imag(tau) > 0 || throw(ArgumentError("Invalid `tau`."))
    if derivative != 1
      weier = _wpFromTau(z, tau)
      if derivative == 0
        return weier
      end
      if derivative == 2
        g2 = _g2_from_omega1_and_tau(omega1, tau)
        return 6 * weier * weier .- g2/2
      end
    end
    omega1 = 0.5
  end
  if !ismissing(omega)
    omega1, omega2 = omega
    tau = omega2/omega1
    imag(tau) > 0 || throw(ArgumentError("Invalid `omega`."))
    if derivative != 1
      weier = _wpFromTau(z/omega1/2, tau) / omega1 / omega1 / 4
      if derivative == 0
        return weier
      end
      if derivative == 2
        g2 = _g2_from_omega1_and_tau(omega1, tau)
        return 6 * weier * weier .- g2/2
      end
    end
  end
  if !ismissing(g)
    omega1, tau = _omega1_and_tau(g)
    if derivative != 1
      weier = _wpFromTau(z/omega1/2, tau) / omega1 / omega1 / 4
      if derivative == 0
        return weier
      end
      if derivative == 2
        g2, g3 = g
        return 6 * weier * weier .- g2/2
      end
    end
  end
  weierPrime = _wpDerivative(z, omega1, tau)
  if derivative == 1
    return weierPrime
  end
  return 12 * weier * weierPrime # derivative = 3
end

"""
    wsigma(z; tau, omega, g)

Weierstrass sigma-function. One and only one of the parameters `tau`, `omega` or `g` must be given.

# Arguments
- `z`: real or complex number or array of numbers
- `tau`: half-periods ratio, complex number with non negative imaginary part
- `omega`: half-periods, a pair (tuple) of complex numbers
- `g`: elliptic invariants, a pair (tuple) of complex numbers
"""
function wsigma(z; tau::Union{Missing,Complex}=missing, omega::Union{Missing,Tuple{Number,Number}}=missing, g::Union{Missing,Tuple{Number,Number}}=missing)
  local omega1
  nmissing = ismissing(tau) + ismissing(omega) + ismissing(g)
  nmissing == 2 || throw(ArgumentError("You must supply either `tau`, `omega` or `g`."))
  if !ismissing(tau)
    imag(tau) > 0 || throw(ArgumentError("Invalid `tau`."))
    omega1 = 0.5
  elseif !ismissing(omega)
    omega1 = omega[1]
    tau = omega[2]/omega1
    imag(tau) > 0 || throw(ArgumentError("Invalid `omega`."))
  elseif !ismissing(g)
    omega1, tau = _omega1_and_tau(g)
  end
  w1 = -2 * omega1 / pi
  j1 = _jtheta1.(z/w1, tau)
  f = _jtheta1dash0(tau)
  h = -pi/6/w1 * _jtheta1dashdashdash0(tau) / f
  return w1 * exp.(h * z .* z / w1 / pi) .* j1 / f
end

"""
    wzeta(z; tau, omega, g)

Weierstrass zeta-function. One and only one of the parameters `tau`, `omega` or `g` must be given.

# Arguments
- `z`: real or complex number or array of numbers
- `tau`: half-periods ratio, complex number with non negative imaginary part
- `omega`: half-periods, a pair of complex numbers
- `g`: elliptic invariants, a pair of complex numbers
"""
function wzeta(z; tau::Union{Missing,Complex}=missing, omega::Union{Missing,Tuple{Number,Number}}=missing, g::Union{Missing,Tuple{Number,Number}}=missing)
  local omega1, omega2
  nmissing = ismissing(tau) + ismissing(omega) + ismissing(g)
  nmissing == 2 || throw(ArgumentError("You must supply either `tau`, `omega` or `g`."))
  if !ismissing(tau) || !ismissing(omega)
    if !ismissing(tau)
      imag(tau) > 0 || throw(ArgumentError("Invalid `tau`."))
      omega1 = 0.5
      omega2 = tau / 2
    elseif !ismissing(omega)
      omega1 = omega[1]
      omega2 = omega[2]
      tau = omega[2]/omega1
      imag(tau) > 0 || throw(ArgumentError("Invalid `omega`."))
    end
    if omega1 == Inf && omega2 == 1im*Inf # i.e. g2=0 g3=0
      return 1 ./ z
    end
    if omega1 == pi/sqrt(6) && omega2 == 1im*Inf # i.e. g2=3 g3=1
      return z/2 + sqrt(3/2) / tan.(sqrt(3/2)*z)
    end
  end
  if !ismissing(g)
    g2, g3 = g
    if g2 == 0 && g3 == 0
      return 1 ./ z
    end
    if g2 == 3 && g3 == 1
      return z/2 + sqrt(3/2) / tan.(sqrt(3/2)*z)
    end
    omega1, tau = _omega1_and_tau(g)
  end
  w1 = - omega1 / pi
  p = 1.0 / w1 / 2.0
  eta1 = p / 6.0 / w1 * _jtheta1dashdashdash0(tau) / _jtheta1dash0(tau)
  return - eta1 * z + p * _dljtheta1.(p * z, tau)
end

"""
    thetaS(z; tau, m)

Neville S-theta function. Only one of the parameters `tau` or `m` must be supplied.

# Arguments
- `z`: real or complex number or array of numbers
- `tau`: complex number with nonnegative imaginary part
- `m`: real or complex number, square of the elliptic modulus
"""
function thetaS(z; tau::Union{Missing,Complex}=missing, m::Union{Missing,Number}=missing)
  tau = _check_and_get_tau_from_m(tau, m)
  return _thetaS(z, tau)
end

"""
    thetaC(z; tau, m)

Neville C-theta function. Only one of the parameters `tau` or `m` must be supplied.

# Arguments
- `z`: real or complex number or array of numbers
- `tau`: complex number with nonnegative imaginary part
- `m`: real or complex number, square of the elliptic modulus
"""
function thetaC(z; tau::Union{Missing,Complex}=missing, m::Union{Missing,Number}=missing)
  tau = _check_and_get_tau_from_m(tau, m)
  return _thetaC(z, tau)
end

"""
    thetaD(z; tau, m)

Neville D-theta function. Only one of the parameters `tau` or `m` must be supplied.

# Arguments
- `z`: real or complex number or array of numbers
- `tau`: complex number with nonnegative imaginary part
- `m`: real or complex number, square of the elliptic modulus
"""
function thetaD(z; tau::Union{Missing,Complex}=missing, m::Union{Missing,Number}=missing)
  tau = _check_and_get_tau_from_m(tau, m)
  return _thetaD(z, tau)
end

"""
    thetaN(z; tau, m)

Neville N-theta function. Only one of the parameters `tau` or `m` must be supplied.

# Arguments
- `z`: real or complex number or vector/array of complex numbers
- `tau`: complex number with nonnegative imaginary part
- `m`: real or complex number, square of the elliptic modulus
"""
function thetaN(z; tau::Union{Missing,Complex}=missing, m::Union{Missing,Number}=missing)
  tau = _check_and_get_tau_from_m(tau, m)
  return _thetaN(z, tau)
end

"""
    jellip(kind, u; tau, m)

Jacobi elliptic functions. Only one of the parameters `tau` or `m` must be supplied.

# Arguments
- `kind`: a string with two characters among 'c', 'd', 'n' or 's'; this string specifies the function: the two letters respectively denote the basic functions `sn`, `cn`, `dn` and `1`, and the string specifies the ratio of two such functions, e.g. `ns=1/sn` and `cd=cn/dn`
- `u`: a real or complex number or array of numbers
- `tau`: complex number with nonnegative imaginary part
- `m`: real or complex number, square of the elliptic modulus
"""
function jellip(kind::String, u; tau::Union{Missing,Complex}=missing, m::Union{Missing,Number}=missing)
  local num, den
  length(kind) == 2 || throw(ArgumentError("The string `kind` must contain two characters."))
  f1, f2 = kind
  (f1 in "cdns" && f2 in "cdns") || throw(ArgumentError("Invalid string `kind`."))
  tau = _check_and_get_tau_from_m(tau, m)
  if f1 == 'c'
      num = _thetaC(u, tau)
  elseif f1 == 'd'
      num = _thetaD(u, tau)
  elseif f1 == 'n'
      num = _thetaN(u, tau)
  else
      num = _thetaS(u, tau)
  end
  if f2 == 'c'
      den = _thetaC(u, tau)
  elseif f2 == 'd'
      den = _thetaD(u, tau)
  elseif f2 == 'n'
      den = _thetaN(u, tau)
  else
      den = _thetaS(u, tau)
  end
  return num ./ den
end

"""
    am(u, m)

Amplitude function.

# Arguments
- `u`: real or complex number or array of numbers
- `m`: real or complex number, square of the elliptic modulus
"""
function am(u, m::Number)
  w = asin.(jellip("sn", u; m = m))
  k = round.(real.(u) / pi) + round.(real.(w) / pi)
  return (-1).^k .* w + k * pi
end

end  # module Jacobi
