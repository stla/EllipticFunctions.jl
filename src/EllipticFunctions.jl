module EllipticFunctions

export ljtheta1
export jtheta1
export ljtheta2
export jtheta2
export ljtheta3
export jtheta3
export ljtheta4
export jtheta4
export jtheta1dash
export etaDedekind
export lambda
export kleinj
export kleinjinv
export CarlsonRF
export CarlsonRD
export ellipticF
export ellipticK
export ellipticE
export agm
export EisensteinE2
export EisensteinE4
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
    return exp(1im * pi * x)
end

function isqrt(x::Number)
  return sqrt(Complex(x))
end

function areclose(z1::Number, z2::Number)
  isinf(z1) && isinf(z2) && (return true)
  eps2 = eps()^2
  mod2_z2 = abs2(z2)
  maxmod2 = (mod2_z2 < eps2) ? 1.0 : max(abs2(z1), mod2_z2)
  return abs2(z1 - z2) < 4.0 * eps2 * maxmod2
end

################################# Beginning of althetas branch mods ##########################
# Use NIST definition where q = cispi(tau)
function _calctheta1_alt1(z::Number, q::Number)
  n = -1
  series = zero(promote_type(typeof(z), typeof(q)))
  maxiter = 3000
  q² = q * q
  q²ⁿ = one(q)
  qⁿ⁽ⁿ⁺¹⁾ = one(q)
  while n < maxiter
    n += 1
    if n > 0
      q²ⁿ *= q²
      qⁿ⁽ⁿ⁺¹⁾ *= q²ⁿ
    end
    term = qⁿ⁽ⁿ⁺¹⁾ * sin((2n+1)*z)
    isodd(n) && (term = -term)
    nextseries = series + term
    if n ≥ 2 && areclose(nextseries, series)
      return 2 * sqrt(sqrt(q)) * series
    else
      series = nextseries
    end
  end
  error("Reached $maxiter iterations.")
end

# Use Poisson transform of NIST definition:
"""
    _calctheta1_alt2(zopi::Number, topi::Number)
Calculate the Jacobian elliptic theta function θ₁ using the Poisson summation
formula.  Most useful for 0 < Im(tau) ≤ 1.3.  
# Input Arguments:
- `zopi`: z/π, where z is the first argument of the theta function.
- `topi`: t/π, where t = -i*τ, and τ is the second argument of the theta function.
"""
function _calctheta1_alt2(zopi::Number, topi::Number)
  nminus = round(Int, 0.5 - real(zopi)) + 1
  nplus = nminus - 1
  series = zero(promote_type(typeof(zopi), typeof(topi)))
  maxterms = 3000
  while nplus - nminus < maxterms
    nplus += 1
    nminus -= 1
    termm = exp(-(nminus - 0.5 + zopi)^2 / topi)
    termp = exp(-(nplus - 0.5 + zopi)^2 / topi)
    if isodd(nplus) 
      termp = -termp
    else
      termm = -termm
    end
    nextseries = series + (termp + termm)
    if nplus - nminus > 2 && areclose(nextseries, series)
      return sqrt(1/(pi * topi)) * series
    else
      series = nextseries
    end
  end
  error("Reached $maxterms terms.")
end

function alpha(z::Number, tau::Complex)
  return sqrt(-im * tau) * exp(im / tau * z * z / pi)
end

function _jtheta1(z::Number, tau::Complex)
  if imag(tau) > 1.3 # Chosen empirically
    # Large imag(tau) case: compute in terms of q
    q = xcispi(tau)
    if isreal(q) 
      if isreal(z)
        # Both inputs are real
        outr = _calctheta1_alt1(real(z), real(q))
        out = complex(outr, zero(outr))
      else
        # q is real, but z isn't
        out = _calctheta1_alt1(z, real(q))
      end
    else
      # q is not real
      out = im * _calctheta1_alt2(z / pi / tau, im / pi / tau) / alpha(z, tau)
    end
  else
    # Small imag(tau) case: compute in terms of t/pi where t = -im * tau
    topi = -im * (tau/pi)
    if isreal(topi)
      if isreal(z)
        # both z and t are real
        outr = _calctheta1_alt2(real(z)/pi, real(topi))
        out = complex(outr, 0)
      else
        # t is real but z isn't
        out = _calctheta1_alt2(z/pi, real(topi))
      end
    else
      # t is not real. No point in special casing real z here
      out = im * _calctheta1_alt1(z / tau, exp(-im * pi / tau)) / alpha(z, tau)
    end
  end
  return out
end

function _jtheta2(z::Number, tau::Complex) 
  return _jtheta1(z + pi/2, tau)
end

function expM(z::Number, tau::Complex) 
  return exp(im * (z + tau * pi/4))
end

function _jtheta3(z::Number, tau::Complex) 
  return _jtheta2(z - pi/2 * tau, tau) * expM(-z, tau)
end

function _jtheta4(z::Number, tau::Complex) 
  return _jtheta3(z + pi/2, tau);
end

function _ljtheta1(z::Number, tau::Complex)
  return log(_jtheta1(z, tau))
end

function _ljtheta2(z::Number, tau::Complex)
  return log(_jtheta2(z, tau))
end

function _ljtheta3(z::Number, tau::Complex)
  return log(_jtheta3(z, tau))
end

function _ljtheta4(z::Number, tau::Complex)
  return log(_jtheta4(z, tau))
end

function _jtheta1_raw(z::Number, tau::Complex)
  return _jtheta1(z * pi, tau)
end

function _jtheta2_raw(z::Number, tau::Complex)
  return _jtheta2(z * pi, tau)
end

function _jtheta3_raw(z::Number, tau::Complex)
  return _jtheta3(z * pi, tau)
end

function _jtheta4_raw(z::Number, tau::Complex)
  return _jtheta4(z * pi, tau)
end

function _jtheta1dash(z::Number, tau::Complex)
  q = xcispi(tau)
  out = complex(0.0, 0.0)
  alt = -1.0
  for n = 0:3000
    alt = -alt
    k = 2.0 * n + 1.0
    outnew = out + alt * q^(n * (n + 1)) * k * cos(k * z)
    if areclose(out, outnew)
      return 2 * sqrt(isqrt(q)) * out
    end
    out = outnew
  end
  error("Reached 3000 iterations.")
end

function _etaDedekind(tau::Complex)
  return xcispi(-1 / tau / 12.0) *
    _jtheta3_raw((-1 / tau + 1.0) / 2.0, -3.0 / tau) /
    sqrt(-1im * tau)
end

function isvector(x)
  return length(size(x)) == 1
end

function _EisensteinE2(tau::Complex)
  j3 = _jtheta3(0, tau)
  lbd = (_jtheta2(0, tau) / j3)^4
  j3sq = j3^2
  return 6.0/pi * ellipticE(lbd) * j3sq - j3sq^2 - _jtheta4(0, tau)^4
end

function _jtheta1dash0(tau::Complex)
  return exp(_ljtheta2(0.0, tau) + _ljtheta3(0.0, tau) + _ljtheta4(0.0, tau))
end

function _jtheta1dashdashdash0(tau::Complex)
  return -2.0 * _etaDedekind(tau)^3 * _EisensteinE2(tau)
end

function _dljtheta1(z::Number, tau::Complex)
  if z == 0
    return _jtheta1dash0(tau) / _jtheta1(0.0, tau)
  end
  return _jtheta1dash(z, tau) / _jtheta1(z, tau)
end

function _E4(tau::Complex)
  return (_jtheta2(0, tau)^8 + _jtheta3(0, tau)^8 + _jtheta4(0, tau)^8) / 2.0
end

function _omega1_and_tau(g)
  g2, g3 = g
  g2cube = g2*g2*g2
  j = 1728 * g2cube / (g2cube - 27*g3*g3)
  if isinf(j)
    return (-1im*pi/2/sqrt(3), complex(Inf, Inf))
  end
  tau = kleinjinv(j)
  omega1 = 1im * pi * sqrt(sqrt(1.0 / g2 / 12 * _E4(tau)))
  return (omega1, tau)
end

function _g2_from_omega1_and_tau(omega1::Number, tau::Complex)
  j2 = _jtheta2_raw(0, tau)
  j3 = _jtheta3_raw(0, tau)
  return 4/3 * (pi/2/omega1)^4 * (j2^8 - (j2*j3)^4 + j3^8)
end

function _wpFromTau(z, tau::Complex)
  j2 = _jtheta2(0, tau)
  j3 = _jtheta3(0, tau)
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
  j3sq = _jtheta3_raw(0, tau)^2
  zprime = z / j3sq / pi
  return j3sq * _jtheta1_raw.(zprime, tau) / _jtheta1dash0(tau)
end

function _thetaC(z, tau::Complex)
  zprime = z / _jtheta3_raw(0, tau)^2 / pi
  return _jtheta2_raw.(zprime, tau) / _jtheta2_raw(0, tau)
end

function _thetaN(z, tau::Complex)
  zprime = z / _jtheta3_raw(0, tau)^2 / pi
  return _jtheta4_raw.(zprime, tau) / _jtheta4_raw(0, tau)
end

function _thetaD(z, tau::Complex)
  j3 = _jtheta3_raw(0, tau)
  zprime = z / j3^2 / pi
  return _jtheta3_raw.(zprime, tau) / j3
end

function _tau_from_m(m::Number)
  1im * ellipticK(1.0 - m) / ellipticK(m)
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
    ljtheta1(z, tau)

Logarithm of the first Jacobi theta function.

# Arguments
- `z`: complex number or vector/array of complex numbers
- `tau`: complex number with nonnegative imaginary part
"""
function ljtheta1(z, tau::Complex)
  @assert imag(tau) > 0 ArgumentError("Invalid `tau`.")
  return _ljtheta1.(z, tau)
end

"""
    jtheta1(z, tau)

First Jacobi theta function.

# Arguments
- `z`: complex number or vector/array of complex numbers
- `tau`: complex number with nonnegative imaginary part
"""
function jtheta1(z, tau::Complex)
  @assert imag(tau) > 0 "Invalid `tau`."
  return _jtheta1.(z, tau)
end

"""
    ljtheta2(z, tau)

Logarithm of the second Jacobi theta function.

# Arguments
- `z`: complex number or vector/array of complex numbers
- `tau`: complex number with nonnegative imaginary part
"""
function ljtheta2(z, tau::Complex)
  @assert imag(tau) > 0 "Invalid `tau`."
  return _ljtheta2.(z, tau)
end

"""
    jtheta2(z, tau)

Second Jacobi theta function.

# Arguments
- `z`: complex number or vector/array of complex numbers
- `tau`: complex number with nonnegative imaginary part
"""
function jtheta2(z, tau::Complex)
  @assert imag(tau) > 0 "Invalid `tau`."
  return _jtheta2.(z, tau)
end

"""
    ljtheta3(z, tau)

Logarithm of the third Jacobi theta function.

# Arguments
- `z`: complex number or vector/array of complex numbers
- `tau`: complex number with nonnegative imaginary part
"""
function ljtheta3(z, tau::Complex)
  @assert imag(tau) > 0 "Invalid `tau`."
  return _ljtheta3.(z, tau)
end

"""
    jtheta3(z, tau)

Third Jacobi theta function.

# Arguments
- `z`: complex number or vector/array of complex numbers
- `tau`: complex number with nonnegative imaginary part
"""
function jtheta3(z, tau::Complex)
  @assert imag(tau) > 0 "Invalid `tau`."
  return _jtheta3.(z, tau)
end

"""
    ljtheta4(z, tau)

Logarithm of the fourth Jacobi theta function.

# Arguments
- `z`: complex number or vector/array of complex numbers
- `tau`: complex number with nonnegative imaginary part
"""
function ljtheta4(z, tau::Complex)
  @assert imag(tau) > 0 "Invalid `tau`."
  return _ljtheta4.(z, tau)
end

"""
    jtheta4(z, tau)

Fourth Jacobi theta function.

# Arguments
- `z`: complex number or vector/array of complex numbers
- `tau`: complex number with nonnegative imaginary part
"""
function jtheta4(z, tau::Complex)
  @assert imag(tau) > 0 "Invalid `tau`."
  return _jtheta4.(z, tau)
end

"""
    jtheta1dash(z, tau)

Derivative of the first Jacobi theta function.

# Arguments
- `z`: complex number or vector/array of complex numbers
- `tau`: complex number with nonnegative imaginary part
"""
function jtheta1dash(z, tau::Complex)
  @assert imag(tau) > 0 "Invalid `tau`."
  return _jtheta1dash.(z, tau)
end

"""
    etaDedekind(tau)

Dedekind eta function.

# Arguments
- `tau`: complex number with nonnegative imaginary part
"""
function etaDedekind(tau::Complex)
  @assert imag(tau) > 0 "Invalid `tau`."
  return _etaDedekind(tau)
end

"""
    lambda(tau)

Lambda modular function.

# Arguments
- `tau`: complex number with nonnegative imaginary part
"""
function lambda(tau::Complex)
  @assert imag(tau) > 0 "Invalid `tau`."
  return (_jtheta2(0, tau) / _jtheta3(0, tau))^4
end

"""
    kleinj(tau)

Klein j-invariant function.

# Arguments
- `tau`: complex number with nonnegative imaginary part
"""
function kleinj(tau::Complex)
  @assert imag(tau) > 0 "Invalid `tau`."
  lbd = (_jtheta2(0, tau) / _jtheta3(0, tau))^4
  x = lbd * (1.0 - lbd)
  return 256 * (1-x)^3 / x^2
end

"""
    CarlsonRF(x, y, z)

Carlson 'RF' integral.

# Arguments
- `x`,`y`,`z`: complex numbers; at most one of them can be zero
"""
function CarlsonRF(x::Number, y::Number, z::Number)
  local A
  xzero = x == 0
  yzero = y == 0
  zzero = z == 0
  @assert xzero + yzero + zzero <= 1 "At most one of `x`, `y`, `z` can be 0."
  dx = typemax(Float64)
  dy = typemax(Float64)
  dz = typemax(Float64)
  epsilon = 2.0 * eps()^2
  while dx > epsilon || dy > epsilon || dz > epsilon
    lambda = isqrt(x)*isqrt(y) + isqrt(y)*isqrt(z) + isqrt(z)*isqrt(x)
    x = (x + lambda) / 4.0
    y = (y + lambda) / 4.0
    z = (z + lambda) / 4.0
    A = (x + y + z) / 3.0
    dx = abs2(1.0 - x/A)
    dy = abs2(1.0 - y/A)
    dz = abs2(1.0 - z/A)
  end
  E2 = sqrt(dx*dy) + sqrt(dy*dz) + sqrt(dz*dx)
  E3 = sqrt(dy*dx*dz)
  return (1 - E2/10 + E3/14 + E2*E2/24 - 3*E2*E3/44 - 5*E2*E2*E2/208 +
    3*E3*E3/104 + E2*E2*E3/16) / sqrt(A)
end

"""
    CarlsonRD(x, y, z)

Carlson 'RD' integral.

# Arguments
- `x`,`y`,`z`: complex numbers; at most one of them can be zero
"""
function CarlsonRD(x::Number, y::Number, z::Number)
  local A
  xzero = x == 0
  yzero = y == 0
  zzero = z == 0
  @assert xzero + yzero + zzero <= 1 "At most one of `x`, `y`, `z` can be 0."
  dx = typemax(Float64)
  dy = typemax(Float64)
  dz = typemax(Float64)
  epsilon = 2.0 * eps()^2
  s = complex(0.0, 0.0)
  fac = complex(1.0, 0.0)
  while dx > epsilon || dy > epsilon || dz > epsilon
    lambda = isqrt(x)*isqrt(y) + isqrt(y)*isqrt(z) + isqrt(z)*isqrt(x)
    s = s + fac/(isqrt(z) * (z + lambda))
    fac = fac / 4.0
    x = (x + lambda) / 4.0
    y = (y + lambda) / 4.0
    z = (z + lambda) / 4.0
    A = (x + y + 3*z) / 5.0
    dx = abs2(1.0 - x/A)
    dy = abs2(1.0 - y/A)
    dz = abs2(1.0 - z/A)
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
    ellipticF(phi, m)

Incomplete elliptic integral of the first kind.

# Arguments
- `phi`: complex number, the amplitude
- `m`: complex number, the squared modulus
"""
function ellipticF(phi::Number, m::Number)
  local k
  if phi == 0 || m == Inf || m == -Inf
    return complex(0.0, 0.0)
  end
  rphi = real(phi)
  if rphi == 0 && imag(phi) == Inf && imag(m) == 0 && real(m) > 0 && real(m) < 1
    return sign(imag(phi)) *
        (ellipticF(pi/2, m) - ellipticF(pi/2, 1/m) / isqrt(m))
  end
  if abs(rphi) == pi/2 && m == 1
    return complex(NaN, NaN)
  end
  if rphi >= -pi/2 && rphi <= pi/2
    if m == 1 && abs(rphi) < pi/2
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
  if rphi > pi/2
    k = ceil((rphi-pi/2) / pi)
    phi = phi - k * pi
  else
    k = -floor((pi/2-rphi) / pi)
    phi = phi - k * pi
  end
  return 2*k*ellipticF(pi/2, m) + ellipticF(phi, m)
end

"""
    ellipticK(m)

Complete elliptic integral of the first kind.

# Arguments
- `m`: complex number, the squared modulus
"""
function ellipticK(m::Number)
  return ellipticF(pi/2, m)
end

"""
    ellipticE(phi, m)

Incomplete elliptic integral of the second kind.

# Arguments
- `phi`: complex number, the amplitude
- `m`: complex number, the squared modulus
"""
function ellipticE(phi::Number, m::Number)
  local k
  if phi == 0
    return complex(0.0, 0.0)
  end
  if real(m) == Inf && imag(m) == 0
    return complex(NaN, NaN)
  end
  rphi = real(phi)
  if rphi >= -pi/2 && rphi <= pi/2
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
  if rphi > pi/2
    k = ceil((rphi-pi/2) / pi)
    phi = phi - k * pi
  else
    k = -floor((pi/2-rphi) / pi)
    phi = phi - k * pi
  end
  return 2 * k * ellipticE(pi/2, m) + ellipticE(phi, m)
end

"""
    ellipticE(m)

Complete elliptic integral of the second kind.

# Arguments
- `m`: complex number, the squared modulus
"""
function ellipticE(m::Number)
  return ellipticE(pi/2, m)
end

"""
    agm(x, y)

Arithmetic-geometric mean.

# Arguments
- `x`,`y`: complex numbers
"""
function agm(x::Number, y::Number)
  if x + y == 0 || x == 0 || y == 0
    return complex(0.0, 0.0)
  end
  return pi/4 * (x + y) / ellipticK(((x-y)/(x+y))^2)
end

"""
    kleinjinv(j)

Inverse of the Klein j-invariant function.

# Arguments
- `j`: complex number
"""
function kleinjinv(j::Number)
  local x
  if isinf(j)
    x = complex(0.0, 0.0)
  else
    j2 = j * j
    j3 = j2 * j
    t = (-j3 + 2304 * j2 + 12288 *
          isqrt(3 * (1728 * j2 - j3)) - 884736 * j)^(1/3)
    x = 1/768 * t - (1536 * j - j2) / (768 * t) + (1 - j/768)
  end
  lbd = -(-1 - sqrt(1 - 4*x)) / 2
  return 1im * agm(1, sqrt(1-lbd)) / agm(1, sqrt(lbd))
end

"""
    EisensteinE2(q)

Eisenstein E-series of weight 2.

# Arguments
- `q`: nome, complex number; it must not be a negative real number and its modulus must be strictly smaller than 1
"""
function EisensteinE2(q::Number)
  @assert abs(q) < 1 "Invalid `q`."
  @assert imag(q) != 0 || real(q) > 0 "Invalid `q`."
  tau = -1im * log(q) / pi / 2.0
  return _EisensteinE2(tau)
end

"""
    EisensteinE4(q)

Eisenstein E-series of weight 4.

# Arguments
- `q`: nome, complex number; it must not be a negative real number and its modulus must be strictly smaller than 1
"""
function EisensteinE4(q::Number)
  @assert abs(q) < 1 "Invalid `q`."
  @assert imag(q) != 0 || real(q) > 0 "Invalid `q`."
  tau = -1im * log(q) / pi / 2.0
  return (_jtheta2(0, tau)^8 + _jtheta3(0, tau)^8 + _jtheta4(0, tau)^8) / 2.0
end

"""
    wp(z; tau, omega)

Weierstrass p-function. One and only one of the parameters `tau`, `omega` or `g` must be given.

# Arguments
- `z`: complex number or vector/array of complex numbers
- `tau`: half-periods ratio, complex number with non negative imaginary part
- `omega`: half-periods, a pair of complex numbers
- `g`: elliptic invariants, a pair of complex numbers
- `derivative`: order of differentiation, 0, 1, 2 or 3
"""
function wp(z; tau::Union{Missing,Number}=missing, omega::Union{Missing,Tuple{Number,Number}}=missing, g::Union{Missing,Tuple{Number,Number}}=missing, derivative::Int64=0)
  local omega1, weier, weierPrime
  @assert derivative >= 0 && derivative <= 3 ArgumentError("`derivative` must be beetween 0 and 3.")
  nmissing = ismissing(tau) + ismissing(omega) + ismissing(g)
  @assert nmissing == 2 ArgumentError("You must supply either `tau`, `omega` or `g`.")
  if !ismissing(tau)
    @assert imag(tau) > 0 "Invalid `tau`."
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
    @assert imag(tau) > 0 "Invalid `omega`."
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
    wsigma(z; tau, omega)

Weierstrass sigma-function. One and only one of the parameters `tau`, `omega` or `g` must be given.

# Arguments
- `z`: complex number or vector/array of complex numbers
- `tau`: half-periods ratio, complex number with non negative imaginary part
- `omega`: half-periods, a pair of complex numbers
- `g`: elliptic invariants, a pair of complex numbers
"""
function wsigma(z; tau::Union{Missing,Number}=missing, omega::Union{Missing,Tuple{Number,Number}}=missing, g::Union{Missing,Tuple{Number,Number}}=missing)
  local omega1
  nmissing = ismissing(tau) + ismissing(omega) + ismissing(g)
  @assert nmissing == 2 ArgumentError("You must supply either `tau`, `omega` or `g`.")
  if !ismissing(tau)
    @assert imag(tau) > 0 "Invalid `tau`."
    omega1 = 0.5
  elseif !ismissing(omega)
    omega1 = omega[1]
    tau = omega[2]/omega1
    @assert imag(tau) > 0 "Invalid `omega`."
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
    wzeta(z; tau, omega)

Weierstrass zeta-function. One and only one of the parameters `tau`, `omega` or `g` must be given.

# Arguments
- `z`: complex number or vector/array of complex numbers
- `tau`: half-periods ratio, complex number with non negative imaginary part
- `omega`: half-periods, a pair of complex numbers
- `g`: elliptic invariants, a pair of complex numbers
"""
function wzeta(z; tau::Union{Missing,Number}=missing, omega::Union{Missing,Tuple{Number,Number}}=missing, g::Union{Missing,Tuple{Number,Number}}=missing)
  local omega1, omega2
  nmissing = ismissing(tau) + ismissing(omega) + ismissing(g)
  @assert nmissing == 2 ArgumentError("You must supply either `tau`, `omega` or `g`.")
  if !ismissing(tau) || !ismissing(omega)
    if !ismissing(tau)
      @assert imag(tau) > 0 "Invalid `tau`."
      omega1 = 0.5
      omega2 = tau / 2
    elseif !ismissing(omega)
      omega1 = omega[1]
      omega2 = omega[2]
      tau = omega[2]/omega1
      @assert imag(tau) > 0 "Invalid `omega`."
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
    thetaS(z, tau)

Neville S-theta function. Only one of the parameters `tau` or `m` must be supplied.

# Arguments
- `z`: complex number or vector/array of complex numbers
- `tau`: complex number with nonnegative imaginary part
- `m`: complex number, square of the elliptic modulus
"""
function thetaS(z; tau::Union{Missing,Number}=missing, m::Union{Missing,Number}=missing)
  tau = _check_and_get_tau_from_m(tau, m)
  return _thetaS(z, tau)
end

"""
    thetaC(z, tau)

Neville C-theta function. Only one of the parameters `tau` or `m` must be supplied.

# Arguments
- `z`: complex number or vector/array of complex numbers
- `tau`: complex number with nonnegative imaginary part
- `m`: complex number, square of the elliptic modulus
"""
function thetaC(z; tau::Union{Missing,Number}=missing, m::Union{Missing,Number}=missing)
  tau = _check_and_get_tau_from_m(tau, m)
  return _thetaC(z, tau)
end

"""
    thetaD(z, tau)

Neville D-theta function. Only one of the parameters `tau` or `m` must be supplied.

# Arguments
- `z`: complex number or vector/array of complex numbers
- `tau`: complex number with nonnegative imaginary part
- `m`: complex number, square of the elliptic modulus
"""
function thetaD(z; tau::Union{Missing,Number}=missing, m::Union{Missing,Number}=missing)
  tau = _check_and_get_tau_from_m(tau, m)
  return _thetaD(z, tau)
end

"""
    thetaN(z, tau)

Neville N-theta function. Only one of the parameters `tau` or `m` must be supplied.

# Arguments
- `z`: complex number or vector/array of complex numbers
- `tau`: complex number with nonnegative imaginary part
- `m`: complex number, square of the elliptic modulus
"""
function thetaN(z; tau::Union{Missing,Number}=missing, m::Union{Missing,Number}=missing)
  tau = _check_and_get_tau_from_m(tau, m)
  return _thetaN(z, tau)
end

"""
    jellip(kind, u; tau, m)

Jacobi elliptic functions. Only one of the parameters `tau` or `m` must be supplied.

# Arguments
- `kind`: a string with two characters among 'c', 'd', 'n' or 's'; this string specifies the function: the two letters respectively denote the basic functions `sn`, `cn`, `dn` and `1`, and the string specifies the ratio of two such functions, e.g. `ns=1/sn` and `cd=cn/dn`
- `u`: complex number or vector/array of complex numbers
- `tau`: complex number with nonnegative imaginary part
- `m`: complex number, square of the elliptic modulus
"""
function jellip(kind::String, u; tau::Union{Missing,Number}=missing, m::Union{Missing,Number}=missing)
  local num, den
  @assert length(kind) == 2 ArgumentError("The string `kind` must contain two characters.")
  f1 = kind[1]
  f2 = kind[2]
  @assert f1 == 'c' || f1 == 'd'  || f1 == 'n'  || f1 == 's' ArgumentError("Invalid string `kind`.")
  @assert f2 == 'c' || f2 == 'd'  || f2 == 'n'  || f2 == 's' ArgumentError("Invalid string `kind`.")
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
- `u`: complex number or vector/array of complex numbers
- `m`: complex number, square of the elliptic modulus
"""
function am(u, m::Number)
  w = asin.(jellip("sn", u; m = m))
  k = round.(real.(u) / pi) + round.(real.(w) / pi)
  return (-1).^k .* w + k * pi
end

end  # module Jacobi
