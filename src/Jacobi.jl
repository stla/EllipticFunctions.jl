module Jacobi

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
export CarlsonRF
export CarlsonRD
export ellipticF
export ellipticK
export ellipticE
export agm

function areclose(z1::Number, z2::Number)
  mod_z2 = abs(z2)
  maxmod = (mod_z2 < eps()) ? 1.0 : max(abs(z1), mod_z2)
  return abs(z1 - z2) < 2.0 * eps() * maxmod
end

function modulo(a::Real, p::Real)
  i = a > 0 ? floor(a / p) : ceil(a / p)
  return a - i * p
end

function calctheta3(z::Number, tau::Number, passes::Int64)
  out = complex(1.0, 0.0)
  n = 0
  while n < 2000
    n += 1
    qweight =
      exp(n * 1im * pi * (n * tau + 2 * z)) +
      exp(n * 1im * pi * (n * tau - 2 * z))
    out += qweight
    if n >= 3 && areclose(out + qweight, out)
      return log(out)
    end
  end
  error("Reached 2000 iterations.")
end

function argtheta3(z::Number, tau::Number, passin::Int64)
  local out
  passes = passin + 1
  if passes > 2000
    error("Reached 2000 iterations.")
  end
  zimg = imag(z)
  h = imag(tau) / 2
  zuse = complex(modulo(real(z), 1), zimg)
  if zimg < -h
    out = argtheta3(-zuse, tau, passes)
  elseif zimg >= h
    zmin = zuse - tau
    out = -2 * pi * 1im * zmin + argtheta3(zmin, tau, passes) - 1im * pi * tau
  else
    out = calctheta3(zuse, tau, passes)
  end
  return out
end

function dologtheta3(z::Number, tau::Number, passin::Int64)
  local tau2, out
  passes = passin + 1
  if passes > 2000
    error("Reached 2000 iterations.")
  end
  rl = real(tau)
  if rl >= 0
    tau2 = modulo(rl + 1.0, 2.0) - 1.0 + 1im * imag(tau)
  else
    tau2 = modulo(rl - 1.0, 2.0) + 1.0 + 1im * imag(tau)
  end
  if abs(tau2) < 0.98 && imag(tau2) < 0.98
    tauprime = -1.0 / tau2
    out =
      1im * pi * tauprime * z * z +
      dologtheta3(-z * tauprime, tauprime, passes) - log(sqrt(-1im * tau2))
  elseif rl >= 0.6
    out = dologtheta3(z + 0.5, tau2 - 1.0, passes)
  elseif rl <= -0.6
    out = dologtheta3(z + 0.5, tau2 + 1.0, passes)
  else
    out = argtheta3(z, tau2, passes)
  end
  return out
end

function M(z::Number, tau::Number)
  return 1im * (z + pi * tau / 4.0)
end

function _ljtheta2(z::Number, tau::Number)
  return M(z, tau) + dologtheta3(z / pi + 0.5 * tau, tau, 0)
end

function _jtheta2(z::Number, tau::Number)
  return exp(_ljtheta2(z, tau))
end

function _ljtheta1(z::Number, tau::Number)
  return _ljtheta2(z - pi / 2, tau)
end

function _jtheta1(z::Number, tau::Number)
  return exp(_ljtheta1(z, tau))
end

function _ljtheta3(z::Number, tau::Number)
  return dologtheta3(z / pi, tau, 0)
end

function _jtheta3(z::Number, tau::Number)
  return exp(_ljtheta3(z, tau))
end

function _ljtheta4(z::Number, tau::Number)
  return _ljtheta3(z + pi / 2, tau)
end

function _jtheta4(z::Number, tau::Number)
  return exp(_ljtheta4(z, tau))
end

function _jtheta1dash(z::Number, tau::Number)
  q = exp(1im * pi * tau)
  out = complex(0.0, 0.0)
  alt = -1.0
  for n = 0:2000
    alt = -alt
    k = 2.0 * n + 1.0
    outnew = out + alt * q^(n * (n + 1)) * k * cos(k * z)
    if areclose(out, outnew)
      return 2 * q^0.25 * out
    end
    out = outnew
  end
  error("Reached 2000 iterations.")
end

# exports ####

"""
    ljtheta1(z, tau)

Logarithm of the first Jacobi theta function.

# Arguments
- `z`: complex number
- `tau`: complex number with nonnegative imaginary part
"""
function ljtheta1(z::Number, tau::Number)
  if imag(tau) <= 0
    ArgumentError("Invalid `tau`.")
  end
  return _ljtheta1(z, tau)
end

"""
    jtheta1(z, tau)

First Jacobi theta function.

# Arguments
- `z`: complex number
- `tau`: complex number with nonnegative imaginary part
"""
function jtheta1(z::Number, tau::Number)
  if imag(tau) <= 0
    ArgumentError("Invalid `tau`.")
  end
  return _jtheta1(z, tau)
end

"""
    ljtheta2(z, tau)

Logarithm of the second Jacobi theta function.

# Arguments
- `z`: complex number
- `tau`: complex number with nonnegative imaginary part
"""
function ljtheta2(z::Number, tau::Number)
  if imag(tau) <= 0
    ArgumentError("Invalid `tau`.")
  end
  return _ljtheta2(z, tau)
end

"""
    jtheta2(z, tau)

Second Jacobi theta function.

# Arguments
- `z`: complex number
- `tau`: complex number with nonnegative imaginary part
"""
function jtheta2(z::Number, tau::Number)
  if imag(tau) <= 0
    ArgumentError("Invalid `tau`.")
  end
  return _jtheta2(z, tau)
end

"""
    ljtheta3(z, tau)

Logarithm of the third Jacobi theta function.

# Arguments
- `z`: complex number
- `tau`: complex number with nonnegative imaginary part
"""
function ljtheta3(z::Number, tau::Number)
  if imag(tau) <= 0
    ArgumentError("Invalid `tau`.")
  end
  return _ljtheta3(z, tau)
end

"""
    jtheta3(z, tau)

Third Jacobi theta function.

# Arguments
- `z`: complex number
- `tau`: complex number with nonnegative imaginary part
"""
function jtheta3(z::Number, tau::Number)
  if imag(tau) <= 0
    ArgumentError("Invalid `tau`.")
  end
  return _jtheta3(z, tau)
end

"""
    ljtheta4(z, tau)

Logarithm of the fourth Jacobi theta function.

# Arguments
- `z`: complex number
- `tau`: complex number with nonnegative imaginary part
"""
function ljtheta4(z::Number, tau::Number)
  if imag(tau) <= 0
    ArgumentError("Invalid `tau`.")
  end
  return _ljtheta4(z, tau)
end

"""
    jtheta4(z, tau)

Fourth Jacobi theta function.

# Arguments
- `z`: complex number
- `tau`: complex number with nonnegative imaginary part
"""
function jtheta4(z::Number, tau::Number)
  if imag(tau) <= 0
    ArgumentError("Invalid `tau`.")
  end
  return _jtheta4(z, tau)
end

"""
    jtheta1dash(z, tau)

Derivative of the first Jacobi theta function.

# Arguments
- `z`: complex number
- `tau`: complex number with nonnegative imaginary part
"""
function jtheta1dash(z::Number, tau::Number)
  if imag(tau) <= 0
    ArgumentError("Invalid `tau`.")
  end
  return _jtheta1dash(z, tau)
end

"""
    etaDedekind(tau)

Dedekind eta function.

# Arguments
- `tau`: complex number with nonnegative imaginary part
"""
function etaDedekind(tau::Number)
  if imag(tau) <= 0
    ArgumentError("Invalid `tau`.")
  end
  return exp(
    1im * pi * tau / 12.0 + dologtheta3((tau + 1.0) / 2.0, 3.0 * tau, 0)
  )
end

"""
    lambda(tau)

Lambda modular function.

# Arguments
- `tau`: complex number with nonnegative imaginary part
"""
function lambda(tau::Number)
  if imag(tau) <= 0
    ArgumentError("Invalid `tau`.")
  end
  return (_jtheta2(0, tau) / _jtheta3(0, tau))^4
end

"""
    kleinj(tau)

Klein j-invariant function.

# Arguments
- `tau`: complex number with nonnegative imaginary part
"""
function kleinj(tau::Number)
  if imag(tau) <= 0
    ArgumentError("Invalid `tau`.")
  end
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
  if xzero + yzero + zzero >= 2
    ArgumentError("At most one of `x`, `y`, `z` can be 0.")
  end
  dx = typemax(Float64)
  dy = typemax(Float64)
  dz = typemax(Float64)
  epsilon = 2.0 * eps()
  while dx > epsilon || dy > epsilon || dz > epsilon
    lambda = sqrt(x)*sqrt(y) + sqrt(y)*sqrt(z) + sqrt(z)*sqrt(x)
    x = (x + lambda) / 4.0
    y = (y + lambda) / 4.0
    z = (z + lambda) / 4.0
    A = (x + y + z) / 3.0
    dx = abs(1.0 - x/A)
    dy = abs(1.0 - y/A)
    dz = abs(1.0 - z/A)
  end
  E2 = dx*dy + dy*dz + dz*dx
  E3 = dy*dx*dz
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
  if xzero + yzero + zzero >= 2
    ArgumentError("At most one of `x`, `y`, `z` can be 0.")
  end
  dx = typemax(Float64)
  dy = typemax(Float64)
  dz = typemax(Float64)
  epsilon = 2.0 * eps()
  s = complex(0.0, 0.0)
  fac = complex(1.0, 0.0)
  while dx > epsilon || dy > epsilon || dz > epsilon
    lambda = sqrt(x)*sqrt(y) + sqrt(y)*sqrt(z) + sqrt(z)*sqrt(x)
    s = s + fac/(sqrt(z) * (z + lambda))
    fac = fac / 4.0
    x = (x + lambda) / 4.0
    y = (y + lambda) / 4.0
    z = (z + lambda) / 4.0
    A = (x + y + 3*z) / 5.0
    dx = abs(1.0 - x/A)
    dy = abs(1.0 - y/A)
    dz = abs(1.0 - z/A)
  end
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
  if real(phi) == 0 && imag(phi) == Inf && imag(m) == 0 && real(m) > 0 && real(m) < 1
    return sign(imag(phi)) *
        (ellipticF(pi/2, m) - ellipticF(pi/2, 1/m) / sqrt(m))
  end
  if abs(real(phi)) == pi/2 && m == 1
    return complex(NaN, NaN)
  end
  if real(phi) >= -pi/2 && real(phi) <= pi/2
    if m == 1 && abs(real(phi)) < pi/2
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
  if real(phi) > pi/2
    k = ceil((real(phi)-pi/2) / pi)
    phi = phi - k * pi
  else
    k = -floor((pi/2-real(phi)) / pi)
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
  if real(phi) >= -pi/2 && real(phi) <= pi/2
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
  if real(phi) > pi/2
    k = ceil((real(phi)-pi/2) / pi)
    phi = phi - k * pi
  else
    k = -floor((pi/2-real(phi)) / pi)
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

end  # module Jacobi
