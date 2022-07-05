module Jacobi

export ljtheta1
export jtheta1
export ljtheta2
export jtheta2
export ljtheta3
export jtheta3
export ljtheta4
export jtheta4

function close(z1::Number, z2::Number)
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
    if n >= 3 && close(out + qweight, out)
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
      1i * pi * tauprime * z * z +
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
  return M(z, tau) + dologtheta3(z/pi + 0.5 * tau, tau, 0)
end

function _jtheta2(z::Number, tau::Number)
  return exp(_ljtheta2(z, tau))
end

function _ljtheta1(z::Number, tau::Number)
    return _ljtheta2(z - pi/2, tau)
end

function _jtheta1(z::Number, tau::Number)
  return exp(_ljtheta1(z, tau))
end

function _ljtheta3(z::Number, tau::Number)
  return dologtheta3(z/pi, tau, 0)
end

function _jtheta3(z::Number, tau::Number)
  return exp(_ljtheta3(z, tau))
end

function _ljtheta4(z::Number, tau::Number)
  return _ljtheta3(z + pi/2, tau)
end

function _jtheta4(z::Number, tau::Number)
  return exp(_ljtheta4(z, tau))
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
        ArgumentError("Invalid tau.")
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
        ArgumentError("Invalid tau.")
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
        ArgumentError("Invalid tau.")
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
        ArgumentError("Invalid tau.")
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
        ArgumentError("Invalid tau.")
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
        ArgumentError("Invalid tau.")
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
        ArgumentError("Invalid tau.")
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
        ArgumentError("Invalid tau.")
    end
    return _jtheta4(z, tau)
end

end  # module Jacobi
