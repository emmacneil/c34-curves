"""
  Summary of costs for tripling a divisor
  
  deg | I |  M | S |  A | 
  ----+---+----+---+----+
    1 | 1 | 35 | 2 | 43 | Typically
    1 | 1 | 26 | 1 | 33 | If type(3D) = 31 and tangent line at point is vert. or hor.
    1 | 1 | 25 | 0 | 33 | If type(3D) = 32
    1 | 0 | 23 | 0 | 32 | If type(3D) = 33
  ----+---+----+---+----+

  Summary of costs for tripling, including reduction step afterwards
  
  deg | I |  M | S |  A | 
  ----+---+----+---+----+
    1 | 1 | 35 | 2 | 43 | Typically
    1 | 1 | 26 | 1 | 33 | If type(3D) = 31 tangent line at point is vert. or hor.
    1 | 1 | 33 | 0 | 44 | If type(3D) = 32
    1 | 0 | 23 | 0 | 32 | If type(3D) = 33
  ----+---+----+---+----+
"""

def triple(D) :
  if D.degree == 0 :
    return C34CurveDivisor(D.C, [[D.K.one()], [], []])
  elif D.degree == 1 :
    return triple1(D)
  else :
    raise NotImplementedError("Tripling of divisors of degree {} not implemented.".format(D.degree))

def triple1(D) :
  K = D.K
  a, b = -D.f[0], -D.g[0]
  c0, c1, c2, c3, c4, c5, c6, c7, c8 = D.C.coefficients()

  # Precompute some powers of a and b
  aa  = a  * a
  ab  = a  * b
  bb  = b  * b
  aaa = aa * a

  # Precompute some multiples of things
  _3aa = aa + aa + aa
  _6aa = _3aa + _3aa
  _2aaa = aaa + aaa
  _4aaa = _2aaa + _2aaa
  _2a = a + a
  _3a = _2a + a
  _2b = b + b
  _3b = _2b + b
  _2bb = bb + bb
  _3bb = _2bb + bb
  _2ab = ab + ab
  
  d1 = _4aaa + c8*bb + c7*_2ab + c6*_3aa + c4*b + c3*_2a + c1
  d2 = _3bb + c8*_2ab + c7*aa + c5*_2b + c4*a + c2
  d3 = _6aa + c7*b + c6*_3a + c3
  d4 = c8*_2b + c7*_2a + c4
  d5 = _3b + c8*a + c5
  alpha = d1*(d1*d5 - d2*d4) + d2*d2*d3
  # Subtotal : 23M 32A
  # d1 and d2 are not both 0, since P is not a singular point.
  
  if (alpha != 0) :
    # XXX : I claim that P is not an inflection point.
    alpha = 1/alpha
    
    if (d1 == 0) :
      # Tangent line at P is horizontal.
      r0 = alpha * d2^2
      r2 = r0 * d2

      u0 = aa - r2*b
      u1 = - _2a
      u2 = r2
      v0 = ab
      v1 = -b
      v2 = -a
      w0 = bb
      w2 = - _2b
      # Subtotal : 1I 3M 1S 1A
      # Total    : 1I 26M 1S 33A
      # 3D is type 31
      return C34CurveDivisor(D.C, [[u0, u1, u2, 1], [v0, v1, v2, 0, 1], [w0, 0, w2, 0, 0, 1]],
                             degree = 3, typ = 31, typical = (u2 != 0), reduced = True)

    if (d2 == 0) :
      # Tangent line at P is vertical
      t0 = alpha * d1^2
      t1 = t0 * d1
      u0 = aa
      u1 = - _2a
      v0 = ab
      v1 = -b
      v2 = -a
      w0 = bb - t1*a
      w1 = t1
      w2 = - _2b
      # Subtotal : 1I 3M 1S 2A
      # Total    : 1I 26M 1S 33A
      # 3D is type 31
      return C34CurveDivisor(D.C, [[u0, u1, 0, 1], [v0, v1, v2, 0, 1], [w0, w1, w2, 0, 0, 1]],
                             degree = 3, typ = 31, typical = False, reduced = True)
    else :
      # Tangent line at P is neither horizontal nor vertical
      r0 = alpha * d2^2
      r1 = r0 * d1
      r2 = r0 * d2
      t0 = alpha * d1^2
      t1 = t0 * d1
      t2 = t0 * d2
      s1 = -t2
      s2 = -r1
      u0 = aa - r2*b - r1*a
      u1 = - _2a + r1
      u2 = r2
      v0 = ab - s2*b - s1*a
      v1 = s1 - b
      v2 = s2 - a
      w0 = bb - t2*b - t1*a
      w1 = t1
      w2 = -_2b + t2
      # Subtotal : 1I 12M 2S 11A
      # Total    : 1I 35M 2S 43A
      # 3D is type 31
      return C34CurveDivisor(D.C, [[u0, u1, u2, 1], [v0, v1, v2, 0, 1], [w0, w1, w2, 0, 0, 1]],
                             degree = 3, typ = 31, typical = (u2 != 0), reduced = True)

  else :
    # P is an inflection point.
    if (d2 != 0) :
      # Tangent line at P is not vertical
      beta = 1/d2
      z = beta*d1
      u0 = -b - a*z
      u1 = z
      v0 = -aaa
      v1 = _3aa
      v3 = - _3a
      # Subtotal : 1I 2M 1A
      # Total :    1I 25M 33A
      # 3D is type 32
      return C34CurveDivisor(D.C, [[u0, u1, 1], [v0, v1, 0, v3, 0, 0, 1], []],
                             degree = 3, typ = 32, typical = False, reduced = False)
    else :
      # Total : 23M 32A
      # 3D is type 33 (principal)
      return C34CurveDivisor(D.C, [[-a, 1], [], []],
                             degree = 3, typ = 33, typical = False, reduced = False)

