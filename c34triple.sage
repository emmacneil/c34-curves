"""
  Summary of costs for doubling a divisor
  
  
  deg |   I |   M |
  ----+-----+-----+
    0 |   0 |   0 |
  ----+-----+-----+
    1 |   0 |   8 | 
"""

def triple(D) :
  if D.degree == 0 :
    return C34CrvDiv(D.C, [[D.K.one()], [], []])
  elif D.degree == 1 :
    return triple1(D)
  else :
    raise NotImplementedError("Tripling of divisors of degree {} not implemented.".format(D.degree))

def triple1(D) :
  K = D.K
  a, b = -D.f[0], -D.g[0]
  c = D.C.coefficients()
  new_f = []
  new_g = []
  new_h = []

  # Precompute some powers of a and b
  aa  = a  * a
  ab  = a  * b
  bb  = b  * b
  aaa = aa * a
  
  d1 = 4*aaa + c[8]*bb + 2*c[7]*ab + 3*c[6]*aa + c[4]*b + 2*c[3]*a + c[1]
  d2 = 3*bb + 2*c[8]*ab + c[7]*aa + 2*c[5]*b + c[4]*a + c[2]
  d3 = 6*aa + c[7]*b + 3*c[6]*a + c[3]
  d4 = 2*c[8]*b + 2*c[7]*a + c[4]
  d5 = 3*b + c[8]*a + c[5]
  alpha = d1*(d1*d5 - d2*d4) + d2*d2*d3
  # Subtotal : 9M 26C
  # d1 and d2 are not both 0, since P is not a singular point.
  
  if (alpha != 0) :
    # XXX : I claim that P is not an inflection point.
    alpha = 1/alpha
    
    if (d1 == 0) :
      # Tangent line at P is horizontal.
      r0 = alpha * d2^2
      r1 = r0 * d1
      r2 = r0 * d2
      new_f = [ aa - r2*b - r1*a, -2*a + r1, r2, K.one() ]
      new_g = [ ab, -b, -a, K.zero(), K.one() ]
      new_h = [ bb, K.zero(), -2*b, K.zero(), K.zero(), K.one() ]
      # Total : 1I 14M 1SQ 28C
    elif (d2 == 0) :
      # Tangent line at P is vertical
      t0 = alpha * d1^2
      t1 = t0 * d1
      t2 = t0 * d2
      new_f = [ aa, -2*a, K.zero(), K.one() ]
      new_g = [ ab, -b, -a, K.zero(), K.one() ]
      new_h = [ bb - t2*b - t1*a, t1, -2*b + t2, K.zero(), K.zero(), K.one() ]
      # Total : 1I 14M 1SQ 28C
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
      new_f = [ aa - r2*b - r1*a, -2*a + r1, r2, K.one() ]
      new_g = [ ab - s2*b - s1*a, s1 - b, s2 - a, K.zero(), K.one() ]
      new_h = [ bb - t2*b - t1*a, t1, -2*b + t2, K.zero(), K.zero(), K.one() ]
      # Total : 1I 21M 2SQ 28C

  else :
    # XXX : I claim that P is an inflection point.
    if (d2 != 0) :
      # Tangent line at P is not vertical
      beta = 1/d2
      z = beta*d1
      new_f = [ - b - z*a, z, K.one() ]
      new_g = [ -aaa, 3*aa, K.zero(), -3*a, K.zero(), K.zero(), K.one() ]
      # Subtotal : 1I 11M 28C
    else :
      # Tangent line at P is vertical
      new_f = [ -a, K.one() ]
      # Total : 9M 26C

  return C34CrvDiv(D.C, [new_f, new_g, new_h])

