"""
  Summary of costs for flipping a divisor
  
  
  Deg | Type | I |  M | S |  A |
  ----+------+---+----+---+----+
    0 |    0 | 0 |  0 | 0 |  0 |
  ----+------+---+----+---+----+
    1 |   11 | 0 |  4 | 0 |  5 |
  ----+------+---+----+---+----+
    2 |   21 | 0 |  7 | 0 | 12 |
      |   22 | 0 |  1 | 0 |  2 |
  ----+------+---+----+---+----+
    3 |   31 | 1 | 15 | 0 | 18 | Typical
      |   31 | 0 | 12 | 0 | 20 | Semi-typical
      |   32 | 0 |  3 | 0 |  6 | 
      |   33 | 0 |  0 | 0 |  0 | 
  ----+------+---+----+---+----+
    4 |   41 | 1 | 23 | 1 | 26 | Typical
      |   41 | 1 | 32 | 1 | 38 | Semi-typical (Can be improved?)
      |   41 | 0 | 28 | 2 | 46 | Case where <f, g> =/= <f, g, h> =/= <f, h>
      |   42 | 0 |  8 | 0 |  9 | 
      |   43 | 0 |  5 | 0 |  7 | 
      |   44 | 0 |  0 | 0 |  0 | 
  ----+------+---+----+---+----+
    5 |   51 | 1 | 24 | 0 | 28 | Typical
      |   51 | 1 | 21 | 0 | 31 | Semi-typical
      |   51 | 0 | 33 | 0 | 52 | Case where <f, g> =/= <f, g, h> =/= <f, h>
      |   52 | 0 | 22 | 1 | 27 | 
      |   53 | 0 | 25 | 1 | 31 | 
      |   54 | 0 |  3 | 1 |  6 | 
  ----+------+---+----+---+----+
    6 |   61 | 1 | 24 | 0 | 41 | Typical
      |   61 | 0 | 25 | 0 | 45 | Semi-typical
      |   61 | 0 | 44 | 0 | 80 | Case where <f, g> =/= <f, g, h> =/= <f, h>
      |   62 | 0 |  4 | 0 |  6 | 
      |   63 | 0 | 18 | 0 | 28 | 
      |   64 | 0 |  3 | 0 |  9 |
      |   65 | 0 |  0 | 0 |  0 | 
"""


def flip(D) :
  """
    Compute the 'flip' of a divisor.
    
    A typical (or semi-typical) divisor D is represented by two (or three) polynomials, f and g
    (or f, g, h). Then the divisor of f is
    
      div f = D + A
    
    The divisor A is the 'flip' of D.
    It is simply negation in the Jacobian of the curve.
    
    Input : A typical or semi-typical C34CurveDivisor D.
    Output : A typical or semi-typical C34CurveDivisor A, the flip of D.
  """
  if (D.type == 0) :
    return flip_0(D)
  elif (D.type == 11) :
    return flip_11(D)
  elif (D.type == 21) :
    return flip_21(D)
  elif (D.type == 22) :
    return flip_22(D)
  elif (D.type == 31) :
    return flip_31(D)
  elif (D.type == 32) :
    return flip_32(D)
  elif (D.type == 33) :
    return flip_33(D)
  elif (D.type == 41) :
    return flip_41(D)
  elif (D.type == 42) :
    return flip_42(D)
  elif (D.type == 43) :
    return flip_43(D)
  elif (D.type == 44) :
    return flip_44(D)
  elif (D.type == 51) :
    return flip_51(D)
  elif (D.type == 52) :
    return flip_52(D)
  elif (D.type == 53) :
    return flip_53(D)
  elif (D.type == 54) :
    return flip_54(D)
  elif (D.type == 61) :
    return flip_61(D)
  elif (D.type == 62) :
    return flip_62(D)
  elif (D.type == 63) :
    return flip_63(D)
  elif (D.type == 64) :
    return flip_64(D)
  elif (D.type == 65) :
    return flip_65(D)
  else :
    raise NotImplementedError("Flipping of divisors of type {} not implemented.\nD = {}.".format(D.type, D))



def flip_0(D) :
  # A is of type 0
  # Total : 0I 0M
  return D.C.zero_divisor()



def flip_11(D) :
  c0, c1, c2, c3, c4, c5, c6, c7, c8 = D.C.coefficients()
  f0, g0 = D.f[0], D.g[0]
  c = D.C.coefficients()

  a = (c4 - c7*f0)*f0 - c2
  b = c8*f0 - c5 + g0
  gg0 = b*g0 - a
  gg2 = -b
  # A is type 22
  # Total : 0I 4M 5A
  return C34CurveDivisor(D.C, [[f0, 1], [gg0, 0, gg2, 0, 0, 1], []],
                         degree = 2, typ = 22, reduced = True, typical = False)



def flip_21(D) :
  c0, c1, c2, c3, c4, c5, c6, c7, c8 = D.C.coefficients()
  f0, f1 = D.f[0:2]
  g0, g1 = D.g[0:2]
  
  a = c7 + f1*(f1 - c8)
  b = g0 - c3 + f0*a + f1*(c4 + f0*(f1 - c8) + f1*(f0 - c5))
  c = g1 - c6 + f1*a
  gg0 = c*g1 - b
  gg1 = -c
  # A is type 21
  # Total : 0I 7M 12A
  return C34CurveDivisor(D.C, [[f0, f1, 1], [gg0, gg1, 0, 1], []],
                         degree = 2, typ = 21, reduced = True, typical = False)



def flip_22(D) :
  c0, c1, c2, c3, c4, c5, c6, c7, c8 = D.C.coefficients()
  f0 = D.f[0]
  g2 = D.g[2]
  gg0 = c5 - g2 - c8*f0
  # A is type 11
  # Total : 0I 1M 2A
  return C34CurveDivisor(D.C, [[f0, 1], [gg0, 0, 1], []],
                         degree = 1, typ = 11, reduced = True, typical = False)



def flip_31(D) :
  c0, c1, c2, c3, c4, c5, c6, c7, c8 = D.C.coefficients()
  f0, f1, f2 = D.f[0:3]
  g0, g1, g2 = D.g[0:3]
  h0, h1, h2 = D.h[0:3]
  
  if D.typical : 
    s0 = f1 - g2
    
    z0 = c7 - f2
    z1 = c6 - f1
    a1 = h1 - c4 + z0*(g2 + s0) + f2*z1
    a2 = h2 - c5 + f2*z0
    
    v0 = f2*(a1 + g1*c8)
    v1 =  - (a2 + g2*c8)
    v2 =  - f2*c8
    
    v0 = v0 + v1*s0
    v2 = v2 + s0

    u0, u1, u2 = f0, f1, f2

    alpha = 1 / u2
    v2mu1 = v2 - u1
    w0 = alpha*(u0*v1 + v0*v2mu1)
    w1 = alpha*(v1*v2 - v0)
    w2 = v1 + alpha*(u0 + v2*v2mu1)
    # A is of type 31, typical
    # Total : 1I 15M 18A
    # (In high characterstic, can be reduced to 12M)
    return C34CurveDivisor(D.C, [[u0, u1, u2, 1], [v0, v1, v2, 0, 1], [w0, w1, w2, 0, 0, 1]],
                           degree = 3, typ = 31, typical = True, reduced = True)
    
  else :
    s0 = f1 - g2
    s2 = -g1 + c5 - c8*s0
    
    r0 = -c8*f0
    r1 = c5 - c8*f1

    t0 = c5 - h2
    t1 = c8
    t2 = r0
    t3 = r1 - h2
    t4 = c2 - c8*r0 - c7*f0 - h0 + c5*(h2 - c5)
    t5 = c4 - c7*f1 - h1 + c8*(h2 - r1 - c5)
    
    u0 = s0*t0 + t2
    u1 = s0*t1 + t3
    u2 = s0
    v0 = s2*t0 + t4
    v1 = s2*t1 + t5
    v2 = s2

    # A is type 31, non-typical
    # Total: 0I 12M 20A
    return C34CurveDivisor(D.C, [[f0, f1, 0, 1], [u0, u1, u2, 0, 1], [v0, v1, v2, 0, 0, 1]],
                           degree = 3, typ = 31, typical = False, reduced = True)



def flip_32(D) :
  c0, c1, c2, c3, c4, c5, c6, c7, c8 = D.C.coefficients()
  f0, f1 = D.f[0:2]
  g3     = D.g[3]
  
  a1 = -f1*(f1 - c8)
  a2 = g3 + f1*(c7 - a1)
  a3 = -f0 + f1*(c6 - a2)

  ff0 = c6 - a2
  gg0 = -a3
  # A is of type 11
  # Total 0I 3M 6A
  return C34CurveDivisor(D.C, [[ff0, 1], [gg0, 0, 1], []],
                         degree = 1, typ = 11, reduced = True, typical = False)



def flip_33(D) :
  # A is of type 0
  # Total 0I 0M
  return D.C.zero_divisor()



def flip_41(D) :
  c0, c1, c2, c3, c4, c5, c6, c7, c8 = D.C.coefficients()
  f0, f1, f2, f3 = D.f[0:4]
  g0, g1, g2, g3 = D.g[0:4]
  h0, h1, h2, h3 = D.h[0:4]

  if D.typical :
    s0 = f2
    t0 = - f3^2 - g3
    
    c7f3 = c7*f3
    z0 = c8*f3
    z1 = h3 - f2 + c7f3
    a1 = c5 - f2*c8
    a2 = h2 + c5*f3 - f2*z0
    a3 = c6 - h3 - c8*(g3 + t0) - c7f3
    a4 = f1 + f3*(z1 - c6) + z0*(g3 + t0)
    
    f3g3 = f3*g3
    u0 = t0*(g2 - a1)
    u1 = a3 + f3g3
    u2 = - t0
    v0 = t0*(a2 - g2*f3)
    v1 = a4 - f3*f3g3
    v2 = t0*f3
    
    u0 = u0 + u1*s0
    u1 = u1 + s0
    
    v0 = v0 + v1*s0
    v2 = v2 + s0

    alpha = 1 / u2
    v2mu1 = v2 - u1
    w0 = alpha*(u0*v1 + v0*v2mu1)
    w1 = alpha*(v1*v2 - v0)
    w2 = v1 + alpha*(u0 + v2*v2mu1)

    # A is of type 31, typical
    # Total : 1I 23M 1S 26A
    # (In high characteristic, is is only 18M)
    return C34CurveDivisor(D.C, [[u0, u1, u2, 1], [v0, v1, v2, 0, 1], [w0, w1, w2, 0, 0, 1]],
                           degree = 3, typ = 31, typical = True, reduced = True)

  elif h2 != 0 :
    e0 = f1 - f2*f3
    e1 = - f2^2
    e2 = g1 - f3*g2
    e3 = c3 + c6*f2 - e2*(c8 - f3) - c7*f1 - g3*(c5 + c8*f2 - f1) - c4*f3
    e4 = c6 + f2 - g3*(c8 - f3) - c7*f3
    e7 = - f2*e1
    e8 = - g3*e3 - g2*e0 + g0
    e9 = - g3*e4 + e2
    
    t0 = h1 - c8*h2
    a2 = h0 - c5*h2 - e7 - h3*e1 - f2*t0
    a3 = - e0 - h3*f3
    a4 = h1 - e3 - h2*f3
    a6 = - e8 - h3*e2 - g3*t0 + h2*(c7*f3 - c6)
    a7 = h3 - e4
    a9 = - h2 - e9 - h3*g3

    gamma = 1/h2
    u1 = -a7
    u0 = -(a4 + h3*u1)
    v2 = f2
    v1 = f1 - f3*u1
    v0 = f0 - f3*u0
    w2 = -gamma*a2
    w1 = -(a9 - f3*w2)
    w0 = -(a6 + a3*w2 + h3*w1)
    
    # A is of type 31, semi-typical
    # Total : 1I 32M 1S 38A
    return C34CurveDivisor(D.C, [[u0, u1, 0, 1], [v0, v1, v2, 0, 1], [w0, w1, w2, 0, 0, 1]],
                           degree = 3, typ = 31, typical = False, reduced = True)

  else :
    c8f2 = c8*f2
    e0 = f1 - f2*f3
    e1 = - f2^2
    e2 = g1 - f3*g2
    e3 = c3 + c6*f2 - e2*(c8 - f3) - c7*f1 - g3*(c5 + c8f2 - f1) - c4*f3
    e4 = c6 + f2 - g3*(c8 - f3) - c7*f3
    e6 = e0 - f3*e4

    s0 = g0 - c2 - c5*(g2 - c5)
    s3 = g1 - c4 - c8*(g2 - c5 - c5)
    s4 = c6*c8 - g2 + c5
    s5 = - c6 + c7*c8
    s6 = g3 - c7 + c8^2

    a1 = g2 - c5 + c8f2
    a2 = s0 - s6*e1 - s3*f2
    a3 = s4 + e6 - c8*e4 - s6*g3 - s5*f3
    
    p0 = f2
    q2 = a3
    q0 = -(a2 + a1*a3)
    r0 = e0
    r1 = f3
    s0 = h3*(h3 - e4) + e3 - h1
    s1 = e4 - h3
    
    t0 = r0 + r1*(p0 - s1)
    t1 = r1*(r1*s1 + q2 - r0 - r0)
    u0 = s0
    u1 = s1
    v0 = t0*p0
    v1 = t0
    v2 = p0
    w0 = t1*p0 + q0
    w1 = t1
    w2 = q2
    
    # A is of type 31, semi-typical
    # Total : 0I 28M 2S 46A
    return C34CurveDivisor(D.C, [[u0, u1, 0, 1], [v0, v1, v2, 0, 1], [w0, w1, w2, 0, 0, 1]],
                           degree = 3, typ = 31, typical = False, reduced = True)



def flip_42(D) :
  K = D.K
  f, g, h = D.f, D.g, D.h
  c = D.C.coefficients()
  
  s0 =      - c[2] - c[4]*g[2]
  #s1 = g[0]        - c[5]*g[2]
  s2 =      - c[4] - c[7]*g[2]
  s3 = g[1] - c[5] - c[8]*g[2]
  
  a1 = s0 + c[7]*f[0] - f[1]*(s2 + c[7]*f[1])
  #a2 = s1 + c[8]*f[0]
  a3 = s3 + c[8]*f[1]
  
  u2 = -a3
  u0 = -(a1 + g[1]*u2)
  
  new_f = [f[1] - g[2], K.one()]
  new_g = [u0, K.zero(), u2, K.zero(), K.zero(), K.one()]
 
  # A is of type 22
  # Total 0I 8M 9A
  return C34CurveDivisor(D.C, [new_f, new_g, []])



def flip_43(D) :
  K = D.K
  f, g, h = D.f, D.g, D.h
  c = D.C.coefficients()
  
  s0 = c[5] + f[2]*(f[2] - c[7])
  
  a1 = - f[2]*g[4]
  a2 = g[2] - s0
  a3 = g[4] - c[8]
  
  u1 = -a3
  u0 = -(a2 + a1*u1)
  v1 = f[1] - f[2]*u1
  v0 = f[0] - f[2]*u0
  
  new_f = [u0, u1, K.one()]
  new_g = [v0, v1, K.zero(), K.one()]
  # A is of type 21
  # Total 0I 5M 7A
  return C34CurveDivisor(D.C, [new_f, new_g, []])



def flip_44(D) :
  return C34CurveDivisor(D.C, [[D.K.one()], [], []])



def flip_51(D) :
  K = D.K
  f, g, h = D.f, D.g, D.h
  c = D.C.coefficients()
  new_f, new_g, new_h = [], [], []
  
  if D.typical :
    r0 = c[5] - f[2]
    r1 = c[8] - f[4]
    s0 = c[6] - g[3] - f[3]*r1
    t0 = f[4]*r1 + f[3] - c[7] + g[4]

    # D is typical if and only if t0 =/= 0
    alpha = 1 / t0
    a1 = - f[3]*g[4]
    a2 = h[3]*t0
    a3 = - t0*(f[1] - f[3]*(g[3] + s0) + f[3]*h[4])
    a4 = alpha*(g[3] - f[4]*g[4])
    a5 = h[4]
    a6 = h[3] - f[2] - f[3]*(t0 - g[4]) - f[4]*h[4]
    a7 = - t0*f[4]

    u2 = - t0
    v2 = - a7
    u1 = -(a5 + a4*u2)
    v1 = -(a6 + a4*v2)
    u0 = -(a2 + a1*u2)
    v0 = -(a3 + a1*v2)
    
    # Change of basis
    u0, u1 = u1*s0 + u0, u1 + s0
    v0, v2 = v1*s0 + v0, v2 + s0
    
    # Compute third polynomial, h
    w2 = f[2] - f[3]*u2 - f[4]*v2 
    w1 = f[1] - f[3]*u1 - f[4]*v1 
    w0 = f[0] - f[3]*u0 - f[4]*v0    

    new_f = [u0, u1, u2, K.one()]
    new_g = [v0, v1, v2, K.zero(), K.one()]
    new_h = [w0, w1, w2, K.zero(), K.zero(), K.one()]
    
    # A is of type 31, typical
    # Total : 1I 24M 28A
  else : 
    r0 = c[5] - f[2]
    r1 = c[8] - f[4]
    s0 = c[6] - g[3] - f[3]*r1
    R0 = f[3]*r0 + h[2]
    R1 = f[3]*r1 + h[4]
    R2 = f[3]
    S0 = f[3]*(h[4] - s0 - g[3]) - f[4]*h[3] + f[1]
    T0 = f[2] - h[3] - f[3]*g[4]
    T1 = f[4]
    # Subtotal : 0I 6M 12A

    if S0 != 0 :
      l = g[4]*(h[4] - s0) - g[2]
      a2 = h[2] - g[1] - g[3]*s0 + g[4]*(T0 + T1*h[4]) - f[4]*l
      a3 = h[4] - g[3] - s0 + g[4]*T1

      alpha = - 1 / S0
      u1 = - a3
      u0 = -(a2 + h[4]*u1)
      v2 = s0
      v1 = T0 + T1*(s0 - u1)
      v0 = s0*T0 - T1*u0
      
      # Compute third polynomial, h
      w2 = f[2] - f[4]*v2  
      w1 = f[1] - f[3]*u1 - f[4]*v1 
      w0 = f[0] - f[3]*u0 - f[4]*v0    

      new_f = [u0, u1, K.zero(), K.one()]
      new_g = [v0, v1, v2, K.zero(), K.one()]
      new_h = [w0, w1, w2, K.zero(), K.zero(), K.one()]
      # A is of type 31, non-typical
      # Total : 1I 21M 31A
    else :
      # <f, g> represents a type 62 divisor and
      # <f, h> represents a type 63 divisor.
      # Compute the intersection of (f : g) and (f : h)
      D1 = C34CurveDivisor(D.C, [D.f, D.g, []])
      D2 = C34CurveDivisor(D.C, [D.f, D.h, []])
      A1 = flip_62(D1) # Costs 4M 6A
      A2 = flip_63(D2) # Costs 18M 28A
      # Subtotal :      0I 22M 34A
      # Running total : 0I 28M 46A
      p0 = A1.f[0]
      q0, q2 = A1.g[0], A1.g[2]
      r0, r1 = A2.f[0], A2.f[1]
      s0, s1 = A2.g[0], A2.g[1]
      t0 = r0 + r1*(p0 - s1)
      t1 = r1*(r1*s1 + q2 - r0 - r0)
      new_f = [s0, s1, K.zero(), K.one()]
      new_g = [t0*p0, t0, p0, K.zero(), K.one()]
      new_h = [q0 + t1*p0, t1, q2, K.zero(), K.zero(), K.one()]
      # A is of type 31, non-typical
      # Subtotal :      0I  5M  6A
      # Running total : 0I 33M 52A
  
  return C34CurveDivisor(D.C, [new_f, new_g, new_h])
  
  
  
def flip_52(D) :
  K = D.K
  f, g, h = D.f, D.g, D.h
  c = D.C.coefficients()
  
  #r3 = c[3]*f[3]
  #r4 = c[4]*f[3]
  r5 = c[5]*f[3]
  r6 = c[6]*f[3]
  r7 = c[7]*f[3]
  r8 = c[8]*f[3]
  r9 = f[1] - r6

  t0 =      - c[2]      - c[5]*r9
  #t1 = f[0]        - r3 - c[6]*r9
  #t2 =      - c[3] - r4 - c[7]*r9
  t3 =      - c[4] - r5 - c[8]*r9
  t4 =      - c[5]          -  r9
  t5 = f[2] - c[6] - r7
  t6 =      - c[7] - r8
  t7 =      - c[8] - f[3]

  s0  = g[3] # = -f[3]*f[3]
  s1  = - r5 + r8*f[2]
  #s2  = r9 + r8*s0 - f[3]*(f[2] - r7)
  s3  = -f[3]*s1 - f[2]^2
  #s4  = -f[3]*s2 - f[2]*s0 - f[1]*f[3]
  s5  = f[0] - f[3]*s3 - f[1]*f[2]
  #s6  = -f[3]*s4 - f[1]*s0
  s7  = f[2] + f[3]*s0
  s8  = - t0 + t7*s5 + t6*s3 + t5*s1 + t3*f[2]
  #s9  = - t1 + t7*s6 + t6*s4 + t5*s2 + t3*s0 + t2*f[3]
  s10 = - t4 + t7*s7 - t6*s0 - t5*f[3]

  a1 = g[2]
  a2 = g[0] - s8 - g[3]*s3 - f[2]*g[1]
  #a3 = - f[3]*g[3]
  #a4 = - s9 - g[3]*s4 - g[1]*s0
  a5 = g[2] - s10 + g[3]*s0

  u0 = f[2]
  v1 = -a5
  v0 = -(a2 + a1*v1)
  new_f = [ u0, K.one() ]
  new_g = [ v0, K.zero(), v1, K.zero(), K.zero(), K.one() ]
  
  # A has type 22
  # Total : 0I 22M 1SQ 27A
  return C34CurveDivisor(D.C, [new_f, new_g, []])



def flip_53(D) :
  K = D.K
  f, g, h = D.f, D.g, D.h
  c = D.C.coefficients()
  new_f, new_g, new_h = [], [], []

  s0 =      - c[3] - c[6]*f[2]
  s1 =      - c[4] - c[7]*f[2]
  s2 = f[1] - c[5] - c[8]*f[2]
  s3 =      - c[6] -      f[2]
  s5 = f[3] - c[8]
  
  t0 =      - c[6]*g[5]
  t1 = g[3] - c[7]*g[5]
  t2 =      - c[8]*g[5]
  
  j0 = -f[3]^2
  
  k0 =      - f[1]*f[3]
  k1 = f[1] - f[2]*f[3]
  k2 = f[2]
  k3 = j0
  
  l0 = s0
  l1 = s1 - f[1]*s5
  l2 = s2 - f[2]*s5
  l3 = s3
  l4 = -c[7] - f[3]*s5
  
  z0 = g[2] + l1 - g[5]*k1 - l4*f[2]
  z1 =        l2 - g[5]*k2
  
  a1 = g[3] - g[5]*j0
  a2 = t0 - f[1] - t2*j0 - (t1 - f[2])*f[3]
  a3 = g[1] + l0 - g[5]*k0 - l4*f[1] - z1*j0 - z0*f[3]
  a4 = -g[5] - f[3]
  a5 = g[3] + l3 - g[5]*k3 - l4*f[3]
  
  u1 = -a4
  v1 = -a5
  u0 = -(a2 + a1*u1)
  v0 = -(a3 + a1*v1)
  
  new_f = [u0, u1, K.one()]
  new_g = [v0, v1, K.zero(), K.one()]
  
  # A is of type 21
  # Total 0I 25M 1S 31A
  return C34CurveDivisor(D.C, [new_f, new_g, new_h])



def flip_54(D) :
  K = D.K
  f, g, h = D.f, D.g, D.h
  c = D.C.coefficients()

  s0 = f[1] - c[8]*f[2]
  s1 = -f[2]^2
  a1 = g[5] - s0
  a2 = g[4] - c[5] + s1 - c[8]*a1 + c[7]*f[2]

  new_f = [-a1, K.one()]
  new_g = [-a2, K.zero(), K.one()]
  
  # A is of type 21
  # Total 0I 3M 1SQ 6A
  return C34CurveDivisor(D.C, [new_f, new_g, []])



def flip_61(D) :
  K = D.K
  f, g, h = D.f, D.g, D.h
  c = D.C.coefficients()
  new_f, new_g, new_h = [], [], []
  
  if D.typical :
    s0 = f[3] - g[4] + f[5]*(f[4] - c[7])
    t0 = f[4] - g[5] + f[5]*(f[5] - c[8])
    
    z0 = c[7] - f[4]
    z1 = h[5] + c[6] - f[3]
    
    a1 = g[3] + g[5]*(f[4] - c[7])
    a2 = c[4] + c[7]*h[5] - f[2] - h[3] - z0*(g[4] + s0) - z1*f[4]
    a3 = g[4] + g[5]*(f[5] - c[8])
    a4 = c[5] + c[8]*h[5] - h[4] - z0*(g[5] + t0) - z1*f[5]
    
    u0 = t0*(h[4] - a1)
    u1 = a3 - h[5]
    u2 = t0
    v0 = t0*(a1*(c[8] - f[5]) - a2)
    v1 = a4 - a3*(c[8] - f[5])
    v2 = -t0*(c[8] - f[5])
    
    u0, u1 = u1*s0 + u0, u1 + s0
    v0, v2 = v1*s0 + v0, v2 + s0

    alpha = 1 / u2
    w0 = alpha*(v1*u0 + v0*(v2 - u1))
    w1 = alpha*(v1*v2 - v0)
    w2 = v1 + alpha*(u0 + v2*(v2 - u1))
    
    new_f = [u0, u1, u2, K.one()]
    new_g = [v0, v1, v2, K.zero(), K.one()]
    new_h = [w0, w1, w2, K.zero(), K.zero(), K.one()]
    # A is of type 31, typical
    # Total : 1I 24M 41A
    # XXX : Some additions can still be saved
  else :
    s0 = f[3] - g[4] + f[5]*(f[4] - c[7])
    
    T1 = c[8] - f[5]
    R2 = c[7] - f[4]
    K0 = f[5]*R2 + h[5]
    R1 = K0 + c[6] - f[3]
    T0 = R2*(c[8]*f[5] - f[4]) + f[5]*(h[5] - R1) + c[5] - h[4]
    S0 = f[3]*R2 + f[4]*R1 + f[2] + h[4]*T1 + h[3] - c[7]*K0 - c[4]
    # Subtotal : 9M 19A

    if S0 != 0 :
      t0 = g[5]*R2
      z0 = h[4] - g[3] + t0
      z1 = h[3] + g[5]*(c[6] + s0 - f[3])
      a2 = g[3] - t0
      a3 = - g[2] - s0*g[4] + g[5]*(c[5] + c[8]*s0 - z0) - z1*f[5]
      a4 = g[4] - g[5]*T1
      a5 = h[5] - s0 - a4
      
      u0 = h[5]*a5 - a3
      u1 = - a5
      v2 = s0
      w0 = S0*(h[5] - a4)
      w1 = - S0
      w2 = h[4] - a2
      
      # Change of basis
      v0 = T0*v2 - T1*u0
      v1 = T0 + T1*(v2 - u1)
      w0 = w0 + T0*w2 - T1*v0
      w1 = w1 + T1*(w2 - v1)
      w2 = w2 + T0 - T1*v2
      # Subtotal : 16M 26A
      # Total : 25M 45A
      return C34CurveDivisor(D.C, [[u0, u1, 0, 1], [v0, v1, v2, 0, 1], [w0, w1, w2, 0, 0, 1]],
                             degree = 3, typ = 31, typical = False, reduced = True)

    else :
      # Compute (f : g) = <x + p0, y^2 + q2*y + q0>
      z0 = g[4] - h[5] + g[5]*(T1 - c[8])
      z1 = - h[3] - T1*(h[4] - g[3]) + g[5]*(g[4] + s0 - c[6])
      z2 = g[2] + T1*(h[3] + g[5]*s0) - g[5]*(c[5] + g[3]) - z0*h[5]
      a1 = g[3] + g[5]*(f[4] - c[7])
      a2 = - h[1] - T0*h[4] + T1*(c[7]*h[3] - h[2] + g[1] + g[3]*s0) + g[5]*(g[2] + g[4]*s0 - c[3] - c[7]*g[3]) - z0*(h[3] + T1*h[4]) - z1*(g[4] + s0) + z2*(f[4] - c[7])
      a5 = g[3] - h[4] - T0 + T1*(g[4] + s0 - h[5] - z0) + g[5]*(g[5] - c[7])
      p0 = s0
      q0 = a1*a5 - a2
      q2 = - a5
      # Subtotal :      22M 39A
      # Running total : 31M 58A

      # Compute (f : h) = <y + r1*x + r0, x^2 + s1*x + s0>
      z0 = h[4] - g[3] + g[5]*(c[7] - f[4])
      z1 = h[3] + g[5]*(c[6] - f[3] + s0)
      a2 = - g[2] - g[4]*s0 + g[5]*(c[5] + c[8]*s0 - z0) - z1*f[5]
      a3 = h[5] - g[4] - s0 + g[5]*(c[8] - f[5])
      r0 = T0
      r1 = T1
      s0 = h[5]*a3 - a2
      s1 = - a3
      # Subtotal :       8M 16A
      # Running total : 39M 74A

      # Compute intersection of (f : g) and (f : h)
      t0 = r0 + r1*(p0 - s1)
      t1 = r1*(r1*s1 + q2 - r0 - r0)
      u0 = s0
      u1 = s1
      v0 = t0*p0
      v1 = t0
      v2 = p0
      w0 = q0 + t1*p0
      w1 = t1
      w2 = q2
      # Subtotal :       5M  6A
      # Running total : 44M 80A
      return C34CurveDivisor(D.C, [[u0, u1, 0, 1], [v0, v1, v2, 0, 1], [w0, w1, w2, 0, 0, 1]],
                             degree = 3, typ = 31, typical = False, reduced = True)
    
  
  return C34CurveDivisor(D.C, [new_f, new_g, new_h])



def flip_62(D) :
  K = D.K
  f, g, h = D.f, D.g, D.h
  c = D.C.coefficients()
  new_f, new_g, new_h = [], [], []
  
  u0 = c[6] + f[3]*(f[4] - c[8]) - g[3]
  v0 = (f[3]*u0 - f[1])*u0 + f[0]
  v1 = f[2] - f[4]*u0
  new_f = [u0, K.one()]
  new_g = [v0, K.zero(), v1, K.zero(), K.zero(), K.one()]
  # A is of type 22
  # Total : 0I 4M 6A
  return C34CurveDivisor(D.C, [new_f, new_g, []])



def flip_63(D) :
  K = D.K
  f, g, h = D.f, D.g, D.h
  c = D.C.coefficients()
  
  q0 = f[4] - c[8]
  q1 = f[2] - c[5]
  r0 = f[1]*q0
  r1 = f[2]*q0
  r2 = f[3]*q0
  r3 = f[4]*q0
  r5 = f[3]*q1
  r6 = f[4]*q1
  
  s0 = c[4] - f[1] + r1
  s1 = c[5] - f[2]
  s2 = c[3] + r0 + r5
  s3 = s0 + r6
  s4 = c[6] + r2
  s5 = c[7] + r3 - f[3]

  a1 = g[4] - g[6]*s5
  a2 = g[3] - f[2] + f[3]*s5 - f[4]*g[4]
  a4 = g[6] - f[4]
  a3 = g[2] - s2 - g[6]*s0 + s5*(f[2] - g[3] + g[6]*s4 - s5*f[3]) + f[4]*(s3 + g[6]*s1)
  a5 = a1 - s4 + s5*f[4]

  u1 = -a4
  v1 = -a5
  u0 = -(a2 + a1*u1)
  v0 = -(a3 + a1*v1)

  new_f = [u0, u1, K.one()]
  new_g = [v0, v1, K.zero(), K.one()]
  # A is of type 21
  # Total : 0I 18M 28A
  return C34CurveDivisor(D.C, [new_f, new_g, []])



def flip_64(D) :
  K = D.K
  f, g, h = D.f, D.g, D.h
  c = D.C.coefficients()
  
  s0 = c[6] + f[2] - f[3]*(c[7] - f[3]*(c[8] - f[3]))
  a1 = g[6] - s0
  a2 = -g[5] - f[1] + f[3]*(s0 + f[2] - g[6])
  new_f = [-a1, K.one()]
  new_g = [-a2, K.zero(), K.one()]
  # A is of type 11
  # Total : 0I 3M 9A
  return C34CurveDivisor(D.C, [new_f, new_g, []])



def flip_65(D) :
  return D.C.zero_divisor()

