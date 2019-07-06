"""
  Summary of costs for flipping a divisor
  
  
  Deg | Type |   I |   M |
  ----+------+-----+-----+
    0 |    0 |   0 |   0 |
  ----+------+-----+-----+
    1 |   11 |   0 |   4 |
  ----+------+-----+-----+
    2 |   21 |   0 |   7 | 
      |   22 |   0 |   1 | 
  ----+------+-----+-----+
    3 |   31 |   1 |  16?| Typical
      |   31 |   0 |  12 | Semi-typical
      |   32 |   0 |   3 | 
      |   33 |   0 |   0 | 
  ----+------+-----+-----+
    4 |   41 |   1 |  24 | Typical
      |   41 |   1 |  31 | Semi-typical (Can be improved)
      |   42 |   0 |  10 | 
      |   43 |   0 |   5 | 
      |   44 |   0 |   0 | 
  ----+------+-----+-----+
    5 |   51 |   1 |  24 | Typical
      |   51 |   1 |  21 | Semi-typical
      |   51 |   0 |  27 | Case where <f, g> =/= <f, g, h> =/= <f, h>
      |   52 |   0 |  23 | 
      |   53 |   0 |  26 | 
      |   54 |   0 |   3 | 
  ----+------+-----+-----+
    6 |   61 |   1 |  24 | Typical
      |   61 |     |     | Semi-typical, NOT IMPLEMENTED
      |   62 |   0 |   4 | 
      |   63 |   0 |  18 | 
      |   64 |   0 |   3 |
      |   65 |   0 |   0 | 
"""


def flip(D) :
  """
    Compute the 'flip' of a divisor.
    
    A typical (or semi-typical) divisor D is represented by two (or three) polynomials, f and g
    (or f, g, h). Then the divisor of f is
    
      div f = D + A
    
    The divisor A is the 'flip' of D.
    It is simply negation in the Jacobian of the curve.
    
    Input : A typical or semi-typical C34CrvDiv D.
    Output : A typical or semi-typical C34CrvDiv A, the flip of D.
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
  return C34CrvDiv(D.C, [[D.K.one()], [], []])



def flip_11(D) :
  K = D.K
  f, g = D.f, D.g
  c = D.C.coefficients()

  a = (c[4] - c[7]*f[0])*f[0] - c[2]
  b = c[8]*f[0] - c[5] + g[0]
  new_f = [f[0], K.one()]
  new_g = [b*g[0] - a, K.zero(), -b, K.zero(), K.zero(), K.one()]
  new_h = []
  # A is type 22
  # Total : 0I 4M
  return C34CrvDiv(D.C, [new_f, new_g, new_h])



def flip_21(D) :
  K = D.K
  f, g, h = D.f, D.g, D.h
  c = D.C.coefficients()
  a = c[7] + f[1]*(f[1] - c[8])
  b = g[0] - c[3] + f[0]*a + f[1]*(c[4] + f[0]*(f[1] - c[8]) + f[1]*(f[0] - c[5]))
  c = g[1] - c[6] + f[1]*a
  new_f = [f[0], f[1], K.one()]
  new_g = [c*g[1] - b, -c, K.zero(), K.one()]
  # A is type 21
  # Total : 0I 7M
  return C34CrvDiv(D.C, [new_f, new_g, []])



def flip_22(D) :
  K = D.K
  f, g, h = D.f, D.g, D.h
  c = D.C.coefficients()
  a = g[2] - c[5] + c[8]*f[0]
  new_f = [f[0], K.one()]
  new_g = [-a, K.zero(), K.one()]
  # A is type 11
  # Total : 0I 1M
  return C34CrvDiv(D.C, [new_f, new_g, []])



def flip_31(D) :
  K = D.K
  f, g, h = D.f, D.g, D.h
  c = D.C.coefficients()
  new_f, new_g, new_h = [], [], []
  
  if D.typical : 
    s0 = f[1] - g[2]
    
    z0 = c[7] - f[2]
    z1 = c[6] - f[1]
    a1 = h[1] - c[4] + z0*(g[2] + s0) + f[2]*z1
    a2 = h[2] - c[5] + f[2]*z0
    
    v0 = f[2]*(a1 + g[1]*c[8])
    v1 =  - (a2 + g[2]*c[8])
    v2 =  - f[2]*c[8]
    
    v0 = v0 + v1*s0
    v2 = v2 + s0

    u0, u1, u2 = f[0:3]

    alpha = 1 / u2
    v2mu1 = v2 - u1
    w0 = alpha*(u0*v1 + v0*v2mu1)
    w1 = alpha*(v1*v2 - v0)
    w2 = v1 + alpha*(u0 + v2*v2mu1)
    new_f = [u0, u1, u2, K.one()]
    new_g = [v0, v1, v2, K.zero(), K.one()]
    new_h = [w0, w1, w2, K.zero(), K.zero(), K.one() ]
    # A is of type 31, typical
    # Total : 1I 15M 18A
    # (In high characterstic, is is only 12M)
    
  else :
    s0 = f[1] - g[2]
    s2 = -g[1] + c[5] - c[8]*s0
    
    r0 = -c[8]*f[0]
    r1 = c[5] - c[8]*f[1]

    t0 = c[5] - h[2]
    t1 = c[8]
    t2 = r0
    t3 = r1 - h[2]
    t4 = c[2] - c[8]*r0 - c[7]*f[0] - h[0] + c[5]*(h[2] - c[5])
    t5 = c[4] - c[7]*f[1] - h[1] + c[8]*(h[2] - r1 - c[5])
    
    u0 = s0*t0 + t2
    u1 = s0*t1 + t3
    u2 = s0
    v0 = s2*t0 + t4
    v1 = s2*t1 + t5
    v2 = s2

    new_f = [f[0], f[1], K.zero(), K.one() ]
    new_g = [u0, u1, u2, K.zero(), K.one() ]
    new_h = [v0, v1, v2, K.zero(), K.zero(), K.one() ]
    # A is type 31, non-typical
    # Total: 0I 12M
  
  return C34CrvDiv(D.C, [new_f, new_g, new_h])



def flip_32(D) :
  K = D.K
  f, g, h = D.f, D.g, D.h
  c = D.C.coefficients()
  
  a1 = -f[1]*(f[1] - c[8])
  a2 = g[3] + f[1]*(c[7] - a1)
  a3 = -f[0] + f[1]*(c[6] - a2)
  new_f = [c[6] - a2, K.one()]
  new_g = [-a3, K.zero(), K.one()]
  # A is of type 11
  # Total 0I 3M
  return C34CrvDiv(D.C, [new_f, new_g, []])



def flip_33(D) :
  # A is of type 0
  # Total 0I 0M
  return C34CrvDiv(D.C, [[D.K.one()], [], []])



def flip_41(D) :
  K = D.K
  f, g, h = D.f, D.g, D.h
  c = D.C.coefficients()
  new_f, new_g, new_h = [], [], []

  #f3f3 = f[3]^2
  #alpha = g[3] + f3f3
  if D.typical :
    s0 = f[2]
    t0 = - f[3]*f[3] - g[3]
    
    c7f3 = c[7]*f[3]
    z0 = c[8]*f[3]
    z1 = h[3] - f[2] + c7f3
    a1 = c[5] - f[2]*c[8]
    a2 = h[2] + c[5]*f[3] - f[2]*z0
    a3 = c[6] - h[3] - c[8]*(g[3] + t0) - c7f3
    a4 = f[1] + f[3]*(z1 - c[6]) + z0*(g[3] + t0)
    
    f3g3 = f[3]*g[3]
    u0 = t0*(g[2] - a1)
    u1 = a3 + f3g3
    u2 = - t0
    v0 = t0*(a2 - g[2]*f[3])
    v1 = a4 - f[3]*f3g3
    v2 = t0*f[3]
    
    u0 = u0 + u1*s0
    u1 = u1 + s0
    
    v0 = v0 + v1*s0
    v2 = v2 + s0

    alpha = 1 / u2
    v2mu1 = v2 - u1
    w0 = alpha*(u0*v1 + v0*v2mu1)
    w1 = alpha*(v1*v2 - v0)
    w2 = v1 + alpha*(u0 + v2*v2mu1)

    new_f = [u0, u1, u2, K.one()]
    new_g = [v0, v1, v2, K.zero(), K.one()]
    new_h = [w0, w1, w2, K.zero(), K.zero(), K.one() ]
    # A is of type 41, typical
    # Total : 1I 24M 26A
    # (In high characteristic, is is only 18M)
  elif h[2] != 0 :
    e0 = f[1] - f[2]*f[3]
    e1 = - f[2]*f[2]
    e2 = g[1] - f[3]*g[2]
    e3 = c[3] + c[6]*f[2] - e2*(c[8] - f[3]) - c[7]*f[1] - g[3]*(c[5] + c[8]*f[2] - f[1]) - c[4]*f[3]
    e4 = c[6] + f[2] - g[3]*(c[8] - f[3]) - c[7]*f[3]
    e7 = - f[2]*e1
    e8 = - g[3]*e3 - g[2]*e0 + g[0]
    e9 = - g[3]*e4 + e2
    
    a2 = h[0] - c[5]*h[2] - e7 - h[3]*e1 - f[2]*(h[1] - c[8]*h[2])
    a3 = - e0 - h[3]*f[3]
    a4 = h[1] - e3 - h[2]*f[3]
    a6 = - c[6]*h[2] - e8 - h[3]*e2 - g[3]*(h[1] - c[8]*h[2]) + c[7]*h[2]*f[3]
    a7 = h[3] - e4
    a9 = - h[2] - e9 - h[3]*g[3]

    gamma = 1/h[2]
    u1 = -a7
    u0 = -(a4 + h[3]*u1)
    v2 = f[2]
    v1 = f[1] - f[3]*u1
    v0 = f[0] - f[3]*u0
    w2 = -gamma*a2
    w1 = -(a9 - f[3]*w2)
    w0 = -(a6 + a3*w2 + h[3]*w1)
    
    new_f = [u0, u1, K.zero(), K.one()]
    new_g = [v0, v1, v2, K.zero(), K.one()]
    new_h = [w0, w1, w2, K.zero(), K.zero(), K.one()]
    # A is of type 31, semi-typical
    # Total : 1I 31M
  else :
    e0 = f[1] - f[2]*f[3]
    e1 = - f[2]*f[2]
    e2 = g[1] - f[3]*g[2]
    e3 = c[3] + c[6]*f[2] - e2*(c[8] - f[3]) - c[7]*f[1] - g[3]*(c[5] + c[8]*f[2] - f[1]) - c[4]*f[3]
    e4 = c[6] + f[2] - g[3]*(c[8] - f[3]) - c[7]*f[3]
    e6 = e0 - f[3]*e4
    #e7 = - f[2]*e1
    #e8 = - g[3]*e3 - g[2]*e0 + g[0]
    #e9 = - g[3]*e4 + e2

    s0 = g[0] - c[2] - c[5]*(g[2] - c[5])
    #s1 = c[3]*c[8] - c[6]*(g[2] - c[5])
    #s2 = - c[3] + c[4]*c[8] - c[7]*(g[2] - c[5])
    s3 = g[1] - c[4] + c[5]*c[8] - c[8]*(g[2] - c[5])
    s4 = c[6]*c[8] - g[2] + c[5]
    s5 = - c[6] + c[7]*c[8]
    s6 = g[3] - c[7] + c[8]*c[8]

    a1 = g[2] - c[5] + c[8]*f[2]
    a2 = s0 - s6*e1 - s3*f[2]
    a3 = s4 + e6 - c[8]*e4 - s6*g[3] - s5*f[3]
    
    p0 = f[2]
    q2 = a3
    q0 = -(a2 + a1*a3)
    r0 = e0
    r1 = f[3]
    s0 = h[3]*(h[3] - e4) + e3 - h[1]
    s1 = e4 - h[3]
    
    #FG = C34CrvDiv(D.C, [[p0, 1], [q0, 0, q2, 0, 0, 1], []])
    #FH = C34CrvDiv(D.C, [[r0, r1, 1], [s0, s1, 0, 1], []])
    #print "(f : g) = "
    #print FG
    #print "(f : h) = "
    #print FH
    t0 = r0 + r1*(p0 - s1)
    t1 = r1*(r1*s1 + q2 - 2*r0)
    u0 = s0
    u1 = s1
    v0 = t0*p0
    v1 = t0
    v2 = p0
    w0 = t1*p0 + q0
    w1 = t1
    w2 = q2
    
    new_f = [u0, u1, K.zero(), K.one()]
    new_g = [v0, v1, v2, K.zero(), K.one()]
    new_h = [w0, w1, w2, K.zero(), K.zero(), K.one()]
    # A is of type 31, semi-typical
    # Total : 0I 31M 1C

  return C34CrvDiv(D.C, [new_f, new_g, new_h])



def flip_42(D) :
  K = D.K
  f, g, h = D.f, D.g, D.h
  c = D.C.coefficients()
  
  s0 =      - c[2] - c[4]*g[2]
  s1 = g[0]        - c[5]*g[2]
  s2 =      - c[4] - c[7]*g[2]
  s3 = g[1] - c[5] - c[8]*g[2]
  
  a1 = s0 + c[7]*f[0] - f[1]*(s2 + c[7]*f[1])
  a2 = s1 + c[8]*f[0]
  a3 = s3 + c[8]*f[1]
  
  u1 = -a3
  u0 = -(a1 + g[1]*u1)
  
  new_f = [f[1] - g[2], K.one()]
  new_g = [u0, K.zero(), u1, K.zero(), K.zero(), K.one()]
 
  # A is of type 22
  # Total 0I 10M
  return C34CrvDiv(D.C, [new_f, new_g, []])



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
  # Total 0I 5M
  return C34CrvDiv(D.C, [new_f, new_g, []])



def flip_44(D) :
  return C34CrvDiv(D.C, [[D.K.one()], [], []])



def flip_4(D) :
  K = D.K
  f, g, h = D.f, D.g, D.h
  c = D.C.coefficients()
  new_f, new_g, new_h = [], [], []
  T = (len(f), len(g), len(h))
  
  if T == (5, 6, 7) : # Most typical case, type 41
    alpha = g[3] + f[3]*f[3]
    if alpha == 0 :
      raise ValueError("Matrix does not have full rank")
    a1 = g[3] - c[7] + c[8]*f[3] 
    a2 = g[2] - c[5] + c[8]*f[2] 
    a3 = - c[6] - a1*f[3] 
    a4 = g[1] + f[1]*f[3] 
    a5 = g[2] - f[1] + 2*f[2]*f[3] 
    a6 = alpha 
    a7 = f[2]*f[2] 
    a8 = a4 - a5*f[3] 
    a9 = a3 - f[2] + f[3]*(a6 - g[3]) 
    a10 = g[2] - c[5] + f[1] - f[2]*f[3] 
    a11 = g[1] - c[4] + f[1]*(c[8] - f[3]) + f[2]*(a6 - a1 - g[3] - c[7]) - f[3]*a10 
    a12 = f[0] - f[2]*(a10 + c[5]) 
    a13 = -c[3] - c[6]*f[2] + f[1]*(a6 - a1 - g[3]) - f[3]*a11

    a9  = f[2] - f[3]*a6
    a10 = g[2] - f[1] + f[2]*(c[8] + f[3])
    a11 = a4 + f[2]*(c[7] - a6) - f[3]*a10
    a12 = f[2]*c[6] - f[1]*a6 - f[3]*a11
    a13 = - f[0] + f[2]*(c[5] - a10)
    
    #M = Matrix(K, 3, 5, [
    #  [1, -f[2], a2, a7, a13],
    #  [0, a6, a3, a8, a12],
    #  [0, 0, -1, a6, a9]])
    #print M
    #print

    alpha = 1/alpha 
    b1 = alpha*a8 + a3 # = alpha*(a8 + a3*a6) 
    b2 = alpha*(a12 + a3*a9) 
    v0 = -(a7 + f[2]*b1 + a2*a6) 
    v1 = -b1 
    v2 = a6 
    w0 = -(a13 + f[2]*b2 + a2*a9) 
    w1 = -b2 
    w2 = a9
    new_f = [v0, v1, v2, K.one()]
    new_g = [w0, w1, w2, K.zero(), K.one()]
    # Total : 1I 24M
    
    # TODO: Compute h polynomial :
  elif (T == (4, 5, 0)):
    # In this case we have
    # f = x^2 + f1*x + f0
    # g = xy + g2*y + g1*x + g0
    # XXX : Not sure that this is true in general, but for now assume that
    # f = (x - x0)^2
    # g = (x - x0)*(y - y0)
    # I.e. D = 2*(x0 : y0 : 1) + (x0 : y1 : 1) + (x0 : y2 : 1) for some x0, y0, y1, y2.
    # Assuming also y0, y1, y2 are distinct, then D reduced is
    # D = (x0 : y0 : 1)
    # Moreover, xy + g2*y + g1*x + g0 = (x - x0)*(y - y0) = xy - x0*y - y0*x - x0*y0
    # implies that x0 = -g2 and y0 = -g1.
    # Under these assumptions, this quick hack...
    return flip(D.C.divisor([D.C.point(-g[2], -g[1])]))
    
  else :
    raise NotImplementedError("Flipping of degree 4 divisors only implemented for typical divisors of type <xy, y^2>.\nD = {}".format(D))
  return C34CrvDiv(D.C, [new_f, new_g, new_h])



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
      D1 = C34CrvDiv(D.C, [D.f, D.g, []])
      D2 = C34CrvDiv(D.C, [D.f, D.h, []])
      A1 = flip_62(D1) # Costs 4M
      A2 = flip_63(D2) # Costs 18M
      p0 = A1.f[0]
      q0, q2 = A1.g[0], A1.g[2]
      r0, r1 = A2.f[0], A2.f[1]
      s0, s1 = A2.g[0], A2.g[1]
      t0 = r0 + r1*(p0 - s1)
      t1 = r1*(r1*s1 + q2 - 2*r0)
      new_f = [s0, s1, K.zero(), K.one()]
      new_g = [t0*p0, t0, p0, K.zero(), K.one()]
      new_h = [q0 + t1*p0, t1, q2, K.zero(), K.zero(), K.one()]
      # A is of type 31, non-typical
      # Total : 0I 27M 1CM ??A
  
  return C34CrvDiv(D.C, [new_f, new_g, new_h])
  
  
  
def old_flip_51(D) :
  K = D.K
  f, g, h = D.f, D.g, D.h
  c = D.C.coefficients()
  new_f, new_g, new_h = [], [], []
  
  if D.typical :
    d0 = c[3] + f[1]*(f[4] - c[8]) + f[3]*(f[2] - c[5])
    d1 = c[4] - f[1] + f[2]*(f[4] - c[8]) + f[4]*(f[2] - c[5])
    e0 = c[6] + f[3]*(f[4] - c[8])
    e1 = c[7] - f[3] + f[4]*(f[4] - c[8])
    e2 = f[1] - f[3]*e0
    e3 = f[2] - f[3]*e1
    e4 = d0 - e0*e0
    e5 = d1 - e0*e1
    e6 = - e1*e2 - d1*f[3]
    e7 = d0 - e1*e3 - d1*f[4]
    e8 = e0 - e1*f[4]
    
    a1 = g[3] - e0
    a2 = - g[4]*f[3]
    a3 = g[1] - e4 - g[3]*e0
    a4 = - e6 - g[4]*e2 - g[2]*f[3]
    a5 = g[4] - e1
    a6 = g[3] - g[4]*f[4]
    a7 = g[2] - e5 - g[3]*e1
    a8 = g[1] - e7 - g[4]*e3 - g[2]*f[4]
    a9 = g[3] - e8 - g[4]*f[4]
    
    alpha = 1/a5
    u2 = -a5
    v2 = -a9
    u1 = -alpha*(a7 + a6*u2)
    v1 = -alpha*(a8 + a6*v2)
    u0 = -(a3 + a2*u2 + a1*u1)
    v0 = -(a4 + a2*v2 + a1*v1)
    w0 = alpha*(v0*(u1 - v2) - u0*v1)
    w1 = alpha*(v0 - v1*v2)
    w2 = v1 + alpha*(v2*(u1 - v2) - u0)
    
    new_f = [u0, u1, u2, K.one()]
    new_g = [v0, v1, v2, K.zero(), K.one()]
    new_h = [w0, w1, w2, K.zero(), K.zero(), K.one()]
    
    # A is of type 31, typical
    # Total : 1I 39M 58A
  elif False : # Ignore this old code for now
    q1 = f[1] - c[4] 
    q2 = f[2] - c[5] 
    q3 = f[3] - c[7] 
    q4 = f[4] - c[8] 
    
    s0 = f[1] - c[6]*f[3] 
    s1 = f[2] + f[3]*q3 
    s2 = f[3]*q4 
    s3 = f[0] - c[6]*s0 - c[3]*f[3] 
    s4 = - c[7]*s0 + f[3]*q1 - c[3]*f[4] 
    s5 = - c[8]*s0 + f[3]*q2 + f[4]*q1 
    s6 =  - s0 + f[4]*q2 
    s7 = s1 - c[6]*f[4] 
    s8 = s2 + f[4]*q3 
    s9 = f[4]*q4
    
    t0 = h[4] - s2 
    t1 = q3 - s9 
    t2 = q2 - h[3] 
    t3 = q1 - c[8]*h[3] - f[2]*q4 - s2*t1 - f[4]*t2 
    t4 = h[4] - s8 + f[4]*s9 
    t5 = h[2] - s5 + f[2]*s9 - s2*t4 + f[4]*s6 
    
    a1 = - s0 - f[3]*t0 
    a2 = h[1] - c[6]*h[3] - s0*t1 - f[3]*t3 
    a3 = - s3 - s0*t4 - f[3]*t5 
    a4 = h[3] - s1 - f[4]*t0 
    a5 = h[2] - c[3] - c[7]*h[3] - f[1]*q4 - s1*t1 - f[3]*t2 - f[4]*t3 
    a6 = h[1] - s4 + f[1]*s9 - s1*t4 + f[3]*s6 - f[4]*t5 
    a7 = t0 - c[6] - f[4]*t1 
    a8 = h[3] - s7 + f[3]*s9 - f[4]*t4 
    
    assert a1 + f[4]*h[3] != 0, "Flipping of super-non-typical divisors not implemented. D = {}".format(D)    
    alpha = 1 / (a1 + f[4]*h[3]) 
    u2 = -alpha*(a2 - a7*h[3]) 
    v2 = -alpha*(a3 - a8*h[3]) 
    u1 = -(a7 + f[4]*u2) 
    v1 = -(a8 + f[4]*v2) 
    u0 = -(a5 + h[4]*u1 + a4*u2) 
    v0 = -(a6 + h[4]*v1 + a4*v2) 
    w2 = f[2] - f[3]*u2 - f[4]*v2 
    w1 = f[1] - f[3]*u1 - f[4]*v1 
    w0 = f[0] - f[3]*u0 - f[4]*v0    
    new_f = [u0, u1, u2, K.one()]
    new_g = [v0, v1, v2, K.zero(), K.one()]
    new_h = [w0, w1, w2, K.zero(), K.zero(), K.one()]
  else :
    d0 = c[3] + f[1]*(f[4] - c[8]) + f[3]*(f[2] - c[5])
    d1 = c[4] - f[1] + f[2]*(f[4] - c[8]) + f[4]*(f[2] - c[5])
    e0 = c[6] + f[3]*(f[4] - c[8])
    e1 = c[7] - f[3] + f[4]*(f[4] - c[8])
    e2 = f[1] - f[3]*e0
    e3 = f[2] - f[3]*e1
    e4 = d0 - e0*e0
    e5 = d1 - e0*e1
    e6 = - e1*e2 - d1*f[3]
    e7 = d0 - e1*e3 - d1*f[4]
    e8 = e0 - e1*f[4]
    e9  = f[0] - f[4]*e6 - f[3]*e4 - f[1]*e0
    e10 = - f[4]*e7 - f[3]*e5 - f[1]*e1
    e11 = f[2] - f[4]*e8 - f[3]*e1

    a1 = - f[3]*h[4]
    a2 = h[1] - e6 - h[3]*e0
    a3 = - e9 - h[4]*e2 - f[3]*h[2]
    a4 = h[3] - f[4]*h[4]
    a5 = h[2] - e7 - h[3]*e1
    a6 = h[1] - e10 - h[4]*e3 - f[4]*h[2]
    a7 = h[4] - e8
    a8 = h[3] - e11 - h[4]*f[4]

    assert a1 + f[4]*h[3] != 0, "Flipping of super-non-typical divisors not implemented. D = {}".format(D)    
    alpha = 1 / (a1 + f[4]*h[3])
    u2 = -alpha*(a2 - h[3]*a7)
    v2 = -alpha*(a3 - h[3]*a8)
    u1 = -(a7 - f[4]*u2)
    v1 = -(a8 - f[4]*v2)
    u0 = -(a5 + h[4]*u1 + a4*u2)
    v0 = -(a6 + h[4]*v1 + a4*v2)
    w2 = f[2] - f[3]*u2 - f[4]*v2 
    w1 = f[1] - f[3]*u1 - f[4]*v1 
    w0 = f[0] - f[3]*u0 - f[4]*v0    
    new_f = [u0, u1, u2, K.one()]
    new_g = [v0, v1, v2, K.zero(), K.one()]
    new_h = [w0, w1, w2, K.zero(), K.zero(), K.one()]
    # A is of type 31, atypical(?)
    # Total : 1I 49M
  
  return C34CrvDiv(D.C, [new_f, new_g, new_h])



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
  s3  = -f[3]*s1 - f[2]*f[2]
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
  # Total : 0I 23M
  return C34CrvDiv(D.C, [new_f, new_g, []])



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
  
  j0 = -f[3]*f[3]
  
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
  # Total 0I 26M
  return C34CrvDiv(D.C, [new_f, new_g, new_h])



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
  return C34CrvDiv(D.C, [new_f, new_g, []])



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

    if S0 != 0 :
      z0 = h[4] - g[3] + g[5]*R2
      z1 = h[3] + g[5]*(c[6] + s0 - f[3])
      a2 = g[3] - g[5]*R2
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

      new_f = [u0, u1, K.zero(), K.one()]
      new_g = [v0, v1, v2, K.zero(), K.one()]
      new_h = [w0, w1, w2, K.zero(), K.zero(), K.one()]
      # A is of type 31, non-typical
      # Total : 0I 25M
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

      # Compute (f : h) = <y + r1*x + r0, x^2 + s1*x + s0>
      z0 = h[4] - g[3] + g[5]*(c[7] - f[4])
      z1 = h[3] + g[5]*(c[6] - f[3] + s0)
      a2 = - g[2] - g[4]*s0 + g[5]*(c[5] + c[8]*s0 - z0) - z1*f[5]
      a3 = h[5] - g[4] - s0 + g[5]*(c[8] - f[5])
      r0 = T0
      r1 = T1
      s0 = h[5]*a3 - a2
      s1 = - a3

      # Compute intersection of (f : g) and (f : h)
      t0 = r0 + r1*(p0 - s1)
      t1 = r1*(r1*s1 + q2 - 2*r0)
      new_f = [s0, s1, K.zero(), K.one()]
      new_g = [t0*p0, t0, p0, K.zero(), K.one()]
      new_h = [q0 + t1*p0, t1, q2, K.zero(), K.zero(), K.one()]
  
  return C34CrvDiv(D.C, [new_f, new_g, new_h])



def old_flip_61(D) :
  K = D.K
  f, g, h = D.f, D.g, D.h
  c = D.C.coefficients()
  new_f, new_g, new_h = [], [], []
  
  if D.typical :
    d1 = c[3] - f[1]
    d2 = c[4] - f[2]
    d3 = c[6] - f[3]
    d4 = c[7] - f[4]
    d5 = c[8] - f[5]
    e1 = -f[5]*d1
    e2 = f[1] - f[5]*d2
    e3 = f[2] - f[5]*c[5]
    e4 = -f[5]*d3
    e5 = f[3] - f[5]*d4
    e6 = f[4] - f[5]*d5
    e7 = d2 - f[4]*d3
    e8 = c[5] + e4
    
    a1 = g[4] - e5
    a5 = g[5] - e6
    a2 = g[3] - g[5]*d4
    a6 = g[4] - g[5]*d5
    a3 = g[2] - f[1] + f[5]*e7 - a1*e5 - f[4]*g[3]
    a7 = - f[2] - a1*e6 - f[5]*(g[3] - e8)
    t  = g[3] - e4 - a5*d4
    a4 = g[1] - e1 - a5*e7 - t*e5 + e3*d4
    a8 = g[2] - e2 - t*e6 + e3*d5
    a9 = a1 - a5*d5

    assert a5 != 0, "a5 should not be 0 if D is typical. D = {}".format(D)
    # Subtotal : 21M 32A

    alpha = 1 / a5
    u2 = -a5
    u1 = -(alpha*a7 - a6)
    u0 = -(a3 + a2*u2 + a1*u1)
    v2 = -a9
    v1 = -(alpha*(a8 + a6*v2) - e8)
    v0 = -(a4 + a2*v2 + a1*v1)
    w2 = v1 - alpha*(v2*(v2 - u1) + u0)
    w1 = -alpha*(v1*v2 - v0)
    w0 = -alpha*(u0*v1 + v0*(v2 - u1))
    new_f = [u0, u1, u2, K.one()]
    new_g = [v0, v1, v2, K.zero(), K.one()]
    new_h = [w0, w1, w2, K.zero(), K.zero(), K.one()]
    # Total : 1I 35M 45A
  else :
    raise NotImplementedError("Flipping of atypical type 61 divisors not implemented. D = {}".format(D))
  return C34CrvDiv(D.C, [new_f, new_g, new_h])

def old_old_flip_61(D) :
  K = D.K
  f, g, h = D.f, D.g, D.h
  c = D.C.coefficients()
  new_f, new_g, new_h = [], [], []
  
  # If D is typical and of the form <x^3, x^2y>

  F31 = f[2] - c[4] + c[7]*(c[6] - f[3])
  F32 = - c[5] + c[8]*(c[6] - f[3])
  
  a1 = f[3] + f[5]*(f[4] - c[7])
  a2 = f[4] + f[5]*(f[5] - c[8])
  
  b1 = f[1] - c[3] + f[5]*F31
  b2 = f[2] - c[4] + f[5]*F32
  b3 =      - c[5] + f[5]*(c[6] - f[3])
  b4 = a1 - c[6]
  b5 = a2 - c[7]
  
  c1 =        f[5]*(f[1] - c[3])
  c2 = f[1] + f[5]*(f[2] - c[4])
  c3 = f[2] - f[5]*c[5]
  c4 = - b3 - c[5]

  j1 = g[2] - c[7]*g[3] - c[3] - b1
  j2 =      - c[8]*g[3] - c[4] - b2
  j3 =           - g[3] - c[5] - b3
  j4 =             g[4] - c[6] - b4
  j5 =             g[5] - c[7] - b5

  k1 = g[1] - c1
  k2 = g[2] - c2
  k3 =      - c3
  k4 = g[3] - c4
  k5 = g[4] - a1
  k6 = g[5] - a2
  
  l1 = k1 + k6*(F31)
  l2 = k2 + k6*(F32)
  l3 = k3 + k6*(c[6] - f[3])
  l4 = k4 + k6*(f[4] - c[7])
  l5 = k5 + k6*(f[5] - c[8])

  # One multiplication can be saved when computing m3, m4, m7, m8
  # with Strassen's technique
  m1 = g[4] - a1
  m2 = g[3] + g[5]*(f[4] - c[7])
  m3 = j1 - j4*a1 + j3*(f[4] - c[7])
  m4 = l1 - l4*a1 + l3*(f[4] - c[7])
  m5 = g[5] - a2
  m6 = g[4] + g[5]*(f[5] - c[8])
  m7 = j2 - j4*a2 + j3*(f[5] - c[8])
  m8 = l2 - l4*a2 + l3*(f[5] - c[8])

  if m5 == 0 :
    raise ValueError("Leftmost 3x3 submatrix is not invertible.")

  mu = 1/m5
  v2 = -j5
  v1 = -mu*(m7 + m6*v2)
  v0 = -(m3 + m1*v1 + m2*v2)
  w2 = -l5
  w1 = -mu*(m8 + m6*w2)
  w0 = -(m4 + m1*w1 + m2*w2)
  
  new_f = [v0, v1, v2, K.one()]
  new_g = [w0, w1, w2, K.zero(), K.one()]
  # Total 1I 35M

  return C34CrvDiv(D.C, [new_f, new_g, new_h])



def flip_62(D) :
  K = D.K
  f, g, h = D.f, D.g, D.h
  c = D.C.coefficients()
  new_f, new_g, new_h = [], [], []
  
  # TODO : Not sure if this is correct at all!
  u0 = c[6] + f[3]*(f[4] - c[8]) - g[3]
  v0 = (f[3]*u0 - f[1])*u0 + f[0]
  v1 = f[2] - f[4]*u0
  new_f = [u0, K.one()]
  new_g = [v0, K.zero(), v1, K.zero(), K.zero(), K.one()]
  # A is of type 22
  # Total : 0I 4M
  return C34CrvDiv(D.C, [new_f, new_g, []])



def old_flip_63(D) :
  # XXX : This passes unit tests
  K = D.K
  f, g, h = D.f, D.g, D.h
  c = D.C.coefficients()

  s0 = f[2] + f[3]*(f[3] - c[7])
  s1 = f[3]*(f[4] - c[8])
  s2 = g[3] - c[6]*g[6]
  
  t0 = f[4] - c[8] - g[6]
  t1 = f[3] - c[7] - c[8]*g[6] - f[4]*t0
  t2 = f[2] - c[5] - s2
  t3 = f[1] - c[4] - c[5]*g[6] - c[8]*s2 - f[2]*t0 - s1*t1 - f[4]*t2
  
  a1 = g[4] + g[6]*(f[3] - c[7] - f[4]*(f[4] - c[8]))
  a2 = g[3] - s0 - f[4]*(g[4] - s1)
  a4 = g[6] - f[4]
  a3 = g[2] - c[3] - c[4]*g[6] - c[7]*s2 - f[1]*t0 - s0*t1 - f[3]*t2 - f[4]*t3
  a5 = g[4] - c[6] - c[7]*g[6] - f[3]*t0 - f[4]*t1
  
  u1 = - a4
  v1 = - a5
  u0 = -(a2 + a1*u1)
  v0 = -(a3 + a1*v1)

  new_f = [u0, u1, K.one()]
  new_g = [v0, v1, K.zero(), K.one()]
  # A is of type 21
  # Total : 0I 24M
  
  return C34CrvDiv(D.C, [new_f, new_g, []])

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
  # Total : 0I 18M
  return C34CrvDiv(D.C, [new_f, new_g, []])



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
  return C34CrvDiv(D.C, [new_f, new_g, []])



def flip_65(D) :
  return C34CrvDiv(D.C, [[D.K.one()], [], []])

