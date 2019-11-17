"""
  Summary of costs for doubling a divisor
  
  
  deg |   I |   M |
  ----+-----+-----+
    0 |   0 |   0 |
  ----+-----+-----+
    1 |   1 |  16 | When dy/dx =/= 0
      |   0 |   9 | 
  ----+-----+-----+
    2 |   1 |  63 | Distinct x-coordinates
      |   1 |  24 | Same x-coordinate (this can be brought down to 17M)
  ----+-----+-----+
    3 |   3 | 136 | Typical
    3 |   5 | 177 | Semi-typical
"""

def double(D) :
  """
    Double a divisor D.
    Assumes D is a typical divisor of degree 3.
    
    Input: A typical C34CurveDivisor D.
    Output: The C34CurveDivisor E equivalent to D + D.
            May be typical or semi-typical (or neither?).
  """
  if not D.reduced :
    raise ValueError("Divisor must be reduced. D = {}".format(D))
  
  if D.type == 0 :
    DD = double_0(D)
  elif D.type == 11 :
    DD = double_11(D)
  elif D.type == 21 :
    DD = double_21(D)
  elif D.type == 22 :
    DD = double_22(D)
  elif D.type == 31 :
    DD = double_31(D)
  else :
    raise ValueError("Divisor is of unexpected type. D = {}".format(D))
  if DD.reduced :
    return DD
  return flip(flip(DD))



def thesis_double_31(D) :
  C = D.C
  f0, f1, f2 = D.f[0:3]
  g0, g1, g2 = D.g[0:3]
  h0, h1, h2 = D.h[0:3]
  c0, c1, c2, c3, c4, c5, c6, c7, c8 = C.coefficients()
  print_matrices = false

  if (D.inv == 0) :
    D.inv = 1/f2
  f2_inv = D.inv
  
  if (D.type != 31) :
    raise ValueError("Divisor is not of type 31.\nD = {}".format(D))

  if (D.typical) :
    # Compute polynomials G and H of the forms
    #
    #   G = xy + G2*y + G1*x + G0
    #   H = y^2 + H3*x^2 + H2*y + H1*x + H0
    #
    # satisfying fH = gG (mod C)
    H3 = f2
    G2 = -c8*f2 + f1 - g2
    H2 = c5 + (G2*g2 - f0)*f2_inv 
    G1 = f2*(f2 - c7) + H2 - g1
    H1 = f2*(c6 - f1)
    G0 = H2*f1 + H1*f2 - c4*f2 - G2*g1 - G1*g2 - g0
    H0 = -f0*f2 - H1*f1 + c3*f2 + G1*g1

    a1  = G0 - g0
    a6  = G1 - g1
    a11 = G2 - g2
    a2  = H0 - h0 - f0*f2
    a7  = H1 - h1 - f1*f2
    a12 = H2 - h2 - f2*f2
    a4  =    - f0*a6 - g0*a11
    a9  = a1 - f1*a6 - g1*a11
    a14 =    - f2*a6 - g2*a11
    a5  =    - f0*a7 - g0*a12
    a10 = a2 - f1*a7 - g1*a12
    a15 =    - f2*a7 - g2*a12
    a3  = (-g0*a6 - h0*a11 + g1*a1 - a5  - (f1 - g2)*a2)*f2_inv
    a8  = (       - h1*a11         - a10 - (f1 - g2)*a7)*f2_inv
    a13 = ( a1 - g2*a6 - (h2 - g1)*a11 - a15 - (f1 - g2)*a12)*f2_inv

  elif (something == 0) :
    raise NotImplementedError("Do this later")
  else :
    # Do slow doubling
    raise ValueError("Funny case")
  
  if (print_matrices) :
    print("M = ")
    print(Matrix(C.K, [
      [a1, a2, a3, a4, a5],
      [a6, a7, a8, a9, a10],
      [a11, a12, a13, a14, a15]]))
    print
  
  if (a1 == 0) and (a6 == 0) and (a11 == 0) :
    raise ValueError("Double of divisor is not typical.")
  
  if (a1 == 0) :
    if (a6 != 0) :
      a1, a2, a3, a4, a5, a6, a7, a8, a9, a10 = \
          a6, a7, a8, a9, a10, a1, a2, a3, a4, a5
    else :
      a1, a2, a3, a4, a5, a11, a12, a13, a14, a15 = \
          a11, a12, a13, a14, a15, a1, a2, a3, a4, a5

  # Row reduce M' to row echelon form
  #
  #        [ a1  a2  a3  a4  a5 ]
  #   M' = [ 0   b1  b2  b3  b4 ]
  #        [ 0   0   c1  c2  c3 ]
  d1 = a1*a12 - a2*a11
  d2 = a6*a12 - a7*a11
  b1 = a1*a7  - a2*a6
  b2 = a1*a8  - a3*a6
  b3 = a1*a9  - a4*a6
  b4 = a1*a10 - a5*a6
  c1 = b1*a13 - d1*a8  + d2*a3
  c2 = b1*a14 - d1*a9  + d2*a4
  c3 = b1*a15 - d1*a10 + d2*a5
  # Subtotal : 0I 21M 12A
  
  if (b1 == 0) or (c1 == 0) :
    raise ValueError("Sum is not typical.")
  
  # Reduce M even more via back-substitution
  #
  #         [ Z  0  0  A1  A2 ]
  #   M'' = [ 0  Z  0  B1  B2 ]
  #         [ 0  0  Z  C1  C2 ]
  e1 = b3*c1 - b2*c2
  e2 = b4*c1 - b2*c3
  AB = a1*b1
  Z  = AB*c1
  C1 = AB*c2
  C2 = AB*c3
  B1 = a1*e1
  B2 = a1*e2
  A1 = b1*(a4*c1 - c2*a3) - a2*e1
  A2 = b1*(a5*c1 - c3*a3) - a2*e2
  # Subtotal : 0I 18M 6A
  
  # Compute
  #
  #   u = Z*x*f - C1*h - B1*g - A1*f
  #   v = Z*x*g - C2*h - B2*g - A2*f
  # u0 =      - C1*h0 - B1*g0 - A1*f0
  u1 = Z*f0 - C1*h1 - B1*g1 - A1*f1
  u2 =      - C1*h2 - B1*g2 - A1*f2
  u3 = Z*f1 - A1
  u4 = Z*f2 - B1
  u5 =      - C1
  # v0 =      - C2*h0 - B2*g0 - A2*f0
  v1 = Z*g0 - C2*h1 - B2*g1 - A2*f1
  v2 =      - C2*h2 - B2*g2 - A2*f2
  v3 = Z*g1 - A2
  v4 = Z*g2 - B2
  v5 =      - C2
  # Subtotal : 0I 18M 14A
  
  # Compute some inverses
  c0, c1, c2, c3, c4, c5, c6, c7, c8 = C.coefficients()
  ZZt0      = u5^2 - Z*(u5*c8 - u4 + v5)
  if (ZZt0 == 0) :
    raise ValueError("Sum of divisors is non-typical.")
  ZZZt0     = Z*ZZt0
  ZZZt0_inv = 1/ZZZt0
  ZZt0_inv  = Z*ZZZt0_inv
  zeta      = ZZt0*ZZZt0_inv # 1/Z
  tau       = Z*Z*ZZt0_inv   # 1/t0
  # Subtotal : 1I 7M 3A
  
  # Rescale u and v polynomials by 1/Z
  u1 = zeta*u1
  u2 = zeta*u2
  u3 = zeta*u3
  u4 = zeta*u4
  u5 = zeta*u5
  v1 = zeta*v1
  v2 = zeta*v2
  v3 = zeta*v3
  v4 = zeta*v4
  v5 = zeta*v5
  # Subtotal : 0I 10M 0A
  
  # Compute ff, gg such that gg*u = ff*v (mod C)
  gg3 = u5
  ff2 = u5*(u5 - c8) + u4 - v5
  gg2 = v4 + v5*(u5 - c8) + tau*(u5*(u5*(u3 - c6) + v5*(u4 - c7) + c5 - v3) + v5*(u3 - v4) - u2)
  e3 = ff2*v5 - gg2*u5
  ff1 = u5*(u4 - c7) + gg2 + u3 - v4
  gg1 = u5*(c6 - u3) - e3 + v3
  ff0 = c7*e3 + u5*(u2 - c4) + gg2*u3 + gg1*u4 - ff2*v3 - ff1*v4 + u1 - v2
  gg0 = -c6*e3 + u5*(c3 - u1) - gg1*u3 + ff1*v3 + v1
  # Subtotal : 0I 21M 35A
  
  # Reduce gg modulo ff
  gg2 = gg2 - gg3*ff2
  gg1 = gg1 - gg3*ff1
  gg0 = gg0 - gg3*ff0 # Save 1M by combining this with calculation of gg0 above
  # Subtotal : 0I 3M 3A
  
  # Compute third polynomial ...
  hh0 = tau*(ff0*gg1 + gg0*(gg2 - ff1))
  hh1 = tau*(gg1*gg2 - gg0)
  hh2 = gg1 + tau*(gg2*(gg2 - ff1) + ff0)
  # Subtotal : 0I 7M 6A
  
  ret = C34CurveDivisor(C, [[ff0, ff1, ff2, 1],
                       [gg0, gg1, gg2, 0, 1],
                       [hh0, hh1, hh2, 0, 0, 1]],
                       degree = 3, typ = 31, typical = True, reduced = True)
  ret.inv = tau
  return ret



def fast_double_31_high_char(D) :
  # In this version, we assume that the curve polynomial has coefficients
  # c5, c6, c8 = 0.
  C = D.C
  f0, f1, f2 = D.f[0:3]
  g0, g1, g2 = D.g[0:3]
  h0, h1, h2 = D.h[0:3]
  c0, c1, c2, c3, c4, c5, c6, c7, c8 = C.coefficients()
  
  print_matrices = False
  
  if (D.type != 31) :
    raise ValueError("Divisor is not of type 31.\nD = {}".format(D))
  if (c5 != 0) or (c6 != 0) or (c8 != 0):
    raise ValueError("Curve equation is not in short form. C = {}".format(C))
  if (f2 == 0) :
    raise ValueError("Divisor is not typical.\nD = {}".format(D))
    
  # Compute two solutions
  #
  #   rf + sg + th = 0
  #   Rf + Sg + Th = C
  #
  # where
  # 
  #   r = y + r0
  #   s = -(x + s0)
  #   t = t0 = -f2
  #   R = x^2 + R2y + R1x + R0
  #   S = S0
  #   T = y + T1x + T0
  r0 = g1
  s0 = f1 - g2
  t0 = -f2

  # T1 = 0
  # R1 = - f1
  R2 = c7 - f2
  T0 = - h2 - f2*R2
  S0 = c4 - h1 + f1*(f2 - R2)
  R0 = c3 + f1*f1 - f0
  # Subtotal : 0I 2M 1SQ 8A

  # Compute
  #
  #   df = St - sT (mod f, g, h)
  #   dg = Tr - tR (mod f, g, h)
  #   dh = Rs - rSe2g, h)
  e1 = -f1 - g2
  e2 = R2 - f2
  df2 = s0 - g2
  df1 = T0 - g1
  df0 = T0*s0 + S0*t0 - g0
  dg2 = T0 - h2 + r0 - t0*e2
  dg1 = t0*(f1 + f1) - h1 # XXX : r0 - g1 = 0 Change this in other version
  dg0 = T0*r0 + t0*(f0 - R0) - h0
  dh2 = f2*(e1 - g2) + R2*(g2 - s0) - S0
  dh1 = f1*e1 + g1*e2 + f1*s0 - R0 + f0
  dh0 = f0*e1 + g0*e2 - S0*r0 - R0*s0
  # Subtotal : 0I 15M 25A
  # Running total : 0I 17M 1SQ 33A

  # Compute M
  #
  #       [ a1   a2   a3   a4   a5  ]
  #   M = [ a6   a7   a8   a9   a10 ]
  #       [ a11  a12  a13  a14  a15 ]
  #
  # Where the columns represent the reductions of G, H, yG, xG, xH, modulo f, g
  a1,  a2,  a3  = df0, dg0, dh0
  a6,  a7,  a8  = df1, dg1, dh1
  a11, a12, a13 = df2, dg2, dh2
  a4  =    - f0*a6 - g0*a11
  a9  = a1 - f1*a6 - g1*a11
  a14 =    - f2*a6 - g2*a11
  a5  =    - f0*a7 - g0*a12
  a10 = a2 - f1*a7 - g1*a12
  a15 =    - f2*a7 - g2*a12
  # Subtotal : 0I 12M 8A
  # Running total : 0I 29M 1SQ 41A
  
  if (print_matrices) :
    print("M = ")
    print(Matrix(C.K, [
      [a1, a2, a3, a4, a5],
      [a6, a7, a8, a9, a10],
      [a11, a12, a13, a14, a15]]))
    print

  if (a1 == 0) and (a6 == 0) and (a11 == 0) :
    raise ValueError("Sum is not typical.")

  if (a1 == 0) :
    if (a6 != 0) :
      a1, a2, a3, a4, a5, a6, a7, a8, a9, a10 = \
          a6, a7, a8, a9, a10, a1, a2, a3, a4, a5
    else :
      a1, a2, a3, a4, a5, a11, a12, a13, a14, a15 = \
          a11, a12, a13, a14, a15, a1, a2, a3, a4, a5
  
  # Row reduce M' to row echelon form
  #
  #        [ a1  a2  a3  a4  a5 ]
  #   M' = [ 0   b1  b2  b3  b4 ]
  #        [ 0   0   c1  c2  c3 ]
  d1 = a1*a12 - a2*a11
  d2 = a6*a12 - a7*a11
  b1 = a1*a7  - a2*a6
  b2 = a1*a8  - a3*a6
  b3 = a1*a9  - a4*a6
  b4 = a1*a10 - a5*a6
  c1 = b1*a13 - d1*a8  + d2*a3
  c2 = b1*a14 - d1*a9  + d2*a4
  c3 = b1*a15 - d1*a10 + d2*a5
  # Subtotal : 0I 21M 12A
  # Running total : 0I 50M 1SQ 53A
  
  if (b1 == 0) or (c1 == 0) :
    raise ValueError("Sum is not typical.")

  # Reduce M even more via back-substitution
  #
  #         [ Z  0  0  A1  A2 ]
  #   M'' = [ 0  Z  0  B1  B2 ]
  #         [ 0  0  Z  C1  C2 ]
  e1 = b3*c1 - b2*c2
  e2 = b4*c1 - b2*c3
  AB = a1*b1
  Z  = AB*c1
  C1 = AB*c2
  C2 = AB*c3
  B1 = a1*e1
  B2 = a1*e2
  A1 = b1*(a4*c1 - c2*a3) - a2*e1
  A2 = b1*(a5*c1 - c3*a3) - a2*e2
  # Subtotal : 0I 18M 6A
  # Running total : 0I 68M 1SQ 59A

  # Compute
  #
  #   u = Z*x*f - C1*h - B1*g - A1*f
  #   v = Z*x*g - C2*h - B2*g - A2*f
  # u0 =      - C1*h0 - B1*g0 - A1*f0
  u1 = Z*f0 - C1*h1 - B1*g1 - A1*f1
  u2 =      - C1*h2 - B1*g2 - A1*f2
  u3 = Z*f1 - A1
  u4 = Z*f2 - B1
  u5 =      - C1
  # v0 =      - C2*h0 - B2*g0 - A2*f0
  v1 = Z*g0 - C2*h1 - B2*g1 - A2*f1
  v2 =      - C2*h2 - B2*g2 - A2*f2
  v3 = Z*g1 - A2
  v4 = Z*g2 - B2
  v5 =      - C2
  # Subtotal : 0I 18M 14A
  # Running total : 0I 86M 1SQ 73A
  
  # Compute some inverses
  c0, c1, c2, c3, c4, c5, c6, c7, c8 = C.coefficients()
  u5u5 = u5^2
  ZZt0      = Z*(u4 - v5)
  if (ZZt0 == 0) :
    raise ValueError("Sum of divisors is non-typical.")
  ZZZt0     = Z*ZZt0
  ZZZt0_inv = 1/ZZZt0
  ZZt0_inv  = Z*ZZZt0_inv
  zeta      = ZZt0*ZZZt0_inv # 1/Z
  tau       = Z*Z*ZZt0_inv   # 1/t0
  # Subtotal : 1I 6M 1SQ 1A
  # Running total : 1I 92M 2SQ 74A
  
  # Rescale u and v polynomials by 1/Z
  u1 = zeta*u1
  u2 = zeta*u2
  u3 = zeta*u3
  u4 = zeta*u4
  u5 = zeta*u5
  v1 = zeta*v1
  v2 = zeta*v2
  v3 = zeta*v3
  v4 = zeta*v4
  v5 = zeta*v5
  # Subtotal : 0I 10M 0A
  # Running total : 1I 102M 2SQ 74A
  
  # Compute ff, gg such that gg*u = ff*v (mod C)
  u3u5 = u3*u5
  gg3 = u5
  ff2 = u5u5 + u4 - v5
  gg2 = v4 + u5*v5 + tau*(u5*(u3u5 + v5*(u4 - c7) + c5 - v3) + v5*(u3 - v4) - u2)
  e3 = ff2*v5 - gg2*u5
  ff1 = u5*(u4 - c7) + gg2 + u3 - v4
  gg1 = u3u5 - e3 + v3
  ff0 = c7*e3 + u5*(u2 - c4) + gg2*u3 + gg1*u4 - ff2*v3 - ff1*v4 + u1 - v2
  gg0 = u5*(c3 - u1) - gg1*u3 + ff1*v3 + v1
  # Subtotal : 0I 18M 30A
  # Running total : 1I 120M 2SQ 104A

  # Reduce gg modulo ff
  gg2 = gg2 - gg3*ff2
  gg1 = gg1 - gg3*ff1
  gg0 = gg0 - gg3*ff0 # Save 1M by combining this with calculation of gg0 above
  # Subtotal : 0I 3M 3A
  # Running total : 1I 123M 2SQ 107A

  # Compute third polynomial ...
  hh0 = tau*(ff0*gg1 + gg0*(gg2 - ff1))
  hh1 = tau*(gg1*gg2 - gg0)
  hh2 = gg1 + tau*(gg2*(gg2 - ff1) + ff0)
  # Subtotal : 0I 7M 6A
  # Running total : 1I 130M 2SQ 113A
  
  return C34CurveDivisor(C, [[ff0, ff1, ff2, 1],
                       [gg0, gg1, gg2, 0, 1],
                       [hh0, hh1, hh2, 0, 0, 1]],
                       degree = 3, typ = 31, typical = True, reduced = True)



def fast_double_31(D) :
  # TODO : Change this so that it works over char K = 2 or 3
  C = D.C
  f0, f1, f2 = D.f[0:3]
  g0, g1, g2 = D.g[0:3]
  h0, h1, h2 = D.h[0:3]
  c0, c1, c2, c3, c4, c5, c6, c7, c8 = C.coefficients()
  
  print_matrices = False
  
  if (f2 == 0) :
    raise ValueError("Divisor is not typical.\nD = {}".format(D))
    
  # Compute two solutions
  #
  #   rf + sg + th = 0
  #   Rf + Sg + Th = C
  #
  # where
  # 
  #   r = y + r0
  #   s = -(x + s0)
  #   t = t0 = -f2
  #   R = x^2 + R2y + R1x + R0
  #   S = S0
  #   T = y + T1x + T0
  r0 = g1
  s0 = f1 - g2
  t0 = -f2

  T1 = c8
  R2 = c7 - f2
  R1 = c6 - f1
  T0 = c5 - h2 - f2*R2
  S0 = c4 - h2*T1 - h1 - f2*R1 - f1*R2
  R0 = c3 - h1*T1 - f1*R1 - f0
  # Subtotal : 0I 6M 12A (in high characteristic, reduces to 2M 1SQ 7A)

  # Compute
  #
  #   df = St - sT (mod f, g, h)
  #   dg = Tr - tR (mod f, g, h)
  #   dh = Rs - rS (mod f, g, h)
  df2 = s0 - g2 - T1*f2
  df1 = T0 - g1 + T1*(s0 - f1)
  df0 = T0*s0 + S0*t0 - T1*f0 - g0
  dg2 = T0 - h2 + r0 - T1*g2 + t0*(f2 - R2)
  dg1 = T1*(r0 - g1) + t0*(f1 - R1) - h1
  dg0 = T0*r0 - T1*g0 + t0*(f0 - R0) - h0
  dh2 = f2*(R1 - f1 - g2 + s0) + R2*(g2 - s0) - S0
  dh1 = f1*(R1 - f1 + s0) + g1*(R2 - f2) - R1*s0 - R0 + f0
  dh0 = f0*(R1 - f1 + s0) + g0*(R2 - f2) - S0*r0 - R0*s0
  # Subtotal : 0I 21M 40A
  # Running total : 0I 27M 52A

  # Compute M
  #
  #       [ a1   a2   a3   a4   a5  ]
  #   M = [ a6   a7   a8   a9   a10 ]
  #       [ a11  a12  a13  a14  a15 ]
  #
  # Where the columns represent the reductions of G, H, yG, xG, xH, modulo f, g
  a1,  a2,  a3  = df0, dg0, dh0
  a6,  a7,  a8  = df1, dg1, dh1
  a11, a12, a13 = df2, dg2, dh2
  a4  =    - f0*a6 - g0*a11
  a9  = a1 - f1*a6 - g1*a11
  a14 =    - f2*a6 - g2*a11
  a5  =    - f0*a7 - g0*a12
  a10 = a2 - f1*a7 - g1*a12
  a15 =    - f2*a7 - g2*a12
  # Subtotal : 0I 12M 8A
  # Running total : 0I 39M 60A (Several multiplications saved in high characteristic)
  
  if (print_matrices) :
    print("M = ")
    print(Matrix(C.K, [
      [a1, a2, a3, a4, a5],
      [a6, a7, a8, a9, a10],
      [a11, a12, a13, a14, a15]]))
    print

  if (a1 == 0) and (a6 == 0) and (a11 == 0) :
    raise ValueError("Sum is not typical.")

  if (a1 == 0) :
    if (a6 != 0) :
      a1, a2, a3, a4, a5, a6, a7, a8, a9, a10 = \
          a6, a7, a8, a9, a10, a1, a2, a3, a4, a5
    else :
      a1, a2, a3, a4, a5, a11, a12, a13, a14, a15 = \
          a11, a12, a13, a14, a15, a1, a2, a3, a4, a5
  
  # Row reduce M' to row echelon form
  #
  #        [ a1  a2  a3  a4  a5 ]
  #   M' = [ 0   b1  b2  b3  b4 ]
  #        [ 0   0   c1  c2  c3 ]
  d1 = a1*a12 - a2*a11
  d2 = a6*a12 - a7*a11
  b1 = a1*a7  - a2*a6
  b2 = a1*a8  - a3*a6
  b3 = a1*a9  - a4*a6
  b4 = a1*a10 - a5*a6
  c1 = b1*a13 - d1*a8  + d2*a3
  c2 = b1*a14 - d1*a9  + d2*a4
  c3 = b1*a15 - d1*a10 + d2*a5
  # Subtotal : 0I 21M 12A
  # Running total : 0I 60M 72A
  
  if (b1 == 0) or (c1 == 0) :
    raise ValueError("Sum is not typical.")

  # Reduce M even more via back-substitution
  #
  #         [ Z  0  0  A1  A2 ]
  #   M'' = [ 0  Z  0  B1  B2 ]
  #         [ 0  0  Z  C1  C2 ]
  e1 = b3*c1 - b2*c2
  e2 = b4*c1 - b2*c3
  AB = a1*b1
  Z  = AB*c1
  C1 = AB*c2
  C2 = AB*c3
  B1 = a1*e1
  B2 = a1*e2
  A1 = b1*(a4*c1 - c2*a3) - a2*e1
  A2 = b1*(a5*c1 - c3*a3) - a2*e2
  # Subtotal : 0I 18M 6A
  # Running total : 0I 78M 78A

  # Compute
  #
  #   u = Z*x*f - C1*h - B1*g - A1*f
  #   v = Z*x*g - C2*h - B2*g - A2*f
  # u0 =      - C1*h0 - B1*g0 - A1*f0
  u1 = Z*f0 - C1*h1 - B1*g1 - A1*f1
  u2 =      - C1*h2 - B1*g2 - A1*f2
  u3 = Z*f1 - A1
  u4 = Z*f2 - B1
  u5 =      - C1
  # v0 =      - C2*h0 - B2*g0 - A2*f0
  v1 = Z*g0 - C2*h1 - B2*g1 - A2*f1
  v2 =      - C2*h2 - B2*g2 - A2*f2
  v3 = Z*g1 - A2
  v4 = Z*g2 - B2
  v5 =      - C2
  # Subtotal : 0I 18M 14A
  # Running total : 0I 96M 92A
  
  # Compute some inverses
  c0, c1, c2, c3, c4, c5, c6, c7, c8 = C.coefficients()
  ZZt0      = u5^2 - Z*(u5*c8 - u4 + v5)
  if (ZZt0 == 0) :
    raise ValueError("Sum of divisors is non-typical.")
  ZZZt0     = Z*ZZt0
  ZZZt0_inv = 1/ZZZt0
  ZZt0_inv  = Z*ZZZt0_inv
  zeta      = ZZt0*ZZZt0_inv # 1/Z
  tau       = Z*Z*ZZt0_inv   # 1/t0
  # Subtotal : 1I 7M 3A
  # Running total : 1I 103M 95A
  
  # Rescale u and v polynomials by 1/Z
  u1 = zeta*u1
  u2 = zeta*u2
  u3 = zeta*u3
  u4 = zeta*u4
  u5 = zeta*u5
  v1 = zeta*v1
  v2 = zeta*v2
  v3 = zeta*v3
  v4 = zeta*v4
  v5 = zeta*v5
  # Subtotal : 0I 10M 0A
  # Running total : 1I 113M 95A
  
  # Compute ff, gg such that gg*u = ff*v (mod C)
  gg3 = u5
  ff2 = u5*(u5 - c8) + u4 - v5
  gg2 = v4 + v5*(u5 - c8) + tau*(u5*(u5*(u3 - c6) + v5*(u4 - c7) + c5 - v3) + v5*(u3 - v4) - u2)
  e3 = ff2*v5 - gg2*u5
  ff1 = u5*(u4 - c7) + gg2 + u3 - v4
  gg1 = u5*(c6 - u3) - e3 + v3
  ff0 = c7*e3 + u5*(u2 - c4) + gg2*u3 + gg1*u4 - ff2*v3 - ff1*v4 + u1 - v2
  gg0 = -c6*e3 + u5*(c3 - u1) - gg1*u3 + ff1*v3 + v1
  # Subtotal : 0I 21M 35A
  # Running total : 1I 134M 130A

  # Reduce gg modulo ff
  gg2 = gg2 - gg3*ff2
  gg1 = gg1 - gg3*ff1
  gg0 = gg0 - gg3*ff0 # Save 1M by combining this with calculation of gg0 above
  # Subtotal : 0I 3M 3A
  # Running total : 1I 137M 133A

  # Compute third polynomial ...
  hh0 = tau*(ff0*gg1 + gg0*(gg2 - ff1))
  hh1 = tau*(gg1*gg2 - gg0)
  hh2 = gg1 + tau*(gg2*(gg2 - ff1) + ff0)
  # Subtotal : 0I 7M 6A
  # Running total : 1I 144M 139A
  
  return C34CurveDivisor(C, [[ff0, ff1, ff2, 1],
                       [gg0, gg1, gg2, 0, 1],
                       [hh0, hh1, hh2, 0, 0, 1]],
                       degree = 3, typ = 31, typical = True, reduced = True)



def km_double_31(D) :
  C = D.C
  f0, f1, f2 = D.f[0:3]
  g0, g1, g2 = D.g[0:3]
  c0, c1, c2, c3, c4, c5, c6, c7, c8 = C.coefficients()
  
  print_matrices = False
  strassen = True
  toom_cook = True
  karatsuba = True

  if (D.type != 31) :
    raise ValueError("Divisor is not of type 31.\nD = {}".format(D))
  if (f2 == 0) :
    raise ValueError("Divisor is not typical.\nD = {}".format(D))
  if (C.K.characteristic() <= 3) :
    raise ValueError("Curve's base field is of characteristic 3 or less.")
  if (c5 != 0) or (c6 != 0) or (c8 != 0) :
    raise ValueError("Curve equation is not in short form.\nC = {}".format(C))
  half = C.K(1/2) # Assumed to be a "free" computation in Kamal's model
  if (D.inv == 0) :
    D.inv = 1/f2
  f2_inv = D.inv

  # Find polynomials
  #
  #   G = xy + G2*y + G1*x + G0
  #   H = y^2 + H3*x^2 + H2*y + H1*x + H0
  #
  # Satisfying fH - gG = 0 (mod C)
  l = f0 + g2*(g2 - f1)
  m = f2*(c7 - f2) + g1
  l_over_f2 = l*f2_inv
  f1f2 = f1*f2
  G2 = f1 - g2
  G1 = -l_over_f2 - m
  G0 = g2*m + (l_over_f2 + g1)*(g2 - f1) - f2*(f1f2 + c4) - g0
  H3 = f2
  H2 = -l_over_f2
  H1 = -f1f2
  H0 = -(l_over_f2 + m)*g1 + f2*(f1^2 - f0 + c3)
  # Subtotal : 0I 9M 1SQ 16A
  # Running total : 0I 9M 1SQ 16A

  # Compute the matrix Ty
  #
  #        [ 0  -g0  -h0 ]
  #   Ty = [ 0  -g1  -h1 ]
  #        [ 1  -g2  -h2 ]
  #
  # (This is equivalent to computing the third polynomial in the reduced Groebner basis
  # for D, h = y^2 + h2*y + h1*x + h0).
  h2 = l_over_f2 + g1
  h1 = f2_inv*(g1*g2 - g0)
  h0 = f2_inv*(g1*f0 + g0*(g2 - f1))
  # print("h = {}".format(y^2 + h2*y + h1*x + h0))
  # Subtotal : 0I 5M 4A
  # Running total : 0I 14M 1SQ 20A

  # Compute M
  #
  #       [ a1   a2   a3   a4   a5  ]
  #   M = [ a6   a7   a8   a9   a10 ]
  #       [ a11  a12  a13  a14  a15 ]
  #
  # Where the columns represent the reductions of G, xG, yG, H, xH, modulo f, g
  a1  = G0 - g0
  a6  = G1 - g1
  a11 = G2 - g2
  a4  = H0 - h0 - f0*f2
  a9  = H1 - h1 - f1f2
  a14 = H2 - h2 - f2^2
  a3  =    - a6*g0 - a11*h0
  a8  =    - a6*g1 - a11*h1
  a13 = a1 - a6*g2 - a11*h2
  
  if (strassen) :
    m1 = (-f1 - g2)*(a6 + a14)
    m2 = (-f2 - g2)*a6
    m3 = -f1*(a9 - a14)
    m4 = -g2*(a11 - a6)
    m5 = (-f1 - g1)*a14
    m6 = (-f2 + f1)*(a6 + a9)
    m7 = (-g1 + g2)*(a11 + a14)
    a2  = -f0*a6 - g0*a11
    a7  = a1 + m1 + m4 - m5 + m7
    a12 = m2 + m4
    a5  = -f0*a9 - g0*a14
    a10 = a4 + m3 + m5
    a15 = m1 - m2 + m3 + m6
  else :
    a2  =    - f0*a6 - g0*a11
    a7  = a1 - f1*a6 - g1*a11
    a12 =    - f2*a6 - g2*a11
    a5  =    - f0*a9 - g0*a14
    a10 = a4 - f1*a9 - g1*a14
    a15 =    - f2*a9 - g2*a14
  
  if (print_matrices) :
    print("M = ")
    print(Matrix(C.K, [
      [a1, a2, a3, a4, a5],
      [a6, a7, a8, a9, a10],
      [a11, a12, a13, a14, a15]]))
    print
  # Subtotal : 0I 18M 1SQ 35A (assuming Strassed used)
  # Running total : 0I 32M 2SQ 55A

  if (a1 == 0) and (a6 == 0) and (a11 == 0) :
    raise ValueError("Sum is not typical".format(D2))
  if (a1 == 0) :
    if (a6 != 0) :
      a1, a2, a3, a4, a5, a6, a7, a8, a9, a10 = \
          a6, a7, a8, a9, a10, a1, a2, a3, a4, a5
    else :
      a1, a2, a3, a4, a5, a11, a12, a13, a14, a15 = \
          a11, a12, a13, a14, a15, a1, a2, a3, a4, a5

  # Construct the modified matrix M'
  #
  #        [ A1   A2   A3   A4   A5  ]
  #   M' = [ A6   A7   A8   A9   A10 ]
  #        [ A11  A12  A13  A14  A15 ]
  A1,  A2 , A3,  A4,  A5  = a1,  a4,  a3  - a5,  a2,  a5
  A6,  A7,  A8,  A9,  A10 = a6,  a9,  a8  - a10, a7,  a10
  A11, A12, A13, A14, A15 = a11, a14, a13 - a15, a12, a15
  # Subtotal : 0I 0M 5A
  # Running total : 0I 32M 2SQ 60A

  if (print_matrices) :
    print("M' = ")
    print(Matrix(C.K, [
      [A1, A2, A3, A4, A5],
      [A6, A7, A8, A9, A10],
      [A11, A12, A13, A14, A15]]))
    print

  # Find a basis for ker M'
  # Begin by row reducing M' to echelon form
  #
  #            [ A1  A2  A3  A4  A5 ]
  #   M'_ref = [ 0   B1  B2  B3  B4 ]
  #            [ 0   0   C1  C2  C3 ]
  D1 = A1*A12 - A2*A11
  D2 = A6*A12 - A7*A11
  B1 = A1*A7  - A2*A6
  B2 = A1*A8  - A3*A6
  B3 = A1*A9  - A4*A6
  B4 = A1*A10 - A5*A6
  C1 = B1*A13 - D1*A8  + D2*A3
  C2 = B1*A14 - D1*A9  + D2*A4
  C3 = B1*A15 - D1*A10 + D2*A5
  # Subtotal : 0I 21M 12A
  # Running total : 0I 53M 2SQ 72A

  if (B1 == 0) and (D1 == 0) :
    raise ValueError("Sum is not typical".format(D2))
  if (B1 == 0) :
    # M' is row-reducible to 
    #
    #   [ A1   A2   A3   A4   A5  ]
    #   [ A11  A12  A13  A14  A15 ]
    #   [ 0    0    B2   B3   B4  ]
    #
    # M' may still be full rank if D1 is non-zero.
    # We do some re-labeling and attempt to reduce M' more.
    C1, C2, C3 = B2, B3, B4
    B1 = D1
    B2 = A1*A13 - A3*A11
    B3 = A1*A14 - A4*A11
    B4 = A1*A15 - A5*A11
    # After relabelling, M' is full rank iff C1 != 0
  if (C1 == 0) :
    raise ValueError("Sum is not typical".format(D2))
  
  if (print_matrices) :
    print("M'_ref = ")
    print(Matrix(C.K, [
      [ A1, A2, A3, A4, A5 ],
      [ 0, B1, B2, B3, B4 ],
      [ 0, 0, C1, C2, C3 ]]))
    print

  # Compute inverses of A1, B1, C1
  A1B1       = A1 * B1
  A1B1C1     = A1B1 * C1
  A1B1C1_inv = 1 / A1B1C1
  C1_inv     = A1B1 * A1B1C1_inv
  A1B1_inv   = C1 * A1B1C1_inv
  B1_inv     = A1 * A1B1_inv
  A1_inv     = B1 * A1B1_inv
  # Subtotal : 1I 6M 0A
  # Running total : 1I 59M 2SQ 72A
  
  # Compute RREF of M'
  #
  #             [ 1  0  0  -r0  -s0 ]
  #   M'_rref = [ 0  1  0  -r1  -s1 ]
  #             [ 0  0  1  -r2  -s2 ]
  r2 = - C1_inv * C2
  r1 = - B1_inv * (B3 + B2*r2)
  r0 = - A1_inv * (A4 + A3*r2 + A2*r1)
  s2 = - C1_inv * C3
  s1 = - B1_inv * (B4 + B2*s2)
  s0 = - A1_inv * (A5 + A3*s2 + A2*s1)
  # Subtotal : 0I 12M 6A
  # Running total : 1I 71M 2SQ 78A

  if (print_matrices) :
    print("M'_rref = ")
    print(Matrix(C.K, [
      [1, 0, 0, -r0, -s0],
      [0, 1, 0, -r1, -s1],
      [0, 0, 1, -r2, -s2]]))
    print

  # Find polynomials u, v generating the ideal of D1 + D2
  #
  #   u = xf + r2(yf - xg) + r1*g + r0*f
  #   v = xg + s2(yf - xg) + s1*g + s0*f
  #
  # Computing u, e.g., requires computing r0*f0, r0*f2 + r2*f0, r2*f2
  # as well as r1*g0, r1*g1 - r2*g0, - r2*g1. Karatsuba multiplication may be applied
  # t0 save 2 multiplications here.
  if (karatsuba) :
    r0f0 = r0*f0
    r2f2 = r2*f2
    r1g0 = r1*g0
    r2g1 = r2*g1
    u0 = r0f0 + r1g0
    u1 = r0*f1 + f0 + (r1 - r2)*(g0 + g1) - r1g0 + r2g1
    u2 = r1*g2 + (r0 + r2)*(f0 + f2) - r0f0 - r2f2
    u3 = r0 - r2g1 + f1
    u4 = r1 + r2*(f1 - g2) + f2
    u5 = r2f2
    
    s0f0 = s0*f0
    s2f2 = s2*f2
    s1g0 = s1*g0
    s2g1 = s2*g1
    v0 = s0f0 + s1g0
    v1 = s0*f1 + g0 + (s1 - s2)*(g0 + g1) - s1g0 + s2g1
    v2 = s1*g2 + (s0 + s2)*(f0 + f2) - s0f0 - s2f2
    v3 = s0 - s2g1 + g1
    v4 = s1 + s2*(f1 - g2) + g2
    v5 = s2f2
  else :
    u0 = r0*f0 + r1*g0
    u1 = r0*f1 + r1*g1 - r2*g0 + f0
    u2 = r0*f2 + r1*g2 + r2*f0
    u3 = r0 - r2*g1 + f1
    u4 = r1 + r2*(f1 - g2) + f2
    u5 = r2*f2
    v0 = s0*f0 + r1*g0
    v1 = s0*f1 + r1*g1 - r2*g0 + g0
    v2 = s0*f2 + r1*g2 + r2*f0
    v3 = s0 - r2*g1 + g1
    v4 = s1 + r2*(f1 - g2) + g2
    v5 = s2*f2
  # Subtotal : 0I 18M 34A (assuming Karatsuba used)
  # Running total : 1I 89M 2SQ 112A

  # Now reduce the divisor div(u, v) following the improvment in [KM2018]
  # The values l1, l2, l3 here correspond to l1, l2, l3 in [KM2018]
  # The values l4, .., l7 here correspond to m0, .., m3 in [KM2018]
  # The values t1, .., t6 here save 2M in computing l4, .., l7, as per the advice of KM.
  l1 = v5 - u4 - u5*u5
  if (l1 == 0) :
    raise ValueError("Sum is not typical".format(D2))
  l2 = v4 - u3 - u5*(u4 - c7)
  l3 = - v3 + u3*u5

  if (toom_cook) :
    t1 = l1*v5
    t2 = l2*v3
    t3 = (v5 + v4 + v3)*(l1 + l2)
    t4 = (v5 - v4 + v3)*(l1 - l2)
    t5 = t3 - t1 - t2
    t6 = t4 - t1 + t2

    l4 = l3 - t1
    l5 = - u2 - half*(t5 - t6) + l4*u5
    l6 = v2 - u1 - u5*(u2 - c4) + c7*l3 - half*(t5 + t6) + l4*(u4 - c7)
    l7 = v1 - u5*(u1 - c3) - t2 + l4*u3
  else :
    l1 = v5 - u4 - u5^2
    l2 = v4 - u3 - u5*(u4 - c7)
    l3 = - v3 + u3*u5
    l4 = l3 - l1*v5
    l5 = - u2 - (l1*v4 + l2*v5) + l4*u5
    l6 = v2 - u1 - u5*(u2 - c4) + c7*l3 - v3*l1 - l2*v4 + l4*(u4 - c7)
    l7 = v1 - u5*(u1 - c3) - l2*v3 + l4*u3
  # Subtotal : 0I 11M 1SQ 1CC 2CM 32A (assuming Toom-Cook used)
  # Running total : 1I 100M 3SQ 1CC 2CM 144A

  # Compute polynomials f, g generating the reduction of D1 + D2
  l1_inv = 1/l1
  l5_over_l1 = l5*l1_inv
  new_f2 = - l1
  new_f1 = - l5_over_l1 - l2
  new_f0 =   l5_over_l1*l2 - l6
  new_g2 = - l5_over_l1 - u5*new_f2
  new_g1 = - u5*(l5_over_l1 + new_f1) - l4
  new_g0 = l7 + l3*l5_over_l1 - u5*new_f0
  # Subtotal : 1I 6M 7A
  # Running total : 2I 106M 3SQ 1CC 2CM 151A

  ret = C34CurveDivisor(C, [[new_f0, new_f1, new_f2, 1],
                      [new_g0, new_g1, new_g2, 0, 1], []],
                      degree = 3, typ = 31, reduced = True, typical = True, inv = -l1_inv)
  return ret



def double_0(D):
  return C34CurveDivisor(D.C, [[D.K.one()], [], []])



def double_11(D):
  K = D.K
  f, g = D.f, D.g
  c = D.C.coefficients()
  new_f, new_g, new_h = [], [], []
  
  # Need to find G and H such that fH = gG (mod c)
  # G is just the polynomial in flip(D) = A = <F,G>
  A = flip(D) # Costs 0I 4M
  G = A.g
  t = f[0]*(f[0] - c[6]) + c[3]
  H = [ c[1] -f[0]*t,
        t,
        c[4] - c[7]*f[0],
        c[6] - f[0],
        c[7],
        c[8],
        K.one() ]
  # I claim that :
  #   * fH + gG = c
  #   * G(-f0,-g0) = (dc/dy)(-f0,-g0)
  #   * H(-f0,-g0) = (dc/dx)(-f0,-g0)

  # Compute alpha = G(-f[0], -g[0])
  alpha = (g[0] - G[2])*g[0] + G[0]
  # Subtotal : 0I 8M
  
  if alpha == 0 :
    new_f = [f[0], K.one()]
    new_g = [g[0]*g[0], K.zero(), g[0] + g[0], K.zero(), K.zero(), K.one()]
    # Total : 0I 9M
  
  else :
    # Compute beta = H(-f[0], -g[0])
    beta = ((H[3] - f[0])*f[0] + H[4]*g[0] - H[1])*f[0] + (H[5]*g[0] - H[2])*g[0] + H[0] # beta = h(-f0,-g0)
    a = beta*(1/alpha)
    new_f = [g[0] + f[0]*a, a, K.one()]
    new_g = [f[0]*f[0], f[0] + f[0], K.zero(), K.one()]
    # Total : 1I 16M
  
  return C34CurveDivisor(D.C, [new_f, new_g, new_h])



def double_21(D) :
  # In this case, x-coordinates are distinct.
  # XXX: This is wrong... the double of a point can be of type 2a.
  #      2a is more accurately when the ideal is of form <y + ..., x^2 + ...>
  # Maybe we can consider derivatives with respect only to x.
  # Is it faster to use the partial derivates of the curve equation rather than G1, H1?
  # Note in this case, partial der. of D_y(C) and h' agree at the divisor's points. 
  K = D.K
  f, g = D.f, D.g
  c = D.C.coefficients()
  new_f, new_g, new_h = [], [], []
  
  # Find polynomials G, H satisfying fH + gG = C
  A = flip(D) # Costs 0I 7M
  G = A.g

  if G == g :
    # If A.g = D.g,
    # then A = D and 2D = D + A = 0.
    # (More precisely, 2D = <D.f> is principal)
    return C34CurveDivisor(D.C, [[f[0], f[1], K.one()], [], []])
    #return C34CurveDivisor(D.C, [[K.one()], [], []])
  
  t1 = c[8] - f[1]
  t2 = c[5] - f[0]
  H = [ c[2] - f[0]*t2,
        c[4] - f[1]*t2 - f[0]*t1,
        t2,
        c[7] - f[1]*t1,
        t1,
        K.one() ]
  # Subtotal : 0I 11M
  
  # Construct matrix M
  # M = [ a1  a2  a3  a4  a5  ]
  #     [ a6  a7  a8  a9  a10 ]
  #
  # H = t2*x^2 + t1*x + t0
  #   = (t1 - g1*t2)*x + (t0 - g0*t2)
  #
  # T_x = [ 0  -g0 ]  T_y = [ -f0  t3 ]
  #       [ 1  -g1 ]        [ -f1  t4 ]
  t0 = -f[0]*(-f[0] + H[2]) + H[0]
  t1 = -f[1]*(-f[0] + H[2]) - f[0]*(-f[1] + H[4]) + H[1]
  t2 = -f[1]*(-f[1] + H[4]) + H[3]
  t3 = f[1]*g[0]
  t4 = f[1]*g[1] - f[0]
  a1 = G[0] - g[0]
  a6 = G[1] - g[1]
  a2 = g[0]*t2 - t0
  a7 = g[1]*t2 - t1
  a3 = -g[0]*a6
  a8 = a1 - g[1]*a6
  a4 = -f[0]*a1 + t3*a6
  a9 = -f[1]*a1 + t4*a6
  a5 = -g[0]*a7
  a10 = a2 - g[1]*a7
  # Subtotal : 0I 16M

  # If A =/= D, then either a1 =/= 0 or a6 =/= 0
  # If a1 == 0, then swap rows.
  if a1 == 0 :
    a1, a6  = a6,  a1
    a2, a7  = a7,  a2
    a3, a8  = a8,  a3
    a4, a9  = a9,  a4
    a5, a10 = a10, a5

  # Get M to row echelon form
  # M_ref = [ a1  a2  a3  a4  a5 ]
  #         [  0  b1  b2  b3  b4 ]
  b1 = a1*a7  - a2*a6
  b2 = a1*a8  - a3*a6
  b3 = a1*a9  - a4*a6
  b4 = a1*a10 - a5*a6
  # Subtotal : 1I 8M
  
  if b1 != 0 :
    # Reduce matrix M
    # M_rref = [ 1  0  -u0  -v0  -w0 ]
    #          [ 0  1  -u1  -v1  -w1 ]
    # Find kernel of M
    #         [ u0 ]  [ v0 ]  [ w0 ]
    #         [ u1 ]  [ v1 ]  [ w1 ]
    # ker M = [  1 ], [  0 ], [  0 ]
    #         [  0 ]  [  1 ]  [  0 ]
    #         [  0 ]  [  0 ]  [  1 ]
    gamma = 1/(a1*b1)
    alpha = gamma*b1
    beta = gamma*a1
    u1 = -beta*b2
    v1 = -beta*b3
    w1 = -beta*b4
    u0 = -alpha*(a3 + a2*u1)
    v0 = -alpha*(a4 + a2*v1)
    w0 = -alpha*(a5 + a2*w1)
    # Subtotal : 1I 12M

    # Find polynomials forming an ideal generating set for 2D
    new_f = [ f[0]*u0 + g[0]*u1,
              f[1]*u0 + g[1]*u1 + f[0],
              u0,
              u1 + f[1],
              K.one() ]
    new_g = [ f[0]*v0 + g[0]*v1 - f[1]*new_f[0],
              f[1]*v0 + g[1]*v1 - f[1]*new_f[1],
              v0 + f[0] - f[1]*new_f[2],
              v1 - f[1]*new_f[3],
              K.zero(),
              K.one() ]
    new_h = [ f[0]*w0 + g[0]*w1,
              f[1]*w0 + g[1]*w1 + g[0],
              w0,
              w1 + g[1],
              K.zero(),
              K.zero(),
              K.one() ]
    # Subtotal : 0I 16M
    # Total : 1I 63M
    # Notes : This is a lot of multiplications.
    #         If D was computed as the flip of another divisor, then can we save 7M when finding A, above,
    #         since it was already computed earlier. When 2D is typical, the polynomial new_h is redundant.
    #         Not computing it saves 9M. This brings us down to 1I 47M.
    # 2D is of type 41
  elif b2 != 0 :
    # Reduce matrix M
    # M_rref = [ 1  -u0  0  -v0  -w0 ]
    #          [ 0    0  1  -v1  -w1 ]
    # Find kernel of M
    #         [ u0 ]  [ v0 ]  [ w0 ]
    #         [ u1 ]  [ v1 ]  [ w1 ]
    # ker M = [  1 ], [  0 ], [  0 ]
    #         [  0 ]  [  1 ]  [  0 ]
    #         [  0 ]  [  0 ]  [  1 ]
    gamma = 1/(a1*b2)
    alpha = gamma*b2
    beta = gamma*a1
    v1 = -beta*b3
    w1 = -beta*b4
    u0 = -alpha*a2
    v0 = -alpha*(a4 + a3*v1)
    w0 = -alpha*(a5 + a3*w1)
    # Subtotal : 1I 10M

    # Find polynomials forming an ideal generating set for 2D
    # f'' = x^2 + ...
    # g'' = y^2 + ...
    # h'' = x^3 + ...
    # The leading term of f'' divides that of h'', so assume h'' is redundant.
    new_f = [ f[0]*u0 + g[0],
              f[1]*u0 + g[1],
              u0,
              K.one() ]
    t = f[1]*v1
    new_g = [ f[0]*v0 - t*new_f[0],
              f[1]*v0 + f[0]*v1 - t*new_f[1],
              v0 + f[0] - t*new_f[2],
              K.zero(),
              v1 + f[1],
              K.one() ]
    # Subtotal : 0I 9M
    # Total : 1I 54M
    # 2D is of type 43
  elif b3 != 0 :
    # TODO : Verify this case is possible
    # Reduce matrix M
    # M_rref = [ 1  -u0  -v0  0  -w0 ]
    #          [ 0    0    0  1  -w1 ]
    # Find kernel of M
    #         [ u0 ]  [ v0 ]  [ w0 ]
    #         [  1 ]  [  0 ]  [  0 ]
    # ker M = [  0 ], [  1 ], [  0 ]
    #         [  0 ]  [  0 ]  [ w1 ]
    #         [  0 ]  [  0 ]  [  1 ]
    gamma = 1/(a1*b3)
    alpha = gamma*b3
    beta = gamma*a1
    w1 = -beta*b4
    u0 = -alpha*a2
    v0 = -alpha*a3
    w0 = -alpha*(a5 + a4*w1)
    # Subtotal : 1I 8M

    # Find polynomials forming an ideal generating set for 2D
    # f'' = x^2 + ...
    # g'' = xy  + ...
    # h'' = x^3 + ...
    # The leading term of f'' divides that of h'', so assume h'' is redundant.
    new_f = [ f[0]*u0 + g[0],
              f[1]*u0 + g[1],
              u[0],
              K.one() ]
    new_g = [ f[0]*v0        - f[1]*new_f[0],
              f[1]*v0 + f[0] - f[1]*new_f[1],
              v0             - f[1]*new_f[2],
              K.zero(),
              K.one() ]
    # Subtotal : 0I 7M
    # Total : 1I 50M
    # 2D is of type 42
  else :
    # XXX : This would occur if b1 = b2 = b3 = 0, but I don't think that's possible
    raise ValueError("Cannot reduce matrix.")
  return C34CurveDivisor(D.C, [new_f, new_g, new_h])



def double_22(D) :
  K = D.K
  f, g = D.f, D.g
  c = D.C.coefficients()
  new_f, new_g, new_h = [], [], []
  # Find polynomials G, H satisfying fH + gG = C
  A = flip(D) # Costs 0I 1M
  G = A.g
  H = [0, 0, 0, c[6] - f[0], c[7], c[8], K.one()]
  H[2] = c[4] - f[0]*H[4]
  H[1] = c[3] - f[0]*H[3]
  H[0] = c[1] - f[0]*H[1]
  # Subtotal : 0I 4M
  
  # Construct matrix M
  # M = [ a1  a2  a3  a4  a5 ]
  #     [  1  a6  a7  a8  a9 ]
  a1 = G[0]
  a2 = -f[0]*a1
  a3 = -g[0]
  a4 = -H[0] + H[5]*g[0] + f[0]*(H[1] - f[0]*(H[3] - f[0]))
  a5 = -f[0]*a2
  a6 = -f[0]
  a7 = a1 - g[2]
  a8 = -H[2] + H[5]*g[2] + H[4]*f[0]
  a9 = -f[0]*a6
  # Subtotal : 0I 8M

  # Get matrix M into row echelon form.
  # M_ref = [ 1  a6  a7  a8  a9 ]
  #         [ 0   0  b1  b2  b3 ]
  b1 = a3 - a1*a7
  b2 = a4 - a1*a8
  b3 = a5 - a1*a9
  # Subtotal : 0I 3M
 
  if b1 != 0 :
    # Reduce matrix M
    # M_rref = [ 1  -u0  0  -v0  -w0 ]
    #          [ 0    0  1  -v1  -w1 ]
    # Find kernel of M
    #         [ u0 ]  [ v0 ]  [ w0 ]
    #         [  1 ]  [  0 ]  [  0 ]
    # ker M = [  0 ], [ v1 ], [ w1 ]
    #         [  0 ]  [  1 ]  [  0 ]
    #         [  0 ]  [  0 ]  [  1 ]
    beta = 1/b1
    v1 = -beta*b2
    w1 = -beta*b3
    u0 = -a6
    v0 = -(a8 + a7*v1)
    w0 = -(a9 + a7*w1)
    # Subtotal : 1I 4M
    
    # Find polynomials forming an ideal generating set for 2D
    # f'' = x^2 + ...
    # g'' = y^2 + ...
    # h'' = x^3 + ...
    new_f = [ f[0]*u0, f[0] + u0, K.zero(), K.one() ]
    new_g = [ f[0]*v0 + g[0], v0, f[0]*v1 + g[2], K.zero(), v1, K.one() ]
    # Subtotal : 0I 5M
    # Total : 1I 27M
  elif b2 != 0 :
    # Reduce matrix M
    # M_rref = [ 1  -u0  -v0  0  -w0 ]
    #          [ 0    0    0  1  -w1 ]
    # Find kernel of M
    #         [ u0 ]  [ v0 ]  [ w0 ]
    #         [  1 ]  [  0 ]  [  0 ]
    # ker M = [  0 ], [  1 ], [  0 ]
    #         [  0 ]  [  0 ]  [ w1 ]
    #         [  0 ]  [  0 ]  [  1 ]
    beta = 1/b2
    w1 = -beta*b3
    u0 = -a6
    v0 = -a7
    w0 = -(a9 + a8*w1)
    # Subtotal : 1I 2M
    
    # Find polynomials forming an ideal generating set for 2D
    # f'' = x^2 + ...
    # g'' =  xy + ...
    # h'' = x^3 + ...
    new_f = [ f[0]*u0, f[0] + u0, K.zero(), K.one() ]
    new_g = [ f[0]*v0, v0, f[0], K.zero(), K.one() ]
    # new_h = [ f[0]*w0 + g[0]*w1, w0, g[2]*w1, f[0], K.zero(), w1, K.one() ]
    # Subtotal : 0I 2M
    # Total : 1I 19M
  else :
    # Occurs if b1 = b2 = 0, but I don't think this is possible.
    raise ValueError("Cannot reduce matrix.")
  return C34CurveDivisor(D.C, [new_f, new_g, new_h])



def double_31(D) :
  C = D.C
  f0, f1, f2 = D.f[0:3]
  g0, g1, g2 = D.g[0:3]
  h0, h1, h2 = D.h[0:3]
  c0, c1, c2, c3, c4, c5, c6, c7, c8 = C.coefficients()

  if (D.typical) :
    # Compute two solutions
    #
    #   rf + sg + th = 0
    #   Rf + Sg + Th = C
    #
    # where
    # 
    #   r = y + r0
    #   s = -(x + s0)
    #   t = t0 = -f2
    #   R = x^2 + R2y + R1x + R0
    #   S = S0
    #   T = y + T1x + T0
    r0 = g1
    s0 = f1 - g2
    t0 = -f2

    T1 = c8
    R2 = c7 - f2
    R1 = c6 - f1
    T0 = c5 - h2 - f2*R2
    S0 = c4 - h2*T1 - h1 - f2*R1 - f1*R2
    R0 = c3 - h1*T1 - f1*R1 - f0
    # Subtotal : 0I 6M 12A

    # Compute
    #
    #   df = St - sT (mod f, g, h)
    #   dg = Tr - tR (mod f, g, h)
    #   dh = Rs - rS (mod f, g, h)
    df2 = s0 - g2 - T1*f2
    df1 = T0 - g1 + T1*(s0 - f1)
    df0 = T0*s0 + S0*t0 - T1*f0 - g0
    dg2 = T0 - h2 + r0 - T1*g2 + t0*(f2 - R2)
    dg1 = T1*(r0 - g1) + t0*(f1 - R1) - h1
    dg0 = T0*r0 - T1*g0 + t0*(f0 - R0) - h0
    dh2 = f2*(R1 - f1 - g2 + s0) + R2*(g2 - s0) - S0
    dh1 = f1*(R1 - f1 + s0) + g1*(R2 - f2) - R1*s0 - R0 + f0
    dh0 = f0*(R1 - f1 + s0) + g0*(R2 - f2) - S0*r0 - R0*s0
    # Subtotal : 0I 21M 40A
  
  elif (h1 != 0) :
    # Compute two solutions
    #
    #   rf + sg - th = 0
    #   Rf + Sg + Th = C
    #
    # where
    # 
    #   r =     r0 = h1
    #   s = y + s0 = y + h2 - g1
    #   t = x + t0 = x + g2
    #   R = x^2 + R2y + R1x + R0
    #   S = S0
    #   T = y + T1x + T0
    r0 = h1
    s0 = h2 - g1
    t0 = g2

    T1 = c8
    R2 = c7
    R1 = c6 - f1
    T0 = c5 - h2
    S0 = c4 - h2*T1 - h1 - f1*R2
    R0 = c3 - h1*T1 - f1*R1 - f0
    # Subtotal : 0I 4M 9A
    
    # Compute
    #
    #   df = -(sT + tS) (mod f, g, h)
    #   dg = rT + tR    (mod f, g, h)
    #   dh = sR - rS    (mod f, g, h)
    # ...
    # Subtotal : 0I 0M 0A
    df2 = T1*g2 - T0 + h2 - s0
    df1 = T1*(g1 - s0) - S0 + h1
    df0 = T1*g0 - T0*s0 - S0*t0 + h0
    dg3 = R1 - f1 + t0
    dg2 = R2*(t0 - g2) + r0
    dg1 = -R2*g1 + T1*r0 + R1*t0 + R0 - f0 - f1*dg3
    dg0 = -R2*g0 + T0*r0 + R0*t0 - f0*dg3
    dh4 = R1 - g2
    dh3 = s0 - g1
    dh2 = -R2*g1 + R0 - dh4*g2
    dh1 = -R2*h1 + R1*s0 - g0 - dh4*g1 - dh3*f1
    dh0 = -R2*h0 - S0*r0 + R0*s0 - dh4*g0 - dh3*f0

  else :
    # D = <f, g, h>, but
    # <f, g> =/= <f, g, h>
    # <g, h> =/= <f, g, h>
    # Fall back on slowed method.
    # Compute
    # df = fx*Cy - fy*Cx
    # dg = gx*Cy - gy*Cx
    # dh = hx*Cy - hy*Cx
    df4 = 2*c8*f1 + 4*c5 - 6*h2 - 4*c8*g2
    df3 = c7*f1 + 2*c4 - 4*c8*g1 - 2*c7*f1
    df2 = 2*c5*f1 - 3*f1*h2 - df4*g2
    df1 = c4*f1 + 2*c2 - 6*h0 - 4*c8*g0 - 2*c7*f0 - df4*g1 - df3*f1
    df0 = c2*f1 - 3*f1*h0 - df4*g0 - df3*f0
    
    dg6 = -6*c6 - 4*g2 + 7*f1
    dg5 = -c8*g2 - c5 + 3*g1
    dg4 = 2*c8*g1 - 2*c7*g2 - 3*c4 + 2*c8*h2 + 4*c7*g2
    dg3 = c7*g1 - 3*c6*g2 - 5*c3 + 7*f0 + 4*c7*g1 - dg6*f1
    dg2 = 2*c5*g1 - c4*g2 - 2*c2 - dg5*h2 - dg4*g2
    dg1 = c4*g1 - 2*c3*g2 - 4*c1 + 2*c8*h0 + 4*c7*g0 - dg6*f0 - dg5*h1 - dg4*g1 - dg3*f1
    dg0 = c2*g1 - c1*g2 - 3*c0 - dg5*h0 - dg4*g0 - dg3*f0

    dh8 = 2*c8^2 - 4*c7
    dh7 = 2*c7*c8 - 6*c6 + 8*f1
    dh6 = 2*c6*c8 - 4*h2 - 2*c8*f1
    dh5 = 2*c5*c8 - c8*h2 - 2*c4
    dh4 = 2*c4*c8 - 2*c7*h2 - 4*c3 + 8*f0 - dh8*h2 - dh7*g2
    dh3 = 2*c3*c8 - 3*c6*h2 - 2*c8*f0 - dh7*g1 - dh6*f1
    dh2 = 2*c2*c8 - c4*h2 - 2*c1 - dh5*h2 - dh4*g2
    dh1 = 2*c1*c8 - 2*c3*h2 - dh8*h0 - dh7*g0 - dh6*f0 - dh4*g1 - dh3*f1
    dh0 = 2*c0*c8 - c1*h2 - dh5*h0 - dh4*g0 - dh3*f0

  #       [ a1   a2   a3   a4   a5   a6   a7  ]
  #   M = [ a8   a9   a10  a11  a12  a13  a14 ]
  #       [ a15  a16  a17  a18  a19  a20  a21 ]
  #
  #   [ a4   a5   a6  ]   [ 0  -f0  -g0 ] [ a1   a2   a3  ]
  #   [ a11  a12  a13 ] = [ 1  -f1  -g1 ]*[ a8   a9   a10 ]
  #   [ a18  a19  a20 ]   [ 0  -f2  -g2 ] [ a15  a16  a17 ]
  #
  #   [ a7  ]   [ 0  -f0  -g0 ] [ a4  ]
  #   [ a14 ] = [ 1  -f1  -g1 ]*[ a11 ]
  #   [ a21 ]   [ 0  -f2  -g2 ] [ a18 ]
  a1  = df0
  a2  = dg0
  a3  = dh0
  a8  = df1
  a9  = dg1
  a10 = dh1
  a15 = df2
  a16 = dg2
  a17 = dh2

  a4  =    - f0*a8  - g0*a15
  a5  =    - f0*a9  - g0*a16
  a6  =    - f0*a10 - g0*a17
  a11 = a1 - f1*a8  - g1*a15
  a12 = a2 - f1*a9  - g1*a16
  a13 = a3 - f1*a10 - g1*a17
  a18 =    - f2*a8  - g2*a15
  a19 =    - f2*a9  - g2*a16
  a20 =    - f2*a10 - g2*a17

  a7  =    - f0*a11 - g0*a18
  a14 = a4 - f1*a11 - g1*a18
  a21 =    - f2*a11 - g2*a18
  # Subtotal : 0I 24M 16A

  if (a1 != 0) or (a8 != 0) or (a15 != 0) :
    if (a1 == 0) :
      if (a8 != 0) :
        a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14 = \
            a8, a9, a10, a11, a12, a13, a14, a1, a2, a3, a4, a5, a6, a7
      else :
        a1, a2, a3, a4, a5, a6, a7, a15, a16, a17, a18, a19, a20, a21 = \
            a15, a16, a17, a18, a19, a20, a21, a1, a2, a3, a4, a5, a6, a7
    #        [ a1  a2  a3  a4  a5   a6   a7  ]
    #   M' = [ 0   b1  b2  b3  b4   b5   b6  ]
    #        [ 0   b7  b8  b9  b10  b11  b12 ]
    b1  = a1*a9  - a2*a8
    b2  = a1*a10 - a3*a8
    b3  = a1*a11 - a4*a8
    b4  = a1*a12 - a5*a8
    b5  = a1*a13 - a6*a8
    b6  = a1*a14 - a7*a8
    b7  = a1*a16 - a2*a15
    b8  = a1*a17 - a3*a15
    b9  = a1*a18 - a4*a15
    b10 = a1*a19 - a5*a15
    b11 = a1*a20 - a6*a15
    b12 = a1*a21 - a7*a15
    # Subtotal : 0I 24M 12A

    if (b1 != 0) or (b7 != 0) :
      if (b1 == 0) :
        b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12 = \
            b7, b8, b9, b10, b11, b12, b1, b2, b3, b4, b5, b6
      #           [ a1  a2  a3  a4  a5  a6  a7 ]
      #   M_ref = [ 0   b1  b2  b3  b4  b5  b6 ]
      #           [ 0   0   c1  c2  c3  c4  c5 ]
      c1 = b1*b8  - b2*b7
      c2 = b1*b9  - b3*b7
      c3 = b1*b10 - b4*b7
      c4 = b1*b11 - b5*b7
      c5 = b1*b12 - b6*b7
      # Subtotal : 0I 10M 5A

      if (c1 != 0) :
        #            [ 1  0  0  -r0  -s0  -t0  * ]
        #   M_rref = [ 0  1  0  -r1  -s1  -t1  * ]
        #            [ 0  0  1  -r2  -s2  -t2  * ]
        ab = a1*b1
        abc = ab*c1
        delta = 1/abc
        alpha = delta*b1*c1
        beta = delta*a1*c1
        gamma = delta*ab

        r2 = -gamma*c2
        s2 = -gamma*c3
        t2 = -gamma*c4
        r1 = -beta*(b3 + b2*r2)
        s1 = -beta*(b4 + b2*s2)
        t1 = -beta*(b5 + b2*t2)
        r0 = -alpha*(a4 + a3*r2 + a2*r1)
        s0 = -alpha*(a5 + a3*s2 + a2*s1)
        t0 = -alpha*(a6 + a3*t2 + a2*t1)
        # Subtotal : 1I 25M 9A

        # 2D = <u, v, w>, where
        # 
        #   u = r0*f + r1*g + r2*h + xf
        #   v = s0*f + s1*g + s2*h + xg
        #   w = t0*f + t1*g + t2*h + xh
        u0 = f0*r0 + g0*r1 + h0*r2
        u1 = f1*r0 + g1*r1 + h1*r2 + f0
        u2 = f2*r0 + g2*r1 + h2*r2
        u3 = r0 + f1
        u4 = r1 + f2
        u5 = r2
        v0 = f0*s0 + g0*s1 + h0*s2
        v1 = f1*s0 + g1*s1 + h1*s2 + g0
        v2 = f2*s0 + g2*s1 + h2*s2
        v3 = s0 + g1
        v4 = s1 + g2
        v5 = s2
        w0 = f0*t0 + g0*t1 + h0*t2
        w1 = f1*t0 + g1*t1 + h1*t2 + h0
        w2 = f2*t0 + g2*t1 + h2*t2
        w3 = t0 + h1
        w4 = t1 + h2
        w5 = t2
        # Subtotal : 0I 27M 27A

        # 2D is of type 61
        # Total : 0I 137M 121A
        # Approx 40M unneeded
        return C34CurveDivisor(C, [[u0, u1, u2, u3, u4, u5, 1],
                             [v0, v1, v2, v3, v4, v5, 0, 1],
                             [w0, w1, w2, w3, w4, w5, 0, 0, 1]])
      elif (c2 != 0) :
        #           [ a1  a2  a3  a4  a5  a6  a7 ]
        #   M_ref = [ 0   b1  b2  b3  b4  b5  b6 ]
        #           [ 0   0   0   c2  c3  c4  c5 ]
        #
        #            [ 1  0  -r0  0  -s0  *  * ]
        #   M_rref = [ 0  1  -r1  0  -s1  *  * ]
        #            [ 0  0   0   1  -s2  *  * ]
        ab = a1*b1
        abc = ab*c2
        delta = 1/abc
        alpha = delta*b1*c2
        beta = delta*a1*c2
        gamma = delta*ab
        s2 = -gamma*c3
        r1 = -beta*b2
        s1 = -beta*(b4 + b3*s2)
        r0 = -alpha*(a3 + a2*r1)
        s0 = -alpha*(a5 + a4*s2 + a2*s1)

        #   u = r0*f + r1*g + h
        #   v = s0*f + s1*g + s2*xf + xg
        u0 = f0*r0 + g0*r1 + h0
        u1 = f1*r0 + g1*r1 + h1
        u2 = f2*r0 + g2*r1 + h2
        u3 = r0
        u4 = r1
        v0 = f0*s0 + g0*s1
        v1 = f1*s0 + g1*s1 + f0*s2 + g0
        v2 = f2*s0 + g2*s1
        v3 = s0 + f1*s2 + g1
        v4 = s1 + f2*s2 + g2
        v6 = s2

        # 2D is of type 63
        return C34CurveDivisor(C, [[u0, u1, u2, u3, u4, 1],
                             [v0, v1, v2, v3, v4, 0, v6, 1]])
      
      elif (c3 != 0) :
        #           [ a1  a2  a3  a4  a5  a6  a7 ]
        #   M_ref = [ 0   b1  b2  b3  b4  b5  b6 ]
        #           [ 0   0   0   0   c3  c4  c5 ]
        #
        #            [ 1  0  -r0  -s0  0  *  * ]
        #   M_rref = [ 0  1  -r1  -s1  0  *  * ]
        #            [ 0  0   0    0   1  *  * ]
        gamma = 1/(a1*b1)
        alpha = gamma*b1
        beta = gamma*a1
        r1 = -beta*b2
        s1 = -beta*b3
        r0 = -alpha*(a3 + a2*r1)
        s0 = -alpha*(a4 + a2*s1)

        #   u = r0*f + r1*g + h
        #   v = s0*f + s1*g + xf
        u0 = f0*r0 + g0*r1 + h0
        u1 = f1*r0 + g1*r1 + h1
        u2 = f2*r0 + g2*r1 + h2
        u3 = r0
        u4 = r1
        v0 = f0*s0 + g0*s1
        v1 = f1*s0 + g1*s1 + f0
        v2 = f2*s0 + g2*s1
        v3 = s0 + f1
        v4 = s1 + f2

        # 2D is of type 62
        return C34CurveDivisor(C, [[u0, u1, u2, u3, u4, 1],
                             [v0, v1, v2, v3, v4, 0, 1]])
      else :
        #           [ a1  a2  a3  a4  a5  a6  a7 ]
        #   M_ref = [ 0   b1  b2  b3  b4  b5  b6 ]
        #           [ 0   0   0   0   0   0   0  ]
        # 
        #            [ 1  0  -r0  -s0  -t0  *  * ]
        #   M_rref = [ 0  1  -r1  -s1  -t1  *  * ]
        #            [ 0  0   0    0    0   0  0 ]
        gamma = 1/(a1*b1)
        alpha = gamma*b1
        beta = gamma*a1
        r1 = -beta*b2
        s1 = -beta*b3
        t1 = -beta*b4
        r0 = -alpha*(a3 + a2*r1)
        s0 = -alpha*(a4 + a2*s1)
        t0 = -alpha*(a5 + a2*t1)

        # TODO : Experimentally, it seems like this situation arises only when f2 = 0.
        #   u = r0*f + r1*g + h
        #   v = s0*f + s1*g + xf
        u0 = f0*r0 + g0*r1 + h0
        u1 = f1*r0 + g1*r1 + h1
        u2 = f2*r0 + g2*r1 + h2
        u3 = r0
        u4 = r1
        v0 = f0*s0 + g0*s1
        v1 = f1*s0 + g1*s1 + f0
        v2 = f2*s0 + g2*s1
        v3 = s0 + f1
        v4 = s1 + f2
        w0 = f0*t0 + g0*t1
        w1 = f1*t0 + g1*t1 + g0
        w2 = f2*t0 + g2*t1
        w3 = t0 + g1
        w4 = t1 + g2

        L = C34CurveDivisor(C, [[u0, u1, u2, u3, u4, 1],
                          [v0, v1, v2, v3, v4, 0, 1],
                          [w0, w1, w2, w3, w4, 0, 0, 1]])
        
        # G is type 11.
        # A basis for G is given by the 1st and 2nd columns of M
        #
        # [ df0  dg0 ]     [ n2 m3 ]     [ p0  q0 ]
        # [ df1  dg1 ] ==> [ n1 m2 ] ==> [ 1   0  ]
        # [ df2  dg2 ]     [ 0  m1 ]     [ 0   1  ]
        if (dg2 != 0) :
          m1, m2, m3 = dg2, dg1, dg0
          n1 = dg2*df1 - dg1*df2
          n2 = dg2*df0 - dg0*df2
        else :
          m1, m2, m3 = df2, df1, df0
          n1 = df2*dg1 - df1*dg2
          n2 = df2*dg0 - df0*dg2
        mu = 1/(m1*n1)
        p0 = mu*m1*n2
        q0 = mu*n1*(m3 - m2*p0)
        G = C34CurveDivisor(C, [[p0, 1], [q0, 0, 1]])
        
        return flip(flip(L)) + G

    elif (b2 != 0) or (b8 != 0) :
      if (b2 == 0) :
        b2, b3, b4, b5, b6, b8, b9, b10, b11, b12 = \
            b8, b9, b10, b11, b12, b2, b3, b4, b5, b6
      #        [ a1  a2  a3  a4  a5   a6   a7  ]
      #   M' = [ 0   0   b2  b3  b4   b5   b6  ]
      #        [ 0   0   b8  b9  b10  b11  b12 ]
      #
      #           [ a1  a2  a3  a4  a5  a6  a7 ]
      #   M_ref = [ 0   0   b2  b3  b4  b5  b6 ]
      #           [ 0   0   0   c1  c2  c3  c4 ]
      c1 = b2*b9  - b3*b8
      c2 = b2*b10 - b4*b8
      c3 = b2*b11 - b5*b8
      c4 = b2*b12 - b6*b8

      if (c1 != 0) :
        #            [ 1  -r0  0  0  *  *  -s0 ]
        #   M_rref = [ 0   0   1  0  *  *  -s1 ]
        #            [ 0   0   0  1  *  *  -s2 ]
        ab = a1*b2
        abc = ab*c1
        delta = 1/abc
        gamma = delta*ab
        beta = delta*a1*c1
        alpha = delta*b2*c1
        s2 = -gamma*c4
        s1 = -beta*(b6 + b3*s2)
        s0 = -alpha*(a7 + a3*s1 + a4*s2)
        r0 = -alpha*a2

        #   u = r0*f + g
        #   v = s0*f + s1*h + s2*xf + x^2f - f2*xu - f2*(s2 - u2)*u
        u0 = f0*r0 + g0
        u1 = f1*r0 + g1
        u2 = f2*r0 + g2
        u3 = r0
        
        z  = f2*(s2 - u2)
        v0 = f0*s0 + h0*s1 - z*u0
        v1 = f1*s0 + h1*s1 + f0*s2 - f2*u0 - z*u1
        v2 = f2*s0 + h2*s1 - z*u2
        v3 = s0 + f1*s2 + f0 - f2*u1 - z*u3
        v5 = s1
        v6 = s2 + f1 - f2*u3

        return C34CurveDivisor(C, [[u0, u1, u2, u3, 1],
                             [v0, v1, v2, v3, 0, v5, v6, 0, 0, 1]])
      
      else :
        #            [ 1  -r0  0  -s0  *  *  * ]
        #   M_rref = [ 0   0   1  -s1  *  *  * ]
        #            [ 0   0   0   0   0  0  0 ]
        gamma = 1/a1*b2
        beta = gamma*a1
        alpha = gamma*b2
        s1 = -beta*b3
        s0 = -alpha*(a4 + a3*s1)
        r0 = -alpha*a2

        # Compute G = <u,v>, where
        # 
        #   u = r0*f + g
        #   v = s0*f + s1*h + xf
        #
        # TODO: Experimentally, it seems this case only occurs when f2 = 0
        u0 = f0*r0 + g0
        u1 = f1*r0 + g1
        u2 = f2*r0 + g2
        u3 = r0
        v0 = f0*s0 + h0*s1 - f2*u0
        v1 = f1*s0 + h1*s1 + f0 - f2*u1
        v2 = f2*(s0 - u2) + h2*s1
        v3 = s0 + f1 - f2*u3
        v5 = s1
        L = C34CurveDivisor(C, [[u0, u1, u2, u3, 1],
                          [v0, v1, v2, v3, 0, v5, 1]])

        # G is type 11.
        # A basis for G is given by the 1st and 2nd columns of M
        #
        # [ df0  dh0 ]     [ n2 m3 ]     [ p0  q0 ]
        # [ df1  dh1 ] ==> [ n1 m2 ] ==> [ 1   0  ]
        # [ df2  dh2 ]     [ 0  m1 ]     [ 0   1  ]
        if (dh2 != 0) :
          m1, m2, m3 = dh2, dh1, dh0
          n1 = dh2*df1 - dh1*df2
          n2 = dh2*df0 - dh0*df2
        else :
          m1, m2, m3 = df2, df1, df0
          n1 = df2*dh1 - df1*dh2
          n2 = df2*dh0 - df0*dh2
        mu = 1/(m1*n1)
        p0 = mu*m1*n2
        q0 = mu*n1*(m3 - m2*p0)
        G = C34CurveDivisor(C, [[p0, 1], [q0, 0, 1]])

        return flip(flip(L)) + G

    elif (b3 != 0) or (b9 != 0) :
      #           [ a1  a2  a3  a4  a5  a6  a7 ]
      #   M_ref = [ 0   0   0   b3  b4  b5  b6 ]
      #           [ 0   0   0   0   0   0   0  ]
      #
      #            [ 1  -r0  -s0  0  *  *  * ]
      #   M_rref = [ 0   0    0   1  *  *  * ]
      #            [ 0   0    0   0  0  0  0 ]
      alpha = 1/a1
      r0 = -alpha*a2
      s0 = -alpha*a3
      u0 = f0*r0 + g0
      u1 = f1*r0 + g1
      u2 = f2*r0 + g2
      u3 = r0
      v0 = f0*s0 + h0
      v1 = f1*s0 + h1
      v2 = f2*s0 + h2
      v3 = s0
      L = C34CurveDivisor(C, [[u0, u1, u2, u3, 1], [v0, v1, v2, v3, 0, 1], []])

      # G is type 11.
      # A basis for G is given by the 1st and 4th columns of M
      #
      # [ df0  dxf0 ]     [ n2 m3 ]     [ p0  q0 ]
      # [ df1  dxf1 ] ==> [ n1 m2 ] ==> [ 1   0  ]
      # [ df2  dxf2 ]     [ 0  m1 ]     [ 0   1  ]
      if (dxf2 != 0) :
        m1, m2, m3 = dxf2, dxf1, dxf0
        n1 = dxf2*df1 - dxf1*df2
        n2 = dxf2*df0 - dxf0*df2
      else :
        m1, m2, m3 = df2, df1, df0
        n1 = df2*dxf1 - df1*dxf2
        n2 = df2*dxf0 - df0*dxf2
      mu = 1/(m1*n1)
      p0 = mu*m1*n2
      q0 = mu*n1*(m3 - m2*p0)
      G = C34CurveDivisor(C, [[p0, 1], [q0, 0, 1]])

      return flip(flip(L)) + G

    else :
      #           [ 1  -r0  -s0  -t0  *  *  * ]
      #   M_ref = [ 0   0    0    0   0  0  0 ]
      #           [ 0   0    0    0   0  0  0 ]
      alpha = 1/a1
      r0 = -alpha*a2
      s0 = -alpha*a3
      t0 = -alpha*a4
      u0 = f0*r0 + g0
      u1 = f1*r0 + g1
      u2 = f2*r0 + g2
      u3 = r0
      v0 = f0*s0 + h0
      v1 = f1*s0 + h1
      v2 = f2*s0 + h2
      v3 = s0
      w0 = f0*t0 - f2*u0
      w1 = f1*t0 + f0 - f2*u1
      w2 = f2*(t0 - u2)
      w3 = t0 + f1 - f2*u3
      
      L = C34CurveDivisor(C, [[u0, u1, u2, u3, 1], [v0, v1, v2, v3, 0, 1], [w0, w1, w2, w3, 0, 0, 1]])

      if (df2 != 0) :
        mu = 1/df2
        p0 = mu*df0
        p1 = mu*df1
        G = C34CurveDivisor(C, [[p0, p1, 1], copy(D.f), []])
      else :
        p0 = df0/df1
        G = C34CurveDivisor(C, [[p0, 1], [h0 - h1*p0, 0, h2, 0, 0, 1], []])

      return flip(flip(L)) + G

  elif (a2 != 0) or (a9 != 0) or (a16 != 0) :
    if (a2 == 0) :
      if (a9 != 0) :
        a2, a3, a5, a6, a9, a10, a12, a13 = \
            a9, a10, a12, a13, a2, a3, a5, a6
      else :
        a2, a3, a5, a6, a16, a17, a19, a20 = \
            a16, a17, a19, a20, a2, a3, a5, a6
    #       [ 0  a2   a3   0  a5   a6   0 ]
    #   M = [ 0  a9   a10  0  a12  a13  0 ]
    #       [ 0  a16  a17  0  a19  a20  0 ]
    #
    # Compute
    #
    #        [ 0  a2  a3  0  a5  a6  0 ]
    #   M' = [ 0  0   b1  0  b2  b3  0 ]
    #        [ 0  0   b4  0  b5  b6  0 ]
    b1 = a2*a10 - a3*a9
    b2 = a2*a12 - a5*a9
    b3 = a2*a13 - a6*a9
    b4 = a2*a17 - a3*a16
    b5 = a2*a19 - a5*a16
    b6 = a2*a20 - a6*a16
    
    if (b1 != 0) or (b4 != 0) :
      if (b1 == 0) :
        b1, b2, b3, b4, b5, b6 = b4, b5, b6, b1, b2, b3

      #           [ 0  a2  a3  0  a5  a6  0 ]
      #   M_ref = [ 0  0   b1  0  b2  b3  0 ]
      #           [ 0  0   0   0  c1  c2  0 ]
      c1 = b1*b5 - b2*b4
      c2 = b1*b6 - b3*b4

      assert c1 == 0, "This case is impossible!"
      if (c2 != 0) :
        #            [ 0  1  0  0  *  0  0 ]
        #   M_rref = [ 0  0  1  0  *  0  0 ]
        #            [ 0  0  0  0  0  1  0 ]
        return C34CurveDivisor(C, [copy(D.f), [], []])
      else :
        #            [ 0  1  0  0  *  -r0  0 ]
        #   M_rref = [ 0  0  1  0  *  -r1  0 ]
        #            [ 0  0  0  0  0   0   0 ]
        gamma = 1/(a2*b1)
        beta = gamma*a2
        alpha = gamma*b1
        r1 = -beta*b3
        r0 = -alpha*(a6 + a3*r1)

        # Compute L = <f, v>, where
        #
        # v = r0*g + r1*h + xh
        v0 = g0*r0 + h0*r1 - h1*f0
        v1 = g1*r0 + h1*(r1 - f1) + h0
        v2 = g2*r0 + h2*r1 - h1*f2
        v4 = r0 + h2
        v5 = r1
        L = C34CurveDivisor(C, [copy(D.f), [v0, v1, v2, 0, v4, v5, 0, 0, 1], []])

        # G is type 11.
        # A basis for G is given by the 1st and 2nd columns of M
        #
        # [ dg0  dh0 ]     [ n2 m3 ]     [ p0  q0 ]
        # [ dg1  dh1 ] ==> [ n1 m2 ] ==> [ 1   0  ]
        # [ dg2  dh2 ]     [ 0  m1 ]     [ 0   1  ]
        if (dh2 != 0) :
          m1, m2, m3 = dh2, dh1, dh0
          n1 = dh2*dg1 - dh1*dg2
          n2 = dh2*dg0 - dh0*dg2
        else :
          m1, m2, m3 = dg2, dg1, dg0
          n1 = dg2*dh1 - dg1*dh2
          n2 = dg2*dh0 - dg0*dh2
        mu = 1/(m1*n1)
        p0 = mu*m1*n2
        q0 = mu*n1*(m3 - m2*p0)
        G = C34CurveDivisor(C, [[p0, 1], [q0, 0, 1]])

        return flip(flip(L)) + G
    else :
      #        [ 0  a2  a3  0  a5  a6  0 ]
      #   M' = [ 0  0   0   0  b2  b3  0 ]
      #        [ 0  0   0   0  b5  b6  0 ]
      #
      #            [ 0  1  -r0  0  *  *  0 ]
      #   M_rref = [ 0  0   0   0  0  0  0 ]
      #            [ 0  0   0   0  0  0  0 ]
      r0 = -a3/a2
      v0 = g0*r0 + h0
      v1 = g1*r0 + h1
      v2 = g2*r0 + h2
      v4 = r0
      L = C34CurveDivisor(C, [copy(D.f), [v0, v1, v2, 0, v4, 1], []])

    if (dg2 != 0) :
      mu = 1/dg2
      p0 = mu*dg0
      p1 = mu*dg1
      G = C34CurveDivisor(C, [[p0, p1, 1], copy(D.f), []])
    else :
      p0 = dg0/dg1
      G = C34CurveDivisor(C, [[p0, 1], [h0 - h1*p0, 0, h2, 0, 0, 1], []])
    return flip(flip(L)) + G
  
  else :
    #       [ 0  0  a3   0  0  a6   0 ]
    #   M = [ 0  0  a10  0  0  a13  0 ]
    #       [ 0  0  a17  0  0  a20  0 ]
    L = C34CurveDivisor(C, [copy(D.f), copy(D.g), []])

    if (dh2 != 0) :
      mu = 1/dh2
      p0 = mu*dh0
      p1 = mu*dh1
      G = C34CurveDivisor(C, [[p0, p1, 1], copy(D.f), []])
    else :
      p0 = dh0/dh1
      G = C34CurveDivisor(C, [[p0, 1], [h0 - h1*p0, 0, h2, 0, 0, 1], []])
    return flip(flip(L)) + G



def old_double_31(D):
  K = D.K
  f, g, h = D.f, D.g, D.h
  c = D.C.coefficients()
  new_f, new_g, new_h = [], [], []

  if D.typical :
    # Find polynomials G, H satisfying fH + gG = C
    # Find also 1/f[2]
    A = flip(D) # Costs 1I 17M
    G = A.g
    if G == g :
      # If A.g = D.g,
      # then A = D and 2D = D + A = 0.
      # (More precisely, 2D = <D.f> is principal)
      return C34CurveDivisor(D.C, [[f[0], f[1], K.one()], [], []])
      #return C34CurveDivisor(D.C, [[K.one()], [], []])
    H = [ 0,
          f[2]*(f[1] - c[6]),
          f[2]*(f[2] - c[7]) - g[1] - G[1],
          -f[2],
          K.zero(), -K.one() ]
    H[0] = f[2]*(f[0] - c[3]) - f[1]*H[1] - g[1]*G[1]
    ainv = 1/f[2] # f[2] is guaranteed to be non-zero.
    # Subtotal : 2I 22M
    # Note : We are inverting the same number twice.
    #        We can save 1I by recording 1/f[2] when computing flip(D)
    
    # Construct matrix M
    # M = [ a1   a2   a3   a4   a5   a6   a7  ]
    #     [ a8   a9   a10  a11  a12  a13  a14 ]
    #     [ a15  a16  a17  a18  a19  a20  a21 ]
    a1  = G[0] - g[0]
    a8  = G[1] - g[1]
    a15 = G[2] - g[2]
    a2  = -H[0] - h[0] + H[3]*f[0]
    a9  = -H[1] - h[1] + H[3]*f[1]
    a16 = -H[2] - h[2] + H[3]*f[2]
    a4  =    - f[0]*a8 - g[0]*a15
    a11 = a1 - f[1]*a8 - g[1]*a15
    a18 =    - f[2]*a8 - g[2]*a15
    a5  =    - f[0]*a9 - g[0]*a16
    a12 = a2 - f[1]*a9 - g[1]*a16
    a19 =    - f[2]*a9 - g[2]*a16
    # a3  = ainv*(    - g[0]*a8 - h[0]*a15 + g[1]*a1  - a5  - (f[1] - g[2])*a2  )
    # a10 = ainv*(    - g[1]*a8 - h[1]*a15 + g[1]*a8  - a12 - (f[1] - g[2])*a9  )
    # a17 = ainv*( a1 - g[2]*a8 - h[2]*a15 + g[1]*a15 - a19 - (f[1] - g[2])*a16 )
    a3  = ainv*(    - g[0]*a8 - h[0]*a15 + g[1]*a1  - a5  - (f[1] - g[2])*a2  )
    a10 = ainv*(              - h[1]*a15            - a12 - (f[1] - g[2])*a9  )
    a17 = ainv*( a1 - g[2]*a8 + (g[1] - h[2])*a15   - a19 - (f[1] - g[2])*a16 )
    a6  =    - f[0]*a10 - g[0]*a17
    a13 = a3 - f[1]*a10 - g[1]*a17
    a20 =    - f[2]*a10 - g[2]*a17
    a7  =    - f[0]*a11 - g[0]*a18
    a14 = a4 - f[1]*a11 - g[1]*a18
    a21 =    - f[2]*a11 - g[2]*a18
    # Subtotal : 0I 33M
    
    # Since G != g, at least one of a1, a8, a15 is non-zero
    # Swap rows so that the top-left element of M is non-zero.
    if a1 == 0 :
      if a8 != 0 :
        a1, a2, a3, a4, a5, a6, a7, a8,  a9,  a10, a11, a12, a13, a14 = a8,  a9,  a10, a11, a12, a13, a14, a1, a2, a3, a4, a5, a6, a7
      else : # a13 != 0
        a1, a2, a3, a4, a5, a6, a7, a15, a16, a17, a18, a19, a20, a21 = a15, a16, a17, a18, a19, a20, a21, a1, a2, a3, a4, a5, a6, a7
    
    # Partially reduce M to the form
    #   [ a1  a2  a3  a4  a5  a6  a7  ]
    #   [  0  b1  b2  b3  b4  b5  b6  ]
    #   [  0  b7  b8  b9  b10 b11 b12 ]
    b1 = a1*a9  - a2*a8
    b2 = a1*a10 - a3*a8
    b3 = a1*a11 - a4*a8
    b4 = a1*a12 - a5*a8
    b5 = a1*a13 - a6*a8
    b6 = a1*a14 - a7*a8
    b7 = a1*a16 - a2*a15
    b8 = a1*a17 - a3*a15
    b9 = a1*a18 - a4*a15
    b10 = a1*a19 - a5*a15
    b11 = a1*a20 - a6*a15
    b12 = a1*a21 - a7*a15

    if (b1 == 0) and (b7 == 0) :
      # TODO : Handle this
      # XXX : I expect that one of b2 or b8 is non-zero
      if b2 == 0 :
        b2, b3, b4, b5, b6, b8, b9, b10, b11, b12 = b8, b9, b10, b11, b12, b2, b3, b4, b5, b6
      # Reduce M to row echelon form
      #           [ a1  a2  a3  a4  *  *  a7 ]
      #   M_ref = [  0   0  b2  b3  *  *  b6 ]
      #           [  0   0   0  c2  *  *  c5 ]
      c2 = b2*b9  - b3*b8
      c5 = b2*b12 - b6*b8

      # Reduce to RREF
      #            [ 1  -r0  0  0  *  *  -s0 ]
      #   M_rref = [ 0   0   1  0  *  *  -s1 ]
      #            [ 0   0   0  1  *  *  -s2 ]
      ab = a1*b2
      abc = ab*c2
      delta = 1/(abc)
      alpha = delta*b2*c2
      beta  = delta*a1*c2
      gamma = delta*ab
      s2 = -gamma*c5
      s1 = -beta*(b6 + b3*s2)
      s0 = -alpha*(a7 + a4*s2 + a3*s1)
      r0 = -alpha*a2

      # u = r0*f + g
      # v = s0*f + s1*h + s2*xf + x^2f - f[2]*xu - f[2]*(s2 - u2)*u
      u0 = f[0]*r0 + g[0]
      u1 = f[1]*r0 + g[1]
      u2 = f[2]*r0 + g[2]
      u3 = r0
      
      z = f[2]*(s2 - u2)
      v0 = f[0]*s0 + h[0]*s1 - z*u0
      v1 = f[1]*s0 + h[1]*s1 + f[0]*s2 - f[2]*u0 - z*u1
      v2 = f[2]*s0 + h[2]*s1 - z*u2
      v3 = s0 + f[1]*s2 + f[0] - f[2]*u1 - z*u3
      v5 = s1
      v6 = s2 + f[1] - f[2]*u3

      return C34CurveDivisor(D.C, [[u0, u1, u2, u3, 1], [v0, v1, v2, v3, 0, v5, v6, 0, 0, 1], []])
    
    if b1 == 0 :
      b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12 = b7, b8, b9, b10, b11, b12, b1, b2, b3, b4, b5, b6

    # Reduce M to its reduced row echelon form
    #           [ a1  a2  a3  a4  a5  a6  a7 ]
    #   M_ref = [  0  b1  b2  b3  b4  b5  b6 ]
    #           [  0   0  c1  c2  c3  c4  c5 ]
    c1 = b1*b8  - b2*b7
    c2 = b1*b9  - b3*b7
    c3 = b1*b10 - b4*b7
    c4 = b1*b11 - b5*b7
    c5 = b1*b12 - b6*b7

    if c1 != 0 :
      # Reduce matrix M
      # (Omit last column -- it is not needed in this case)
      # M_rref = [ 1  0  0  -u0  -v0  -w0 ]
      #          [ 0  1  0  -u1  -v1  -w1 ]
      #          [ 0  0  1  -u2  -v2  -w2 ]
      # Find kernel of M
      #         [ u0 ]  [ v0 ]  [ w0 ]
      #         [ u1 ]  [ v1 ]  [ w1 ]
      # ker M = [ u2 ], [ v2 ], [ w2 ]
      #         [  1 ], [  0 ], [  0 ]
      #         [  0 ]  [  1 ]  [  0 ]
      #         [  0 ]  [  0 ]  [  1 ]
      delta = 1/(a1*b1*c1)
      alpha = delta*b1*c1
      beta  = delta*a1*c1
      gamma = delta*a1*b1
      u2 = -gamma*c2
      v2 = -gamma*c3
      w2 = -gamma*c4
      u1 = -beta*(b3 + b2*u2)
      v1 = -beta*(b4 + b2*v2)
      w1 = -beta*(b5 + b2*w2)
      u0 = -alpha*(a4 + a3*u2 + a2*u1)
      v0 = -alpha*(a5 + a3*v2 + a2*v1)
      w0 = -alpha*(a6 + a3*w2 + a2*w1)
      # Subtotal : 1I 54M

      # Find polynomials forming an ideal generating set for 2D
      new_f = [f[0]*u0 + g[0]*u1 + h[0]*u2,
               f[1]*u0 + g[1]*u1 + h[1]*u2 + f[0],
               f[2]*u0 + g[2]*u1 + h[2]*u2,
               u0 + f[1],
               u1 + f[2],
               u2,
               K.one()]
      new_g = [f[0]*v0 + g[0]*v1 + h[0]*v2,
               f[1]*v0 + g[1]*v1 + h[1]*v2 + g[0],
               f[2]*v0 + g[2]*v1 + h[2]*v2,
               v0 + g[1],
               v1 + g[2],
               v2,
               K.zero(), K.one()]
      new_h = [f[0]*w0 + g[0]*w1 + h[0]*w2,
               f[1]*w0 + g[1]*w1 + h[1]*w2 + h[0],
               f[2]*w0 + g[2]*w1 + h[2]*w2,
               w0 + h[1],
               w1 + h[2],
               w2,
               K.zero(), K.zero(), K.one()]
      # Subtotal : 0I 27M
      # Total : 3I 136M
      # Note : If 1/f[2] is precomputed, this saves 2I.
    elif c2 != 0 :
      # Reduce matrix M
      # M_rref = [ 1  0  -u0  0  -v0  -w0 ]
      #          [ 0  1  -u1  0  -v1  -w1 ]
      #          [ 0  0    0  1  -v2  -w2 ]
      # Find kernel of M
      #         [ u0 ]  [ v0 ]  [ w0 ]
      #         [ u1 ]  [ v1 ]  [ w1 ]
      # ker M = [  1 ], [  0 ], [  0 ]
      #         [  0 ], [ v2 ], [ w2 ]
      #         [  0 ]  [  1 ]  [  0 ]
      #         [  0 ]  [  0 ]  [  1 ]
      delta = 1/(a1*b1*c2)
      alpha = delta*b1*c2
      beta  = delta*a1*c2
      gamma = delta*a1*b1
      v2 = -gamma*c3
      w2 = -gamma*c4
      u1 = -beta*b2
      v1 = -beta*(b4 + b3*v2)
      w1 = -beta*(b5 + b3*w2)
      u0 = -alpha*(a3 + a2*u1)
      v0 = -alpha*(a5 + a4*v2 + a2*v1)
      w0 = -alpha*(a6 + a4*w2 + a2*w1)
      # Subtotal : ?I ??M
      
      new_f = [ f[0]*u0 + g[0]*u1 + h[0],
                f[1]*u0 + g[1]*u1 + h[1],
                f[2]*u0 + g[2]*u1 + h[2],
                u0,
                u1,
                K.one() ]
      new_g = [ f[0]*v0 + g[0]*v1,
                f[1]*v0 + g[1]*v1 + f[0]*v2 + g[0],
                f[2]*v0 + g[2]*v1,
                v0 + f[1]*v2 + g[1],
                v1 + f[2]*v2 + g[2],
                K.zero(), v2, K.one() ]
      # new_h = [ f[0]*w0 + g[0]*w1,
      #           f[1]*w0 + g[1]*w1 + f[0]*w2 + h[0],
      #           f[2]*w0 + g[2]*w1,
      #           w0 + f[1]*w2 + h[1],
      #           w1 + f[2]*w2 + h[2],
      #           K.zero(), w2, K.zero(), K.one() ]
      # Subtotal : ?I ??M
      # Total : ?I ??M
    elif c3 != 0 :
      # Reduce matrix M
      # M_rref = [ 1  0  -u0  -v0  0  -w0 ]
      #          [ 0  1  -u1  -v1  0  -w1 ]
      #          [ 0  0    0    0  1  -w2 ]
      # Find kernel of M
      #         [ u0 ]  [ v0 ]  [ w0 ]
      #         [ u1 ]  [ v1 ]  [ w1 ]
      # ker M = [  1 ], [  0 ], [  0 ]
      #         [  0 ], [  1 ], [  0 ]
      #         [  0 ]  [  0 ]  [ w2 ]
      #         [  0 ]  [  0 ]  [  1 ]
      delta = 1/(a1*b1*c3)
      alpha = delta*b1*c3
      beta  = delta*a1*c3
      gamma = delta*a1*b1
      w2 = -gamma*c4
      u1 = -beta*b2
      v1 = -beta*b3
      w1 = -beta*(b5 + b4*w2)
      u0 = -alpha*(a3 + a2*u1)
      v0 = -alpha*(a4 + a2*v1)
      w0 = -alpha*(a6 + a5*w2 + a2*w1)
      # Subtotal : ?I ??M
      
      new_f = [ f[0]*u0 + g[0]*u1 + h[0],
                f[1]*u0 + g[1]*u1 + h[1],
                f[2]*u0 + g[2]*u1 + h[2],
                u0,
                u1,
                K.one() ]
      new_g = [ f[0]*v0 + g[0]*v1,
                f[1]*v0 + g[1]*v1 + f[0],
                f[2]*v0 + g[2]*v1,
                v0 + f[1],
                v1 + f[2],
                K.zero(), K.one() ]
      # new_h = [ f[0]*w0 + g[0]*w1,
      #           f[1]*w0 + g[1]*w1 + f[0]*w2 + h[0],
      #           f[2]*w0 + g[2]*w1,
      #           w0 + f[1]*w2 + h[1],
      #           w1 + f[2]*w2 + h[2],
      #           K.zero(), K.zero(), w2, K.one() ]
      # Subtotal : ?I ??M
      # Total : ?I ??M
      
  else : 
    # Split D into two divisors, A and B, of degrees 1 and 2, respectively.
    # XXX : This way is much more costly than expected!
    x2 = -g[2]
    y1 = -g[1]
    x1 = g[2] - f[1]
    
    fA = [-x1, K.one()]
    gA = [-y1, K.zero(), K.one()]
    fB = [-x2, K.one()]
    gB = [h[0] + h[1]*x2, K.zero(), h[2], K.zero(), K.zero(), K.one()]
    A = C34CurveDivisor(D.C, [fA, gA, []])
    B = C34CurveDivisor(D.C, [fB, gB, []])
    # Subtotal : 0I 1M

    AA = double(A)     # Cost : 1I 16M
    BB = double(B)     # Cost : 1I 24M
    fAA = flip(AA)     # Cost : 0I  7M (worst case) (usually 0I 1M)
    fBB = flip(BB)     # Cost : 1I 24M (if BB is typical) (probably not)
    fAAfBB = fAA + fBB # Cost : 1I 41M (worst case)
    DD = flip(fAAfBB)  # Cost : 1I 64M (Can be improved)
    # Subtotal : 5I 176M

    new_f, new_g, new_h = DD.f, DD.g, DD.h
    # Total : 5I 177M
    # Note : This is way more multiplications that I expected.
    # XXX : This count may not be accurate

  return C34CurveDivisor(D.C, [new_f, new_g, new_h])

