"""
  Summary of costs for doubling a divisor
  
  
  deg | I |   M | S |   A |
  ----+---+-----+---+-----+
    0 | 0 |   0 |   |     |
  ----+---+-----+---+-----+
    1 | 1 |  15 | 1 |  15 | When dy/dx =/= 0
      | 0 |   8 | 1 |   8 | 
  ----+---+-----+---+-----+
    2 | 1 |  63 |   |     | Distinct x-coordinates
      | 1 |  24 |   |     | Same x-coordinate (this can be brought down to 17M)
  ----+---+-----+---+-----+
    3 | 3 | 136 |   |     | Typical
    3 | 5 | 177 |   |     | Semi-typical

  Type
    D | 2D | I |   M | S |   A |
  ----+----+---+-----+---+-----+
   11 | 21 | 1 |  15 | 1 |  20 |
      | 22 | 0 |   8 | 1 |  13 |
  ----+----+---+-----+---+-----+
   21 | 41 | 1 |  63 | 0 |  57 |
      | 42 | 1 |  50 | 0 |  46 |
      | 43 | 1 |  54 | 0 |  49 |
      | 44 | 0 |   7 | 0 |  12 |
  ----+----+---+-----+---+-----+
   22 | 42 | 1 |  22 | 0 |  21 |
      | 43 | 1 |  19 | 0 |  18 |
  ----+----+---+-----+---+-----+
   31 | 61 | 1 | 135 | 3 | 116 | fast_double_31_high_char
      | 61 | 1 | 145 | 2 | 126 | fast_double_31
      | 61 | 1 | 127 | 4 | 112 | faster_double_31_high_char
      | 61 | 1 | 138 | 2 | 130 | faster_double_31
      | 61 | 1 | 124 | 0 | 110 |
      | 62 | 1 |  92 | 0 |  85 |
      | 63 | 1 | 102 | 0 |  90 |
      | 64 | 1 | 107 | 0 |  90 |
      | 65 | 0 |  57 | 0 |  64 |
  ----+----+---+-----+---+-----+
"""

def double(D) :
  """
    Double a divisor D.
    
    Input: A reduced C34CurveDivisor D.
    Output: The unique reduced C34CurveDivisor D' equivalent to D + D.
  """
  if not D.reduced :
    return double(reduce(D))
  
  if D.type == 31 :
    try :
      c = D.C.coefficients()
      if (c[5] == c[6] == c[8] == 0) :
        DD = fast_double_31_high_char(D)
      else :
        DD = fast_double_31(D)
    except ValueError:
      DD = double_31(D)
  elif D.type == 0 :
    DD = double_0(D)
  elif D.type == 11 :
    DD = double_11(D)
  elif D.type == 21 :
    DD = double_21(D)
  elif D.type == 22 :
    DD = double_22(D)
  else :
    raise ValueError("Divisor is of unexpected type. D = {}".format(D))
  if DD.reduced :
    return DD
  return reduce(DD)



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
  R0 = c3 + f1^2 - f0
  # Subtotal : 0I 2M 1SQ 8A
  # Compute
  #
  #   df = St - sT (mod f, g, h)
  #   dg = Tr - tR (mod f, g, h)
  #   dh = Rs - rS (mod f, g, h)
  e1 = -f1 - g2 # e1 + s0 = -2*g2, e1 - s0 = -2*f1
  e2 = R2 - f2
  df2 = s0 - g2
  df1 = T0 - g1
  df0 = T0*s0 + S0*t0 - g0
  dg2 = T0 - h2 + r0 - t0*e2
  dg1 = t0*(f1 + f1) - h1
  dg0 = T0*r0 + t0*(f0 - R0) - h0        # XXX : Can save 1A by saving f0 - R0 for later
  dh2 = f2*(e1 - g2) + R2*(g2 - s0) - S0 # XXX : Can save 1A by R2*(g2 - s0) --> R2*df2
  dh1 = f1*(e1 + s0) + g1*e2 - R0 + f0
  dh0 = f0*e1 + g0*e2 - S0*r0 - R0*s0
  # Subtotal :         14M     25A
  # Running total : 0I 16M 1SQ 33A

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
  # Subtotal :         12M      8A
  # Running total : 0I 28M 1SQ 41A
  
  if (print_matrices) :
    print("M = ")
    print(Matrix(C.K, [
      [a1,  a2,  a3,  a4,  a5 ],
      [a6,  a7,  a8,  a9,  a10],
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
  #        [ 0   0   b5  b6  b7 ]
  d1 = a1*a12 - a2*a11
  d2 = a6*a12 - a7*a11
  b1 = a1*a7  - a2*a6
  b2 = a1*a8  - a3*a6
  b3 = a1*a9  - a4*a6
  b4 = a1*a10 - a5*a6
  b5 = b1*a13 - d1*a8  + d2*a3
  b6 = b1*a14 - d1*a9  + d2*a4
  b7 = b1*a15 - d1*a10 + d2*a5
  # Subtotal :         21M     12A
  # Running total : 0I 49M 1SQ 53A
  
  if (b1 == 0) or (b5 == 0) :
    raise ValueError("Sum is not typical.")
  
  # Reduce M even more via back-substitution
  #
  #         [ Z  0  0  A1  A2 ]
  #   M'' = [ 0  Z  0  B1  B2 ]
  #         [ 0  0  Z  C1  C2 ]
  e1 = b3*b5 - b2*b6
  e2 = b4*b5 - b2*b7
  AB = a1*b1
  Z  = AB*b5
  C1 = AB*b6
  C2 = AB*b7
  B1 = a1*e1
  B2 = a1*e2
  A1 = b1*(a4*b5 - b6*a3) - a2*e1
  A2 = b1*(a5*b5 - b7*a3) - a2*e2
  # Subtotal :         18M      6A
  # Running total : 0I 67M 1SQ 59A
  
  # Compute
  #
  #   U = Z*x*f - C1*h - B1*g - A1*f
  #   V = Z*x*g - C2*h - B2*g - A2*f
  U1 = Z*f0 - C1*h1 - B1*g1 - A1*f1
  U2 =      - C1*h2 - B1*g2 - A1*f2
  U3 = Z*f1 - A1
  U4 = Z*f2 - B1
  U5 =      - C1
  V1 = Z*g0 - C2*h1 - B2*g1 - A2*f1
  V2 =      - C2*h2 - B2*g2 - A2*f2
  V3 = Z*g1 - A2
  V4 = Z*g2 - B2
  V5 =      - C2
  # Subtotal :         18M     14A
  # Running total : 0I 85M 1SQ 73A
  
  # Compute some inverses
  ZZt0      = U5^2 + Z*(U4 - V5)
  if (ZZt0 == 0) :
    raise ValueError("Sum of divisors is non-typical.")
  ZZZt0     = Z*ZZt0
  ZZZt0_inv = 1/ZZZt0
  ZZt0_inv  = Z*ZZZt0_inv
  zeta      = ZZt0*ZZZt0_inv # 1/Z
  tau       = (Z^2)*ZZt0_inv   # 1/t0
  # Subtotal :      1I  5M 2SQ  2A
  # Running total : 1I 90M 3SQ 75A
  
  # Rescale U and V polynomials by 1/Z
  u1 = zeta*U1
  u2 = zeta*U2
  u3 = zeta*U3
  u4 = zeta*U4
  u5 = zeta*U5
  v1 = zeta*V1
  v2 = zeta*V2
  v3 = zeta*V3
  v4 = zeta*V4
  v5 = zeta*V5
  # Subtotal :          10M
  # Running total : 1I 100M 3SQ 75A

  ff2 = u5^2 + u4 - v5
  r0  = u5*(ff2 + u4 - c7) + u3 - v4
  r1  = ff2*(ff2 - u4)
  gg1 = r1 - u5*(u3 + r0) + v3
  gg2 = -u4*u5 + v4 - r0 + (u4*r0 - u5*gg1 - u2)*tau
  ff1 = r0 + gg2
  ff0 = -c7*(r1 + gg2*u5) + u5*(ff2*u3 + ff1*u4 - c4 + u2) + gg2*u3 + gg1*u4 - ff2*v3 - ff1*v4 + u1 - v2
  gg0 = u5*(c3 - ff0 - u1 - ff1*u3) - gg1*u3 + ff1*v3 + v1
  hh0 = tau*(ff0*gg1 - gg0*r0)
  hh1 = tau*(gg1*gg2 - gg0)
  hh2 = gg1 + tau*(ff0 - gg2*r0)
  # Subtotal :          27M 1SQ  37A
  # Running total : 1I 127M 4SQ 112A
  
  ret = C34CurveDivisor(C, [[ff0, ff1, ff2, 1],
                       [gg0, gg1, gg2, 0, 1],
                       [hh0, hh1, hh2, 0, 0, 1]],
                       degree = 3, typ = 31, typical = True, reduced = True)
  ret.inv = tau
  return ret



def fast_double_31(D) :
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
  f2R2 = f2*R2
  f1R1 = f1*R1
  T0 = c5 - h2 - f2R2
  S0 = c4 - h2*T1 - h1 + f1R1 + f2R2 - (f1 + f2)*(R1 + R2)
  R0 = c3 - h1*T1 - f1R1 - f0
  # Subtotal : 0I 5M 15A

  # Compute
  #
  #   df = St - sT (mod f, g, h)
  #   dg = Tr - tR (mod f, g, h)
  #   dh = Rs - rS (mod f, g, h)
  R2mf2 = R2 - f2
  R1mf1 = R1 - f1
  R1mf1ps0 = R1mf1 + s0
  df2 = s0 - g2 - T1*f2
  df1 = T0 - g1 + T1*(s0 - f1)
  df0 = T0*s0 + S0*t0 - T1*f0 - g0
  dg2 = T0 - h2 + r0 - T1*g2 - t0*(R2mf2)
  dg1 = -t0*(R1mf1) - h1
  dg0 = T0*r0 - T1*g0 + t0*(f0 - R0) - h0
  dh2 = f2*(R1mf1ps0 - g2) + R2*(g2 - s0) - S0
  dh1 = f1*(R1mf1ps0) + g1*(R2mf2) - R1*s0 - R0 + f0
  dh0 = f0*(R1mf1ps0) + g0*(R2mf2) - S0*r0 - R0*s0
  # Subtotal :      0I 20M 30A
  # Running total : 0I 25M 45A

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
  # Subtotal :      0I 12M 8A
  # Running total : 0I 37M 53A
  
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
  #        [ 0   0   b5  b6  b7 ]
  d1 = a1*a12 - a2*a11
  d2 = a6*a12 - a7*a11
  b1 = a1*a7  - a2*a6
  b2 = a1*a8  - a3*a6
  b3 = a1*a9  - a4*a6
  b4 = a1*a10 - a5*a6
  b5 = b1*a13 - d1*a8  + d2*a3
  b6 = b1*a14 - d1*a9  + d2*a4
  b7 = b1*a15 - d1*a10 + d2*a5
  # Subtotal :      0I 21M 12A
  # Running total : 0I 58M 67A
  
  if (b1 == 0) or (b5 == 0) :
    raise ValueError("Sum is not typical.")
  
  # Reduce M even more via back-substitution
  #
  #         [ Z  0  0  A1  A2 ]
  #   M'' = [ 0  Z  0  B1  B2 ]
  #         [ 0  0  Z  C1  C2 ]
  e1 = b3*b5 - b2*b6
  e2 = b4*b5 - b2*b7
  AB = a1*b1
  Z  = AB*b5
  C1 = AB*b6
  C2 = AB*b7
  B1 = a1*e1
  B2 = a1*e2
  A1 = b1*(a4*b5 - b6*a3) - a2*e1
  A2 = b1*(a5*b5 - b7*a3) - a2*e2
  # Subtotal :      0I 18M  6A
  # Running total : 0I 74M 73A

  if (print_matrices) :
    print("Z*M_rref = ")
    print(Matrix(C.K, [
      [Z, 0, 0, A1, A2],
      [0, Z, 0, B1, B2],
      [0, 0, Z, C1, C2]]))
    print

  # Compute
  #
  #   U = Z*x*f - C1*h - B1*g - A1*f
  #   V = Z*x*g - C2*h - B2*g - A2*f
  U1 = Z*f0 - C1*h1 - B1*g1 - A1*f1
  U2 =      - C1*h2 - B1*g2 - A1*f2
  U3 = Z*f1 - A1
  U4 = Z*f2 - B1
  U5 =      - C1
  V1 = Z*g0 - C2*h1 - B2*g1 - A2*f1
  V2 =      - C2*h2 - B2*g2 - A2*f2
  V3 = Z*g1 - A2
  V4 = Z*g2 - B2
  V5 =      - C2
  # Subtotal :      0I 18M 14A
  # Running total : 0I 92M 87A

  # Compute some inverses
  ZZt0      = U5^2 - Z*(U5*c8 - U4 + V5)
  if (ZZt0 == 0) :
    raise ValueError("Sum of divisors is non-typical.")
  ZZZt0     = Z*ZZt0
  ZZZt0_inv = 1/ZZZt0
  ZZt0_inv  = Z*ZZZt0_inv
  zeta      = ZZt0*ZZZt0_inv # 1/Z
  tau       = (Z^2)*ZZt0_inv # 1/t0
  # Subtotal :      1I  6M 2S  3A
  # Running total : 1I 98M 2S 90A

  # Rescale U and V polynomials by 1/Z
  u1 = zeta*U1
  u2 = zeta*U2
  u3 = zeta*U3
  u4 = zeta*U4
  u5 = zeta*U5
  v1 = zeta*V1
  v2 = zeta*V2
  v3 = zeta*V3
  v4 = zeta*V4
  v5 = zeta*V5
  # Subtotal :      0I  10M 0S  0A
  # Running total : 1I 108M 2S 90A
  
  # Compute ff, gg such that gg*u = ff*v (mod C)
  gg3 = u5
  ff2 = u5*(u5 - c8) + u4 - v5
  r0 = u5*(u3 - c6) - v3
  gg2 = v4 + v5*(u5 - c8) + tau*(u5*(r0 + v5*(u4 - c7) + c5) + v5*(u3 - v4) - u2)
  r1 = ff2*v5 - gg2*u5
  ff1 = u5*(u4 - c7) + gg2 + u3 - v4
  gg1 = -r0 - r1
  ff0 = c7*r1 + u5*(u2 - c4) + gg2*u3 + gg1*u4 - ff2*v3 - ff1*v4 + u1 - v2
  gg0 = -c6*r1 + u5*(c3 - u1) - gg1*u3 + ff1*v3 + v1
  # Subtotal :          20M     32A
  # Running total : 1I 128M 2S 122A

  # Reduce gg modulo ff
  gg2 = gg2 - gg3*ff2
  gg1 = gg1 - gg3*ff1
  gg0 = gg0 - gg3*ff0
  # Subtotal :           3M      3A
  # Running total : 1I 131M 2S 125A

  # Compute third polynomial ...
  r2 = gg2 - ff1
  hh0 = tau*(ff0*gg1 + gg0*r2)
  hh1 = tau*(gg1*gg2 - gg0)
  hh2 = gg1 + tau*(gg2*r2 + ff0)
  # Subtotal :           7M      5A
  # Running total : 1I 138M 2S 130A

  ret = C34CurveDivisor(C, [[ff0, ff1, ff2, 1],
                        [gg0, gg1, gg2, 0, 1],
                        [hh0, hh1, hh2, 0, 0, 1]],
                        degree = 3, typ = 31, typical = True, reduced = True)
  ret.inv = tau
  return ret



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
  H = [ c[1] - f[0]*t,
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
  # Subtotal : 0I 8M 7A
  
  if alpha == 0 :
    new_f = [f[0], K.one()]
    new_g = [g[0]^2, K.zero(), g[0] + g[0], K.zero(), K.zero(), K.one()]
    # Total : 0I 8M 1S 8A
  
  else :
    # Compute beta = H(-f[0], -g[0])
    beta = ((H[3] - f[0])*f[0] + H[4]*g[0] - H[1])*f[0] + (H[5]*g[0] - H[2])*g[0] + H[0] # beta = h(-f0,-g0)
    a = beta*(1/alpha)
    new_f = [g[0] + f[0]*a, a, K.one()]
    new_g = [f[0]^2, f[0] + f[0], K.zero(), K.one()]
    # Total : 1I 15M 1S 15A
  
  return C34CurveDivisor(D.C, [new_f, new_g, new_h])



def double_21(D) :
  K = D.K
  f, g = D.f, D.g
  c = D.C.coefficients()
  new_f, new_g, new_h = [], [], []
  
  # Find polynomials G, H satisfying fH + gG = C
  A = flip(D) # Costs 0I 7M 12A
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
  # Subtotal : 0I 11M 18A
  
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
  # Subtotal : 0I 16M 17A
  # Running total : 0I 27M 35A

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
  # Subtotal :      0I  8M  4A
  # Running total : 0I 35M 39A
  
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
    # Subtotal :      1I 12M  3A
    # Running total : 1I 47M 42A

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
    # Subtotal :      0I 16M 15A
    # Running total : 1I 63M 57A
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
    # Subtotal :      1I 10M  2A
    # Running total : 1I 45M 41A

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
    # Subtotal :      1I  9M  8A
    # Running total : 1I 54M 49A
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
    # Subtotal :      1I  8M  1A
    # Running total : 1I 43M 40A

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
    # Subtotal :      0I  7M  6A
    # Running total : 1I 50M 46A
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
  A = flip(D) # Costs 0I 1M 2A
  G = A.g
  H = [0, 0, 0, c[6] - f[0], c[7], c[8], K.one()]
  H[2] = c[4] - f[0]*H[4]
  H[1] = c[3] - f[0]*H[3]
  H[0] = c[1] - f[0]*H[1]
  # Subtotal : 0I 4M 6A
  
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
  # Subtotal :      0I  8M  7A
  # Running total : 0I 12M 13A

  # Get matrix M into row echelon form.
  # M_ref = [ 1  a6  a7  a8  a9 ]
  #         [ 0   0  b1  b2  b3 ]
  b1 = a3 - a1*a7
  b2 = a4 - a1*a8
  b3 = a5 - a1*a9
  # Subtotal :      0I  3M  3A
  # Running total : 0I 15M 16A
 
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
    # Subtotal :      1I  4M  2A
    # Running total : 1I 19M 18A
    
    # Find polynomials forming an ideal generating set for 2D
    # f'' = x^2 + ...
    # g'' = y^2 + ...
    # h'' = x^3 + ...
    new_f = [ f[0]*u0, f[0] + u0, K.zero(), K.one() ]
    new_g = [ f[0]*v0 + g[0], v0, f[0]*v1 + g[2], K.zero(), v1, K.one() ]
    # Subtotal :      0I  3M  3A
    # Running total : 1I 22M 21A
    # 2D is type 42
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
    # Subtotal :      1I  2M  1A
    # Running total : 1I 17M 17A
    
    # Find polynomials forming an ideal generating set for 2D
    # f'' = x^2 + ...
    # g'' =  xy + ...
    # h'' = x^3 + ...
    new_f = [ f[0]*u0, f[0] + u0, K.zero(), K.one() ]
    new_g = [ f[0]*v0, v0, f[0], K.zero(), K.one() ]
    # new_h = [ f[0]*w0 + g[0]*w1, w0, g[2]*w1, f[0], K.zero(), w1, K.one() ]
    # Subtotal :      0I  2M  1A
    # Running total : 1I 19M 18A
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

  r0 = g1
  s0 = f1 - g2
  t0 = -f2

  T1 = c8
  R2 = c7 - f2
  R1 = c6 - f1
  f2R2 = f2*R2
  f1R1 = f1*R1
  T0 = c5 - h2 - f2R2
  S0 = c4 - h2*T1 - h1 + f1R1 + f2R2 - (f1 + f2)*(R1 + R2)
  R0 = c3 - h1*T1 - f1R1 - f0
  # Subtotal : 0I 5M 15A

  # Compute
  #
  #   df = St - sT (mod f, g, h)
  #   dg = Tr - tR (mod f, g, h)
  #   dh = Rs - rS (mod f, g, h)
  #
  # (The divisor A = div(df, dg, dh) is degree 4 and equivalent to D)
  R2mf2 = R2 - f2
  R1mf1 = R1 - f1
  R1mf1ps0 = R1mf1 + s0
  df2 = s0 - g2 - T1*f2
  df1 = T0 - g1 + T1*(s0 - f1)
  df0 = T0*s0 + S0*t0 - T1*f0 - g0
  dg2 = T0 - h2 + r0 - T1*g2 - t0*(R2mf2)
  dg1 = -t0*(R1mf1) - h1
  dg0 = T0*r0 - T1*g0 + t0*(f0 - R0) - h0
  dh2 = f2*(R1mf1ps0 - g2) + R2*(g2 - s0) - S0
  dh1 = f1*(R1mf1ps0) + g1*(R2mf2) - R1*s0 - R0 + f0
  dh0 = f0*(R1mf1ps0) + g0*(R2mf2) - S0*r0 - R0*s0
  # Subtotal :      0I 20M 30A
  # Running total : 0I 25M 45A

  # Compute M
  #
  #       [ a1   a2   a3   a4   a5   a6   a7  ]
  #   M = [ a8   a9   a10  a11  a12  a13  a14 ]
  #       [ a15  a16  a17  a18  a19  a20  a21 ]
  #
  # Where the columns represent the reductions of G, H, yG, xG, xH, modulo f, g.
  # a7, 14, a21 are computed later only when actually needed
  a1,  a2,  a3  = df0, dg0, dh0
  a8,  a9,  a10 = df1, dg1, dh1
  a15, a16, a17 = df2, dg2, dh2
  a4  =    - f0*a8 - g0*a15
  a11 = a1 - f1*a8 - g1*a15
  a18 =    - f2*a8 - g2*a15
  a5  =    - f0*a9 - g0*a16
  a12 = a2 - f1*a9 - g1*a16
  a19 =    - f2*a9 - g2*a16
  a6  =    - f0*a10 - g0*a17
  a13 = a3 - f1*a10 - g1*a17
  a20 =    - f2*a10 - g2*a17
  # Subtotal :      0I 18M 12A
  # Running total : 0I 43M 57A

  # If we swap the first row of the matrix with another,
  # aswap records with which row it was swapped.
  aswap = 0 

  if (a1 != 0 or a8 != 0 or a15 != 0) :
    # Swap rows, if necessary, to ensure a1 != 0
    if (a1 == 0) :
      if (a8 != 0) :
        a1, a2, a3,  a4,  a5,  a6,  a8, a9, a10, a11, a12, a13 = \
        a8, a9, a10, a11, a12, a13, a1, a2, a3,  a4,  a5,  a6
        aswap = 1
      else :
        a1,  a2,  a3,  a4,  a5,  a6,  a15, a16, a17, a18, a19, a20 = \
        a15, a16, a17, a18, a19, a20, a1,  a2,  a3,  a4,  a5,  a6
        aswap = 2

    # Partially reduce M.
    #
    #        [ a1  a2  a3  a4  a5   a6   a7  ]
    #   M' = [ 0   b1  b2  b3  b4   b5   b6  ]
    #        [ 0   b7  b8  b9  b10  b11  b12 ]
    #
    # Since we did not compute a7, a14, a21, we do not compute b6, b12 unless needed later
    
    b1  = a1*a9  - a2*a8
    b2  = a1*a10 - a3*a8
    b3  = a1*a11 - a4*a8
    b4  = a1*a12 - a5*a8
    b5  = a1*a13 - a6*a8
    b7  = a1*a16 - a2*a15
    b8  = a1*a17 - a3*a15
    b9  = a1*a18 - a4*a15
    b10 = a1*a19 - a5*a15
    b11 = a1*a20 - a6*a15
    # Subtotal :      0I 20M 10A
    # Running total : 0I 63M 67A
    
    if (b1 != 0 or b7 != 0) :
      # Before reducing any more, ensure b1 != 0 by swapping rows, if necessary
      if (b1 == 0) :
        b1, b2, b3, b4, b5, b7, b8, b9, b10, b11 = b7, b8, b9, b10, b11, b1, b2, b3, b4, b5
      # Continue reducing M.
      #
      #           [ a1  a2  a3  a4  a5  a6  a7 ]
      #   M_ref = [ 0   b1  b2  b3  b4  b5  b6 ]
      #           [ 0   0   c1  c2  c3  c4  c5 ]
      # 
      # Since we did not compute b6, b12, we do not cmopute c5.

      c1 = b1*b8  - b2*b7
      c2 = b1*b9  - b3*b7
      c3 = b1*b10 - b4*b7
      c4 = b1*b11 - b5*b7
      # Subtotal :      0I  8M  4A
      # Running total : 0I 71M 71A
      
      if (c1 != 0) :
        # Now compute reduced row echelon form of (the first six columns of) M
        #
        #            [ 1  0  0  -r0  -s0  -t0  * ]
        #   M_rref = [ 0  1  0  -r1  -s1  -t1  * ]
        #            [ 0  0  1  -r2  -s2  -t2  * ]
        
        ab = a1*b1
        abc = ab*c1
        delta = 1/abc
        alpha = delta*b1*c1 # = 1/a1
        beta  = delta*a1*c1 # = 1/b1
        gamma = delta*ab    # = 1/c1
        r2 = -gamma*c2
        s2 = -gamma*c3
        t2 = -gamma*c4
        r1 = -beta*(b2*r2 + b3)
        s1 = -beta*(b2*s2 + b4)
        t1 = -beta*(b2*t2 + b5)
        r0 = -alpha*(a2*r1 + a3*r2 + a4)
        s0 = -alpha*(a2*s1 + a3*s2 + a5)
        t0 = -alpha*(a2*t1 + a3*t2 + a6)
        # Subtotal :      1I 25M  9A
        # Running total : 1I 96M 80A

        # Compute D1 + D2 = <u, v, w>, where
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
        is_typical = (u5*(u5 - c8) + u4 - v5 != 0)

        # Subtotal :          28M  30A
        # Running total : 1I 124M 110A
        # 2D is of type 61
        return C34CurveDivisor(C, [[u0, u1, u2, u3, u4, u5, 1],
                                   [v0, v1, v2, v3, v4, v5, 0, 1],
                                   [w0, w1, w2, w3, w4, w5, 0, 0, 1]],
                                   degree = 6, typ = 61, typical = is_typical, reduced = False)
      elif (c2 != 0) :
        # M_ref is the matrix
        #
        #       [ a1  a2  a3  a4  a5  a6  a7 ]
        #   M = [ 0   b1  b2  b3  b4  b5  b6 ]
        #       [ 0   0   0   c2  c3  c4  c5 ]
        #
        # Reduce it to
        #
        #            [ 1  0  -r0  0  -s0  *  * ]
        #   M_rref = [ 0  1  -r1  0  -s1  *  * ]
        #            [ 0  0   0   1  -s2  *  * ]
        ab = a1*b1
        abc = ab*c2
        delta = 1/abc
        alpha = delta*b1*c2
        beta  = delta*a1*c2
        gamma = delta*ab
        s2 = -gamma*c3
        r1 = -beta*b2
        s1 = -beta*(b3*s2 + b4)
        r0 = -alpha*(a2*r1 + a3)
        s0 = -alpha*(a2*s1 + a4*s2 + a5)
        # Subtotal :      1I 16M  4A
        # Running total : 1I 87M 75A
        
        # Compute D1 + D2 = <u, v>, where
        #
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
        v5 = s2
        # Subtotal :          15M 15A
        # Running total : 1I 102M 90A
        # 2D is of type 63.
        return C34CurveDivisor(C, [[u0, u1, u2, u3, u4, 1],
                             [v0, v1, v2, v3, v4, 0, v5, 1], []],
                             degree = 6, typ = 63, reduced = False, typical = False)

      elif (c3 != 0) :
        # M_ref is the matrix
        #
        #       [ a1  a2  a3  a4  a5  a6  a7 ]
        #   M = [ 0   b1  b2  b3  b4  b5  b6 ]
        #       [ 0   0   0   0   c3  c4  c5 ]
        #
        # Reduce it to
        #
        #            [ 1  0  -r0  -s0  0  *  * ]
        #   M_rref = [ 0  1  -r1  -s1  0  *  * ]
        #            [ 0  0   0    0   1  *  * ]
        ab = a1*b1
        gamma = 1/ab
        alpha = gamma*b1
        beta  = gamma*a1
        r1 = -beta*b2
        s1 = -beta*b3
        r0 = -alpha*(a2*r1 + a3)
        s0 = -alpha*(a2*s1 + a4)
        # Subtotal :      1I  9M  2A
        # Running total : 1I 80M 73A
        # Compute D1 + D2 = <u, v>, where
        #
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
        # Subtotal :      0I 12M 12A
        # Running total : 1I 92M 85A
        # 2D is type 62
        return C34CurveDivisor(C, [[u0, u1, u2, u3, u4, 1], [v0, v1, v2, v3, v4, 0, 1], []],
                               degree = 6, typ = 62, reduced = False, typical = False)

      else :
        assert c4 == 0
        # D1 and D2 are non-disjoint.
        # We compute D1 + D2 = lcm(D1, D2) + gcd(D1, D2)
      	# lcm(D1, D2) will be of type 51
        #
        # Reduce
        #
        # [ a1  a2  a3  a4  a5 ]
        # [ 0   b1  b2  b3  b4 ]
        # [ 0   0   0   0   0  ]
        # 
        # to its reduced row echelon form.
        #
        # [ 1  0  -r0  -s0  -t0 ]
        # [ 0  1  -r1  -s1  -t1 ]
        # [ 0  0   0    0    0  ]
        gamma = 1/(a1*b1)
        alpha = gamma*b1
        beta  = gamma*a1
        r1 = -beta*b2
        s1 = -beta*b3
        t1 = -beta*b4
        r0 = -alpha*(a3 + a2*r1)
        s0 = -alpha*(a4 + a2*s1)
        t0 = -alpha*(a5 + a2*t1)
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
        L = C34CurveDivisor(C, [[u0, u1, u2, u3, u4, 1], [v0, v1, v2, v3, v4, 0, 1], [w0, w1, w2, w3, w4, 0, 0, 1]])

        # GCD(A, D) is type 11.
        # A basis for the GCD is given by the 1st and 2nd columns of M
        #
        # [ m6  m3 ]     [ n2 m3 ]     [ p0  q0 ]
        # [ m5  m2 ] ==> [ n1 m2 ] ==> [ 1   0  ]
        # [ m4  m1 ]     [ 0  m1 ]     [ 0   1  ]
        if (aswap == 0) :
          m1, m2, m3, m4, m5, m6 = a16, a9, a2, a15, a8, a1
        elif (aswap == 1) :
          m1, m2, m3, m4, m5, m6 = a16, a2, a9, a15, a1, a8
        else :
          m1, m2, m3, m4, m5, m6 = a2, a9, a16, a1, a8, a15
        if (m1 == 0) :
          m1, m2, m3, m4, m5, m6 = m4, m5, m6, m1, m2, m3
        n1 = m5*m1 - m2*m4
        n2 = m6*m1 - m3*m4
        om = 1/(m1*n1)
        mu = om*n1
        nu = om*m1
        p0 = nu*n2
        q0 = mu*(m3 - m2*p0)
        G = C34CurveDivisor(C, [[p0, 1], [q0, 0, 1], []])
        return reduce(L) + G

    elif (b2 != 0 or b8 != 0) :
      if (aswap == 0) :
        a7  =    - f0*a11 - g0*a18
        a14 = a4 - f1*a11 - g1*a18
        a21 =    - f2*a11 - g2*a18
      elif (aswap == 1) :
        a7  = a11 - f1*a4 - g1*a18
        a14 =     - f0*a4 - g0*a18
        a21 =     - f2*a4 - g2*a18
      else :
        a7  =     - f2*a11 - g2*a4
        a14 = a18 - f1*a11 - g1*a4
        a21 =     - f0*a11 - g0*a4
      b6  = a1*a14 - a7*a8
      b12 = a1*a21 - a7*a15
      # Subtotal :         10M  6A
      # Running total : 0I 73M 73A

      # Swap rows if necessary so that b2 != 0
      if (b2 == 0) :
        b2, b3, b4, b5, b6, b8, b9, b10, b11, b12 = b8, b9, b10, b11, b12, b2, b3, b4, b5, b6
      
      # We have the matrix
      #
      #        [ a1  a2  a3  a4  a5   a6   a7  ]
      #   M' = [ 0   0   b2  b3  b4   b5   b6  ]
      #        [ 0   0   b8  b9  b10  b11  b12 ]
      #
      # We wish to reduce it to
      #
      #           [ a1  a2  a3  a4  a5  a6  a7 ]
      #   M_ref = [ 0   0   b2  b3  b4  b5  b6 ]
      #           [ 0   0   0   c1  c2  c3  c4 ]
      c1 = b2*b9  - b3*b8
      c2 = b2*b10 - b4*b8
      c3 = b2*b11 - b5*b8
      c4 = b2*b12 - b6*b8
      # Subtotal :      0I  8M  4A
      # Running total : 0I 81M 77A
      
      if (c1 != 0) :
        # Compute M_rref
        #
        #          [ 1  -r0  0  0  *  *  -s0 ]
        # M_rref = [ 0   0   1  0  *  *  -s1 ]
        #          [ 0   0   0  1  *  *  -s2 ]
        ab = a1*b2
        abc = ab*c1
        delta = 1/abc
        alpha = delta*b2*c1
        beta  = delta*a1*c1
        gamma = delta*ab
        r0 = -alpha*a2
        s2 = -gamma*c4
        s1 = -beta*(b3*s2 + b6)
        s0 = -alpha*(a3*s1 + a4*s2 + a7)
        # Subtotal :      1I 14M  3A
        # Running total : 1I 95M 80A
        
        # Compute D1 + D2 = <u, v>, where
        #
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
        # Subtotal :      0I  12M 10A
        # Running total : 1I 107M 90A
        # 2D is of type 64.
        return C34CurveDivisor(C, [[u0, u1, u2, u3, 1], [v0, v1, v2, v3, 0, v5, v6, 0, 0, 1], []],
                               degree = 6, typ = 64, typical = False, reduced = False)
        
      else :
        assert c2 == c3 == c4 == 0
        # D1 and D2 are non-disjoint.
        # We compute D1 + D2 = lcm(D1, D2) + gcd(D1, D2)
        #
        # M_ref is the matrix
        #
        #           [ a1  a2  a3  a4  a5  a6  a7 ]
        #   M_ref = [ 0   0   b2  b3  b4  b5  b6 ]
        #           [ 0   0   0   0   0   0   0  ]
        # 
        # Compute its reduced row echelon form
        #
        #            [ 1  -r0  0  -s0  *  *  * ]
        #   M_rref = [ 0   0   1  -s1  *  *  * ]
        #            [ 0   0   0   0   0  0  0 ]
        gamma = 1/(a1*b2)
        alpha = gamma*b2
        beta = gamma*a1
        s1 = -beta*b3
        r0 = -alpha*a2
        s0 = -alpha*(a4 + a3*s1)
        # Subtotal :      1I  7M  1A
        # Running total : 0I 88M 78A
        
        # LCM(D1, D2) is of type 53, generated by <u, v>, where
        #
        #   u = r0*f + g
        #   v = s0*f + s1*h + xf
        u0 = f0*r0 + g0
        u1 = f1*r0 + g1
        u2 = f2*r0 + g2
        u3 = r0
        v0 = f0*s0 + h0*s1 - f2*u0
        v1 = f1*s0 + h1*s1 + f0 - f2*u1
        v2 = f2*(s0 - u2) + h2*s1
        v3 = s0 + f1 - f2*u3
        v5 = s1
        L = C34CurveDivisor(C, [[u0, u1, u2, u3, 1], [v0, v1, v2, v3, 0, v5, 1], []])
        # Subtotal : 0I 12M 12A

        # GCD(D1, D2) is type 11, generated by polynomials <p, q> where
        #
        #   p = x + p0
        #   q = y + q0
        #
        # p0 and q0 are determined by columns 1 and 3 of M, and may be
        # computed by performing column operations on M :
        # 
        #   [ a1   a3  ]    [ m6  m3 ]     [ n2 m3 ]     [ p0  q0 ]
        #   [ a8   a10 ] =: [ m5  m2 ] ==> [ n1 m2 ] ==> [ 1   0  ]
        #   [ a15  a17 ]    [ m4  m1 ]     [ 0  m1 ]     [ 0   1  ]
        #
        # We must account for whether the rows of M were swapped.
        if (aswap == 0) :
          m1, m2, m3, m4, m5, m6 = a17, a10, a3, a15, a8, a1
        elif (aswap == 1) :
          m1, m2, m3, m4, m5, m6 = a17, a3, a10, a15, a1, a8
        else :
          m1, m2, m3, m4, m5, m6 = a3, a10, a17, a1, a8, a15
        if (m1 == 0) :
          m1, m2, m3, m4, m5, m6 = m4, m5, m6, m1, m2, m3
        n1 = m5*m1 - m2*m4
        n2 = m6*m1 - m3*m4
        om = 1/(m1*n1)
        mu = om*n1
        nu = om*m1
        p0 = nu*n2
        q0 = mu*(m3 - m2*p0)
        G = C34CurveDivisor(C, [[p0, 1], [q0, 0, 1], []])
        # Subtotal : 1I 10M 3A
        
        ret = reduce(L) + G # add_21_11
        return ret

    elif (b3 != 0) or (b9 != 0) :
      # M_ref is the matrix
      #
      #           [ a1  a2  a3  a4  a5  a6  a7 ]
      #   M_ref = [ 0   0   0   b3  b4  b5  b6 ]
      #           [ 0   0   0   0   0   0   0  ]
      #
      # Compute M_rref
      #
      #            [ 1  -r0  -s0  0  *  *  * ]
      #   M_rref = [ 0   0    0   1  *  *  * ]
      #            [ 0   0    0   0  0  0  0 ]
      alpha = 1/a1
      r0 = -alpha*a2
      s0 = -alpha*a3
      # Subtotal : 1I 2M 0A

      # LCM(D1, D2) is of type 52, generated by <u, v>, where
      #
      #   u = r0*f + g
      #   v = s0*f + h
      u0 = f0*r0 + g0
      u1 = f1*r0 + g1
      u2 = f2*r0 + g2
      u3 = r0
      v0 = f0*s0 + h0
      v1 = f1*s0 + h1
      v2 = f2*s0 + h2
      v3 = s0
      L = C34CurveDivisor(C, [[u0, u1, u2, u3, 1], [v0, v1, v2, v3, 0, 1], []])
      # Subtotal : 0I 6M 6A
      
      # GCD(D1, D2) is type 11, generated by polynomials <p, q> where
      #
      #   p = x + p0
      #   q = y + q0
      #
      # p0 and q0 are determined by columns 1 and 4 of M, and may be
      # computed by performing column operations on M :
      # 
      #   [ a1   a4  ]    [ m6  m3 ]     [ n2 m3 ]     [ p0  q0 ]
      #   [ a8   a11 ] =: [ m5  m2 ] ==> [ n1 m2 ] ==> [ 1   0  ]
      #   [ a15  a18 ]    [ m4  m1 ]     [ 0  m1 ]     [ 0   1  ]
      #
      # We must account for whether the rows of M were swapped.
      if (aswap == 0) :
        m1, m2, m3, m4, m5, m6 = a18, a11, a4, a15, a8, a1
      elif (aswap == 1) :
        m1, m2, m3, m4, m5, m6 = a18, a4, a11, a15, a1, a8
      else :
        m1, m2, m3, m4, m5, m6 = a4, a11, a18, a1, a8, a15
      if (m1 == 0) :
        m1, m2, m3, m4, m5, m6 = m4, m5, m6, m1, m2, m3
      n1 = m5*m1 - m2*m4
      n2 = m6*m1 - m3*m4
      om = 1/(m1*n1)
      mu = om*n1
      nu = om*m1
      p0 = nu*n2
      q0 = mu*(m3 - m2*p0)
      G = C34CurveDivisor(C, [[p0, 1], [q0, 0, 1], []])
      # Subtotal : 1I 10M 3A

      ret = reduce(L) + G # add_11_11
      return ret

    else :
      # M_ref is the matrix
      #
      #           [ a1  a2  a3  a4  a5  a6  a7 ]
      #   M_ref = [ 0   0   0   0   0   0   0  ]
      #           [ 0   0   0   0   0   0   0  ]
      assert False, "This case (deg(gcd(A, D)) = 2) should be impossible!"

  elif (a2 != 0 or a9 != 0 or a16 != 0) :
    if (a2 == 0) :
      if (a9 != 0) :
        a2, a3, a5, a6, a9, a10, a12, a13 = a9, a10, a12, a13, a2, a3, a5, a6
        aswap = 1
      else :
        a2, a3, a5, a6, a16, a17, a19, a20 = a16, a17, a19, a20, a2, a3, a5, a6
        aswap = 2

    # M is the matrix
    #
    #       [ 0  a2   a3   0  a5   a6   0 ]
    #   M = [ 0  a9   a10  0  a12  a13  0 ]
    #       [ 0  a16  a17  0  a19  a20  0 ]
    #
    # Partially reduce M to
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
    # Subtotal :      0I 12M  6A
    # Running total : 0I 55M 63A

    if (b1 != 0 or b4 != 0) :
      if (b1 == 0) :
        b1, b2, b3, b4, b5, b6 = b4, b5, b6, b1, b2, b3
      # Compute the row echelon form of M
      #
      #           [ 0  a2  a3  0  a5  a6  0 ]
      #   M_ref = [ 0  0   b1  0  b2  b3  0 ]
      #           [ 0  0   0   0  c1  c2  0 ]
      #
      # Although c1 must be zero, otherwise LCM(D1, D2) would be degree 5
      # despite M having full rank.
      #c1 = b1*b5 - b2*b4
      c2 = b1*b6 - b3*b4
      # Subtotal :      0I  2M  1A
      # Running total : 0I 57M 64A
      
      if (c2 != 0) :
        # D1 + D2 is type 65, principal, generated by <f>
        new_f = [f0, f1, f2, 1]
        return C34CurveDivisor(C, [[f0, f1, f2, 1], [], []],
                               degree = 6, typ = 65, typical = false, reduced = false)
      else :
        # Compute the reduced row echelon form of M
        #
        #            [ 0  1  0  0  *  -r0  0 ]
        #   M_rref = [ 0  0  1  0  *  -r1  0 ]
        #            [ 0  0  0  0  0   0   0 ]
        # Then LCM(D1, D2) is type 54
        gamma = 1/(a2*b1)
        alpha = gamma*b1
        beta = gamma*a2
        r1 = -beta*b3
        r0 = -alpha*(a6 + a3*r1)
        # Subtotal : 1I 6M 1A

        # D1 + D2 is generated by <u, v>, where
        #
        #   u = f
        #   v = r0*g + r1*h + xh - h1*u
        v0 = g0*r0 + h0*r1 - h1*f0
        v1 = g1*r0 + h1*(r1 - f1) + h0
        v2 = g2*r0 + h2*r1 - h1*f2
        v4 = r0 + h2
        v5 = r1
        L = C34CurveDivisor(C, [copy(D1.f), [v0, v1, v2, 0, v4, v5, 0, 0, 1], []])
        # Subtotal : 0I 8M 8A

        # GCD(D1, D2) is type 11, generated by polynomials <p, q> where
        #
        #   p = x + p0
        #   q = y + q0
        #
        # p0 and q0 are determined by columns 2 and 3 of M, and may be
        # computed by performing column operations on M :
        # 
        #   [ a2   a3  ]    [ m6  m3 ]     [ n2 m3 ]     [ p0  q0 ]
        #   [ a9   a10 ] =: [ m5  m2 ] ==> [ n1 m2 ] ==> [ 1   0  ]
        #   [ a16  a17 ]    [ m4  m1 ]     [ 0  m1 ]     [ 0   1  ]
        #
        # We must account for whether the rows of M were swapped.
        if (aswap == 0) :
          m1, m2, m3, m4, m5, m6 = a17, a10, a3, a16, a9, a2
        elif (aswap == 1) :
          m1, m2, m3, m4, m5, m6 = a17, a3, a10, a16, a2, a9
        else :
          m1, m2, m3, m4, m5, m6 = a3, a10, a17, a2, a9, a16
        if (m1 == 0) :
          m1, m2, m3, m4, m5, m6 = m4, m5, m6, m1, m2, m3
        n1 = m5*m1 - m2*m4
        n2 = m6*m1 - m3*m4
        om = 1/(m1*n1)
        mu = om*n1
        nu = om*m1
        p0 = nu*n2
        q0 = mu*(m3 - m2*p0)
        G = C34CurveDivisor(C, [[p0, 1], [q0, 0, 1], []])
        # Subtotal : 1I 8M 3A

        ret = reduce(L) + G # add_22_11
        return ret

    else :
      # We have the matrix
      #
      #     [ 0  a2  a3  0  a5  a6  0 ]
      # M = [ 0  0   0   0  b2  b3  0 ]
      #     [ 0  0   0   0  b5  b6  0 ]
      #
      # LCM(D1, D2) will be type 43, so the b-values are all zero.
      assert b2 == b3 == b5 == b6 == 0, "{}".format(Matrix(K, 3, 7, [0, a2, a3, 0, a5, a6, 0, 0, 0, 0, 0, b2, b3, 0, 0, 0, 0, 0, b5, b6, 0]))
      alpha = 1/a2
      r0 = -alpha*a3
      # u = r0*g + h
      u0 = g0*r0 + h0
      u1 = g1*r0 + h1
      u2 = g2*r0 + h2
      u4 = r0
      L = C34CurveDivisor(C, [[f0, f1, f2, 1], [u0, u1, u2, 0, u4, 1], []])

      # GCD(D1, D2) is degree 2 (type 21 or 22).
      # One generator of the GCD is given by the 2nd column of M
      #
      # [ m3 ]     [ u0 ]    [ u0 ]
      # [ m2 ] ==> [ u1 ] or [ 1  ]
      # [ m1 ]     [ 1  ]    [ 0  ]
      if (aswap == 0) :
        m1, m2, m3 = a16, a9, a2
      elif (aswap == 1) :
        m1, m2, m3 = a16, a2, a9
      else :
        m1, m2, m3 = a2, a9, a16
      if (m1 != 0) :
        mu = 1/m1
        u1 = mu*m2
        u0 = mu*m3
        # Compute v = x^2 + v1*x + v0 by taking f modulo u = y + u1*x + u0
        v0 = f0 - f2*u0
        v1 = f1 - f2*u1
        G = C34CurveDivisor(C, [[u0, u1, 1], [v0, v1, 0, 1], []])
      else :
        assert m2 != 0
        mu = 1/m2
        u0 = mu*m3
        # Compute v = y^2 + v2*y + v0 by taking h modulo u = x + u0
        v0 = h0 - h1*u0
        v2 = h2
        G = C34CurveDivisor(C, [[u0, 1], [v0, 0, v2, 0, 0, 1], []])
      return reduce(L) + G

  else :
    # We have the matrix
    #
    #     [ 0  0  a3   0  0  a6   0 ]
    # M = [ 0  0  a10  0  0  a13  0 ]
    #     [ 0  0  a17  0  0  a20  0 ]
    #
    assert False, "This case (deg(gcd(A, D)) = 2) should be impossible!"

