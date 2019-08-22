"""
  Summary of costs for adding divisors.
  Only applies to when D1 and D2 are disjoint.
  
    D1     D2     D3
  type | type | type |   I |   M |   A |
  -----+------+------+-----+-----+-----+
     * |    0 |    * |   0 |   0 |   0 |
  -----+------+------+-----+-----+-----+
    11 |   11 |   21 |   1 |   3 |   4 |
       |      |   22 |   0 |   1 |   3 |
  -----+------+------+-----+-----+-----+
    21 |   11 |   31 |   1 |  13 |  14 |
       |      |   32 |   0 |   4 |   6 |
  -----+------+------+-----+-----+-----+
    21 |   21 |   41 |   1 |  45 |  30 |
       |      |   42 |   1 |  27 |  18 |
       |      |   43 |   1 |  33 |  21 |
       |      |   44 |   0 |  12 |   9 |
  -----+------+------+-----+-----+-----+
    21 |   22 |   41 |   1 |  17 |  13 |
       |      |   42 |   0 |   2 |   1 |
  -----+------+------+-----+-----+-----+
    22 |   11 |   31 |   1 |   5 |   5 |
       |      |   33 |   0 |   1 |   3 |
  -----+------+------+-----+-----+-----+
    22 |   22 |   43 |   1 |   5 |   6 |
  -----+------+------+-----+-----+-----+

    D1     D2     D3
  type | type | type |   I |   M |   A |
  -----+------+------+-----+-----+-----+
    31 |   11 |   41 |   1 |  20 |  21 |
       |      |   42 |   0 |   6 |   9 |
       |      |   43 |   1 |  10 |  21 |
  -----+------+------+-----+-----+-----+
    31 |   21 |   51 |   1 |  57 |  44 |
       |      |   52 |   1 |  35 |  29 |
       |      |   53 |   1 |  46 |  36 |
       |      |   54 |   1 |  37 |  30 |
  -----+------+------+-----+-----+-----+
    31 |   22 |   51 |   1 |  45 |  32 |
       |      |   52 |   1 |  23 |  17 |
       |      |   53 |   1 |  34 |  24 |
       |      |   54 |   1 |  29 |  19 |
  -----+------+------+-----+-----+-----+
    31 |   31 |   61 |   1 |  98 |  71 |
       |      |   62 |   1 |  67 |  49 |
       |      |   63 |   1 |  77 |  54 |
       |      |   64 |   1 |  80 |  58 |
       |      |   65 |   0 |  32 |  28 |
  -----+------+------+-----+-----+-----+
"""



def add(D1, D2) :
  """
    Add two divisors, D1 and D2.
    
    Will throw an exception if the support of D1 and D2 share a point in common.
    
    Input : Two distinct typical C34CurveDivisors, D1 and D2.
    Output : The reduced C34CurveDivisor D3 equivalent to D1 + D2. May be typical or semi-typical (or
             neither?)
  """
  if (D2.C != D1.C) :
    raise ValueError("Divisors are of different curves.")
  if (D1 == D2) :
    raise ValueError("Divisors are non-distinct.\nD1 = {}\nD2 = {}".format(D1, D2))
  if (not D1.reduced) :
    raise ValueError("Divisor is not reduced.\nD1 = {}".format(D1))
  if (not D2.reduced) :
    raise ValueError("Divisor is not reduced.\nD2 = {}".format(D2))

  # Arrange divisors so that deg(D1) >= deg(D2)
  if (D2.degree > D1.degree) :
    return add(D2, D1)
  
  if D2.degree == 0 :
    return D1

  # If D1 = <f, g, h> is a reduced degree 3 divisor but missing its h polynomial, compute h.
  if (D1.degree == 3) and (len(D1.h) == 0) :
    print("Computing h polynomial for D1 = {}".format(D1))
    if D1.f[2] == 0 :
      raise ValueError("D1 is a degree 4 divisor.\nD1 = {}".format(D1))
    a = 1/D1.f[2]
    K = D1.K

    # This gives h of the form y^2 + ay + bx + c in 1I 7M
    D1.h = [a*(D1.g[1]*D1.f[0] - (D1.f[1] - D1.g[2])*D1.g[0]),
            a*(- D1.g[0] + D1.g[1]*D1.g[2]),
            a*(D1.f[0] - (D1.f[1] - D1.g[2])*D1.g[2]) + D1.g[1],
            K.zero(), K.zero(), K.one()]

  # If D2 = <f, g, h> is a reduced degree 3 divisor but missing its h polynomial, compute h.
  if (D2.degree == 3) and (len(D2.h) == 0) :
    print("Computing h polynomial for D2 = {}".format(D2))
    if D2.f[2] == 0 :
      raise ValueError("D2 is a degree 4 divisor.\nD2 = {}".format(D2))
    a = 1/D2.f[2]
    K = D1.K

    # This gives h of the form y^2 + ay + bx + c in 1I 7M
    D2.h = [a*(D2.g[1]*D2.f[0] - (D2.f[1] - D2.g[2])*D2.g[0]),
            D2.a*(- D2.g[0] + D2.g[1]*D2.g[2]),
            a*(D2.f[0] - (D2.f[1] - D2.g[2])*D2.g[2]) + D2.g[1],
            K.zero(), K.zero(), K.one()]

  # Examine the types of D1 and D2 and call the appropriate function
  T = (D1.type, D2.type)
  if (T == (11, 11)) :
    D3 = add_11_11(D1, D2)
  elif (T == (21, 11)) :
    D3 = add_21_11(D1, D2)
  elif (T == (21, 21)) :
    D3 = add_21_21(D1, D2)
  elif (T == (21, 22)) :
    D3 = add_21_22(D1, D2)
  elif (T == (22, 11)) :
    D3 = add_22_11(D1, D2)
  elif (T == (22, 21)) :
    D3 = add_21_22(D2, D1)
  elif (T == (22, 22)) :
    D3 = add_22_22(D1, D2)
  elif (T == (31, 11)) :
    D3 = add_31_11(D1, D2)
  elif (T == (31, 21)) :
    D3 = add_31_21(D1, D2)
  elif (T == (31, 22)) :
    D3 = add_31_22(D1, D2)
  elif (T == (31, 31)) :
    D3 = fast_add_31_31(D1, D2)
  else :
    raise ValueError("Divisors are of unhandled types.\nD1 = {}\nD2 = {}".format(D1, D2))
  
  if D3.reduced :
    return D3
  return reduce(D3)



def km_add_31_31(D1, D2) :
  """
    An implementation of Kamal Khuri-Makdisi's typical divisor class arithmetic, based on his
    2007/18 papers. Variable names have changed, and some changes are made to account for a
    difference in curve model -- KM uses a curve equation
      
      y^3 - x^4 + ...,
    
    but I use
      
      y^3 + x^4 + ...
    
    The method and inversion/multiplication counts remain the same.
  """
  C = D1.C
  f0, f1, f2 = D1.f[0:3]
  g0, g1, g2 = D1.g[0:3]
  F0, F1, F2 = D2.f[0:3]
  G0, G1, G2 = D2.g[0:3]

  print_matrices = False
  strassen = True
  toom_cook = True
  karatsuba = True

  if (f2 == 0) :
    raise ValueError("Divisor is not typical.\nD1 = {}".format(D1))
  if (F2 == 0) :
    raise ValueError("Divisor is not typical.\nD2 = {}".format(D2))

  c0, c1, c2, c3, c4, c5, c6, c7, c8 = C.coefficients()
  if (C.K.characteristic() <= 3) :
    raise ValueError("Curve's base field is of characteristic 3 or less.")
  if (c5 != 0) or (c6 != 0) or (c8 != 0) :
    raise ValueError("Curve equation is not in short form.")
  half = C.K(1/2) # Assumed to be a "free" computation in Kamal's model
  if (D2.inv == 0) :
    D2.inv = 1/F2
    # print("Computed 1/F2")
  F2_inv = D2.inv

  # Compute M
  #
  #       [ a1   a2   a3   a4   a5  ]
  #   M = [ a6   a7   a8   a9   a10 ]
  #       [ a11  a12  a13  a14  a15 ]

  # Columns 1 and 4 are reductions of f and g, modulo f' and g', respectively.
  a1  = f0 - F0
  a6  = f1 - F1
  a11 = f2 - F2
  a4  = g0 - G0
  a9  = g1 - G1
  a14 = g2 - G2

  # Column 2 and 5 are derived from 1 and 4 via a linear transformation.
  #
  #   [ a2   a5  ]   [ 0  -F0  -G0 ] [ a1   a4  ]
  #   [ a7   a10 ] = [ 1  -F1  -G1 ]*[ a6   a9  ]
  #   [ a12  a15 ]   [ 0  -F2  -G2 ] [ a11  a14 ]
  #
  # Strassen's technique is used to save 1M at the cost of 14A
  #
  #   [ a2   a5  ]   [ 0  -F0  -G0 ] [ a1   a4  ]   [ 0   0    0  ] [ 0    0   ]
  #   [ a7   a10 ] = [ 1   0    0  ]*[ a6   a9  ] + [ 0  -F1  -G1 ]*[ a6   a9  ]
  #   [ a12  a15 ]   [ 0   0    0  ] [ a11  a14 ]   [ 0  -F2  -G2 ] [ a11  a14 ]
  if (strassen) :
    m1 = (-F1 - G2)*(a6 + a14)
    m2 = (-F2 - G2)*a6
    m3 = -F1*(a9 - a14)
    m4 = -G2*(a11 - a6)
    m5 = (-F1 - G1)*a14
    m6 = (-F2 + F1)*(a6 + a9)
    m7 = (-G1 + G2)*(a11 + a14)
    a2  = -F0*a6 - G0*a11
    a7  = a1 + m1 + m4 - m5 + m7
    a12 = m2 + m4
    a5  = -F0*a9 - G0*a14
    a10 = a4 + m3 + m5
    a15 = m1 - m2 + m3 + m6
  else :
    a2  =    - F0*a6 - G0*a11
    a7  = a1 - F1*a6 - G1*a11
    a12 =    - F2*a6 - G2*a11
    a5  =    - F0*a9 - G0*a14
    a10 = a4 - F1*a9 - G1*a14
    a15 =    - F2*a9 - G2*a14

  # Column 3 represents the coefficients of ((y + g1)F - (x + f1 - g2)G)/F2
  # These are computed via Prop. 3.2. in [KM2007]
  n0 = a11*F2_inv
  n1 = n0*G0
  n2 = a1 - n0*F0
  n3 = n0*G1
  n4 = a6 - n0*(F1 - G2)
  a3  =    - n3*F0 - n4*G0
  a8  = n1 - n3*F1 - n4*G1
  a13 = n2 - n3*F2 - n4*G2
  # Subtotal : 0I 22M 36A (Assuming Strassen used)
  # Running total : 0I 22M 36A

  if (print_matrices) :
    print "M = "
    print Matrix(C.K, [
      [a1, a2, a3, a4, a5],
      [a6, a7, a8, a9, a10],
      [a11, a12, a13, a14, a15]])
    print

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
  # Running total : 0I 22M 41A

  if (print_matrices) :
    print "M' = "
    print Matrix(C.K, [
      [A1, A2, A3, A4, A5],
      [A6, A7, A8, A9, A10],
      [A11, A12, A13, A14, A15]])
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
  # Running total : 0I 43M 53A

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
    print "M'_ref = "
    print Matrix(C.K, [
      [ A1, A2, A3, A4, A5 ],
      [ 0, B1, B2, B3, B4 ],
      [ 0, 0, C1, C2, C3 ]])
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
  # Running total : 1I 49M 53A
  
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
  # Subtotal : 1I 12M 6A
  # Running total : 1I 61M 59A

  if (print_matrices) :
    print("M'_rref = ")
    print Matrix(C.K, [
      [1, 0, 0, -r0, -s0],
      [0, 1, 0, -r1, -s1],
      [0, 0, 1, -r2, -s2]])
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
  # Running total : 1I 79M 93A

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
  # Running total : 1I 90M 1SQ 1CC 2CM 125A

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
  # Running total : 1I 96M 1SQ 1CC 2CM 132A

  ret = C34CurveDivisor(C, [[new_f0, new_f1, new_f2, 1],
                      [new_g0, new_g1, new_g2, 0, 1], []],
                      degree = 3, typ = 31, reduced = True, typical = True, inv = l1_inv)
  return ret



def fast_add_31_31(D1, D2) :
  C = D1.C
  f0, f1, f2 = D1.f[0:3]
  g0, g1, g2 = D1.g[0:3]
  h0, h1, h2 = D1.h[0:3]
  F0, F1, F2 = D2.f[0:3]
  G0, G1, G2 = D2.g[0:3]
  H0, H1, H2 = D2.h[0:3]

  # Compute M
  #
  #       [ a1   a2   a3   a4   a5  ]
  #   M = [ a6   a7   a8   a9   a10 ]
  #       [ a11  a12  a13  a14  a15 ]
  #
  # Columns 4 and 5 are computed by
  #
  #   [ a4   a5  ]   [ 0  -F0  -G0 ] [ a1   a2  ] 
  #   [ a9   a10 ] = [ 1  -F1  -G1 ]*[ a6   a7  ]
  #   [ a14  a15 ]   [ 0  -F2  -G2 ] [ a11  a12 ]
  
  a1  = f0 - F0
  a2  = g0 - G0
  a3  = h0 - H0
  a6  = f1 - F1
  a7  = g1 - G1
  a8  = h1 - H1
  a11 = f2 - F2
  a12 = g2 - G2
  a13 = h2 - H2
  
  a4  =    - F0*a6  - G0*a11
  a5  =    - F0*a7  - G0*a12
  a9  = a1 - F1*a6  - G1*a11
  a10 = a2 - F1*a7  - G1*a12
  a14 =    - F2*a6  - G2*a11
  a15 =    - F2*a7  - G2*a12
  # Subtotal : 0I 12M 17A

  # print "M = "
  # print Matrix(C.K, 3, 5, [a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15])
  # print

  if (a1 == 0) and (a6 == 0) and (a11 == 0) :
    raise ValueError("Sum is not typical.")

  if (a1 == 0) :
    if (a6 != 0) :
      a1, a2, a3, a4, a5, a6, a7, a8, a9, a10 = \
          a6, a7, a8, a9, a10, a1, a2, a3, a4, a5
    else :
      a1, a2, a3, a4, a5, a11, a12, a13, a14, a15 = \
          a11, a12, a13, a14, a15, a1, a2, a3, a4, a5

  # Partially reduce M.
  #
  #        [ a1  a2  a3  a4  a5 ]
  #   M' = [ 0   b1  b2  b3  b4 ]
  #        [ 0   b5  b6  b7  b8 ]
  b1 = a1*a7  - a2*a6
  b2 = a1*a8  - a3*a6
  b3 = a1*a9  - a4*a6
  b4 = a1*a10 - a5*a6
  b5 = a1*a12 - a2*a11
  b6 = a1*a13 - a3*a11
  b7 = a1*a14 - a4*a11
  b8 = a1*a15 - a5*a11
  # Subtotal : 0I 16M 8A
  # Running total : 0I 28M 25A

  # print "M' = "
  # print Matrix(C.K, 3, 5, [a1, a2, a3, a4, a5, 0, b1, b2, b3, b4, 0, b5, b6, b7, b8])
  # print

  if (b1 == 0) and (b5 == 0) :
    raise ValueError("Sum is not typical.")
  
  if (b1 == 0) :
    b1, b2, b3, b4, b5, b6, b7, b8 = b5, b6, b7, b8, b1, b2, b3, b4
  
  # Reduce M to row echelon form
  #
  #           [ a1  a2  a3  a4  a5 ]
  #   M_ref = [  0  b1  b2  b3  b4 ]
  #           [  0   0  c1  c2  c3 ]
  c1 = b1*b6 - b2*b5
  c2 = b1*b7 - b3*b5
  c3 = b1*b8 - b4*b5
  # Subtotal : 0I 6M 3A
  # Running total : 0I 34M 28A
  
  # print "M_ref = "
  # print Matrix(C.K, 3, 5, [a1, a2, a3, a4, a5, 0, b1, b2, b3, b4, 0, 0, c1, c2, c3])
  # print


  if (c1 == 0) :
    raise ValueError("Sum is not typical.")

  # Reduce M even more via back-substitution
  #
  #         [ ABC    0    0  A1  A2 ]
  #   M'' = [   0  ABC    0  B1  B2 ]
  #         [   0    0  ABC  C1  C2 ]
  AB = a1*b1
  Z  = AB*c1
  C1 = AB*c2
  C2 = AB*c3
  B1 = a1*(b3*c1 - b2*c2) # TODO : This value in parentheses reused in A1
  B2 = a1*(b4*c1 - b2*c3) #      :  "    "     "  "           "      " A2
  A1 = b1*(a4*c1 - c2*a3) - a2*(b3*c1 - b2*c2)
  A2 = b1*(a5*c1 - c3*a3) - a2*(b4*c1 - b2*c3)
  # Subtotal : 0I 22M 8A
  # Running total : 0I 56M 36A

  # print("{}*M_rref = ".format(Z))
  # print Matrix(C.K, [
  #   [Z, 0, 0, A1, A2],
  #   [0, Z, 0, B1, B2],
  #   [0, 0, Z, C1, C2]])

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
  # Running total : 0I 74M 50A

  # Compute some inverses
  c3, c4, c5, c6, c7, c8 = C.coefficients()[3:]
  ZZt0      = u5^2 - Z*(u5*c8 - u4 + v5)
  if (ZZt0 == 0) :
    raise ValueError("Sum of divisors is non-typical.")
  ZZZt0     = Z*ZZt0
  ZZZt0_inv = 1/ZZZt0
  ZZt0_inv  = Z*ZZZt0_inv
  zeta      = ZZt0*ZZZt0_inv # 1/Z
  tau       = Z*Z*ZZt0_inv   # 1/t0
  # Subtotal : 1I 7M 3A
  # Running total : 1I 81M 53A

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
  # Running total : 1I 91M 53A
  
  # Compute ff, gg such that gg*u = ff*v (mod C)
  gg3 = u5
  ff2 = u5*(u5 - c8) + u4 - v5
  gg2 = v4 + v5*(u5 - c8) + tau*(u5*(u5*(u3 - c6) + v5*(u4 - c7) + c5 - v3) + v5*(u3 - v4) - u2)
  ff1 = u5*(u4 - c7) + gg2 + u3 - v4
  gg1 = u5*(c6 - u3) - (ff2*v5 - gg2*u5) + v3
  ff0 = c7*(ff2*v5 - gg2*u5) + u5*(u2 - c4) + gg2*u3 + gg1*u4 - ff2*v3 - ff1*v4 + u1 - v2
  gg0 = c6*(gg2*u5 - ff2*v5) + u5*(c3 - u1) - gg1*u3 + ff1*v3 + v1
  # Subtotal : 0I 25M 38A
  # Running total : 1I 116M 91A

  # Reduce gg modulo ff
  gg2 = gg2 - gg3*ff2
  gg1 = gg1 - gg3*ff1
  gg0 = gg0 - gg3*ff0
  # Subtotal : 0I 3M 3A
  # Running total : 1I 119M 94A

  # Compute third polynomial ...
  hh0 = tau*(ff0*gg1 + gg0*(gg2 - ff1))
  hh1 = tau*(gg1*gg2 - gg0)
  hh2 = gg1 + tau*(gg2*(gg2 - ff1) + ff0)
  # Subtotal : 0I 7M 6A
  # Running total : 1I 126M 100A

  return C34CurveDivisor(C, [[ff0, ff1, ff2, 1],
                       [gg0, gg1, gg2, 0, 1],
                       [hh0, hh1, hh2, 0, 0, 1]],
                       degree = 3, typ = 31, typical = True, reduced = True)

  
def add_11_11(D1, D2):
  """
    Add two divisors, D1 and D2, each of type 11 (degree 1).
    
    Divisors are assumed to be distinct, otherwise the doubling function should be used.
  """
  C = D1.C
  f, g = D1.f, D1.g
  F, G = D2.f, D2.g

  a1 = f[0] - F[0]
  a2 = g[0] - G[0]

  if (a1 != 0) :
    alpha = 1/a1
    r0 = -alpha*a2
    s0 = F[0]
    u0 = f[0]*r0 + g[0]
    u1 = r0
    v0 = f[0]*s0
    v1 = s0 + f[0]

    # D1 + D2 is of type 21
    # Total : 1I 3M 4A
    return C34CurveDivisor(C, [[u0, u1, 1], [v0, v1, 0, 1], []])

  else :
    r0 = G[0]
    v0 = g[0]*r0
    v2 = g[0] + r0
    
    # D1 + D2 is of type 22
    # Total : 0I 1M 3A
    return C34CurveDivisor(C, [copy(D1.f), [v0, 0, v2, 0, 0, 1], []])



def add_21_11(D1, D2) :
  C = D1.C
  f, g, h = D1.f, D1.g, D1.h
  F, G, H = D2.f, D2.g, D2.h
  new_f, new_g, new_h = [], [], []

  # Construct the matrix M
  #
  #   M = [ a1  a2  a3  a4 ]
  #
  # The colums are, respectively, the reductions of f, g, xf, and yf.
  # The last two values are not explicitly needed.

  a1 = f[0] - G[0] - f[1]*F[0]
  a2 = g[0] - F[0]*(g[1] - F[0])
  
  if a1 != 0 :
    # Compute reduced row echelon form of M
    #
    #   M_rref = [ 1  -r0  -s0  -t0 ]
    #
    # and compute
    #
    #  u = r0*f + g
    #  v = s0*f + x*f
    #  w = t0*f + y*f
    alpha = 1 / a1
    r0 = -alpha*a2
    s0 = F[0]
    t0 = G[0]

    u0 = f[0]*r0 + g[0]
    u1 = f[1]*r0 + g[1]
    u2 = r0
    v0 = f[0]*s0 - f[1]*u0
    v1 = f[1]*(s0 - u1) + f[0]
    v2 = s0 - f[1]*u2
    w0 = f[0]*t0 - f[1]*v0
    w1 = f[1]*(t0 - v1)
    w2 = t0 + f[0] - f[1]*v2

    # D1 + D2 is of type 31
    # Total 1I 13M 14A
    return C34CurveDivisor(C, [[u0, u1, u2, 1], [v0, v1, v2, 0, 1], [w0, w1, w2, 0, 0, 1]])
    
  elif a2 != 0 :
    # D1 + D2 is of type 32, generated by polynomials
    #
    #   u = f
    #   v = g*F
    v0 = g[0]*F[0]
    v1 = g[0] + g[1]*F[0]
    v3 = g[1] + F[0]

    # D1 + D2 is of type 32
    # Total : 0I 4M 6A
    return C34CurveDivisor(C, [copy(D1.f), [v0, v1, 0, v3, 0, 0, 1], []])
  
  else :
    # Divisors are non-disjoint.
    # So we have D1 = P + Q and D2 = P, for some points P and Q.
    # We must determine whether P and Q are distinct and compute 2P + Q or 3P as appropriate.
    # We have P = Q if and only if g[1] - 2F[0] = 0.
    
    if g[1] - 2*F[0] == 0 :
      return triple(D2)
    else :
      u0 = g[1] - F[0]
      v0 = f[0] - f[1]*u0
      Q = C34CurveDivisor(C, [[u0, 1], [v0, 0, 1], []])
      return Q + double(D2)



def add_21_21(D1, D2) :
  C = D1.C
  f, g = D1.f, D1.g
  F, G = D2.f, D2.g

  # Compute the matrix M
  #
  #          f   g  xf  yf   xg
  #          y  xx  xy  yy  xxx
  #   M = [ a1  a2  a3  a4  a5  ]
  #       [ a6  a7  a8  a9  a10 ]
  #
  # The columns are, from left to right, the reductions of f, g, xf, yf, xg modulo F, G.
  # The last three columns may be compute from the first two by
  #
  #   [ a3  a5  ] = [ 0  -G[0] ]*[ a1  a2 ]
  #   [ a8  a10 ]   [ 1  -G[1] ] [ a6  a7 ]
  #
  #   [ a4 ] = [ -F[0]  F[1]*G[0]        ]*[ a1 ]
  #   [ a9 ]   [ -F[1]  F[1]*G[1] - F[0] ] [ a6 ]
  a1 = f[0] - F[0]
  a2 = g[0] - G[0]
  a6 = f[1] - F[1]
  a7 = g[1] - G[1]

  a3  =    - G[0]*a6
  a5  =    - G[0]*a7
  a8  = a1 - G[1]*a6
  a10 = a2 - G[1]*a7

  a4 = - F[0]*a1 + F[1]*G[0]*a6
  a9 = - F[1]*a1 + (F[1]*G[1] - F[0])*a6
  
  aswap = 0
  # Subtotal : 0I 10M 8A

  if (a1 != 0) or (a6 != 0) :
    if (a1 == 0) :
      a1, a2, a3, a4, a5, a6, a7, a8, a9, a10 = a6, a7, a8, a9, a10, a1, a2, a3, a4, a5
      aswap = 1

    # Compute the row echelon form of M
    #
    #   M_ref = [ a1  a2  a3  a4  a5 ]
    #           [  0  b1  b2  b3  b4 ]
    b1 = a1*a7  - a2*a6
    b2 = a1*a8  - a3*a6
    b3 = a1*a9  - a4*a6
    b4 = a1*a10 - a5*a6
    # Subtotal : 0I 8M 4A

    if (b1 != 0) :
      # Compute the reduced row echelon form of M
      #
      #   M_rref = [ 1  0  -r0  -s0  -t0 ]
      #            [ 0  1  -r1  -s1  -t1 ]
      gamma = 1/(a1*b1)
      alpha = gamma*b1
      beta = gamma*a1
      r1 = -beta*b2
      s1 = -beta*b3
      t1 = -beta*b4
      r0 = -alpha*(a3 + a2*r1)
      s0 = -alpha*(a4 + a2*s1)
      t0 = -alpha*(a5 + a2*t1)
      # Subtotal : 1I 12M 3A

      # Compute D1 + D2 = <u, v, w> where
      #
      #   u = r0*f + r1*g + x*f
      #   v = s0*f + s1*g + y*f - f[1]*u
      #   w = t0*f + t1*g + x*g
      u0 = f[0]*r0 + g[0]*r1
      u1 = f[1]*r0 + g[1]*r1 + f[0]
      u2 = r0
      u3 = r1 + f[1]
      v0 = f[0]*s0 + g[0]*s1 - f[1]*u0
      v1 = f[1]*(s0 - u1) + g[1]*s1
      v2 = s0 + f[0] - f[1]*u2
      v3 = s1 - f[1]*u3
      w0 = f[0]*t0 + g[0]*t1
      w1 = f[1]*t0 + g[1]*t1 + g[0]
      w2 = t0
      w3 = t1 + g[1]
      # Subtotal : 0I 15M 15A

      # D1 + D2 is of type 41
      # Total : 1I 45M 30A
      return C34CurveDivisor(C, [[u0, u1, u2, u3, 1], [v0, v1, v2, v3, 0, 1], [w0, w1, w2, w3, 0, 0, 1]])
    
    elif (b2 != 0) :
      # M_ref is the matrix
      #
      #   M_ref = [ a1  a2  a3  a4  a5 ]
      #           [  0   0  b2  b3  b4 ]
      #
      # Compute its reduced row echelon form
      #
      #   M_ref = [ 1  -r0  0  -s0  * ]
      #           [ 0   0   1  -s1  * ]
      gamma = 1/(a1*b2)
      alpha = gamma*b2
      beta = gamma*a1
      s1 = -beta*b3
      r0 = -alpha*a2
      s0 = -alpha*(a4 + a3*s1)
      # Subtotal : 1I 7M 1A

      # Compute D1 + D2 = <u, v> where
      #
      # u = r0*f + g
      # v = s0*f + s1*x*f + y*f - f[1]*s1*u
      z = f[1]*s1
      u0 = f[0]*r0 + g[0]
      u1 = f[1]*r0 + g[1]
      u2 = r0
      v0 = f[0]*s0 - z*u0
      v1 = f[1]*s0 + f[0]*s1 - z*u1
      v2 = s0 + f[0] - z*u2
      v4 = s1 + f[1]
      # Subtotal : 0I 8M 8A

      # D1 + D2 is of type 43
      # Total : 1I 33M 21A
      return C34CurveDivisor(C, [[u0, u1, u2, 1], [v0, v1, v2, 0, v4, 1], []])

    elif (b3 != 0) :
      # M_ref is the matrix
      #
      #   M_ref = [ a1  a2  a3  a4  a5 ]
      #           [  0  0   0   b3  b4 ]
      #
      # Compute its reduced row echelon form
      #
      #   M_ref = [ 1  -r0  -s0  0  * ]
      #           [ 0   0    0   1  * ]
      alpha = 1/a1
      r0 = -alpha*a2
      s0 = -alpha*a3
      # Subtotal : 1I 2M 0A

      # Compute D1 + D2 = <u, v> where
      #
      # u = r0*f + g
      # v = s0*f + x*f - f[1]*u
      z = f[1]*s1
      u0 = f[0]*r0 + g[0]
      u1 = f[1]*r0 + g[1]
      u2 = r0
      v0 = f[0]*s0 - f[1]*u0
      v1 = f[1]*(s0 - u1) + f[0]
      v2 = s0 - f[1]*u2
      # Subtotal : 0I 7M 6A

      # D1 + D2 is of type 42
      # Total : 1I 27M 18A
      return C34CurveDivisor(C, [[u0, u1, u2, 1], [v0, v1, v2, 0, 1], []])

    else : 
      # M_ref is the matrix
      #
      #   M_ref = [ a1  a2  a3  a4  a5 ]
      #           [  0  0   0   0   0  ]
      #
      # Compute its reduced row echelon form
      #
      #   M_ref = [ 1  -r0  -s0  -t0  * ]
      #           [ 0   0    0    0   0 ]
      alpha = 1/a1
      r0 = -alpha*a2
      s0 = -alpha*a3 
      t0 = -alpha*a4

      # LCM(D1, D2) is type 31, generated by <u, v, w>, where
      #
      #   u = r0*f + g
      #   v = s0*f + x*f - f[1]*u
      #   w = t0*f + y*f - f[1]*v
      u0 = f[0]*r0 + g[0]
      u1 = f[1]*r0 + g[1]
      u2 = r0
      v0 = f[0]*s0 - f[1]*u0
      v1 = f[1]*(s0 - u1) + f[0]
      v2 = s0 - f[1]*u2
      w0 = f[0]*t0 - f[1]*v0
      w1 = f[1]*(t0 - v1)
      w2 = t0 + f[0] - f[1]*v2
      L = C34CurveDivisor(C, [[u0, u1, u2, 1], [v0, v1, v2, 0, 1], [w0, w1, w2, 0, 0, 1]])

      # GCD(D1, D2) is type 11, generated by <p, q> , where
      #
      #   p = x + a1/a6
      #   q = f mod p
      #
      # assuming a1 and a6 have not been swapped.
      if (aswap == 0) :
        mu = 1/a6
        p0 = mu*a1
      else :
        mu = 1/a1
        p0 = mu*a6
      q0 = f[0] - f[1]*p0
      G = C34CurveDivisor(C, [[p0, 1], [q0, 0, 1], []])

      return L + G

  else :
    if (a2 == 0) :
      a2, a5, a7, a10 = a7, a10, a2, a5
      aswap = 1
    # M is the matrix
    #
    #   M = [ 0  a2  0  0  a5  ]
    #       [ 0  a7  0  0  a10 ]
    #
    # Reduce it to
    #
    #   M_ref = [ 0  a2  0  0  a5 ]
    #           [ 0  0   0  0  b1 ]
    b1 = a2*a10 - a5*a7
    # Subtotal : 0I 2M 1A

    if (b1 != 0) :
      # D1 + D2 is type 44, principal, generated by f alone.
      # Total : 0I 12M 9A
      return C34CurveDivisor(C, [copy(D1.f), [], []])
    else :
      #   M_ref = [ 0  a2  0  0  a5 ]
      #           [ 0  0   0  0  b1 ]
      #
      # and
      #   M_rref = [ 0  1  0  0  -r0 ]
      #            [ 0  0  0  0   0  ]
      alpha = 1/a2
      r0 = -alpha*a5

      # LCM(D1, D2) is type 32, generated by (u, v) where
      #
      #   u = f
      #   v = r0*g + x*g
      v0 = f[0]*r0
      v1 = f[1]*r0 + g[0]
      v3 = r0 + g[1]

      L = C34CurveDivisor(C, [copy(D1.f), [v0, v1, 0, v3, 0, 0, 1], []])

      # GCD(D1, D2) is type 11, generated by (p, q) where
      #
      #   p = x + a2/a7
      #   q = f mod p
      #
      # assuming a2 and a7 have not been swapped.
      if (aswap == 0) :
        mu = 1/a7
        p0 = mu*a2
      else :
        mu = 1/a2
        p0 = mu*a7
      q0 = f[0] - f[1]*p0
      G = C34CurveDivisor(C, [[p0, 1], [q0, 0, 1], []])

      return flip(flip(L)) + G



def add_22_11(D1, D2) :
  C = D1.C
  f, g = D1.f, D1.g
  F, G = D2.f, D2.g

  # Construct the matrix M
  #
  #   M = [ a1  a2  a3  a4  a5 ]
  #
  # whose columns are the reductions of f, xf, yf, g, x^2f modulo F, G.
  # Only a1 and a4 are needed explicitly.
  a1 = f[0] - F[0]
  a4 = G[0]*(G[0] - g[2]) + g[0]
  # Subtotal : 0I 1M 3A

  if (a1 != 0) :
    # Compute the reduced row echelon form of M,
    #
    # M_rref = [ 1  -r0  -s0  -t0 * ]
    alpha = 1/a1
    r0 = F[0]
    s0 = G[0]
    t0 = -alpha*a4
    # Subtotal : 1I 1M 0A

    # Compute D1 + D2 = <u, v, w>, where
    #
    #   u = r0*f + x*f
    #   v = s0*f + y*f
    #   w = t0*f + g
    u0 = f[0]*r0
    u1 = r0 + f[0]
    v0 = f[0]*s0
    v1 = s0
    v2 = f[0]
    w0 = f[0]*t0 + g[0]
    w1 = t0
    w2 = g[2]
    # Subtotal : 0I 3M 2A

    # D1 + D2 is type 31, atypical
    # Total : 1I 5M 5A
    return C34CurveDivisor(C, [[u0, u1, 0, 1], [v0, v1, v2, 0, 1], [w0, w1, w2, 0, 0, 1]])

  elif (a4 != 0) :
    # M is the matrix
    #
    #   M = [ 0  0  0  a4  0 ]
    #
    # with reduced row echelon form
    #
    #   M_rref = [ 0  0  0  1  0 ]

    # D1 + D2 is of type 33 (principal)
    # Total : 0I 1M 3A
    return C34CurveDivisor(C, [copy(D1.f), [], []])
  
  else :
    # Divisors are non-disjoint
    # D1 = D2 + P for some point P.
    P = C34CurveDivisor(C, [copy(D1.f), [g[2] - G[0], 0, 1], []])
    
    # If D2 = P, then D1 + D2 = 3P
    # Otherwise, D1 + D2 is 2*D2 + P
    if (D2 == P) :
      return triple(P)
    else :
      return double(D2) + P



def add_21_22(D1, D2) :
  C = D1.C
  f, g = D1.f, D1.g
  F, G = D2.f, D2.g
  
  # Construct the matrix M
  #
  #   M = [ a1  a2  a3  a4  a5 ]
  #       [ 1   0   a6  a7  0  ]
  # 
  # whose columns are the reductions of f, g, xf, yf, xg modulo F, G.
  # The last three columns may be computed by
  #
  # [ a3  a5 ] = [ -F[0]     0  ]*[ a1  a2 ]
  # [ a6  0  ]   [    0   -F[0] ] [  1  0  ]
  #
  # [ a4 ] = [ 0  -G[0] ][ a1 ]
  # [ a7 ]   [ 1  -G[2] ][ 1  ]
  a1 = f[0] - f[1]*F[0]
  a2 = F[0]*(F[0] - g[1]) + g[0]

  #a3 = -F[0]*a1
  a5 = -F[0]*a2
  a6 = -F[0]

  a4 =    - G[0]
  a7 = a1 - G[2]
  # Subtotal : 0I 3M 2A

  # Compute the row echelon form of M
  #
  #   M_ref = [ 1  0   a6  a7  0  ]
  #           [ 0  b1  0   b2  b3 ]
  b1 = a2
  b2 = a4 - a7*a1
  b3 = a5
  # Subtotal : 0I 1M 1A
  
  if (b1 != 0) :
    # Compute the reduced row echelon form of M
    #
    #   M_rref = [ 1  0  -r0  -s0   0  ]
    #            [ 0  1   0   -s1  -t1 ]
    beta = 1/b1
    s1 = -beta*b2
    t1 = -beta*b3
    r0 = -a6
    s0 = -a7
    # Subtotal : 1I 2M 0A
    
    # Compute D1 + D2 = <u, v, w>, where
    #
    #   u = r0*f + xf
    #   v = s0*f + s1*g + yf - f[1]*u
    #   w = t1*g + xg
    u0 = f[0]*r0
    u1 = f[1]*r0 + f[0]
    u2 = r0
    u3 = f[1]
    v0 = f[0]*s0 + g[0]*s1 - f[1]*u0
    v1 = f[1]*(s0 - u1) + g[1]*s1
    v2 = s0 + f[0] - f[1]*u2
    v3 = s1 - f[1]*u3
    w0 = g[0]*t1
    w1 = g[1]*t1 + g[0]
    w3 = t1 + g[1]
    # Subtotal : 0I 11M 10A

    # D1 + D2 is type 41
    # Total : 1I 17M 13A
    return C34CurveDivisor(C, [[u0, u1, u2, u3, 1], [v0, v1, v2, v3, 0, 1], [w0, w1, 0, w3, 0, 0, 1]])

  elif (b2 != 0) :
    # M_ref is the matrix
    #
    #   M_ref = [ 1  0  a6  a7  0 ]
    #           [ 0  0  0   b2  0 ]
    #
    # with reduced row echelon form
    #
    #   M_rref = [ 1  0  -r0  0  0 ]
    #            [ 0  0   0   1  0 ]
    r0 = -a6
    
    # Compute D1 + D2 = <u, v>, where
    # 
    #   u = g
    #   v = r0*f + xf - f[1]*u
    v0 = f[0]*r0 - f[1]*g[0]
    v1 = f[1]*(r0 - g[1]) + f[0]
    v2 = r0
    # Subtotal : 0I 2M 1A

    # D1 + D2 is type 42
    # Total : 0I 6M 4A
    return C34CurveDivisor(C, [copy(D1.g), [v0, v1, v2, 0, 1], []])

  else :
    # M_ref is the matrix
    #
    #   M_ref = [ 1  0  a6  a7  0 ]
    #           [ 0  0  0   0   0 ]
    #
    # and M_rref is the matrix
    #
    #   M_rref = [ 1  0  -r0  -s0  0 ]
    #            [ 0  0   0    0   0 ]
    r0 = -a6
    s0 = -a7

    # LCM(D1, D2) is type 31, generated by
    #
    #   u = g
    #   v = r0*f + x*f - f[1]*u
    #   w = s0*f + y*f - f[1]*v
    v0 = f[0]*r0 - f[1]*g[0]
    v1 = f[1]*(r0 - g[1]) + f[0]
    v2 = r0
    w0 = f[0]*s0 - f[1]*v0
    w1 = f[1]*(s0 - v1)
    w2 = s0 + f[0] - f[1]*v2
    L = C34CurveDivisor(C, [copy(D1.g), [v0, v1, v2, 0, 1], [w0, w1, w2, 0, 0, 1]])

    # GCD(D1, D2) is type 11, generated by
    #
    #   p = F
    #   q = y + a1
    G = C34CurveDivisor(C, [copy(D2.f), [a1, 0, 1], []])
    return L + G



def add_22_22(D1, D2) :
  C = D1.C
  f, g = D1.f, D1.g
  F, G = D2.f, D2.g

  # Construct the matrix M
  #
  #   M = [ a1  a2  0   a3  ]
  #       [  0  0   a1  a5  ]
  #
  # where the columns, from left to right, are the reductions of f, xf, yf, g, x^2f modulo F, G.
  # Column 2 is simply -F[0] times columns 1 and does not need to be computed explicitly.
  a1 = f[0] - F[0]
  a3 = g[0] - G[0]
  a5 = g[2] - G[2]
  # Subtotal : 0I 0M 3A

  if (a1 != 0) :
    # Reduce M to reduced row echelon form
    #
    #   M_rref = [ 1  -r0  0  -s0  * ]
    #              0   0   1  -s1  * ]
    alpha = 1/a1
    r0 = F[0]
    s1 = -alpha*a5
    s0 = -alpha*a3
    # Subtotal : 1I 2M 0A

    # D1 + D2 is generated by <u, v>, where
    #
    #   u = r0*f + xf
    #   v = s0*f + s1*yf + g
    u0 = f[0]*r0
    u1 = r0 + f[0]
    v0 = f[0]*s0 + g[0]
    v1 = s0
    v2 = f[0]*s1 + g[2]
    v4 = s1
    # Subtotal : 0I 3M 3A

    # D1 + D2 is of type 43
    # Total : 1I 5M 6A
    return C34CurveDivisor(C, [[u0, u1, 0, 1], [v0, v1, v2, 0, v4, 1], []])

  else :
    # M is the matrix
    #
    #   M = [ 0  0  0  a3  0 ]
    #       [ 0  0  0  a5  0 ]
    #
    # Since g and G are distinct, a3 and a5 are not both zero. So M has RREF
    #
    #   M_rref = [ 0  0  0  1  0 ]
    #            [ 0  0  0  0  0 ]
    #
    # L = LCM(D1, D2) is of type 33 and G = GCD(D1, D2) is of type 11.
    # G is generated by <p, q>, where
    #
    #   p = f
    #   q = y + a3/a5
    #
    # We compute G and return  D1 + D2 = L + G = G.
    q0 = a3*(1/a5)

    # D1 + D2 is of type 11
    # Total : 1I 1M 3A
    return C34CurveDivisor(C, [copy(D1.f), [q0, 0, 1], []])



def add_31_11(D1, D2):
  C = D1.C
  f, g, h = D1.f, D1.g, D1.h
  F, G    = D2.f, D2.g

  # Compute the matrix M,
  #
  # M = [ a1  a2  a3  a4 ]
  #
  # Where the first three elements are the reductions of f, g, h modulo (F, G),
  # i.e. the values of f, g, h at the point defining D2.
  # The last element is the reduction of x*f, which is merely = -F[0]*a1
  # The only time we actually need a4, we end up dividing it by a1, giving -F[0],
  # so there is no need to actually compute a4 itself at all.
  a1 = -F[0]*(f[1] - F[0]) - G[0]*f[2] + f[0]
  a2 = -F[0]*(g[1] - G[0]) - G[0]*g[2] + g[0]
  a3 = -G[0]*(h[2] - G[0]) - F[0]*h[1] + h[0]
  # Subtotal : 0I 6M 9A

  if (a1 != 0) :
    # Compute the reduced row echelon form of M,
    #
    # M_rref = [ 1  -r0  -s0  -t0 ]
    alpha = 1/a1
    r0 = -alpha*a2
    s0 = -alpha*a3
    t0 = F[0]
    # Subtotal : 1I 2M 0A

    u0 = f[0]*r0 + g[0]
    u1 = f[1]*r0 + g[1]
    u2 = f[2]*r0 + g[2]
    u3 = r0
    v0 = f[0]*s0 + h[0]
    v1 = f[1]*s0 + h[1]
    v2 = f[2]*s0 + h[2]
    v3 = s0
    w0 = f[0]*t0 - f[2]*u0
    w1 = f[1]*t0 + f[0] - f[2]*u1
    w2 = f[2]*(t0 - u2)
    w3 = t0 + f[1] - f[2]*u3
    # Subtotal : 0I 12M 12A
    
    # D1 + D2 is of type 41
    # Total : 1I 20M 21A
    return C34CurveDivisor(C, [[u0, u1, u2, u3, 1], [v0, v1, v2, v3, 0, 1], [w0, w1, w2, w3, 0, 0, 1]])

  elif (a2 != 0) :
    # M = [ 0  a2  a3  0 ]
    #
    # Compute the reduced row echelon form of M,
    #
    # M_rref = [ 0  1  -r0  0 ]
    alpha = 1/a2
    r0 = -alpha*a3
    u0 = g[0]*r0 + h[0]
    u1 = g[1]*r0 + h[1]
    u2 = g[2]*r0 + h[2]
    u4 = r0
    # Subtotal : 1I 4M 3A
    
    # D1 + D2 is of type 43
    # Total : 1I 10M 12A
    return C34CurveDivisor(C, [copy(D1.f), [u0, u1, u2, 0, u4, 1], []])

  elif (a3 != 0) :
    # M = [ 0  0  1  0 ]

    # D1 + D2 is of type 42
    # Total : 0I 6M 9A
    return C34CurveDivisor(C, [copy(D1.f), copy(D1.g), []])
  
  else :
    # M = [ 0  0  0  0 ]
    #
    # In this case, D2 = P for some point P in the support of D1.
    # We find a degree 2 divisor A such that D1 = A + P
    # If P is also in the support of A, then we find the point Q such that A = P + Q
    
    # Find (p, q) such that (f, g, h) = (F, G)(p, q)
    p0, p1, q0, q1, q2 = 0, 0, 0, 0, 0
    if (D1.typical) :
      # If D1 is typical, then div(p, q) is type 21, f[2] =/= 0, and the following solution exists
      p1 = - (g[2] - F[0]) / f[2]
      p0 = g[1] + p1*(f[1] - F[0])
      q1 = g[2] + f[1] - F[0]
      q0 = f[1]*g[2] - f[2]*g[1] - F[0]*q1 + f[0]
      A = C34CurveDivisor(C, [[p0, p1, 1], [q0, q1, 0, 1], []])
    else :
      # TODO : This section is painfully inefficient
      A = D1.slow_add(D2.slow_flip())
      # D1 is not typical
      if (4*f[0] != f[1]^2) : # TODO : Replace this check with (f[1] != 2*F[0])?
        # f has two distinct rational roots
        # f = (x + F[0])*(x + f[1] - F[0])
        if (F[0] == g[2]) :
          # A is type 21
          assert A.type == 21, "A is not type 21."
          Q = C34CurveDivisor(C, [[g[2], 1], [h[2] - G[0], 0, 1], []])
          R = C34CurveDivisor(C, [[f[1] - g[2], 1], [g[1], 0, 1], []])
          A = Q + R
        else :
          # A is type 22
          assert A.type == 22, "A is not type 22."
          p0 = f[1] - F[0]
          q0 = h[0] + h[1]*(F[0] - f[1])
          q2 = h[2]
          A = C34CurveDivisor(C, [[p0, 1], [q0, 0, q2, 0, 0, 1], []])
      else :
        # f has a rational double root
        # f = (x + F[0])^2
        if (G[0] == g[1]) :
          # A is type 22
          assert A.type == 22, "A is not type 22."
          p0 = f[1] - F[0]
          q0 = h[0] + h[1]*(F[0] - f[1])
          q2 = h[2]
          A = C34CurveDivisor(C, [[p0, 1], [q0, 0, q2, 0, 0, 1], []])
        else :
          # A is type 21
          assert A.type == 21, "A is not type 21."
          Q = C34CurveDivisor(C, [[g[2], 1], [h[2] - G[0], 0, 1], []])
          R = C34CurveDivisor(C, [[f[1] - F[0], 1], [g[1], 0, 1], []])
          A = Q + R
      
    assert A.slow_add(D2) == D1, "A + D2 =/= D1"

    # Check if P is in the support of A.
    # If not, return D1 + D2 = A + 2P
    # Otherwise, find Q such that A = P + Q
    if (A.type == 21) and (- G[0] - A.f[1]*F[0] + A.f[0] == 0) and (F[0]*(F[0] - A.g[1]) + A.g[0] == 0) :
      r0 = A.g[1] - F[0]
      s0 = A.f[1]*(F[0] - A.g[1]) + A.f[0]
      Q = C34CurveDivisor(C, [[r0, 1], [s0, 0, 1], []])
      # If P =/= Q, then return D1 + D2 = Q + 3P. Otherwise, return 4P
      assert Q.slow_add(D2) == A, "Q + D2 =/= A"
      if (Q != D2) :
        return Q + flip(flip(triple(D2)))
      else :
        return double(double(D2))
    elif (A.type == 22) and (p0 == F[0]) and (G[0]*(G[0] - A.g[2]) + A.g[0] == 0) :
      r0 = A.f[0]
      s0 = A.g[2] - G[0]
      Q = C34CurveDivisor(C, [[r0, 1], [s0, 0, 1], []])
      # If P =/= Q, then return D1 + D2 = Q + 3P. Otherwise, return 4P
      assert Q.slow_add(D2) == A, "Q + D2 =/= A"
      if (Q != D2) :
        return Q + flip(flip(triple(D2)))
      else :
        return double(double(D2))
    else :
      return A + double(D2)



def add_31_21(D1, D2):
  C = D1.C
  K = C.K
  f, g, h = D1.f, D1.g, D1.h
  F, G    = D2.f, D2.g

  # Compute the matrix M,
  #
  # M = [ a1  a2  a3  a4   a5   a6  ]
  #     [ a7  a8  a9  a10  a11  a12 ]
  #
  # Where the first three columns are the reductions of f, g, h modulo (F, G)
  # The last three columns are the reductions of xf, xg, xh.
  # The last three may be computed via
  #
  # [ a4   a5   a6  ] = [ 0  -G[0] ][ a1  a2  a3 ]
  # [ a10  a11  a12 ]   [ 1  -G[1] ][ a7  a8  a9 ]
  
  a1 = f[0] - G[0] - f[2]*F[0]
  a2 = g[0] + F[1]*G[0] - g[2]*F[0]
  a3 = F[0]*(F[0] - h[2]) + h[0] - F[1]*F[1]*G[0]
  a7 = f[1] - G[1] - f[2]*F[1]
  a8 = g[1] - F[0] + F[1]*(G[1] - g[2])
  a9 = F[1]*(2*F[0] - h[2] - F[1]*G[1]) + h[1]
  
  a4  =    - G[0]*a7
  a5  =    - G[0]*a8
  a6  =    - G[0]*a9
  a10 = a1 - G[1]*a7
  a11 = a2 - G[1]*a8
  a12 = a3 - G[1]*a9
  
  aswap = 0 # Update this if we need to swap rows of M
  # Subtotal : 0I 17M 18A
  
  #print "M = "
  #print Matrix(K, 2, 6, [a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12])
  #print
  
  if (a1 != 0) or (a7 != 0) :
    if (a1 == 0) :
      a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12 = a7, a8, a9, a10, a11, a12, a1, a2, a3, a4, a5, a6
      aswap = 1

    # Reduce M to its row echelon form
    #
    # M_ref = [ a1  a2  a3  a4  a5  a6 ]
    #         [ 0   b1  b2  b3  b4  b5 ]
    b1 = a1*a8  - a2*a7
    b2 = a1*a9  - a3*a7
    b3 = a1*a10 - a4*a7
    b4 = a1*a11 - a5*a7
    b5 = a1*a12 - a6*a7
    # Subtotal : 0I 10M 5A
    
    if (b1 != 0) :
      # Reduce M_ref to its reduced row echelon form
      #
      # M_rref = [ 1  0  -r0  -s0  -t0  * ]
      #          [ 0  1  -r1  -s1  -t1  * ]
      gamma = 1/(a1*b1)
      alpha = gamma*b1
      beta = gamma*a1
      r1 = -beta*b2
      s1 = -beta*b3
      t1 = -beta*b4
      r0 = -alpha*(a3 + a2*r1)
      s0 = -alpha*(a4 + a2*s1)
      t0 = -alpha*(a5 + a2*t1)
      # Subtotal : 1I 12M 3A

      u0 = f[0]*r0 + g[0]*r1 + h[0]
      u1 = f[1]*r0 + g[1]*r1 + h[1]
      u2 = f[2]*r0 + g[2]*r1 + h[2]
      u3 = r0
      u4 = r1
      v0 = f[0]*s0 + g[0]*s1
      v1 = f[1]*s0 + g[1]*s1 + f[0]
      v2 = f[2]*s0 + g[2]*s1
      v3 = s0 + f[1]
      v4 = s1 + f[2]
      w0 = f[0]*t0 + g[0]*t1
      w1 = f[1]*t0 + g[1]*t1 + g[0]
      w2 = f[2]*t0 + g[2]*t1
      w3 = t0 + g[1]
      w4 = t1 + g[2]
      # Subtotal : 0I 18M 18A
      
      # D1 + D2 is of type 51
      # Total : 1I 57M 44A
      return C34CurveDivisor(C, [[u0, u1, u2, u3, u4, 1],
                           [v0, v1, v2, v3, v4, 0, 1],
                           [w0, w1, w2, w3, w4, 0, 0, 1]])
      
    elif (b2 != 0) :
      # M_ref = [ a1  a2  a3  a4  a5  a6 ]
      #         [ 0   0   b2  b3  b4  b5 ]
      #
      # Reduce M_ref to its reduced row echelon form
      #
      # M_rref = [ 1  -r0  0  -s0  *  * ]
      #          [ 0   0   1  -s1  *  * ]
      gamma = 1/(a1*b2)
      alpha = gamma*b2
      beta = gamma*a1
      s1 = -beta*b3
      r0 = -alpha*a2
      s0 = -alpha*(a4 + a3*s1)
      # Subtotal : 1I 7M 1A

      u0 = f[0]*r0 + g[0]
      u1 = f[1]*r0 + g[1]
      u2 = f[2]*r0 + g[2]
      u3 = r0
      v0 = f[0]*s0 + h[0]*s1 - f[2]*u0
      v1 = f[1]*s0 + h[1]*s1 + f[0] - f[2]*u1
      v2 = f[2]*(s0 - u2) + h[2]*s1
      v3 = s0 + f[1] - f[2]*u3
      v5 = s1
      # Subtotal : 0I 12M 12A

      # D1 + D2 is of type 53
      # Total : 1I 46M 36A
      return C34CurveDivisor(C, [[u0, u1, u2, u3, 1], [v0, v1, v2, v3, 0, v5, 1], []])
      
    elif (b3 != 0) :
      # M_ref = [ a1  a2  a3  a4  a5  a6 ]
      #         [ 0   0   0   b3  b4  b5 ]
      #
      # Reduce M_ref to its reduced row echelon form
      #
      # M_rref = [ 1  -r0  -s0  0  *  * ]
      #          [ 0   0    0   1  *  * ]
      alpha = 1/a1
      r0 = -alpha*a2
      s0 = -alpha*a3
      # Subtotal : 1I 2M 0A

      u0 = f[0]*r0 + g[0]
      u1 = f[1]*r0 + g[1]
      u2 = f[2]*r0 + g[2]
      u3 = r0
      v0 = f[0]*s0 + h[0]
      v1 = f[1]*s0 + h[1]
      v2 = f[2]*s0 + h[2]
      v3 = s0
      # Subtotal : 0I 6M 6A

      # D1 + D2 is of type 52
      # Total : 1I 35M 29A
      return C34CurveDivisor(C, [[u0, u1, u2, u3, 1], [v0, v1, v2, v3, 0, 1], []])

    else :
      assert b4 == b5 == 0
      # M_ref = [ a1  a2  a3  a4  a5  a6 ]
      #         [ 0   0   0   0   0   0  ]
      #
      # Reduce M_ref to its reduced row echelon form
      #
      # M_rref = [ 1  -r0  -s0  -t0  *  * ]
      #          [ 0   0    0    0   0  0 ]
      #
      # Compute D1 + D2 = GCD(D1, D2) + LCM(D1, D2)
      # LCM(D1, D2) is of type 41, given by (u, v, w) where
      #
      #   u = r0*f + g
      #   v = s0*f + h
      #   w = t0*f + x*f (mod u)
      #
      # GCD(D1, D2) is of type 11, given by (p, q) where
      #
      #  p = x + a1/a7
      #  q = F (mod p)
      #
      # assuming a1 and a7 have not been swapped.
      alpha = 1/a1
      r0 = -alpha*a2
      s0 = -alpha*a3
      t0 = -alpha*a4
      u0 = f[0]*r0 + g[0]
      u1 = f[1]*r0 + g[1]
      u2 = f[2]*r0 + g[2]
      u3 = r0
      v0 = f[0]*s0 + h[0]
      v1 = f[1]*s0 + h[1]
      v2 = f[2]*s0 + h[2]
      v3 = s0
      w0 = f[0]*t0 - f[2]*u0
      w1 = f[1]*t0 + f[0] - f[2]*u1
      w2 = f[2]*(t0 - u2)
      w3 = t0 + f[1] - f[2]*u3
      L = C34CurveDivisor(C, [[u0, u1, u2, u3, 1], [v0, v1, v2, v3, 0, 1], [w0, w1, w2, w3, 0, 0, 1]])

      if (aswap == 0) :
        mu = 1/a7
        p0 = mu*a1
      else :
        mu = 1/a1
        p0 = mu*a7
      q0 = F[0] - F[1]*p0
      G = C34CurveDivisor(C, [[p0, 1], [q0, 0, 1], []])
      return flip(flip(L)) + G
      
  elif (a2 != 0) or (a8 != 0) :
    if (a2 == 0) :
      a2, a3, a5, a6, a8, a9, a11, a12 = a8, a9, a11, a12, a2, a3, a5, a6
      aswap = 1
    # M = [ 0  a2  a3  0  a5   a6  ]
    #     [ 0  a8  a9  0  a11  a12 ]
    #
    # Reduce M to its row echelon form
    #
    # M_ref = [ 0  a2  a3  0  a5  a6 ]
    #         [ 0  0   b1  0  b2  b3 ]
    b1 = a2*a9  - a3*a8
    b2 = a2*a11 - a5*a8
    b3 = a2*a12 - a6*a8
    # Subtotal : 0I 6M 3A

    if (b1 != 0) :
      # Reduce M_ref to its reduced row echelon form
      #
      # M_ref = [ 0  1  0  0  *  -r0 ]
      #         [ 0  0  1  0  *  -r1 ]
      gamma = 1/(a2*b1)
      alpha = gamma*b1
      beta = gamma*a2
      r1 = -beta*b3
      r0 = -alpha*(a6 + a3*r1)
      # Subtotal : 1I 6M 1A

      u0 = g[0]*r0 + h[0]*r1 - h[1]*f[0]
      u1 = g[1]*r0 + h[1]*(r1 - f[1]) + h[0]
      u2 = g[2]*r0 + h[2]*r1 - h[1]*f[2]
      u4 = r0 + h[2]
      u5 = r1
      # Subtotal : 0I 8M 8A

      # D1 + D2 is of type 54
      # Total : 1I 37M 30A
      return C34CurveDivisor(C, [copy(D1.f), [u0, u1, u2, 0, u4, u5, 0, 0, 1], []])
    else :
      assert b2 == b3 == 0
      # M_ref = [ 0  a2  a3  0  a5  a6 ]
      #         [ 0  0   0   0  0   0  ]
      #
      # Reduce M_ref to its reduced row echelon form
      #
      # M_rref = [ 0  1  -r0  0  *  * ]
      #          [ 0  0   0   0  0  0 ]
      #
      # Compute D1 + D2 = GCD(D1, D2) + LCM(D1, D2)
      # LCM(D1, D2) is of type 43, given by (f, u) where
      #
      #   u = r0*g + h
      #
      # GCD(D1, D2) is of type 11, given by (p, q) where
      #
      #  p = x + a2/a8
      #  q = F (mod p)
      #
      # assuming a2 and a8 have not been swapped.
      alpha = 1/a2
      r0 = -alpha*a3
      u0 = g[0]*r0 + h[0]
      u1 = g[1]*r0 + h[1]
      u2 = g[2]*r0 + h[2]
      u4 = r0
      L = C34CurveDivisor(C, [copy(D1.f), [u0, u1, u2, 0, u4, 1], []])

      if (aswap == 0) :
        mu = 1/a8
        p0 = mu*a2
      else :
        mu = 1/a2
        p0 = mu*a8
      q0 = F[0] - F[1]*p0
      G = C34CurveDivisor(C, [[p0, 1], [q0, 0, 1], []])
      return flip(flip(L)) + G
    
  elif (a3 != 0) or (a9 != 0) :
    # M = [ 0  0  a3  0  0  a6  ]
    #     [ 0  0  a9  0  0  a12 ]
    #
    # Compute D1 + D2 = GCD(D1, D2) + LCM(D1, D2)
    # LCM(D1, D2) is of type 42, given by (f, g).
    # GCD(D1, D2) is of type 11, given by (p, q) where
    #
    #  p = x + a3/a9
    #  q = F (mod p)
    L = C34CurveDivisor(C, [copy(D1.f), copy(D1.g), []])

    mu = 1/a9
    p0 = mu*a3
    q0 = F[0] - F[1]*p0
    G = C34CurveDivisor(C, [[p0, 1], [q0, 0, 1], []])
    return flip(flip(L)) + G

  else :
    # M = [ 0  0  0  0  0  0 ]
    #     [ 0  0  0  0  0  0 ]
    #
    # In this case, D1 = D2 + P for some point (i.e. type 11 divisor) P.
    # We compute P by solving for polynomials p and q such that
    #
    #   (f, g, h) = (F, G)*(p, q)
    P = C34CurveDivisor(D1.C, [ [g[2] + f[2]*F[1], K.one()], [h[2] + g[2]*F[1] - F[0], K.zero(), K.one()], [] ])
    return double(D2) + P



def add_31_22(D1, D2) :
  C = D1.C
  f, g, h = D1.f, D1.g, D1.h
  F, G    = D2.f, D2.g

  # Construct the matrix M
  #
  #   M = [ a1  a2  a3  a4   a5   a6  ]
  #       [ a7  a8  a9  a10  a11  a12 ]
  #
  # where the columns are, from left to right, the reductions of
  # f, g, h, xf, xg, xh modulo F, G.
  # The last three columns are computed from the first three via
  #
  #   [ a4   a5   a6  ] = [ -F[0]     0  ]*[ a1  a2  a3 ]
  #   [ a10  a11  a12 ]   [    0   -F[0] ] [ a7  a8  a9 ]
  #
  # In most cases, a6 and a12 are not needed. They are computed only in the cases they are.
  a1 = F[0]*(F[0] - f[1]) + f[0]
  a2 = g[0] - g[1]*F[0]
  a3 = h[0] - G[0] - h[1]*F[0]
  a7 = f[2]
  a8 = g[2] - F[0]
  a9 = h[2] - G[2]

  a4  = -F[0]*a1
  a5  = -F[0]*a2
  # a6  = -F[0]*a3
  a10 = -F[0]*a7
  a11 = -F[0]*a8
  # a12 = -F[0]*a9

  # print "M = "
  # print Matrix(K, 2, 6, [a1, a2, a3, a4, a5, -F[0]*a3, a7, a8, a9, a10, a11, -F[0]*a9])
  # print

  aswap = 0
  # Subtotal : 0I 7M 7A

  if (a1 != 0) or (a7 != 0) :
    if (a1 == 0) :
      a1, a2, a3, a4, a5, a7, a8, a9, a10, a11 = a7, a8, a9, a10, a11, a1, a2, a3, a4, a5
      aswap = 1

    # Compute the row echelon form of M
    #
    #   M_ref = [ a1  a2  a3  a4  a5  * ]
    #           [ 0   b1  b2  b3  b4  * ]
    b1 = a1*a8  - a2*a7
    b2 = a1*a9  - a3*a7
    b3 = a1*a10 - a4*a7
    b4 = a1*a11 - a5*a7
    # Subtotal : 0I 8M 4A

    if (b1 != 0) :
      # Compute the reduced row echelon form of M
      #
      #   M_rref = [ 1  0  -r0  -s0  -t0 ]
      #            [ 0  1  -r1  -s1  -t1 ]
      gamma = 1/(a1*b1)
      alpha = gamma*b1
      beta = gamma*a1
      r1 = -beta*b2
      s1 = -beta*b3
      t1 = -beta*b4
      r0 = -alpha*(a3 + a2*r1)
      s0 = -alpha*(a4 + a2*s1)
      t0 = -alpha*(a5 + a2*t1)
      # Subtotal : 1I 12M 3A

      # Compute D1 + D2 = <u, v, w>, where
      #
      #   u = r0*f + r1*g + h
      #   v = s0*f + s1*g + xf
      #   w = t0*f + t1*g + xg
      u0 = f[0]*r0 + g[0]*r1 + h[0]
      u1 = f[1]*r0 + g[1]*r1 + h[1]
      u2 = f[2]*r0 + g[2]*r1 + h[2]
      u3 = r0
      u4 = r1
      v0 = f[0]*s0 + g[0]*s1
      v1 = f[1]*s0 + g[1]*s1 + f[0]
      v2 = f[2]*s0 + g[2]*s1
      v3 = s0 + f[1]
      v4 = s1 + f[2]
      w0 = f[0]*t0 + g[0]*t1
      w1 = f[1]*t0 + g[1]*t1 + g[0]
      w2 = f[2]*t0 + g[2]*t1
      w3 = t0 + g[1]
      w4 = t1 + g[2]
      # Subtotal : 0I 18M 18A

      # D1 + D2 is of type 51
      # Total : 1I 45M 32A
      return C34CurveDivisor(C, [[u0, u1, u2, u3, u4, 1],
                           [v0, v1, v2, v3, v4, 0, 1],
                           [w0, w1, w2, w3, w4, 0, 0, 1]])

    if (b2 != 0) :
      # M_ref is the matrix
      #
      #   M_ref = [ a1  a2  a3  a4  a5  * ]
      #           [ 0   0   b2  b3  b4  * ]
      #
      # Compute the reduced row echelon form of M
      #
      #   M_rref = [ 1  -r0  0  -s0  *  * ]
      #            [ 0  0    1  -s1  *  * ]
      gamma = 1/(a1*b2)
      alpha = gamma*b2
      beta = gamma*a1
      s1 = -beta*b3
      r0 = -alpha*a2
      s0 = -alpha*(a4 + a3*s1)
      # Subtotal : 1I 7M 1A

      # Compute D1 + D2 = <u, v>, where
      #
      #   u = r0*f + g
      #   v = s0*f + s1*h + xf - f[2]*u
      u0 = f[0]*r0 + g[0]
      u1 = f[1]*r0 + g[1]
      u2 = f[2]*r0 + g[2]
      u3 = r0
      v0 = f[0]*s0 + h[0]*s1 - f[2]*u0
      v1 = f[1]*s0 + h[1]*s1 + f[0] - f[2]*u1
      v2 = f[2]*(s0 - u2) + h[2]*s1
      v3 = s0 + f[1] - f[2]*u3
      v5 = s1
      # Subtotal : 0I 12M 12A

      # D1 + D2 is of type 53
      # Total : 1I 34M 24A
      return C34CurveDivisor(C, [[u0, u1, u2, u3, 1], 
                           [v0, v1, v2, v3, 0, v5, 1]])

    if (b3 != 0) :
      # M_ref is the matrix
      #
      #   M_ref = [ a1  a2  a3  a4  a5  * ]
      #           [ 0   0   0   b3  b4  * ]
      #
      # Compute the reduced row echelon form of M
      #
      #   M_rref = [ 1  -r0  -s0  0  *  * ]
      #            [ 0  0     0   1  *  * ]
      alpha = 1/a1
      r0 = -alpha*a2
      s0 = -alpha*a3
      # Subtotal : 1I 2M 0A

      # Compute D1 + D2 = <u, v>, where
      #
      #   u = r0*f + g
      #   v = s0*f + h
      u0 = f[0]*r0 + g[0]
      u1 = f[1]*r0 + g[1]
      u2 = f[2]*r0 + g[2]
      u3 = r0
      v0 = f[0]*s0 + h[0]
      v1 = f[1]*s0 + h[1]
      v2 = f[2]*s0 + h[2]
      v3 = s0
      # Subtotal : 0I 6M 6A

      # D1 + D2 is of type 52
      # Total : 1I 23M 17A
      return C34CurveDivisor(C, [[u0, u1, u2, u3, 1],
                           [v0, v1, v2, v3, 0, 1]])

    else :
      # M_ref is the matrix
      #
      #   M_ref = [ a1  a2  a3  a4  a5  a6 ]
      #           [ 0   0   0   0   0   0  ]
      #
      # Compute the reduced row echelon form of M
      #
      #   M_rref = [ 1  -r0  -s0  -t0  *  * ]
      #            [ 0  0     0    0   0  0 ]
      alpha = 1/a1
      r0 = -alpha*a2
      s0 = -alpha*a3
      t0 = -alpha*a4
      # Subtotal : 1I 3M

      # LCM(D1, D2) is of type 41, generated by <u, v, w> where
      #
      #   u = r0*f + g
      #   v = s0*f + h
      #   w = t0*f + xf - f[2]*u
      u0 = f[0]*r0 + g[0]
      u1 = f[1]*r0 + g[1]
      u2 = f[2]*r0 + g[2]
      u3 = r0
      v0 = f[0]*s0 + h[0]
      v1 = f[1]*s0 + h[1]
      v2 = f[2]*s0 + h[2]
      v3 = s0
      w0 = f[0]*t0 - f[2]*u0
      w1 = f[1]*t0 + f[0] - f[2]*u1
      w2 = f[2]*(t0 - u2)
      w3 = t0 + f[1] - f[2]*u3
      L = C34CurveDivisor(C, [[u0, u1, u2, u3, 1],
                        [v0, v1, v2, v3, 0, 1],
                        [w0, w1, w2, w3, 0, 0, 1]])
      # Subtotal : 0I 12M 12A

      # GCD(D1, D2) is of type 11, generated by <p, q> where
      #
      #   p = F
      #   q = y + a1/a7
      #
      # assuming a1 and a7 have not been swapped
      if (aswap == 0) :
        mu = 1/a7
        q0 = mu*a1
      else :
        mu = 1/a1
        q0 = mu*a7
      G = C34CurveDivisor(C, [copy(D2.f), [q0, 0, 1], []])
      # Subtotal : 1I 1M 0A

      return flip(flip(L)) + G

  elif (a2 != 0) or (a8 != 0) :
    # M is the matrix
    #
    #   M = [ 0  a2  a3  0  a5   a6  ]
    #       [ 0  a8  a9  0  a11  a12 ]
    #
    # Compute the row echelon form of M
    #
    #   M_ref = [ 0  a2  a3  0  a5  a6 ]
    #           [ 0  0   b1  0  b2  b3 ]
    a6  = -F[0]*a3
    a12 = -F[0]*a9
    if (a2 == 0) :
      a2, a3, a5, a6, a8, a9, a11, a12 = a8, a9, a11, a12, a2, a3, a5, a6
      aswap = 1
    b1 = a2*a9  - a3*a8
    #b2 = a2*a11 - a5*a8
    b3 = a2*a12 - a6*a8
    # Subtotal : 0I 8M 3A

    if (b1 != 0) :
      # Compute the reduced row echelon form of M
      #
      #   M_rref = [ 0  1  0  0  *  -s0 ]
      #            [ 0  0  1  0  *  -s1 ]
      gamma = 1/(a2*b1)
      alpha = gamma*b1
      beta = gamma*a2
      s1 = -beta*b3
      s0 = -alpha*(a6 + a3*s1)
      # Subtotal : 1I 6M 1A

      # Compute D1 + D2 = <u, v> where
      #
      #   u = f
      #   v = s0*g + s1*h + xh - h[1]*u
      v0 = g[0]*s0 + h[0]*s1 - h[1]*f[0]
      v1 = g[1]*s0 + h[1]*(s1 - f[1]) + h[0]
      v2 = g[2]*s0 + h[2]*s1 - h[1]*f[2]
      v4 = s0 + h[2]
      v5 = s1
      # Subtotal : 0I 8M 8A

      # D1 + D2 is of type 54
      # Total : 1I 29M 19A
      return C34CurveDivisor(C, [copy(D1.f), [v0, v1, v2, 0, v4, v5, 0, 0, 1]])
    
    else :
      # M_ref is the matrix
      #
      #   M_ref = [ 0  a2  a3  0  a5  a6 ]
      #           [ 0  0   0   0  0   0  ]
      #
      # Compute its reduced row echelon form
      #
      #   M_rref = [ 0  1  -r0  0  *  * ]
      #            [ 0  0   0   0  0  0 ]
      alpha = 1/a2
      r0 = -alpha*a3
      # Subtotal : 1I 1M 0A

      # LCM(D1, D2) is type 43, generated by <u, v> where
      #
      #   u = f
      #   v = r0*g + h
      v0 = g[0]*r0 + h[0]
      v1 = g[1]*r0 + h[1]
      v2 = g[2]*r0 + h[2]
      v4 = r0
      L = C34CurveDivisor(C, [copy(D1.f), [v0, v1, v2, 0, v4, 1], []])
      # Subtotal : 0I 3M 3A

      # GCD(D1, D2) is of type 11, generated by <p, q> where
      #
      #   p = F
      #   q = y + a2/a8
      #
      # assuming a2 and a8 have not been swapped
      if (aswap == 0) :
        mu = 1/a8
        q0 = mu*a2
      else :
        mu = 1/a2
        q0 = mu*a8
      G = C34CurveDivisor(C, [copy(D2.f), [q0, 0, 1], []])
      # Subtotal : 1I 1M 0A

      return flip(flip(L)) + G

  elif (a3 != 0) or (a9 != 0) :
    # M is the matrix
    #
    #   M = [ 0  0  a3  0  0  a6  ]
    #       [ 0  0  a9  0  0  a12 ]
    #
    # With reduced row echelon form
    #
    #   M_rref = [ 0  0  1  0  0  * ]
    #            [ 0  0  0  0  0  0 ]
    #
    # LCM(D1, D2) is type 42, generated by <f, g>
    # GCD(D1, D2) is type 11, generated by <p, q>, where
    # 
    #   p = F
    #   q = y + a3/a9
    #
    # assuming a3 and a9 have not been swapped
    if (aswap == 0) :
      mu = 1/a9
      q0 = a3*mu
    else :
      mu = 1/a3
      q0 = a9*mu
    L = C34CurveDivisor(C, [copy(D1.f), copy(D1.g), []])
    G = C34CurveDivisor(C, [copy(D2.f), [q0, 0, 1], []])
    # Subtotal : 1I 1M 0A

    return flip(flip(L)) + G

  else :
    # M = [ 0  0  0  0  0  0 ]
    #     [ 0  0  0  0  0  0 ]
    #
    # In this case, D1 = D2 + P for some point (i.e. type 11 divisor) P.
    # Compute P.
    P = C34CurveDivisor(C, [[f[1] - g[2], 1], [g[1], 0, 1], []])
    return double(D2) + P



def add_31_31(D1, D2) :
  C = D1.C
  K = C.K
  f0, f1, f2 = D1.f[0:3]
  g0, g1, g2 = D1.g[0:3]
  h0, h1, h2 = D1.h[0:3]
  F0, F1, F2 = D2.f[0:3]
  G0, G1, G2 = D2.g[0:3]
  H0, H1, H2 = D2.h[0:3]
  new_f, new_g, new_h = [], [], []

  # Compute M
  #
  #       [ a1   a2   a3   a4   a5   a6   a7  ]
  #   M = [ a8   a9   a10  a11  a12  a13  a14 ]
  #       [ a15  a16  a17  a18  a19  a20  a21 ]
  #
  # Columns 4, 5, and 6 are computed by
  #
  #   [ a4   a5   a6  ]   [ 0  -F0  -G0 ] [ a1   a2   a3  ] 
  #   [ a11  a12  a13 ] = [ 1  -F1  -G1 ]*[ a8   a9   a10 ]
  #   [ a18  a19  a20 ]   [ 0  -F2  -G2 ] [ a15  a16  a17 ]
  #
  # The last column similarly by
  #
  #   [ a7  ]   [ 0  -F0  -G0 ] [ a4  ] 
  #   [ a14 ] = [ 1  -F1  -G1 ]*[ a11 ]
  #   [ a21 ]   [ 0  -F2  -G2 ] [ a18 ]
  #
  # The last column is only computed when needed, which is rarely.
  
  a1  = f0 - F0
  a2  = g0 - G0
  a3  = h0 - H0
  a8  = f1 - F1
  a9  = g1 - G1
  a10 = h1 - H1
  a15 = f2 - F2
  a16 = g2 - G2
  a17 = h2 - H2
  
  a4  =    - F0*a8  - G0*a15
  a5  =    - F0*a9  - G0*a16
  a6  =    - F0*a10 - G0*a17
  a11 = a1 - F1*a8  - G1*a15
  a12 = a2 - F1*a9  - G1*a16
  a13 = a3 - F1*a10 - G1*a17
  a18 =    - F2*a8  - G2*a15
  a19 =    - F2*a9  - G2*a16
  a20 =    - F2*a10 - G2*a17
  # Subtotal : 0I 18M 21A

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
    # Since we did not compute a7, a14, a21, we do not compute b6, b12.
    
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
    # Subtotal : 0I 20M 10A
    
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
      # Subtotal : 0I 8M 4A
      
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
        # Subtotal : 1I 25M 9A

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
        # Subtotal : 0I 27M 27A

        # D1 + D2 is of type 61
        # Total : 1I 98M 71A
        return C34CurveDivisor(C, [[u0, u1, u2, u3, u4, u5, 1],
                             [v0, v1, v2, v3, v4, v5, 0, 1],
                             [w0, w1, w2, w3, w4, w5, 0, 0, 1]])
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
        # Subtotal : 1I 16M 4A
        
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
        # Subtotal : 0I 15M 15A

        # D1 + D2 is of type 63.
        # Total : 1I 77M 54A
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
        # Subtotal : 1I 9M 2A

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
        # Subtotal : 0I 12M 12A

        # D1 + D2 is of type 62.
        # Total : 1I 67M 49A
        return C34CurveDivisor(C, [[u0, u1, u2, u3, u4, 1],
                             [v0, v1, v2, v3, v4, 0, 1], []])

      else :
        assert c4 == 0
        # D1 and D2 are non-disjoint.
        # We compute D1 + D2 = lcm(D1, D2) + gcd(D1, D2)
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
        L = C34CurveDivisor(D1.C, [[u0, u1, u2, u3, u4, 1], [v0, v1, v2, v3, v4, 0, 1], [w0, w1, w2, w3, w4, 0, 0, 1]])

        # GCD(D1, D2) is type 11.
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
        return flip(flip(L)) + G

    elif (b2 != 0 or b8 != 0) :
      if (aswap == 0) :
        a7  =    - F0*a11 - G0*a18
        a14 = a4 - F1*a11 - G1*a18
        a21 =    - F2*a11 - G2*a18
      elif (aswap == 1) :
        a7  = a11 - F1*a4 - G1*a18
        a14 =     - F0*a4 - G0*a18
        a21 =     - F2*a4 - G2*a18
      else :
        a7  =     - F2*a11 - G2*a4
        a14 = a18 - F1*a11 - G1*a4
        a21 =     - F0*a11 - G0*a4
      b6  = a1*a14 - a7*a8
      b12 = a1*a21 - a7*a15

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
      # Subtotal : 0I 18M 10A, including computing a7, a14, a21, b6, b12.
      
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
        # Subtotal : 1I 14M 3A
        
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
        # Subtotal : 0I 12M 10A
        
        # D1 + D2 is of type 64.
        # Total : 1I 80M 58A
        return C34CurveDivisor(C, [[u0, u1, u2, u3, 1],
                             [v0, v1, v2, v3, 0, v5, v6, 0, 0, 1], []])
        
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
        # Subtotal : 1I 7M 1A
        
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
        L = C34CurveDivisor(D1.C, [[u0, u1, u2, u3, 1], [v0, v1, v2, v3, 0, v5, 1], []])
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
        
        L = flip(L) # flip_53
        L = flip(L) # flip_21
        ret = L + G # add_21_11
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

      L = flip(L) # flip_52
      L = flip(L) # flip_22
      ret = L + G # add_11_11
      return ret

    else :
      # M_ref is the matrix
      #
      #           [ a1  a2  a3  a4  a5  a6  a7 ]
      #   M_ref = [ 0   0   0   0   0   0   0  ]
      #           [ 0   0   0   0   0   0   0  ]
      #
      # Compute its reduced row echelon form
      #
      #            [ 1  -r0  -s0  -t0  *  * ]
      #   M_rref = [ 0   0    0    0   0  0 ]
      #            [ 0   0    0    0   0  0 ]
      alpha = 1/a1
      r0 = -alpha*a2
      s0 = -alpha*a3
      t0 = -alpha*a4
      # Subtotal : 1I 3M 0A

      # LCM(D1, D2) is of type 41, generated by <u, v, w>, where
      #
      #   u = r0*f + g
      #   v = s0*f + h
      #   w = t0*f + xf
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
      L = C34CurveDivisor(D1.C, [[u0, u1, u2, u3, 1],
                           [v0, v1, v2, v3, 0, 1], 
                           [w0, w1, w2, w3, 0, 0, 1]])
      # Subtotal : 0I 18M 18A

      # GCD(D1, D2) is degree 2, therefore either type 21 or 22. Assuming the
      # rows of M have not been swapped, if a15 != 0, then GCD(D1, D2) is of
      # type 21. Otherwise, it is of type 22. If it is of type 21, then it is
      # generated by <p, q> where
      #
      #   p = y + (a8/a15)*x + (a1/a15)
      #   q = f mod p
      #
      # If it is of type 22, then it is generated by <p, q> where
      #
      #   p = x + (a1/a8)
      #   q = h mod p
      if (aswap == 0) :
        m1, m2, m3 = a15, a8, a1
      elif (aswap == 1) :
        m1, m2, m3 = a15, a1, a8
      else :
        m1, m2, m3 = a1, a8, a15
      if (m1 != 0) :
        mu = 1/m1
        p1 = mu*m2
        p0 = mu*m3
        q0 = f0 - f2*p0
        q1 = f1 - f2*p1
        G = C34CurveDivisor(C, [[p0, p1, 1], [q0, q1, 0, 1], []])
        # Subtotal : 1I 4M 2A

        L = flip(L) # flip_41
        L = flip(L) # flip_21
        ret = L + G # add_21_21
        return ret
      else :
        assert m2 != 0
        mu = 1/m2
        p0 = mu*m3
        q0 = h0 - h1*p0
        q2 = h2
        G = C34CurveDivisor(C, [[p0, 1], [q0, 0, q2, 0, 0, 1], []])
        # Subtotal : 1I 2M 1A
        
        L = flip(L) # flip_41
        L = flip(L) # flip_21
        ret = L + G # add_21_22
        return ret

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
    # Subtotal : 0I 12M 6A
    
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
      # Subtotal : 0I 2M 1A
      
      if (c2 != 0) :
        # D1 + D2 is type 65, principal, generated by <f>
        new_f = [f0, f1, f2, 1]
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

        L = flip(L) # flip_54
        L = flip(L) # flip_11
        ret = L + G # add_22_11
        return ret

    else :
      # We have the matrix
      #
      #     [ 0  a2  a3  0  a5  a6  0 ]
      # M = [ 0  0   0   0  b2  b3  0 ]
      #     [ 0  0   0   0  b5  b6  0 ]
      #
      # LCM(D1, D2) will be type 43, so the b-values are all zero.
      assert b2 == b3 == b5 == b6 == 0
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
      return flip(flip(L)) + G

  else :
    # We have the matrix
    #
    #     [ 0  0  a3   0  0  a6   0 ]
    # M = [ 0  0  a10  0  0  a13  0 ]
    #     [ 0  0  a17  0  0  a20  0 ]
    #
    # LCM(D1, D2) will be of type 42, generated by <f, g>
    L = C34CurveDivisor(C, [copy(D1.f), copy(D1.g), []])

    # GCD(D1, D2) is degree 2, therefore either type 21 or 22. Assuming the
    # rows of M have not been swapped, if a17 != 0, then GCD(D1, D2) is of
    # type 21. Otherwise, it is of type 22. If it is of type 21, then it is
    # generated by <p, q> where
    #
    #   p = y + (a10/a17)*x + (a3/a17)
    #   q = f mod p
    #
    # If it is of type 22, then it is generated by <p, q> where
    #
    #   p = x + (a1/a8)
    #   q = h mod p
    if (aswap == 0) :
      m1, m2, m3 = a17, a10, a3
    elif (aswap == 1) :
      m1, m2, m3 = a17, a3, a10
    else :
      m1, m2, m3 = a3, a10, a17
    if (m1 != 0) :
      mu = 1/m1
      p1 = mu*m2
      p0 = mu*m3
      q0 = f0 - f2*p0
      q1 = f1 - f2*p1
      G = C34CurveDivisor(C, [[p0, p1, 1], [q0, q1, 0, 1], []])
      # Subtotal : 1I 4M 2A

      L = flip(L) # flip_42
      L = flip(L) # flip_22
      ret = G + L # add_21_11
      return ret
    else :
      assert m2 != 0
      mu = 1/m2
      p0 = mu*m3
      q0 = h0 - h1*p0
      q2 = h2
      G = C34CurveDivisor(C, [[p0, 1], [q0, 0, q2, 0, 0, 1], []])
      # Subtotal : 1I 2M 1A
      
      L = flip(L) # flip_42
      L = flip(L) # flip_22
      ret = G + L # add_22_11
      return ret
  
  return C34CurveDivisor(D1.C, [new_f, new_g, new_h])



