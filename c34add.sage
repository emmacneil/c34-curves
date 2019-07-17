"""
  Summary of costs for adding divisors
  
   D1    D2
  deg | deg |   I |   M |
  ----+-----+-----+-----+
  any |   0 |   0 |   0 |
  ----+-----+-----+-----+
    1 |   1 |   1 |   3 | Complete
  ----+-----+-----+-----+
   2a |   1 |   1 |  13 | Complete
   2a |  2a |   1 |  47 |
  ----+-----+-----+-----+
   2b |   1 |   1 |   5 | Complete
   2b |  2a |   1 |  16 |
   2b |  2b |   1 |   5 |
  ----+-----+-----+-----+
    3 |   1 |   1 |  18 |
    3 |  2a |   1 |  41 |
    3 |  2b |   1 |  27 |
    3 |   3 |   1 |  72 |
"""



def add(D1, D2) :
  """
    Add two divisors, D1 and D2.
    
    Will throw an exception if the support of D1 and D2 share a point in common.
    
    Input : Two distinct typical C34CrvDivs, D1 and D2.
    Output : The reduced C34CrvDiv D3 equivalent to D1 + D2. May be typical or semi-typical (or
             neither?)
  """
  C = D1.C
  if (D2.C != C) :
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
  # TODO : Ensure any function that computes a degree 3 divisor returns h, too.
  #        Then remove the calculation below.
  if (D1.degree == 3) and (len(D1.h) == 0) :
    print("Computing h polynomial for D1 = {}".format(D1))
    if D1.f[2] == 0 :
      raise ValueError("D1 is a degree 4 divisor.\nD1 = {}".format(D1))
    a = 1/D1.f[2]

    # This gives h of the form y^2 + ay + bx + c in 1I 7M
    D1.h = [a*(D1.g[1]*D1.f[0] - (D1.f[1] - D1.g[2])*D1.g[0]),
            a*(- D1.g[0] + D1.g[1]*D1.g[2]),
            a*(D1.f[0] - (D1.f[1] - D1.g[2])*D1.g[2]) + D1.g[1],
            K.zero(), K.zero(), K.one()]

  # If D2 = <f, g, h> is a reduced degree 3 divisor but missing its h polynomial, compute h.
  # TODO : Ensure any function that computes a degree 3 divisor returns h, too.
  #        Then remove the calculation below.
  if (D2.degree == 3) and (len(D2.h) == 0) :
    print("Computing h polynomial for D2 = {}".format(D2))
    if D2.f[2] == 0 :
      raise ValueError("D2 is a degree 4 divisor.\nD2 = {}".format(D2))
    a = 1/D2.f[2]

    # This gives h of the form y^2 + ay + bx + c in 1I 7M
    D2.h = [a*(D2.g[1]*D2.f[0] - (D2.f[1] - D2.g[2])*D2.g[0]),
            D2.a*(- D2.g[0] + D2.g[1]*D2.g[2]),
            a*(D2.f[0] - (D2.f[1] - D2.g[2])*D2.g[2]) + D2.g[1],
            K.zero(), K.zero(), K.one()]

  D3 = C.zero_divisor() 
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
    D3 = add_31_31(D1, D2)
  
  if D3.reduced :
    return D3
  return flip(flip(D3))



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
    return C34CrvDiv(C, [[u0, u1, 1], [v0, v1, 0, 1], []])

  else :
    r0 = G[0]
    v0 = g[0]*r0
    v2 = g[0] + r0
    
    # D1 + D2 is of type 22
    # Total : 0I 1M 3A
    return C34CrvDiv(C, [copy(D1.f), [v0, 0, v2, 0, 0, 1], []])



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
    return C34CrvDiv(C, [[u0, u1, u2, 1], [v0, v1, v2, 0, 1], [w0, w1, w2, 0, 0, 1]])
    
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
    return C34CrvDiv(C, [copy(D1.f), [v0, v1, 0, v3, 0, 0, 1], []])
  
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
      Q = C34CrvDiv(C, [[u0, 1], [v0, 0, 1], []])
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
      return C34CrvDiv(C, [[u0, u1, u2, u3, 1], [v0, v1, v2, v3, 0, 1], [w0, w1, w2, w3, 0, 0, 1]])
    
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
      return C34CrvDiv(C, [[u0, u1, u2, 1], [v0, v1, v2, 0, v4, 1], []])

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
      return C34CrvDiv(C, [[u0, u1, u2, 1], [v0, v1, v2, 0, 1], []])

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
      L = C34CrvDiv(C, [[u0, u1, u2, 1], [v0, v1, v2, 0, 1], [w0, w1, w2, 0, 0, 1]])

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
      G = C34CrvDiv(C, [[p0, 1], [q0, 0, 1], []])

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
      return C34CrvDiv(C, [copy(D1.f), [], []])
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

      L = C34CrvDiv(C, [copy(D1.f), [v0, v1, 0, v3, 0, 0, 1], []])

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
      G = C34CrvDiv(C, [[p0, 1], [q0, 0, 1], []])

      return L + G



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
    return C34CrvDiv(C, [[u0, u1, 0, 1], [v0, v1, v2, 0, 1], [w0, w1, w2, 0, 0, 1]])

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
    return C34CrvDiv(C, [copy(D1.f), [], []])
  
  else :
    # Divisors are non-disjoint
    # D1 = D2 + P for some point P.
    P = C34CrvDiv(C, [copy(D1.f), [g[2] - G[0], 0, 1], []])
    
    # If D2 = P, then D1 + D2 = 3P
    # Otherwise, D1 + D2 is 2*D2 + P
    if (D2 == P) :
      return triple(P)
    else :
      return double(D2) + P



def add_22_21(D1, D2) :
  C = D1.C
  f, g = D1.f, D1.g
  F, G = D2.f, D2.g
  
  # Construct the matrix M
  #
  #   M = [ a1  a2  a3  a4  a5 ]
  #       [ 1   a6  a7  a8  a9 ]
  # 
  # whose columns are the reductions of f, xf, yf, g, x^2*f modulo F, G.
  # Columns 2, 3, and 5 may be computed by
  #
  #   [ a2 ] = [ 0  -G[0] ]*[ a1 ]     [ a5 ] = [ 0  -G[0] ]*[ a2 ]
  #   [ a6 ]   [ 1  -G[1] ] [ a5 ]     [ a9 ]   [ 1  -G[1] ]*[ a6 ]
  #
  #   [ a3 ] = [ -F[0]  F[1]*G[0]        ]*[ a1 ]
  #   [ a7 ]   [ -F[1]  F[1]*G[1] - F[0] ] [ a5 ]
  
  a1 = f[0]
  a4 = -F[1]^2*G[0] + F[0]*(F[0] - g[2]) + g[0]
  a8 = -F[1]*(F[1]*G[1] - 2*F[0] + g[2])

  a2 =    - G[0]*a5
  a6 = a1 - G[1]*a5
  a5 =    - G[0]*a6
  a9 = a2 - G[1]*a6

  a3 = -F[0]*a1 + F[1]*G[0]*a5
  a7 = -F[1]*a1 + (F[1]*G[1] - F[0])*a5
  # Subtotal : 0I 15M 1SQ 10A
  # XXX : This is slower already that add_21_22



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
    return C34CrvDiv(C, [[u0, u1, u2, u3, 1], [v0, v1, v2, v3, 0, 1], [w0, w1, 0, w3, 0, 0, 1]])

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
    return C34CrvDiv(C, [copy(D1.g), [v0, v1, v2, 0, 1], []])

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
    L = C34CrvDiv(C, [copy(D1.g), [v0, v1, v2, 0, 1], [w0, w1, w2, 0, 0, 1]])

    # GCD(D1, D2) is type 11, generated by
    #
    #   p = F
    #   q = y + a1
    G = C34CrvDiv(C, [copy(D2.f), [a1, 0, 1], []])
    return L + G



def old_add_22_21(D1, D2):
  K = D1.K
  f, g = D1.f, D1.g
  F, G = D2.f, D2.g
  new_f, new_g, new_h = [], [], []

  # M = [ a1  a2  a3  a4  a5 ]
  #     [  1   0  a6  a7   0 ]
  a1 = F[0] - F[1]*f[0]
  a2 = f[0]*(f[0] - G[1]) + G[0]
  a3 = -f[0]*a1
  a4 = -g[0]
  a5 = -f[0]*a2
  a6 = -f[0]
  a7 = a1 - g[2]

  # M_ref = [ 1   0  a6  a7  0 ]
  #         [ 0  b1   0  b3  b4 ]
  b1 = a2
  b3 = a4 - a7*a1
  b4 = a5

  
  if (b1 != 0) :
    # M_rref = [ 1  0  -u0  -v0    0 ]
    #          [ 0  1    0  -v1  -w1 ]
    beta = 1/b1
    v1 = -beta*b3
    w1 = -beta*b4
    u0 = -a6
    v0 = -a7

    new_f = [F[0]*u0,
             F[1]*u0 + F[0],
             u0,
             F[1],
             K.one()]
    new_g = [F[0]*v0 + G[0]*v1 - F[1]*new_f[0],
             F[1]*v0 + G[1]*v1 - F[1]*new_f[1],
             v0 + F[0] - F[1]*new_f[2],
             v1 - F[1]*new_f[3],
             K.zero(),
             K.one()]
    new_h = [G[0]*w1,
             G[1]*w1 + G[0],
             K.zero(),
             w1 + G[1],
             K.zero(), K.zero(), K.one()]
  elif (b3 != 0) :
    beta = 1/b3
    v1 = -beta*b4
    u0 = -a6
    v0 = a7*v1

    #z = F[1]*v1
    new_f = [ G[0], G[1], K.zero(), K.one() ]
    new_g = [ F[0]*u0        - F[1]*new_f[0],
              F[1]*u0 + F[0] - F[1]*new_f[1],
              u0,
              K.zero(),
              K.one() ]
    #new_h = [ F[0]*v0 - G[1]*new_f[0] - z*new_g[0],
    #          F[1]*v0 + G[0] - G[1]*new_f[1] - z*new_g[1],
    #          v0 + F[0]*v1 - z*new_g[2],
    #          K.zero(),
    #          K.zero(),
    #          v1,
    #          K.one() ]
  else :
    # Divisors are non-disjoint.
    # So we have D1 = P1 + P2 and D2 = P1 + P3.
    # (Possibly P1 = P2 or P1 = P3, but not P2 = P3.)
    x1 = -f[0]
    x2 = f[0] - G[1]
    y1 = -a1
    y2 = a1 - g[2]
    y3 = -F[1]*(x2) - F[0]
    P1 = C34CrvDiv(D1.C, [[-x1, K.one()], [-y1, K.zero(), K.one()], []])
    P2 = C34CrvDiv(D1.C, [[-x1, K.one()], [-y2, K.zero(), K.one()], []])
    P3 = C34CrvDiv(D1.C, [[-x2, K.one()], [-y3, K.zero(), K.one()], []])

    # D1 = P1 + P2
    # D2 = P1 + P3
    # If P1 = P2 or if P1 = P3, we need to triple and add (3*P1 + P3), resp. (3*P1 + P2)
    # Otherwise, we double and add (2*P1) + (P2 + P3)
    if (P1 == P2) :
      return add(triple(P1), P3)
    elif (P1 == P3) :
      return add(triple(P1), P2)
    else :
      return add(double(P1), add(P2, P3))

    raise ValueError("Divisors are non-disjoint.")

  return C34CrvDiv(D1.C, [new_f, new_g, new_h])



def add_22_22(D1, D2):
  K = D1.K
  f, g = D1.f, D1.g
  F, G = D2.f, D2.g
  new_f, new_g, new_h = [], [], []

  a1 = F[0] - f[0]
  a2 = -f[0]*a1
  a3 = G[0] - g[0]
  a4 = -f[0]*a2
  a5 = G[2] - g[2]

  if a1 != 0 :
    alpha = 1/a1
    u0 = -alpha*a3
    u1 = -alpha*a5
    new_f = [ F[0]*f[0], F[0] + f[0], K.zero(), K.one() ]
    new_g = [ F[0]*u0 + G[0], u0, F[0]*u1 + G[2], K.zero(), u1, K.one() ]
  else :
    # Points are colinear with same x-coordinate and necessarily non-disjoint.
    y0 = -a3*(1/a5) # y-coordinate of common point
    P1 = C34CrvDiv(D1.C, [[f[0], K.one()], [-y0, K.zero(), K.one()], []])
    P2 = C34CrvDiv(D1.C, [[f[0], K.one()], [y0 + g[2], K.zero(), K.one()], []])
    P3 = C34CrvDiv(D1.C, [[f[0], K.one()], [y0 + G[2], K.zero(), K.one()], []])

    # D1 = P1 + P2
    # D2 = P1 + P3
    # If P1 = P2 or if P1 = P3, we need to triple and add (3*P1 + P3), resp. (3*P1 + P2)
    # Otherwise, we double and add (2*P1) + (P2 + P3)
    if (P1 == P2) :
      return add(triple(P1), P3)
    elif (P1 == P3) :
      return add(triple(P1), P2)
    else :
      return add(double(P1), add(P2, P3))

  # What is this?????
  # alpha = 1/alpha
  # a = -alpha*(G[2] - g[2])
  # b = -alpha*(G[0] - g[0])
  # new_f = [ f[0]*F[0], f[0] + F[0], K.zero(), K.one() ]
  # new_g = [ b*F[0] + G[0], b, a*F[0] + G[2], K.zero(), a, K.one() ]

  return C34CrvDiv(D1.C, [new_f, new_g, new_h])



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

  if (a1 != 0) :
    # Compute the reduced row echelon form of M,
    #
    # M_rref = [ 1  -r0  -s0  -t0 ]
    alpha = 1/a1
    r0 = -alpha*a2
    s0 = -alpha*a3
    t0 = F[0]

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
    
    # D1 + D2 is of type 41
    return C34CrvDiv(C, [[u0, u1, u2, u3, 1], [v0, v1, v2, v3, 0, 1], [w0, w1, w2, w3, 0, 0, 1]])

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
    
    # D1 + D2 is of type 43
    return C34CrvDiv(C, [copy(D1.f), [u0, u1, u2, 0, u4, 1], []])

  elif (a3 != 0) :
    # M = [ 0  0  1  0 ]

    # D1 + D2 is of type 42
    return C34CrvDiv(C, [copy(D1.f), copy(D1.g), []])
  
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
      A = C34CrvDiv(C, [[p0, p1, 1], [q0, q1, 0, 1], []])
    else :
      # TODO : This section is painfully inefficient
      A = D1.slow_add(D2.slow_flip())
      # D1 is not typical
      if (4*f[0] != f[1]^2) :
        # f has two distinct rational roots
        # f = (x + F[0])*(x + f[1] - F[0])
        if (F[0] == g[2]) :
          # A is type 21
          assert A.type == 21, "A is not type 21."
          Q = C34CrvDiv(C, [[g[2], 1], [h[2] - G[0], 0, 1], []])
          R = C34CrvDiv(C, [[f[1] - g[2], 1], [g[1], 0, 1], []])
          A = Q + R
        else :
          # A is type 22
          assert A.type == 22, "A is not type 22."
          p0 = f[1] - F[0]
          q0 = h[0] + h[1]*(F[0] - f[1])
          q2 = h[2]
          A = C34CrvDiv(C, [[p0, 1], [q0, 0, q2, 0, 0, 1], []])
      else :
        # f has a rational double root
        # f = (x + F[0])^2
        if (G[0] == g[1]) :
          # A is type 22
          assert A.type == 22, "A is not type 22."
          p0 = f[1] - F[0]
          q0 = h[0] + h[1]*(F[0] - f[1])
          q2 = h[2]
          A = C34CrvDiv(C, [[p0, 1], [q0, 0, q2, 0, 0, 1], []])
        else :
          # A is type 21
          assert A.type == 21, "A is not type 21."
          Q = C34CrvDiv(C, [[g[2], 1], [h[2] - G[0], 0, 1], []])
          R = C34CrvDiv(C, [[f[1] - F[0], 1], [g[1], 0, 1], []])
          A = Q + R
      
    assert A.slow_add(D2) == D1, "A + D2 =/= D1"

    # Check if P is in the support of A.
    # If not, return D1 + D2 = A + 2P
    # Otherwise, find Q such that A = P + Q
    if (A.type == 21) and (- G[0] - A.f[1]*F[0] + A.f[0] == 0) and (F[0]*(F[0] - A.g[1]) + A.g[0] == 0) :
      r0 = A.g[1] - F[0]
      s0 = A.f[1]*(F[0] - A.g[1]) + A.f[0]
      Q = C34CrvDiv(C, [[r0, 1], [s0, 0, 1], []])
      # If P =/= Q, then return D1 + D2 = Q + 3P. Otherwise, return 4P
      assert Q.slow_add(D2) == A, "Q + D2 =/= A"
      if (Q != D2) :
        return Q + flip(flip(triple(D2)))
      else :
        return double(double(D2))
    elif (A.type == 22) and (p0 == F[0]) and (G[0]*(G[0] - A.g[2]) + A.g[0] == 0) :
      r0 = A.f[0]
      s0 = A.g[2] - G[0]
      Q = C34CrvDiv(C, [[r0, 1], [s0, 0, 1], []])
      # If P =/= Q, then return D1 + D2 = Q + 3P. Otherwise, return 4P
      assert Q.slow_add(D2) == A, "Q + D2 =/= A"
      if (Q != D2) :
        return Q + flip(flip(triple(D2)))
      else :
        return double(double(D2))
    else :
      return A + double(D2)



def old_add_31_11(D1, D2):
  K = D1.K
  f, g, h = D1.f, D1.g, D1.h
  F, G, H = D2.f, D2.g, D2.h
  new_f, new_g, new_h = [], [], []
  # Compute matrix M
  # M = [ a1  a2  a3  a4 ]
  # (a4 is not needed)
  a1 = -F[0]*(f[1] - F[0]) - G[0]*f[2] + f[0]
  a2 = -F[0]*(g[1] - G[0]) - G[0]*g[2] + g[0]
  a3 = -G[0]*(h[2] - G[0]) - F[0]*h[1] + h[0]
  #a4 = -F[0]*a1

  if a1 != 0 :
    alpha = 1/a1
    u0 = -alpha*a2
    v0 = -alpha*a3
    w0 = F[0]
    new_f = [ f[0]*u0 + g[0],
              f[1]*u0 + g[1],
              f[2]*u0 + g[2],
              u0, K.one() ]
    new_g = [ f[0]*v0 + h[0],
              f[1]*v0 + h[1],
              f[2]*v0 + h[2],
              v0, K.zero(), K.one() ]
    new_h = [ f[0]*w0        - f[2]*new_f[0],
              f[1]*w0 + f[0] - f[2]*new_f[1],
              f[2]*w0        - f[2]*new_f[2],
                   w0 + f[1] - f[2]*new_f[3],
              K.zero(), K.zero(), K.one() ]
    # Total : 1I 21M
  elif a2 != 0 :
    u0 = -a3*(1/a2)
    new_f = [ f[0], f[1], f[2], K.one() ]
    new_g = [ g[0]*u0 + h[0],
              g[1]*u0 + h[1],
              g[2]*u0 + h[2],
              K.zero(), u0, K.one() ]
    # Total : 1I 10M
  elif a3 != 0 :
    new_f = [ f[0], f[1], f[2], K.one() ]
    new_g = [ g[0], g[1], g[2], K.zero(), K.one() ]
    # Total : 0I 6M
  else :
    # Divisors are non-disjoint.
    # We have D1 = P + Q + R and D2 = P, for some (not necessarily distinct) P, Q, R
    # Extract the largest multiple of P we can from D1.
    # I.e., let n be the order of P in (D1 + D2) and find a divisor E such that D1 + D2 = E + nP.
    if f[2] != 0 :
      # D1 is typical.
      
      # Compute polynomials s = y + s1*x + s0 and t = x^2 + t1*x + t0 such that Q + R = <s, t>
      alpha = 1/f[2]
      s1 = alpha*(F[0] - g[2])
      s0 = g[1] + s1*(f[1] - F[0])
      t1 = f[1] - f[2]*s1
      t0 = f[0] - f[2]*s0

      # If s(P) and t(P) are non-zero,      
      if (-G[0] - s1*F[0] + s0 != 0) or (F[0]*F[0] - t1*F[0] + t0 != 0) :
        # E = <s, t> = Q + R is disjoint from D2.
        # Compute (2P) + (Q + R) = E + 2*D2
        E = C34CrvDiv(D1.C, [[s0, s1, K.one()], [t0, t1, K.zero(), K.one()], []])
        return add(E, double(D2))
      
      # Otherwise P = R (w.l.o.g.)
      # Compute u = x + u0 and v = y + v0 such that Q = <u, v>
      u0 = t1 - F[0]
      v0 = s0 - s1*F[0]
      
      # If u(P) and t(P) are non-zero
      if (u0 != F[0]) or (v0 != G[0]) :
        # E = <u, v> = Q is disjoint from D2.
        # Compute (3P) + Q = 3*D2 + E
        E = C34CrvDiv(D1.C, [[u0, K.one()], [v0, K.zero(), K.one()], []])
        return add(triple(D2), E)
      
      else :
        # D1 = 3P and D2 = P
        # Compute D1 + D2 = 2*(2*D2)
        return double(double(D2))
        
    else :
      # D1 is non-typical
      # D1 is of the form P + Q + R, where Q and R have the same x-coordinate and P has a different x-coordinate.
      # Possibly, Q = R.
      # There are two cases, either D2 = P or D2 =/= P
      # We have D2 = P if and only if F[0] =/= g[2]
      
      # If D2 = P,
      if F[0] != g[2] :
        # Find the type 22 divisor E = Q + R
        # In this case, E = <F, h mod F>
        E = C34CrvDiv(D1.C, [[F[0], K.one()], [h[0] - h[1]*F[0], K.zero(), h[2], K.zero(), K.zero(), K.one()], []])
        return add(E, double(D2))
      # Otherwise, if D2 =/= P
      else :
        # Find the type 21 divisor E = P + R
        EP = C34CrvDiv(D1.C, [[f[1] - F[0], K.one()], [g[1], K.zero(), K.one()], []])
        ER = C34CrvDiv(D1.C, [[F[0], K.one()], [h[2] - G[0], K.zero(), K.one()], []])
        E = add(EP, ER)
        return add(E, double(D2))

  return C34CrvDiv(D1.C, [new_f, new_g, new_h])



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
      
      # D1 + D2 is of type 51
      return C34CrvDiv(C, [[u0, u1, u2, u3, u4, 1], [v0, v1, v2, v3, v4, 0, 1], [w0, w1, w2, w3, w4, 0, 0, 1]])
      
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
      u0 = f[0]*r0 + g[0]
      u1 = f[1]*r0 + g[1]
      u2 = f[2]*r0 + g[2]
      u3 = r0
      v0 = f[0]*s0 + h[0]*s1 - f[2]*u0
      v1 = f[1]*s0 + h[1]*s1 + f[0] - f[2]*u1
      v2 = f[2]*(s0 - u2) + h[2]*s1
      v3 = s0 + f[1] - f[2]*u3
      v5 = s1
      # D1 + D2 is of type 53
      return C34CrvDiv(C, [[u0, u1, u2, u3, 1], [v0, v1, v2, v3, 0, v5, 1], []])
      
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
      u0 = f[0]*r0 + g[0]
      u1 = f[1]*r0 + g[1]
      u2 = f[2]*r0 + g[2]
      u3 = r0
      v0 = f[0]*s0 + h[0]
      v1 = f[1]*s0 + h[1]
      v2 = f[2]*s0 + h[2]
      v3 = s0
      # D1 + D2 is of type 52
      return C34CrvDiv(C, [[u0, u1, u2, u3, 1], [v0, v1, v2, v3, 0, 1], []])

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
      L = C34CrvDiv(C, [[u0, u1, u2, u3, 1], [v0, v1, v2, v3, 0, 1], [w0, w1, w2, w3, 0, 0, 1]])

      if (aswap == 0) :
        mu = 1/a7
        p0 = mu*a1
      else :
        mu = 1/a1
        p0 = mu*a7
      q0 = F[0] - F[1]*p0
      G = C34CrvDiv(C, [[p0, 1], [q0, 0, 1], []])
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
      u0 = g[0]*r0 + h[0]*r1 - h[1]*f[0]
      u1 = g[1]*r0 + h[1]*(r1 - f[1]) + h[0]
      u2 = g[2]*r0 + h[2]*r1 - h[1]*f[2]
      u4 = r0 + h[2]
      u5 = r1

      # D1 + D2 is of type 54
      return C34CrvDiv(C, [copy(D1.f), [u0, u1, u2, 0, u4, u5, 0, 0, 1], []])
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
      L = C34CrvDiv(C, [copy(D1.f), [u0, u1, u2, 0, u4, 1], []])

      if (aswap == 0) :
        mu = 1/a8
        p0 = mu*a2
      else :
        mu = 1/a2
        p0 = mu*a8
      q0 = F[0] - F[1]*p0
      G = C34CrvDiv(C, [[p0, 1], [q0, 0, 1], []])
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
    L = C34CrvDiv(C, [copy(D1.f), copy(D1.g), []])

    mu = 1/a9
    p0 = mu*a3
    q0 = F[0] - F[1]*p0
    G = C34CrvDiv(C, [[p0, 1], [q0, 0, 1], []])
    return flip(flip(L)) + G

  else :
    # M = [ 0  0  0  0  0  0 ]
    #     [ 0  0  0  0  0  0 ]
    #
    # In this case, D1 = D2 + P for some point (i.e. type 11 divisor) P.
    # We compute P by solving for polynomials p and q such that
    #
    #   (f, g, h) = (F, G)*(p, q)
    P = C34CrvDiv(D1.C, [ [g[2] + f[2]*F[1], K.one()], [h[2] + g[2]*F[1] - F[0], K.zero(), K.one()], [] ])
    return double(D2) + P



def old_add_31_21(D1, D2):
  K = D1.K
  f, g, h = D1.f, D1.g, D1.h
  F, G    = D2.f, D2.g
  new_f, new_g, new_h = [], [], []

  disjoint = True
  n_pts_in_common = 0
  # Compute matrix M
  # M = [ a1  a2  a3  a4   a5   a6  ]
  #     [ a7  a8  a9  a10  a11  a12 ]
  F1G0 = F[1]*G[0]
  a1 = f[0] - G[0] - f[2]*F[0]
  a7 = f[1] - G[1] - f[2]*F[1]
  a2 = g[0] + F1G0 - g[2]*F[0]
  a8 = g[1] - F[0] + F[1]*(G[1] - g[2])
  a3 = h[0] - F[1]*F1G0 + F[0]*(F[0] - h[2])
  a9 = h[1] + F[1]*(2*F[0] - F[1]*G[1] - h[2])
  
  a4  =  -   G[0]*a7
  a10 = a1 - G[1]*a7
  a5  =  -   G[0]*a8
  a11 = a2 - G[1]*a8
  a6  =    - G[0]*a9
  a12 = a3 - G[1]*a9
  # Subtotal : 0I 15M 1C

  if (a1 != 0) or (a7 != 0) :
    if a1 == 0 :
      a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12 = a7, a8, a9, a10, a11, a12, a1, a2, a3, a4, a5, a6
    # Compute row echelon form Mref
    # Mref = [ a1  a2  a3  a4  a5  a6 ]
    #        [  0  b1  b2  b3  b4  b5 ]
    b1 = a1*a8  - a2*a7
    b2 = a1*a9  - a3*a7
    b3 = a1*a10 - a4*a7
    b4 = a1*a11 - a5*a7
    b5 = a1*a12 - a6*a7
    # Subtotal : 0I 25M IC

    if (b1 != 0) :
      gamma = 1/(a1*b1)
      alpha = gamma*b1
      beta = gamma*a1
      u1 = -beta*b2
      v1 = -beta*b3
      w1 = -beta*b4
      u0 = -alpha*(a3 + a2*u1)
      v0 = -alpha*(a4 + a2*v1)
      w0 = -alpha*(a5 + a2*w1)
      
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
      new_h = [ f[0]*w0 + g[0]*w1,
                f[1]*w0 + g[1]*w1 + g[0],
                f[2]*w0 + g[2]*w1,
                w0 + g[1],
                w1 + g[2],
                K.zero(), K.zero(), K.one() ]
      # Total : 1I 55M 1C
      # Type <yy, xxx, xxy>
      # (13M are needed for computing new_h)
      # (Another 4M not needed at all (a6, a12, b5)).
    elif (b2 != 0) :
      gamma = 1/(a1*b2)
      alpha = gamma*b2
      beta = gamma*a1
      v1 = -beta*b3
      u0 = -alpha*a2
      v0 = -alpha*(a4 + a3*v1)
      
      new_f = [ f[0]*u0 + g[0],
                f[1]*u0 + g[1],
                f[2]*u0 + g[2],
                u0,
                K.one() ]
      new_g = [ f[0]*v0 + h[0]*v1 - f[2]*new_f[0],
                f[1]*v0 + h[1]*v1 + f[0] - f[2]*new_f[1],
                h[2]*v1 + f[2]*(v0 - new_f[2]),
                v0 + f[1] - f[2]*new_f[3],
                K.zero(),
                v1,
                K.one() ]
      # Total : 1I 44M IC
      # Type <xy, xxx>
    elif (b3 != 0) :
      alpha = 1/a1
      u0 = -alpha*a2
      v0 = -alpha*a3

      new_f = [ f[0]*u0 + g[0],
                f[1]*u0 + g[1],
                f[2]*u0 + g[2],
                u0,
                K.one() ]
      new_g = [ f[0]*v0 + h[0],
                f[1]*v0 + h[1],
                f[2]*v0 + h[2],
                v0,
                K.zero(), K.one() ]
      # Total : 1I 33M IC
      # Type <xy, yy>
    else :
      # D1 and D2 have a single point in common.
      disjoint = False
      n_pts_in_common = 1
  elif (a2 != 0) or (a8 != 0) :
    if (a2 == 0) :
      a2, a3, a5, a6, a8, a9, a11, a12 = a8, a9, a11, a12, a2, a3, a5, a6
    # Compute row echelon form Mref
    # Mref = [ 0  a2  a3  0  a5  a6 ]
    #        [ 0   0  b1  0  b2  b3 ]
    b1 = a2*a9  - a3*a8
    b2 = a2*a11 - a5*a8
    b3 = a2*a12 - a6*a8
    # Subtotal : 0I 21M IC
    
    if (b1 != 0) :
      gamma = 1/(a2*b1)
      alpha = gamma*b1
      beta = gamma*a2
      u1 = -beta*b2
      v1 = -beta*b3
      u0 = -alpha*(a5 + a3*u1)
      u1 = -alpha*(a6 + a3*v1)
      new_f = [ f[0], f[1], f[2], K.one() ]
      new_g = [ g[0]*v0 + h[0]*v1 - h[1]*f[0],
                g[1]*v0 + h[1]*(v1 - f[1]) + h[0],
                g[2]*v0 + h[2]*v1 - h[1]*f[2],
                K.zero(),
                v0 + h[2],
                v1,
                K.zero(), K.zero(), K.one() ]
      # Total : 1I 38M 1C
      # Type <xx, xyy>
    elif (b2 != 0) :
      # XXX : Is this case possible?!
      raise NotImplementedError("Could not compute kernel.")
    else :
      # D1 and D2 have a single point in common.
      disjoint = False
      n_pts_in_common = 1
  else :
    # D1 and D2 have one or two points in common
    # One point if a3 or a9 is non-zero, otherwise two.
    disjoint = False
    if (a3 != 0) or (a9 != 0) :
      n_pts_in_common = 1
    else :
      n_pts_in_common = 2
  
  if (not disjoint) :
    if (n_pts_in_common == 1) :
      # D1 = P + Q + R
      # D2 = P + S
      if (f[2] != 0) :
        r = f[1] - G[1] - f[2]*F[1]
        if (r != 0) :
          px = -(f[0] - G[0] - f[2]*F[0])/r
        else :
          px = -(g[0] - g[2]*F[0] + F[1]*G[0])/(g[1] + F[1]*(G[1] - g[2]) - F[0])
      else :
        if (f[1] != G[1]) :
          px = - (f[0] - G[0])/(f[1] - G[1])
        else :
          px = - (F[0]*(F[0] - h[2]) + h[0] - F[1]^2*G[0])/(F[1]*(2*F[0] - h[2] - F[1]*G[1]) + h[1])
      py = -F[1]*px - F[0]
      P = C34CrvDiv(D1.C, [[-px, K.one()], [-py, K.zero(), K.one()], []])
      
      sx = -G[1] - px
      sy = -F[1]*sx - F[0]
      S = C34CrvDiv(D1.C, [[-sx, K.one()], [-sy, K.zero(), K.one()], []])

      if f[2] != 0 :
        s1 = (-g[2] - px)/f[2]
        s0 = g[1] + s1*(px + f[1])
        t1 = f[1] - s1*f[2]
        t0 = f[0] - s0*f[2]
        QR = C34CrvDiv(D1.C, [[s0, s1, K.one()], [t0, t1, K.zero(), K.one()], []])
      elif g[1] != py :
        s1 = -h[1]/(g[1] + py)
        s0 = h[2] + py + s1*g[2]
        t1 = f[1]
        t0 = f[0]
        QR = C34CrvDiv(D1.C, [[s0, s1, K.one()], [t0, t1, K.zero(), K.one()], []])
      else :
        s0 = g[2]
        t2 = h[2]
        t0 = h[0] - h[1]*g[2]
        QR = C34CrvDiv(D1.C, [[s0, K.one()], [t0, K.zero(), t2, K.zero(), K.zero(), K.one()], []])
      
      # Check if S is P      
      if (px == sx) and (py == sy) :
        # P = S and D1 + D2 = 3P + Q + R
        return add(triple(P), QR)
      
      # Check if Q or R is P
      if len(QR.f) == 3 :
        if (py + QR.f[1]*px + QR.f[0] == 0) and (px*(px + QR.g[1]) + QR.g[0]) :
          if (2*px == -QR.g[1]) :
            # P = Q = R
            # Compute (2(2P)) + (S)
            return add(double(double(P)), S)
          else :
            # P = Q != R
            # Compute (3P) + (R + S)
            rx = -QR.g[1] - px
            ry = -QR.f[1]*rx - QR.f[0]
            R = C34CrvDiv(D1.C, [[-rx, K.one()], [-ry, K.zero(), K.one()], []])
            return add(triple(P), add(R, S))
      else :
        if (px == -QR.f[0]) and (py*(py + QR.g[2]) + QR.g[0]) :
          if (2*py == -QR.g[2]) :
            # P = Q = R
            # Compute (2(2P)) + (S)
            return add(double(double(P)), S)
          else :
            # P = Q != R
            # Compute (3P) + (R + S)
            rx = px
            ry = -py - QR.g[2]
            R = C34CrvDiv(D1.C, [[-rx, K.one()], [-ry, K.zero(), K.one()], []])
            return add(triple(P), add(R, S))
      
      # If we reach this line, then Q, R, S are all distinct from P
      # Compute (2P) + (Q + R + S)
      return add(double(P), add(QR, S))
      
    else :
      # D1 and D2 have two points in common
      # D1 + D2 = (2*D2) + P for some point P
      P = C34CrvDiv(D1.C, [ [g[2] + f[2]*F[1], K.one()], [h[2] + g[2]*F[1] - F[0], K.zero(), K.one()], [] ])
      return add(double(D2), P)
  
  return C34CrvDiv(D1.C, [new_f, new_g, new_h])



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
      return C34CrvDiv(C, [[u0, u1, u2, u3, u4, 1],
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
      return C34CrvDiv(C, [[u0, u1, u2, u3, 1], 
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
      return C34CrvDiv(C, [[u0, u1, u2, u3, 1],
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
      L = C34CrvDiv(C, [[u0, u1, u2, u3, 1],
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
      G = C34CrvDiv(C, [copy(D2.f), [q0, 0, 1], []])
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
      return C34CrvDiv(C, [copy(D1.f), [v0, v1, v2, 0, v4, v5, 0, 0, 1]])
    
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
      L = C34CrvDiv(C, [copy(D1.f), [v0, v1, v2, 0, v4, 1], []])
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
      G = C34CrvDiv(C, [copy(D2.f), [q0, 0, 1], []])
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
    L = C34CrvDiv(C, [copy(D1.f), copy(D1.g), []])
    G = C34CrvDiv(C, [copy(D2.f), [q0, 0, 1], []])
    # Subtotal : 1I 1M 0A

    return flip(flip(L)) + G

  else :
    # M = [ 0  0  0  0  0  0 ]
    #     [ 0  0  0  0  0  0 ]
    #
    # In this case, D1 = D2 + P for some point (i.e. type 11 divisor) P.
    # Compute P.
    P = C34CrvDiv(C, [[f[1] - g[2], 1], [g[1], 0, 1], []])
    return double(D2) + P



def old_add_31_22(D1, D2):
  K = D1.K
  f, g, h = D1.f, D1.g, D1.h
  F, G, H = D2.f, D2.g, D2.h
  new_f, new_g, new_h = [], [], []

  # Construct the matrix
  #
  # [ a1  a2  a3  a4   a5   a6  ]
  # [ a7  a8  a9  a10  a11  a12 ]
  #
  # The first three columns are the reductions of f, g, h modulo F, G.
  # The last three columns are the reductions of xf, xg, xh and computed simply by
  #
  # [ a4   a5   a6  ] = [ -F[0]     0  ]*[ a1  a2  a3 ]
  # [ a10  a11  a12 ]   [    0   -F[0] ] [ a7  a8  a9 ]
  #
  # The last column is only computed when necessary, i.e. when D1 + D2 is of type 54.
  
  a1 = F[0]*(F[0] - f[1]) + f[0]
  a2 = g[0] - g[1]*F[0]
  a3 = h[0] - G[0] - h[1]*F[0]
  a7 = f[2]
  a8 = g[2] - F[0]
  a9 = h[2] - G[2]
  
  a4  = - F[0]*a1
  a5  = - F[0]*a2
  a10 = - F[0]*a7
  a11 = - F[0]*a8

  if (a1 != 0) or (a7 != 0) :
    # Before attempting to row reduce, swap rows and relabel if necessary to ensure a1 != 0.
    if (a1 == 0) :
      a1, a2, a3, a4, a5, a7, a8, a9, a10, a11 = a7, a8, a9, a10, a11, a1, a2, a3, a4, a5

    # Reduce to the matrix
    #
    # [ a1  a2  a3  a4  a5  a6 ]
    # [ 0   b1  b2  b3  b4  b5 ]
    b1 = a1*a8  - a2*a7
    b2 = a1*a9  - a3*a7
    b3 = a1*a10 - a4*a7
    b4 = a1*a11 - a5*a7
    
    if (b1 != 0) :
      # Reduce to RREF
      #
      # [ 1  0  -r0  -s0  -t0  * ]
      # [ 0  1  -r1  -s1  -t1  * ]
      gamma = 1/(a1*b1)
      alpha = gamma*b1
      beta  = gamma*a1
      r1 = -beta*b2
      s1 = -beta*b3
      t1 = -beta*b4
      r0 = -alpha*(a3 + a2*r1)
      s0 = -alpha*(a4 + a2*s1)
      t0 = -alpha*(a5 + a2*t1)
      
      # Find D1 + D2 by applying a change of basis to the kernel of the above matrix.
      #
      # [ u0  v0  w0 ]   [ f[0]  g[0]  h[0]    0     0  ]
      # [ u1  v1  w1 ]   [ f[1]  g[1]  h[1]  f[0]  g[0] ] [ r0  s0  t0 ]
      # [ u2  v2  w2 ]   [ f[2]  g[2]  h[2]    0     0  ] [ r1  s1  t1 ]
      # [ u3  v3  w3 ] = [   1     0     0   f[1]  g[1] ]*[ 1   0   0  ]
      # [ u4  v4  w4 ]   [   0     1     0   f[2]  g[2] ] [ 0   1   0  ]
      # [ 1   0   0  ]   [   0     0     1     0     0  ] [ 0   0   1  ]
      # [ 0   1   0  ]   [   0     0     0     1     0  ]
      # [ 0   0   1  ]   [   0     0     0     0     1  ]
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
      
      new_f = [u0, u1, u2, u3, u4, 1]
      new_g = [v0, v1, v2, v3, v4, 0, 1]
      new_h = [w0, w1, w2, w3, w4, 0, 0, 1]
      # D1 + D2 is of type 51
    elif (b2 != 0) :
      # Reduce the matrix
      #
      # [ a1  a2  a3  a4  a5  a6 ]
      # [ 0   0   b2  b3  b4  b5 ]
      #
      # to its RREF
      #
      # [ 1  -r0  0  -s0  *  * ]
      # [ 0   0   1  -s1  *  * ]
      gamma = 1/(a1*b2)
      alpha = gamma*b2
      beta  = gamma*a1
      s1 = -beta*b3
      r0 = -alpha*a2
      s0 = -alpha*(a4 + a3*s1)

      # Find D1 + D2 by applying a change of basis to the kernel of the above matrix.
      #
      # [ u0  v0 ]   [ f[0]  g[0]  h[0]    0  ]
      # [ u1  v1 ]   [ f[1]  g[1]  h[1]  f[0] ] [ r0  s0 ]
      # [ u2  v2 ]   [ f[2]  g[2]  h[2]    0  ] [ 1   0  ]
      # [ u3  v3 ] = [   1     0     0   f[1] ]*[ 0   s1 ]
      # [ 1   v4 ]   [   0     1     0   f[2] ] [ 0   1  ]
      # [ 0   v5 ]   [   0     0     1     0  ]
      # [ 0   1  ]   [   0     0     0     1  ]
      u0 = f[0]*r0 + g[0]
      u1 = f[1]*r0 + g[1]
      u2 = f[2]*r0 + g[2]
      u3 = r0
      v0 = f[0]*s0 + h[0]*s1
      v1 = f[1]*s0 + h[1]*s1 + f[0]
      v2 = f[2]*s0 + h[2]*s1
      v3 = s0 + f[1]
      v4 = f[2]
      v5 = s1
      
      # Kill off the xy term in v = x^3 + v5*y^2 + v4*xy + v3*x^2 + v2*y + v1*x + v0
      # by subtracting v4*(xy + u3*x^2 + u2*y + u1*x + u0)
      # XXX : Doing this step simultaneously with the above saves 1M
      v0 = v0 - v4*u0
      v1 = v1 - v4*u1
      v2 = v2 - v4*u2
      v3 = v3 - v4*u3
      
      new_f = [u0, u1, u2, u3, 1]
      new_g = [v0, v1, v2, v3, 0, v5, 1]
      # D1 + D2 is of type 53
    elif (b3 != 0) : 
      # Reduce the matrix
      #
      # [ a1  a2  a3  a4  a5  a6 ]
      # [ 0   0   0   b3  b4  b5 ]
      #
      # to its RREF
      #
      # [ 1  -r0  -s0  0  *  * ]
      # [ 0   0    0   1  *  * ]
      alpha = 1/a1
      r0 = -alpha*a2
      s0 = -alpha*a3

      # Find D1 + D2 by applying a change of basis to the kernel of the above matrix.
      #
      # [ u0  v0 ]   [ f[0]  g[0]  h[0] ]
      # [ u1  v1 ]   [ f[1]  g[1]  h[1] ] [ r0  s0 ]
      # [ u2  v2 ]   [ f[2]  g[2]  h[2] ]*[ 1   0  ]
      # [ u3  v3 ] = [   1     0     0  ] [ 0   1  ]
      # [ 1   0  ]   [   0     1     0  ]
      # [ 0   1  ]   [   0     0     1  ]
      u0 = f[0]*r0 + g[0]
      u1 = f[1]*r0 + g[1]
      u2 = f[2]*r0 + g[2]
      u3 = r0
      v0 = f[0]*s0 + h[0]
      v1 = f[1]*s0 + h[1]
      v2 = f[2]*s0 + h[2]
      v3 = s0
      
      new_f = [u0, u1, u2, u3, 1]
      new_g = [v0, v1, v2, v3, 0, 1]
      # D1 + D2 is of type 52
    else :
      raise NotImplementedError("Divisors are non-disjoint.")
  elif (a2 != 0) or (a8 != 0) :
    a6  = - F[0]*a3
    a12 = - F[0]*a9
    
    if (a2 == 0) :
      a2, a3, a5, a6, a8, a9, a11, a12 = a8, a9, a11, a12, a2, a3, a5, a6
    
    # We have the matrix
    #
    # [ 0  a2  a3  0  a5   a6  ]
    # [ 0  a8  a9  0  a11  a12 ]
    #
    # Compute its row echelon form
    #
    # [ 0  a2  a3  0  a5  a6 ]
    # [ 0  0   b1  0  b2  b3 ]
    b1 = a2*a9  - a3*a8
    b2 = a2*a11 - a5*a8
    b3 = a2*a12 - a6*a8
    
    #print "M"
    #print Matrix(K, 2, 6, [0, a2, a3, 0, a5, a6, 0, a8, a9, 0, a11, a12])
    #print
    #print "MREF"
    #print Matrix(K, 2, 6, [0, a2, a3, 0, a5, a6, 0, 0, b1, 0, b2, b3])
    #print

    r0, r1 = 0, 0    
    if (b1 != 0) :
      # Compute the reduced row echelon form
      #
      # [ 0  1  0  *  *  -r0 ]
      # [ 0  0  1  *  *  -r1 ]

      # [ 0  1  0   alpha*(a6 + a3*r1) ]
      # [ 0  0  1  -r1 ]
      gamma = 1/(a2*b1)
      alpha = gamma*b1
      beta = gamma*a2
      r1 = -beta*b3
      r0 = -alpha*(a6 + a3*r1)

    elif (b2 != 0) :
      # Compute the reduced row echelon form
      #
      # [ 0  1  *  *  0  -r0 ]
      # [ 0  0  0  0  1  -r1 ]
      gamma = 1/(a2*b2)
      alpha = gamma*b2
      beta = gamma*a2
      r1 = -beta*b3
      r0 = -alpha*(a6 + a5*r1)
    else :
      raise NotImplementedError("Divisors are non-disjoint.")

    # Find D1 + D2 by applying a change of basis to the kernel of the above matrix.
    #
    # [ f[0]  u0 ]   [ f[0]  g[0]  h[0]    0     0     0  ]
    # [ f[1]  u1 ]   [ f[1]  g[1]  h[1]  f[0]  g[0]  h[0] ][ 1  0  ]
    # [ f[2]  u2 ]   [ f[2]  g[2]  h[2]    0     0     0  ][ 0  r0 ]
    # [   1   u3 ] = [   1     0     0   f[1]  g[1]  h[1] ][ 0  r1 ]
    # [   0   u4 ]   [   0     1     0   f[2]  g[2]  h[2] ][ 0  0  ]
    # [   0   u5 ]   [   0     0     1     0     0     0  ][ 0  0  ]
    # [   0   0  ]   [   0     0     0     1     0     0  ][ 0  1  ]
    # [   0   0  ]   [   0     0     0     0     1     0  ]
    # [   0   1  ]   [   0     0     0     0     0     1  ]
    #
    # Also reduce u modulo f
    u0 = g[0]*r0 + h[0]*r1 - h[1]*f[0]
    u1 = g[1]*r0 + h[1]*(r1 - f[1]) + h[0]
    u2 = g[2]*r0 + h[2]*r1 - h[1]*f[2]
    u4 = r0 + h[2]
    u5 = r1
    
    new_f = [f[0], f[1], f[2], 1]
    new_g = [u0, u1, u2, 0, u4, u5, 0, 0, 1]
    # D1 + D2 is of type 54
  else :
    raise NotImplementedError("Divisors are non-disjoint.")

  # XXX : Old count, 1I 27M for just f'', g''.
  return C34CrvDiv(D1.C, [new_f, new_g, new_h])



def add_31_31(D1, D2) :
  C = D1.C
  K = C.K
  f, g, h = D1.f, D1.g, D1.h
  F, G, H = D2.f, D2.g, D2.h
  new_f, new_g, new_h = [], [], []

  # Compute M
  #
  #     [ a1   a2   a3   a4   a5   a6   a7  ]
  # M = [ a8   a9   a10  a11  a12  a13  a14 ]
  #     [ a15  a16  a17  a18  a19  a20  a21 ]
  #
  # Columns 4, 5, and 6 are computed by
  #
  # [ a4   a5   a6  ]   [ 0  -F[0]  -G[0] ] [ a1   a2   a3  ] 
  # [ a11  a12  a13 ] = [ 1  -F[1]  -G[1] ]*[ a8   a9   a10 ]
  # [ a18  a19  a20 ]   [ 0  -F[2]  -G[2] ] [ a15  a16  a17 ]
  #
  # The last column similarly by
  #
  # [ a7  ]   [ 0  -F[0]  -G[0] ] [ a4  ] 
  # [ a14 ] = [ 1  -F[1]  -G[1] ]*[ a11 ]
  # [ a21 ]   [ 0  -F[2]  -G[2] ] [ a18 ]
  #
  # The last column is computed when needed, which is rarely.
  
  a1  = f[0] - F[0]
  a2  = g[0] - G[0]
  a3  = h[0] - H[0]
  a8  = f[1] - F[1]
  a9  = g[1] - G[1]
  a10 = h[1] - H[1]
  a15 = f[2] - F[2]
  a16 = g[2] - G[2]
  a17 = h[2] - H[2]
  
  a4  =    - F[0]*a8  - G[0]*a15
  a5  =    - F[0]*a9  - G[0]*a16
  a6  =    - F[0]*a10 - G[0]*a17
  a11 = a1 - F[1]*a8  - G[1]*a15
  a12 = a2 - F[1]*a9  - G[1]*a16
  a13 = a3 - F[1]*a10 - G[1]*a17
  a18 =    - F[2]*a8  - G[2]*a15
  a19 =    - F[2]*a9  - G[2]*a16
  a20 =    - F[2]*a10 - G[2]*a17
  aswap = 0 # If we swap the first row of the matrix with another, this gets updated.

  #print "M = "
  #print Matrix(K, 3, 6, [a1, a2, a3, a4, a5, a6, a8, a9, a10, a11, a12, a13, a15, a16, a17, a18, a19, a20])
  #print
  
  if (a1 != 0 or a8 != 0 or a15 != 0) :
    # Swap rows, if necessary, to ensure a1 != 0
    if (a1 == 0) :
      if (a8 != 0) :
        a1, a2, a3, a4, a5, a6, a8, a9, a10, a11, a12, a13 = a8, a9, a10, a11, a12, a13, a1, a2, a3, a4, a5, a6
        aswap = 1
      else :
        a1, a2, a3, a4, a5, a6, a15, a16, a17, a18, a19, a20 = a15, a16, a17, a18, a19, a20, a1, a2, a3, a4, a5, a6
        aswap = 2

    # Partially reduce M.
    # Compute the matrix
    #     [ a1  a2  a3  a4  a5   a6   a7  ]
    # M = [ 0   b1  b2  b3  b4   b5   b6  ]
    #     [ 0   b7  b8  b9  b10  b11  b12 ]
    
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
    
    #print "M' = "
    #print Matrix(K, 3, 6, [a1, a2, a3, a4, a5, a6, 0, b1, b2, b3, b4, b5, 0, b7, b8, b9, b10, b11])
    #print
    
    if (b1 != 0 or b7 != 0) :
      # Before reducing any more, ensure b1 != 0 by swapping rows, if necessary
      if (b1 == 0) :
        b1, b2, b3, b4, b5, b7, b8, b9, b10, b11 = b7, b8, b9, b10, b11, b1, b2, b3, b4, b5
      # Continue reducing M. Compute
      #     [ a1  a2  a3  a4  a5  a6  a7 ]
      # M = [ 0   b1  b2  b3  b4  b5  b6 ]
      #     [ 0   0   c1  c2  c3  c4  c5 ]
      c1 = b1*b8  - b2*b7
      c2 = b1*b9  - b3*b7
      c3 = b1*b10 - b4*b7
      c4 = b1*b11 - b5*b7
      
      #print "M_ref = "
      #print Matrix([
      #  [a1, a2, a3, a4, a5, a6],
      #  [0,  b1, b2, b3, b4, b5],
      #  [0,  0,  c1, c2, c3, c4]])
      #print
      
      if (c1 != 0) :
        # Now compute reduced row echelon form of (the first six columns of) M
        #         [ 1  0  0  -r0  -s0  -t0 ]
        # Mrref = [ 0  1  0  -r1  -s1  -t1 ]
        #         [ 0  0  1  -r2  -s2  -t2 ]
        
        # First compute inverses of a1, b1, c1
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
        #print "M_rref = "
        #print Matrix([
        #  [1, 0, 0, -r0, -s0, -t0],
        #  [0, 1, 0, -r1, -s1, -t1],
        #  [0, 0, 1, -r2, -s2, -t2]])
        # The kernel of M is
        #          [ r0  s0  t0 ]
        #          [ r1  s1  t1 ]
        # ker(M) = [ r2  s2  t2 ]
        #          [ 1   0   0  ]
        #          [ 0   1   0  ]
        #          [ 0   0   1  ]
        #
        # Now perform a change of basis to get D1 + D2
        # [ u0  v0  w0 ]   [ f[0]  g[0]  h[0]    0     0     0  ]
        # [ u1  v1  w1 ]   [ f[1]  g[1]  h[1]  f[0]  g[0]  h[0] ]   [ r0  s0  t0 ]
        # [ u2  v2  w2 ]   [ f[2]  g[2]  h[2]    0     0     0  ]   [ r1  s1  t1 ]
        # [ u3  v3  w3 ]   [   1     0     0   f[1]  g[1]  h[1] ]   [ r2  s2  t2 ]
        # [ u4  v4  w4 ] = [   0     1     0   f[2]  g[2]  h[2] ] * [ 1   0   0  ]
        # [ u5  v5  w5 ]   [   0     0     1     0     0     0  ]   [ 0   1   0  ]
        # [ 1   0   0  ]   [   0     0     0     1     0     0  ]   [ 0   0   1  ]
        # [ 0   1   0  ]   [   0     0     0     0     1     0  ]
        # [ 0   0   1  ]   [   0     0     0     0     0     1  ]
        u0 = f[0]*r0 + g[0]*r1 + h[0]*r2
        u1 = f[1]*r0 + g[1]*r1 + h[1]*r2 + f[0]
        u2 = f[2]*r0 + g[2]*r1 + h[2]*r2
        u3 = r0 + f[1]
        u4 = r1 + f[2]
        u5 = r2
        v0 = f[0]*s0 + g[0]*s1 + h[0]*s2
        v1 = f[1]*s0 + g[1]*s1 + h[1]*s2 + g[0]
        v2 = f[2]*s0 + g[2]*s1 + h[2]*s2
        v3 = s0 + g[1]
        v4 = s1 + g[2]
        v5 = s2
        w0 = f[0]*t0 + g[0]*t1 + h[0]*t2
        w1 = f[1]*t0 + g[1]*t1 + h[1]*t2 + h[0]
        w2 = f[2]*t0 + g[2]*t1 + h[2]*t2
        w3 = t0 + h[1]
        w4 = t1 + h[2]
        w5 = t2
        new_f = [u0, u1, u2, u3, u4, u5, 1]
        new_g = [v0, v1, v2, v3, v4, v5, 0, 1]
        new_h = [w0, w1, w2, w3, w4, w5, 0, 0, 1]
        
        # Requires   1I 67M to compute f', g'
        # Additional 0I 27M to compute h'
        # Total :    1I 94M
      elif (c2 != 0) :
        # D1 + D2 will be of type 63.
        # We need only the first 5 columns of M.
        # So we have
        #
        #     [ a1  a2  a3  a4  a5 ]
        # M = [ 0   b1  b2  b3  b4 ]
        #     [ 0   0   0   c2  c3 ]
        #
        # We wish to reduce it to
        #
        #         [ 1  0  -r0  0  -s0 ]
        # Mrref = [ 0  1  -r1  0  -s1 ]
        #         [ 0  0   0   1  -s2 ]
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
        # The kernel of M is
        #          [ r0  s0 ]
        #          [ r1  s1 ]
        # ker(M) = [ 1   0  ]
        #          [ 0   s2 ]
        #          [ 0   1  ]
        #
        # Now perform a change of basis to get D1 + D2
        # [ u0  v0 ]   [ f[0]  g[0]  h[0]    0     0  ]
        # [ u1  v1 ]   [ f[1]  g[1]  h[1]  f[0]  g[0] ]   [ r0  s0 ]
        # [ u2  v2 ]   [ f[2]  g[2]  h[2]    0     0  ]   [ r1  s1 ]
        # [ u3  v3 ]   [   1     0     0   f[1]  g[1] ] * [ 1   0  ]
        # [ u4  v4 ] = [   0     1     0   f[2]  g[2] ]   [ 0   s2 ]
        # [ 1   0  ]   [   0     0     1     0     0  ]   [ 0   1  ]
        # [ 0   v5 ]   [   0     0     0     1     0  ]
        # [ 0   1  ]   [   0     0     0     0     1  ]
        u0 = f[0]*r0 + g[0]*r1 + h[0]
        u1 = f[1]*r0 + g[1]*r1 + h[1]
        u2 = f[2]*r0 + g[2]*r1 + h[2]
        u3 = r0
        u4 = r1
        v0 = f[0]*s0 + g[0]*s1
        v1 = f[1]*s0 + g[1]*s1 + f[0]*s2 + g[0]
        v2 = f[2]*s0 + g[2]*s1
        v3 = s0 + f[1]*s2 + g[1]
        v4 = s1 + f[2]*s2 + g[2]
        v5 = s2
        new_f = [u0, u1, u2, u3, u4, 1]
        new_g = [v0, v1, v2, v3, v4, 0, v5, 1]

      elif (c3 != 0) :
        # D1 + D2 will be of type 62.
        # We need only the first 4 columns of M.
        # So we have
        #
        #     [ a1  a2  a3  a4 ]
        # M = [ 0   b1  b2  b3 ]
        #     [ 0   0   0   0  ]
        #
        # We wish to reduce it to
        #
        #         [ 1  0  -r2  -s0 ]
        # Mrref = [ 0  1  -r1  -s1 ]
        #         [ 0  0   0    0  ]
        ab = a1*b1
        gamma = 1/ab
        alpha = gamma*b1
        beta  = gamma*a1
        
        r1 = -beta*b2
        s1 = -beta*b3
        r0 = -alpha*(a2*r1 + a3)
        s0 = -alpha*(a2*s1 + a4)
        # The kernel of M is
        #          [ r0  s0 ]
        # ker(M) = [ r1  s1 ]
        #          [ 1   0  ]
        #          [ 0   1  ]
        #
        # Now perform a change of basis to get D1 + D2
        # [ u0  v0 ]   [ f[0]  g[0]  h[0]    0  ]
        # [ u1  v1 ]   [ f[1]  g[1]  h[1]  f[0] ]   [ r0  s0 ]
        # [ u2  v2 ]   [ f[2]  g[2]  h[2]    0  ] * [ r1  s1 ]
        # [ u3  v3 ]   [   1     0     0   f[1] ]   [ 1   0  ]
        # [ u4  v4 ] = [   0     1     0   f[2] ]   [ 0   1  ]
        # [ 1   0  ]   [   0     0     1     0  ]
        # [ 0   1  ]   [   0     0     0     1  ]
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
        new_f = [u0, u1, u2, u3, u4, 1]
        new_g = [v0, v1, v2, v3, v4, 0, 1]

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
        L = C34CrvDiv(D1.C, [[u0, u1, u2, u3, u4, 1], [v0, v1, v2, v3, v4, 0, 1], [w0, w1, w2, w3, w4, 0, 0, 1]])

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
        G = C34CrvDiv(C, [[p0, 1], [q0, 0, 1], []])
        return flip(flip(L)) + G

    elif (b2 != 0 or b8 != 0) :
      a7, a14, a21 = 0, 0, 0
      if (aswap == 0) :
        a7  =    - F[0]*a11 - G[0]*a18
        a14 = a4 - F[1]*a11 - G[1]*a18
        a21 =    - F[2]*a11 - G[2]*a18
      elif (aswap == 1) :
        a7  = a11 - F[1]*a4 - G[1]*a18
        a14 =     - F[0]*a4 - G[0]*a18
        a21 =     - F[2]*a4 - G[2]*a18
      else :
        a7  =     - F[2]*a11 - G[2]*a4
        a14 = a18 - F[1]*a11 - G[1]*a4
        a21 =     - F[0]*a11 - G[0]*a4
      b6  = a1*a14 - a7*a8
      b12 = a1*a21 - a7*a15

      # Swap rows if necessary so that b2 != 0
      if (b2 == 0) :
        b2, b3, b4, b5, b6, b8, b9, b10, b11, b12 = b8, b9, b10, b11, b12, b2, b3, b4, b5, b6
      
      # D1 + D2 will be a type 64 divisor.
      #
      # We have the matrix
      #     [ a1  a2  a3  a4  a5   a6   a7  ]
      # M = [ 0   0   b2  b3  b4   b5   b6  ]
      #     [ 0   0   b8  b9  b10  b11  b12 ]
      #
      # We wish to reduce it to
      #
      #     [ a1  a2  a3  a4  a5  a6  a7 ]
      # M = [ 0   0   b2  b3  b4  b5  b6 ]
      #     [ 0   0   0   c1  c2  c3  c4 ]
      c1 = b2*b9  - b3*b8
      c2 = b2*b10 - b4*b8
      c3 = b2*b11 - b5*b8
      c4 = b2*b12 - b6*b8
      
      #print "M'' = "
      #print Matrix(K, 3, 7, [a1, a2, a3, a4, a5, a6, a7, 0, 0, b2, b3, b4, b5, b6, 0, 0, 0, c1, c2, c3, c4])
      #print
      
      r0, s0, s1, s2 = 0, 0, 0, 0
      if (c1 != 0) :
        # We have the matrix
        #
        #     [ a1  a2  a3  a4  *  *  a7 ]
        # M = [ 0   0   b2  b3  *  *  b6 ]
        #     [ 0   0   0   c1  *  *  c4 ]
        #
        # We wish to reduce it to
        #
        #     [ 1  -r0  0  0  *  *  -s0 ]
        # M = [ 0   0   1  0  *  *  -s1 ]
        #     [ 0   0   0  1  *  *  -s2 ]
        #
        # The values in place of the asterisks are not needed.
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
        
      elif (c2 != 0) :
        # We have the matrix
        #
        #     [ a1  a2  a3  a4  a5  a6  a7 ]
        # M = [ 0   0   b2  b3  b4  b5  b6 ]
        #     [ 0   0   0   0   c2  c3  c4 ]
        #
        # We wish to reduce it to
        #
        #     [ 1  -r0  0  *  0  *  -s0 ]
        # M = [ 0   0   1  *  0  *  -s1 ]
        #     [ 0   0   0  0  1  *  -s2 ]
        #
        # The values in place of the asterisks are not needed.
        ab = a1*b2
        abc = ab*c2
        delta = 1/abc
        alpha = delta*b2*c2
        beta  = delta*a1*c2
        gamma = delta*ab
        
        r0 = -alpha*a2
        s2 = -gamma*c4
        s1 = -beta*(b4*s2 + b6)
        s0 = -alpha*(a3*s1 + a5*s2 + a7)
        
      elif (c3 != 0) :
        # We have the matrix
        #
        #     [ a1  a2  a3  a4  a5  a6  a7 ]
        # M = [ 0   0   b2  b3  b4  b5  b6 ]
        #     [ 0   0   0   0   c2  c3  c4 ]
        #
        # We wish to reduce it to
        #
        #     [ 1  -r0  0  *  *  0  -s0 ]
        # M = [ 0   0   1  *  *  0  -s1 ]
        #     [ 0   0   0  0  0  1  -s2 ]
        #
        # The values in place of the asterisks are not needed.
        ab = a1*b2
        abc = ab*c2
        delta = 1/abc
        alpha = delta*b2*c2
        beta  = delta*a1*c2
        gamma = delta*ab
        
        r0 = -alpha*a2
        s2 = -gamma*c4
        s1 = -beta*(b5*s2 + b6)
        s0 = -alpha*(a3*s1 + a6*s2 + a7)
        
      else :
        assert c4 == 0
        # D1 and D2 are non-disjoint.
        # We compute D1 + D2 = lcm(D1, D2) + gcd(D1, D2)
        #
        # Reduce
        #
        # [ a1  a2  a3  a4  *  *  * ]
        # [ 0   0   b2  b3  *  *  * ]
        # [ 0   0   0   0   0  0  0 ]
        # 
        # to its reduced row echelon form.
        #
        # [ 1  -r0  0  -s0  *  *  * ]
        # [ 0   0   1  -s1  *  *  * ]
        # [ 0   0   0   0   0  0  0 ]
        gamma = 1/(a1*b2)
        alpha = gamma*b2
        beta = gamma*a1
        s1 = -beta*b3
        r0 = -alpha*a2
        s0 = -alpha*(a4 + a3*s1)
        
        # [ u0  v0 ]   [ f[0]  g[0]  h[0]    0  ]
        # [ u1  v1 ]   [ f[1]  g[1]  h[1]  f[0] ]   [ r0  s0 ]
        # [ u2  v2 ]   [ f[2]  g[2]  h[2]    0  ] * [ 1   0  ]
        # [ u3  v3 ]   [   1     0     0   f[1] ]   [ 0   s1 ]
        # [ 1   v4 ] = [   0     1     0   f[2] ]   [ 0   1  ]
        # [ 0   v5 ]   [   0     0     1     0  ]
        # [ 0   1  ]   [   0     0     0     1  ]
        u0 = f[0]*r0 + g[0]
        u1 = f[1]*r0 + g[1]
        u2 = f[2]*r0 + g[2]
        u3 = r0
        v0 = f[0]*s0 + h[0]*s1 - f[2]*u0
        v1 = f[1]*s0 + h[1]*s1 + f[0] - f[2]*u1
        v2 = f[2]*(s0 - u2) + h[2]*s1
        v3 = s0 + f[1] - f[2]*u3
        v5 = s1
        L = C34CrvDiv(D1.C, [[u0, u1, u2, u3, 1], [v0, v1, v2, v3, 0, v5, 1], []])

        # GCD(D1, D2) is type 11.
        # A basis for the GCD is given by the 1st and 3rd columns of M
        #
        # [ m6  m3 ]     [ n2 m3 ]     [ p0  q0 ]
        # [ m5  m2 ] ==> [ n1 m2 ] ==> [ 1   0  ]
        # [ m4  m1 ]     [ 0  m1 ]     [ 0   1  ]
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
        G = C34CrvDiv(C, [[p0, 1], [q0, 0, 1], []])
        
        return flip(flip(L)) + G

      # The kernel of M is
      #          [ r0  s0 ]
      #          [ 1   0  ]
      # ker(M) = [ 0   s1 ]
      #          [ 0   s2 ]
      #          [ 0   1  ]
      #
      # Now perform a change of basis to get D1 + D2
      # [ u0  v0 ]   [ f[0]  g[0]  h[0]    0     0  ]
      # [ u1  v1 ]   [ f[1]  g[1]  h[1]  f[0]    0  ]
      # [ u2  v2 ]   [ f[2]  g[2]  h[2]    0     0  ]   [ r0  s0 ]
      # [ u3  v3 ]   [   1     0     0   f[1]  f[0] ]   [ 1   0  ]
      # [ 1   v4 ] = [   0     1     0   f[2]    0  ] * [ 0   s1 ]
      # [ 0   v5 ]   [   0     0     1     0     0  ]   [ 0   s2 ]
      # [ 0   v6 ]   [   0     0     0     1   f[1] ]   [ 0   1  ]
      # [ 0   v7 ]   [   0     0     0     0   f[2] ]
      # [ 0   0  ]   [   0     0     0     0     0  ]
      # [ 0   1  ]   [   0     0     0     0     1  ]
      u0 = f[0]*r0 + g[0]
      u1 = f[1]*r0 + g[1]
      u2 = f[2]*r0 + g[2]
      u3 = r0
      v0 = f[0]*s0 + h[0]*s1
      v1 = f[1]*s0 + h[1]*s1 + f[0]*s2
      v2 = f[2]*s0 + h[2]*s1
      v3 = s0 + f[1]*s2 + f[0]
      v4 = f[2]*s2
      v5 = s1
      v6 = s2 + f[1]
      v7 = f[2]
      
      # g'' is x^4 + v7*x^2*y + v6*x^3 + v5*y^2 + v4*x*y + v3*x^2 + v2*y + v1*x + v0
      # Kill off monomials divisible by xy by subtracting appropriate multiples of f'' = xy + u3*x^2 + u2*y + u1*x + u0.
      z = v4 - v7*u2
      v0 = v0 - z*u0
      v1 = v1 - v7*u0 - z*u1
      v2 = v2 - z*u2
      v3 = v3 - v7*u1 - z*u3
      #v4 = 0
      #v5 = v5
      v6 = v6 - v7*u3
      #v7 = 0
      new_f = [u0, u1, u2, u3, 1]
      new_g = [v0, v1, v2, v3, 0, v5, v6, 0, 0, 1]
    elif (b3 != 0) or (b9 != 0) :
      #     [ a1  a2  a3  a4  a5  a6  a7 ]
      # M = [ 0   0   0   b3  b4  b5  b6 ]
      #     [ 0   0   0   0   0   0   0  ]
      #
      # Compute M_RREF
      #
      # [ 1  -r0  -s0  0  *  *  * ]
      # [ 0   0    0   1  *  *  * ]
      # [ 0   0    0   0  0  0  0 ]
      #
      # LCM(D1, D2) is type 52.
      alpha = 1/a1
      r0 = -alpha*a2
      s0 = -alpha*a3
      u0 = f[0]*r0 + g[0]
      u1 = f[1]*r0 + g[1]
      u2 = f[2]*r0 + g[2]
      u3 = r0
      v0 = f[0]*s0 + h[0]
      v1 = f[1]*s0 + h[1]
      v2 = f[2]*s0 + h[2]
      v3 = s0
      L = C34CrvDiv(C, [[u0, u1, u2, u3, 1], [v0, v1, v2, v3, 0, 1], []])
      
      # GCD(D1, D2) is type 11.
      # A basis for the GCD is given by the 1st and 4th columns of M
      #
      # [ m6  m3 ]     [ n2 m3 ]     [ p0  q0 ]
      # [ m5  m2 ] ==> [ n1 m2 ] ==> [ 1   0  ]
      # [ m4  m1 ]     [ 0  m1 ]     [ 0   1  ]
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
      G = C34CrvDiv(C, [[p0, 1], [q0, 0, 1], []])
      return flip(flip(L)) + G
    else :
      # D1 and D2 are non-disjoint and their intersection is type 41.
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
      L = C34CrvDiv(D1.C, [[u0, u1, u2, u3, 1], [v0, v1, v2, v3, 0, 1], [w0, w1, w2, w3, 0, 0, 1]])

      # GCD(D1, D2) is degree 2 (type 21 or 22).
      # One generator of the GCD is given by the 1st column of M
      #
      # [ m3 ]     [ u0 ]    [ u0 ]
      # [ m2 ] ==> [ u1 ] or [ 1  ]
      # [ m1 ]     [ 1  ]    [ 0  ]
      if (aswap == 0) :
        m1, m2, m3 = a15, a8, a1
      elif (aswap == 1) :
        m1, m2, m3 = a15, a1, a8
      else :
        m1, m2, m3 = a1, a8, a15
      if (m1 != 0) :
        mu = 1/m1
        u1 = mu*m2
        u0 = mu*m3
        # Compute v = x^2 + v1*x + v0 by taking f modulo u = y + u1*x + u0
        v0 = f[0] - f[2]*u0
        v1 = f[1] - f[2]*u1
        G = C34CrvDiv(C, [[u0, u1, 1], [v0, v1, 0, 1], []])
      else :
        assert m2 != 0
        mu = 1/m2
        u0 = mu*m3
        # Compute v = y^2 + v2*y + v0 by taking h modulo u = x + u0
        v0 = h[0] - h[1]*u0
        v2 = h[2]
        G = C34CrvDiv(C, [[u0, 1], [v0, 0, v2, 0, 0, 1], []])
      return flip(flip(L)) + G

  elif (a2 != 0 or a9 != 0 or a16 != 0) :
    if (a2 == 0) :
      if (a9 != 0) :
        a2, a3, a5, a6, a9, a10, a12, a13 = a9, a10, a12, a13, a2, a3, a5, a6
        aswap = 1
      else :
        a2, a3, a5, a6, a16, a17, a19, a20 = a16, a17, a19, a20, a2, a3, a5, a6
        aswap = 2

    # Here, a1, a8, a15 are zero and M is the matrix
    #
    #     [ 0  a2   a3   0  a5   a6   0 ]
    # M = [ 0  a9   a10  0  a12  a13  0 ]
    #     [ 0  a16  a17  0  a19  a20  0 ]
    #
    # Reduce M to the matrix
    #
    #     [ 0  a2  a3  0  a5  a6  0 ]
    # M = [ 0  0   b1  0  b2  b3  0 ]
    #     [ 0  0   b4  0  b5  b6  0 ]
    b1 = a2*a10 - a3*a9
    b2 = a2*a12 - a5*a9
    b3 = a2*a13 - a6*a9
    b4 = a2*a17 - a3*a16
    b5 = a2*a19 - a5*a16
    b6 = a2*a20 - a6*a16
    
    if (b1 != 0 or b4 != 0) :
      if (b1 == 0) :
        b1, b2, b3, b4, b5, b6 = b4, b5, b6, b1, b2, b3
      # Reduce M to the matrix
      #
      #     [ 0  a2  a3  0  a5  a6  0 ]
      # M = [ 0  0   b1  0  b2  b3  0 ]
      #     [ 0  0   0   0  c1  c2  0 ]
      c1 = b1*b5 - b2*b4
      c2 = b1*b6 - b3*b4
      
      assert c1 == 0 # If c1 != 0, then M has full rank, but produces a degree 5 (type 54) divisor, a contradiction.
      if (c2 != 0) :
        new_f = [f[0], f[1], f[2], 1]
      else :
        # Then LCM(D1, D2) is type 54
        gamma = 1/(a2*b1)
        alpha = gamma*b1
        beta = gamma*a2
        r1 = -beta*b3
        r0 = -alpha*(a6 + a3*r1)
        u0 = g[0]*r0 + h[0]*r1 - f[0]*h[1]
        u1 = g[1]*r0 + h[1]*(r1 - f[1]) + h[0]
        u2 = g[2]*r0 + h[2]*r1 - f[2]*h[1]
        u4 = r0 + h[2]
        u5 = r1
        L = C34CrvDiv(C, [[f[0], f[1], f[2], 1], [u0, u1, u2, 0, u4, u5, 0, 0, 1], []])

        # GCD(D1, D2) is type 11.
        # A basis for the GCD is given by the 2nd and 3rd columns of M
        #
        # [ m6  m3 ]     [ n2 m3 ]     [ p0  q0 ]
        # [ m5  m2 ] ==> [ n1 m2 ] ==> [ 1   0  ]
        # [ m4  m1 ]     [ 0  m1 ]     [ 0   1  ]
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
        G = C34CrvDiv(C, [[p0, 1], [q0, 0, 1], []])
        return flip(flip(L)) + G

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
      u0 = g[0]*r0 + h[0]
      u1 = g[1]*r0 + h[1]
      u2 = g[2]*r0 + h[2]
      u4 = r0
      L = C34CrvDiv(C, [[f[0], f[1], f[2], 1], [u0, u1, u2, 0, u4, 1], []])

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
        v0 = f[0] - f[2]*u0
        v1 = f[1] - f[2]*u1
        G = C34CrvDiv(C, [[u0, u1, 1], [v0, v1, 0, 1], []])
      else :
        assert m2 != 0
        mu = 1/m2
        u0 = mu*m3
        # Compute v = y^2 + v2*y + v0 by taking h modulo u = x + u0
        v0 = h[0] - h[1]*u0
        v2 = h[2]
        G = C34CrvDiv(C, [[u0, 1], [v0, 0, v2, 0, 0, 1], []])
      return flip(flip(L)) + G

  else :
    # We have the matrix
    #
    #     [ 0  0  a3   0  0  a6   0 ]
    # M = [ 0  0  a10  0  0  a13  0 ]
    #     [ 0  0  a17  0  0  a20  0 ]
    #
    # LCM(D1, D2) will be type 42, given by polynomials f'' = f' = f and g'' = g' = g
    L = C34CrvDiv(C, [copy(D1.f), copy(D1.g), []])

    # GCD(D1, D2) is degree 2 (type 21 or 22).
    # One generator of the GCD is given by the 3rd column of M
    #
    # [ m3 ]     [ u0 ]    [ u0 ]
    # [ m2 ] ==> [ u1 ] or [ 1  ]
    # [ m1 ]     [ 1  ]    [ 0  ]
    if (aswap == 0) :
      m1, m2, m3 = a17, a10, a3
    elif (aswap == 1) :
      m1, m2, m3 = a17, a3, a10
    else :
      m1, m2, m3 = a3, a10, a17
    if (m1 != 0) :
      mu = 1/m1
      u1 = mu*m2
      u0 = mu*m3
      # Compute v = x^2 + v1*x + v0 by taking f modulo u = y + u1*x + u0
      v0 = f[0] - f[2]*u0
      v1 = f[1] - f[2]*u1
      G = C34CrvDiv(C, [[u0, u1, 1], [v0, v1, 0, 1], []])
    else :
      assert m2 != 0
      mu = 1/m2
      u0 = mu*m3
      # Compute v = y^2 + v2*y + v0 by taking h modulo u = x + u0
      v0 = h[0] - h[1]*u0
      v2 = h[2]
      G = C34CrvDiv(C, [[u0, 1], [v0, 0, v2, 0, 0, 1], []])
    return flip(flip(L)) + G
  
  return C34CrvDiv(D1.C, [new_f, new_g, new_h])



def old_add_31_31(D1, D2):
  K = D1.K
  f, g, h = D1.f, D1.g, D1.h
  F, G, H = D2.f, D2.g, D2.h
  new_f, new_g, new_h = [], [], []

  if h == [] :
    inv = 1/f[2]
    h = [ inv*(       - g[0]*(f[1] - g[2]) + f[0]*g[1]),
          inv*(- g[0] - g[1]*(f[1] - g[2]) + f[1]*g[1]),
          inv*(  f[0] - g[2]*(f[1] - g[2]) + f[2]*g[1]),
          K.zero(), K.zero(), K.one()]
    D1.h = h
  if H == [] :
    inv = 1/F[2]
    H = [ inv*(       - G[0]*(F[1] - G[2]) + F[0]*G[1]),
          inv*(- G[0] - G[1]*(F[1] - G[2]) + F[1]*G[1]),
          inv*(  F[0] - G[2]*(F[1] - G[2]) + F[2]*G[1]),
          K.zero(), K.zero(), K.one()]
    D2.h = H
  
  
  # Compute M
  #     [ a1  a2  a3  a4  a5  ]
  # M = [ a6  a7  a8  a9  a10 ]
  #     [ a11 a12 a13 a14 a15 ]
  a1  = F[0] - f[0]
  a2  = G[0] - g[0]
  a3  = H[0] - h[0]
  a4  = - (F[2] - f[2])*g[0] - (F[1] - f[1])*f[0]
  a5  = - (G[2] - g[2])*g[0] - (G[1] - g[1])*f[0]
  a6  = F[1] - f[1]
  a7  = G[1] - g[1]
  a8  = H[1] - h[1]
  a9  = F[0] - f[0] - (F[2] - f[2])*g[1] - (F[1] - f[1])*f[1]
  a10 = G[0] - g[0] - (G[2] - g[2])*g[1] - (G[1] - g[1])*f[1]
  a11 = F[2] - f[2]
  a12 = G[2] - g[2]
  a13 = H[2] - h[2]
  a14 = - (F[2] - f[2])*g[2] - (F[1] - f[1])*f[2]
  a15 = - (G[2] - g[2])*g[2] - (G[1] - g[1])*f[2]
  # Subtotal : 0I 12M
  
  # Assuming M's left-most 3x3 submatrix is positive definite,
  # Compute row echelon form of M
  #         [ a1 a2 a3 a4 a5 ]
  # M_ref = [  0 b1 b2 b3 b4 ]
  #         [  0  0 c1 c2 c3 ]
  b1 = a1*a7  - a2*a6
  b2 = a1*a8  - a3*a6
  b3 = a1*a9  - a4*a6
  b4 = a1*a10 - a5*a6
  b5 = a1*a12 - a2*a11
  b6 = a1*a13 - a3*a11
  b7 = a1*a14 - a4*a11
  b8 = a1*a15 - a5*a11
  c1 = b1*b6 - b2*b5
  c2 = b1*b7 - b3*b5
  c3 = b1*b8 - b4*b5
  # Subtotal : 0I 22M

  # Compute reduced row echelon form of M
  #          [ 1  0  0  -v0  -w0 ]
  # M_rref = [ 0  1  0  -v1  -w1 ]
  #          [ 0  0  1  -v2  -w2 ]
  alpha = a1*b1*c1
  if (alpha == 0) :
    raise ValueError("M's left-most 3x3 submatrix is not positive definite.")
  alpha = 1/alpha
  beta  = alpha*b1*c1
  gamma = alpha*a1*c1
  delta = alpha*a1*b1
  v2 = -delta*(c2)
  v1 = -gamma*(b3 + b2*v2)
  v0 = -beta*(a4 + a2*v1 + a3*v2)
  w2 = -delta*(c3)
  w1 = -gamma*(b4 + b2*w2)
  w0 = -beta*(a5 + a2*w1 + a3*w2)
  # Subtotal : 1I 20M

  # Compute polynomials f'' and g'' defining the divisor D1 + D2
  s0 = F[0]*v0 + G[0]*v1 + H[0]*v2
  s1 = F[1]*v0 + G[1]*v1 + H[1]*v2 + F[0]
  s2 = F[2]*v0 + G[2]*v1 + H[2]*v2
  s3 = v0 + F[1]
  s4 = v1 + F[2]
  s5 = v2
  t0 = F[0]*w0 + G[0]*w1 + H[0]*w2
  t1 = F[1]*w0 + G[1]*w1 + H[1]*w2 + G[0]
  t2 = F[2]*w0 + G[2]*w1 + H[2]*w2
  t3 = w0 + G[1]
  t4 = w1 + G[2]
  t5 = w2
  new_f = [s0, s1, s2, s3, s4, s5, K.one()]
  new_g = [t0, t1, t2, t3, t4, t5, K.zero(), K.one()]
  # Subtotal : 0I 18M

  # Total : 1I 72M
  return C34CrvDiv(D1.C, [new_f, new_g, new_h])

