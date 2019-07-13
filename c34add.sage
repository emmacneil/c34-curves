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
  if (not D1.reduced) :
    return add(flip(flip(D1)), D2)
  if (not D2.reduced) :
    return add(D1, flip(flip(D2)))

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

  ret = 0
  # Examine the degrees of D1 and D2 and call the appropriate function
  if (D1.degree, D2.degree) == (1, 1) :
    ret = add_11_11(D1, D2)
  elif (D1.degree, D2.degree) == (2, 1) :
    if (len(D1.f) == 3) :
      ret = add_21_11(D1, D2)
    else :
      ret = add_22_11(D1, D2)
  elif (D1.degree, D2.degree) == (2, 2) :
    if   (len(D1.f) == 3) and (len(D2.f) == 3) :
      ret = add_21_21(D1, D2)
    elif (len(D1.f) == 3) and (len(D2.f) == 2) :
      ret = add_22_21(D2, D1)
    elif (len(D1.f) == 2) and (len(D2.f) == 3) :
      ret = add_22_21(D1, D2)
    else :
      ret = add_22_22(D1, D2)
  elif (D1.degree, D2.degree) == (3, 1) :
    ret = add_31_11(D1, D2)
  elif (D1.degree, D2.degree) == (3, 2) :
    if (len(D2.f) == 3) :
      ret = add_31_21(D1, D2)
    else :
      ret = add_31_22(D1, D2)
  elif (D1.degree, D2.degree) == (3, 3) :
    ret = add_31_31(D1, D2)
  else :
    raise NotImplementedError("Addition of divisors of degrees {} and {} not implemented.".format(D1.degree, D2.degree))
  
  if ret.reduced :
    return ret
  return flip(flip(ret))



def add_11_11(D1, D2):
  """
    Add two divisors, D1 and D2, each of degree 1.
    
    Divisors are assumed to be distinct, otherwise the doubling function should be used.
  """

  K = D1.K
  f, g, h = D1.f, D1.g, D1.h
  F, G, H = D2.f, D2.g, D2.h
  new_f, new_g, new_h = [], [], []

  # There are two cases. Either two points in D1 + D2 have distinct x-coordinates, or not.
  a1 = F[0] - f[0]
  alpha = F[0] - f[0]
  if a1 != 0 :
    # x-coordinates are distinct
    alpha = 1/a1
    u0 = -alpha*(g[0] - G[0])
    new_f = [-F[0]*u0 + G[0], -u0, K.one()]
    new_g = [f[0]*F[0], f[0] + F[0], K.zero(), K.one()]
    # Total : 1I 3M
  else :
    # x-coordinates are the same
    new_f = [f[0], K.one()]
    new_g = [g[0]*G[0], K(0), g[0] + G[0], K(0), K(0), K(1)]
    # Total : 0I 1M
  
  return C34CrvDiv(D1.C, [new_f, new_g, new_h])



def add_21_11(D1, D2) :
  K = D1.K
  f, g, h = D1.f, D1.g, D1.h
  F, G, H = D2.f, D2.g, D2.h
  new_f, new_g, new_h = [], [], []

  # Construct the matrix M
  # M = [ a1  a2  a3  a4 ]
  a1 = f[0] - G[0] - f[1]*F[0]
  a2 = g[0] - F[0]*(g[1] - F[0])
  # a3 = -F[0]*a1
  # a4 = -G[0]*a1
  
  if a1 != 0 :
    # Compute reduced row echelon form of M
    # Mrref = [ 1  -u0  -v0  -w0 ]
    alpha = 1 / a1
    u0 = -alpha*a2
    v0 = F[0]
    w0 = G[0]
    
    if u0 != 0 :
      # Sum will be typical
      new_f = [ u0*f[0] + g[0],
                u0*f[1] + g[1],
                u0,
                K.one() ]
      new_g = [ v0*f[0] - f[1]*new_f[0],
                f[0] + f[1]*(v0 - new_f[1]),
                v0 - f[1]*new_f[2],
                K.zero(), K.one() ]
      new_h = [ w0*f[0] - f[1]*new_g[0],
                f[1]*(w0 - new_g[1]),
                w0 + f[0] - f[1]*new_g[2],
                K.zero(), K.zero(), K.one() ]
      # D3 is of type 31, typical
      # Total 1I 13M 14A
    else :
      # Sum will be semi-typical
      new_f = [ g[0], g[1], K.zero(), K.one() ]
      new_g = [ v0*f[0] - f[1]*g[0],
                f[0] + f[1]*(v0 - g[1]),
                v0,
                K.zero(), K.one() ]
      new_h = [ w0*f[0] - f[1]*new_g[0],
                f[1]*(w0 - new_g[1]),
                w0 + f[0] - f[1]*new_g[2],
                K.zero(), K.zero(), K.one() ]
      # D3 is of type 31, non-typical
      # Total 1I 10M 12A
  elif a2 != 0 :
    # Sum is of type 32 with D1 + D2 = < f, g*F >
    new_f = [ f[0], f[1], K.one() ]
    new_g = [ g[0]*F[0], g[0] + g[1]*F[0], K.zero(), g[1] + F[0], K.zero(), K.zero(), K.one() ]
    # D3 is of type 32
    # Total 4M 6A
  else :
    # Divisors are non-disjoint.
    # So we have D1 = P1 + P2 and D2 = Q, where, wlog, Q = P1.
    # Need to determine whether P1 = P2.
    # We have P1 = P2 iff g[1] - 2F[0] = 0.
    
    if g[1] - F[0] - F[0] == 0 :
      # P1 = P2.
      # So D1 = 2Q and D1 + D2 = 3Q.
      return triple(D2)
    else :
      # P1 =/= P2
      # So D1 = P2 + Q and D1 + D2 = P2 + 2Q.
      # Find P2 and return P2 + (Q + Q)
      new_f = [ g[1] - F[0], K.one() ]
      new_g = [ f[0] - f[1]*new_f[0], K.zero(), K.one() ]
      E = C34CrvDiv(D1.C, [new_f, new_g, new_h])
      return add(E, double(D2))
  return C34CrvDiv(D1.C, [new_f, new_g, new_h])

def old_add_21_11(D1, D2) :
  # TODO : Is it faster to reduce modulo D2 rather than D1?
  K = D1.K
  f, g, h = D1.f, D1.g, D1.h
  F, G, H = D2.f, D2.g, D2.h
  new_f, new_g, new_h = [], [], []

  # Construct matrix M
  # M = [ a1  a2  a3  a4  a5   a6  ]
  #     [  1  a7  a8  a9  a10  a11 ]
  a1  = F[0]
  a2  = G[0] - f[0]
  a7  = -f[1]
  a3  =    - g[0]
  a8  = a1 - g[1]
  a4  =    - g[0]*a7
  a9  = a2 - g[1]*a7
  a5  = -f[0]*a2 +          f[1]*g[0]*a7
  a10 = -f[1]*a2 + (f[1]*g[1] - f[0])*a7
  a6  =    - g[0]*a8
  a11 = a3 - g[1]*a8
  # print Matrix(K, 2, 6, [a1,a2,a3,a4,a5,a6,1,a7,a8,a9,a10,a11])
  # print
  # Subtotal : 0I 10M

  # Reduce M to row echelon form
  # M_ref = [ 1  a7  a8  a9  a10  a11 ]
  #         [ 0  b1  b2  b3  b4   b5  ]
  b1 = a1*a7  - a2
  b2 = a1*a8  - a3
  b3 = a1*a9  - a4
  b4 = a1*a10 - a5
  b5 = a1*a11 - a6
  # print Matrix(K, 2, 6, [1, a7, a8, a9, a10, a11, 0, b1, b2, b3, b4, b5])
  # print
  # Subtotal : 1I 5M

  if b1 != 0 :
    # Reduce M to reduced row echelon form
    # M_rref = [ 1  0  -u0  -v0  -w0  -z0 ]
    #          [ 0  1  -u1  -v1  -w1  -z1 ]
    beta = 1/b1
    u1 = -beta*b2
    v1 = -beta*b3
    w1 = -beta*b4
    z1 = -beta*b5
    u0 = -(a8  + a7*u1)
    v0 = -(a9  + a7*v1)
    w0 = -(a10 + a7*w1)
    z0 = -(a11 + a7*z1)
    # print Matrix(K, 2, 6, [1, 0, -u0, -v0, -w0, -z0, 0, 1, -u1, -v1, -w1, -z1])
    # print
    # Subtotal : 1I 8M

    # Find polynomials forming an ideal generating set for D1 + D2 = <f'',g'',h''>
    # f'' = x^2 + ...
    # g'' = xy + ...
    # h'' = y^2 + ...
    new_f = [ F[0]*u0 + G[0]*u1, u0 + F[0], u1, K.one() ]
    new_g = [ F[0]*v0 + G[0]*v1, v0 + G[0], v1, K.zero(), K.one() ]
    new_h = [ F[0]*w0 + G[0]*w1, w0, w1 + G[0], K.zero(), K.zero(), K.one() ]
    # Subtotal : 0I 6M
    # Total : 1I 29M
    # Notes : 8 of these multiplications are used in computing h''
    #         If D + D' is typical, this can be done faster (in 7M?)
    #         5M are wasted computing z0, z1.
  elif b2 != 0 :
    # b1 = 0
    # Reduce M to reduced row echelon form
    # M_rref = [ 1  -u0  0  -v0  -w0  -z0 ]
    #          [ 0    0  1  -v1  -w1  -z1 ]
    beta = 1/b2
    v1 = -beta*b3
    w1 = -beta*b4
    z1 = -beta*b5
    u0 = -a7
    v0 = -(a9  + a8*v1)
    w0 = -(a10 + a8*w1)
    z0 = -(a11 + a8*z1)
    # print Matrix(K, 2, 6, [1, -u0, 0, -v0, -w0, -z0, 0, 0, 1, -v1, -w1, -z1])
    # print
    # Subtotal : 1I 6M

    # Find polynomials forming an ideal generating set for D1 + D2 = <f'',g''>
    # f'' = y + ...
    # g'' = x^3 + ...
    new_f = [ F[0]*u0 + G[0], u0, K.one() ]
    new_g = [ F[0]*z0, z0 + F[0]*z1, K.zero(), z1 + F[0], K.zero(), K.zero(), K.one() ]
    # Subtotal : 0I 3M
    # Total : 1I 24M
    # Note : If new_f = f, this saves 1M.
    #        11M are wasted computing v0, v1, w0, w1.
  else :
    # b1 = b2 = 0
    # Divisors are non-disjoint.
    # So we have D1 = P1 + P2 and D2 = Q, where, wlog, Q = P1.
    # Need to determine whether P1 = P2.
    # We have P1 = P2 iff g[1] - 2F[0] = 0.
    
    if g[1] - F[0] - F[0] == 0 : # TODO : Is this correct?
      # P1 = P2.
      # So D1 = 2Q and D1 + D2 = 3Q.
      return triple(D2)
    else :
      # P1 =/= P2
      # So D1 = P2 + Q and D1 + D2 = P2 + 2Q.
      # Find P2 and return P2 + (Q + Q)
      new_f = [ g[1] - F[0], K.one() ]
      new_g = [ f[0] - f[1]*new_f[0], K.zero(), K.one() ]
      E = C34CrvDiv(C, [new_f, new_g, new_h])
      return add(E, double(D2))

  return C34CrvDiv(D1.C, [new_f, new_g, new_h])



def add_21_21(D1, D2):
  # TODO: Doesn't work for case when D1 and D2 are non-disjoint and D1.f =/= D2.f
  K = D1.K
  f, g = D1.f, D1.g
  F, G = D2.f, D2.g
  new_f, new_g, new_h = [], [], []

  # Possible cases :
  #   Non-disjoint
  #     Have exactly one point in common. (Two common points handled by doubling.)
  #   Disjoint
  #     4 colinear points (principal divisor)
  #     3 colinear points
  # Compute matrix M
  # M = [ a1  a2  a3  a4  a5  ]
  #     [ a6  a7  a8  a9  a10 ]
  # A multiplication can be saved here
  # XXX: On inspection of adding divisors of the form (P + Q) + (Q + R),
  #      row 2 of M is all zero.
  a1 = F[0] - f[0]
  a2 = G[0] - g[0]
  a3 = -(F[1] - f[1])*g[0]
  a4 = (F[1] - f[1])*f[1]*g[0] - (F[0] - f[0])*f[0]
  a5 = -g[0]*(G[1] - g[1])               # Only needed when D1 + D2 is not typical
  a6 = F[1] - f[1]
  a7 = G[1] - g[1]
  a8 = F[0] - f[0] - (F[1] - f[1])*g[1]
  a9 = -(F[1] - f[1])*f[0] + (F[1] - f[1])*f[1]*g[1] - (F[0] - f[0])*f[1]
  a10 = G[0] - g[0] - g[1]*(G[1] - g[1]) # Only needed when D1 + D2 is not typical
  
  #print "M"
  #print Matrix(K, 2, 5, [a1, a2, a3, a4, a5, a6, a7, a8, a9, a10])
  #print

  # TODO: handle case where a1 = a2 = 0 (i.e. D1.f = D2.f)
  # If a1 = a2 = 0, points are all colinear.
  # If D1, D2 are disjoint, then D1 + D2 = <f>
  # If D1, D2 are non-disjoint, then there are a few cases.
  #   D1 + D2 = (P + Q) + (Q + R)
  #   D1 + D2 = (P + Q) + (Q + Q)
  #   D1 + D2 = (P + P) + (P + Q)
  
  # If a1 is zero, swap rows
  if (a1 == 0) :
    a1, a2, a3, a4, a5, a6, a7, a8, a9, a10 = a6, a7, a8, a9, a10, a1, a2, a3, a4, a5

  if (a1 != 0) :
    # Compute row echelon form of M
    # Mref = [ a1  a2  a3  a4  a5 ]
    #        [  0  b1  b2  b3  b4 ]
    # 
    b1 = a1*a7 - a2*a6
    b2 = a1*a8 - a3*a6
    b3 = a1*a9 - a4*a6
    b4 = a1*a10 - a5*a6 # Only needed when D1 + D2 is not typical
    #print "Mref"
    #print Matrix(K, 2, 5, [a1, a2, a3, a4, a5, 0, b1, b2, b3, b4])
    #print

    if (b1 != 0) :
      # Compute reduced row echelon form of M
      # Mrref = [ 1  0  -u0  -v0  -w1 ]
      #         [ 0  1  -u1  -v1  -w1 ]
      gamma = 1/(a1*b1)
      alpha = gamma*b1
      beta  = gamma*a1
      u1 = -beta*b2
      v1 = -beta*b3
      w1 = -beta*b4            # Only needed when D1 + D2 is not typical
      u0 = -alpha*(a3 + a2*u1)
      v0 = -alpha*(a4 + a2*v1)
      w0 = -alpha*(a5 + a2*w1) # Only needed when D1 + D2 is not typical
      
      new_f = [F[0]*u0 + G[0]*u1,
               F[1]*u0 + G[1]*u1 + F[0],
               u0,
               u1 + F[1],
               K.one()]
      new_g = [F[0]*v0 + G[0]*v1 - F[1]*new_f[0],
               F[1]*v0 + G[1]*v1 - F[1]*new_f[1],
               v0 + F[0] - F[1]*new_f[2],
               v1 - F[1]*new_f[3],
               K.zero(), K.one()]
      new_h = [F[0]*w0 + G[0]*w1,
               F[1]*w0 + G[1]*w1 + G[0],
               w0,
               w1 + G[1],
               K.zero(), K.zero(), K.one()]
    elif (b2 != 0) :
      # Compute reduced row echelon form of M
      # Mrref = [ 1  -u0  0  -v0  -w1 ]
      #         [ 0    0  1  -v1  -w1 ]
      gamma = 1/(a1*b2)
      alpha = gamma*b2
      beta  = gamma*a1
      v1 = -beta*b3
      w1 = -beta*b4            # Only needed when D1 + D2 is not typical
      u0 = -alpha*a2
      v0 = -alpha*(a4 + a3*v1)
      w0 = -alpha*(a5 + a3*w1) # Only needed when D1 + D2 is not typical
      z = F[1]*v1
      new_f = [F[0]*u0 + G[0],
               F[1]*u0 + G[1],
               u0,
               K.one()]
      new_g = [F[0]*v0           - z*new_f[0],
               F[1]*v0 + F[0]*v1 - z*new_f[1],
               v0 + F[0]         - z*new_f[2],
               K.zero(),
               v1 + F[1],
               K.one()]
      # D3 is type 43
    else :
      # If b1 and b2 are 0, then so are b3 and b4 and the divisors are non-disjoint.
      # We have D1 = P1 + P2, D2 = P1 + P3
      # (Possibly P1 = P2 or P1 = P3, but never P2 = P3.)
      r1 = -(F[0] - f[0])/(F[1] - f[1]) if g[1] == G[1] else -(G[0] - g[0])/(G[1] - g[1])
      r2 = -g[1] - r1
      r3 = -G[1] - r1
      s1 = -f[1]*r1 - f[0]
      s2 = -f[1]*r2 - f[0]
      s3 = -F[1]*r3 - F[0]
      P1 = C34CrvDiv(D1.C, [[-r1, K.one()], [-s1, K.zero(), K.one()], []])
      P2 = C34CrvDiv(D1.C, [[-r2, K.one()], [-s2, K.zero(), K.one()], []])
      P3 = C34CrvDiv(D1.C, [[-r3, K.one()], [-s3, K.zero(), K.one()], []])
      # If P1 = P2 or if P1 = P3, we need to triple and add (3*P1 + P3), resp. (3*P1 + P2)
      # Otherwise, we double and add (2*P1) + (P2 + P3)
      if (P1 == P2) :
        return add(triple(P1), P3)
      elif (P1 == P3) :
        return add(triple(P1), P2)
      else :
        return add(double(P1), add(P2, P3))
    """
    elif (b3 != 0) :
      # Compute reduced row echelon form of M
      # Mrref = [ 1  -u0  -v0  0  -w0 ]
      #         [ 0    0    0  1  -w1 ]
      gamma = 1/(a1*b3)
      alpha = alpha*b3
      beta  = alpha*a1
      w1 = -beta*b4            # Only needed when D1 + D2 is not typical
      u0 = -alpha*a2
      v0 = -alpha*a3
      w0 = -alpha*(a5 + a4*w1) # Only needed when D1 + D2 is not typical
      new_f = [F[0]*u0 + G[0],
               F[1]*u0 + G[1],
               u0,
               K.one()]
      new_g = [F[0]*v0,
               F[1]*v0 + F[0],
               v0,
               F[1], # <-- should be 0
               K.one()]
      # new_h = [F[0]*w0,
      #          F[1]*w0 + G[0],
      #          w0 + F[0]*w1,
      #          G[1],
      #          F[1]*w1, # <-- should be 0
      #          w1,
      #          K.one()]
      # XXX: Form <x^2, xy, x^3>. new_h is redundant
    elif (b4 != 0) :
      # Compute reduced row echelon form of M
      # Mrref = [ 1  -u0  -v0  -w0  0 ]
      #         [ 0    0    0    0  1 ]
      # XXX : Is this case possible?
      alpha = 1/a1
      u0 = -alpha*a2
      v0 = -alpha*a3
      w0 = -alpha*a4 # Only needed when D1 + D2 is not typical
      new_f = [F[0]*u0 + G[0],
               F[1]*u0 + G[1],
               u0,
               K.one()]
      new_g = [F[0]*v0,
               F[1]*v0 + F[0],
               v0,
               F[1], # <-- should be 0
               K.one()]
      # new_h = [F[0]*w0,
      #          F[1]*w0,
      #          w0 + F[0],
      #          K.zero(),
      #          F[1], # <-- should be 0
      #          K.one()]
      # XXX: Form <x^2, xy, y^2>
      # XXX: new_h is redundant??? new_h = y*F + w0F
      raise ValueError("Kernel has unexpected form. Ker(M) =\n{}".format(Matrix(K, 5, 3, [u0, v0, w0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0])))
    else :
      raise NotImplementedError("Matrix rows are linearly dependent! Mref = \n{}".format(Matrix(K, 2, 5, [a1, a2, a3, a4, a5, 0, b1, b2, b3, b4])))
    """

  else : # if a1 == 0
    # The f = F and M has the form
    # M = [ 0  a2  0  0  a5  ]
    #     [ 0  a7  0  0  a10 ]
    # Since g != G, one of a2 or a7 is non-zero
    # Swap rows, if necessary
    if (a2 == 0) :
      a2, a5, a7, a10 = a7, a10, a2, a7
    
    # Compute row echelon form of M
    # Mref = [ 0 a2 0 0 a5 ]
    #        [ 0  0 0 0 b1 ]
    b1 = a2*a10 - a5*a7
    
    if (b1 != 0) :
      # Then D1 and D2 are disjoint, and D1 + D2 = <f>
      new_f = [f[0], f[1], K.one()]
    else :
      # D1 and D2 are non-disjoint
      r1 = -(G[0] - g[0])/(G[1] - g[1])
      r2 = -g[1] - r1
      r3 = -G[1] - r1
      s1 = -f[1]*r1 - f[0]
      s2 = -f[1]*r2 - f[0]
      s3 = -F[1]*r3 - F[0]
      P1 = C34CrvDiv(D1.C, [[-r1, K.one()], [-s1, K.zero(), K.one()], []])
      P2 = C34CrvDiv(D1.C, [[-r2, K.one()], [-s2, K.zero(), K.one()], []])
      P3 = C34CrvDiv(D1.C, [[-r3, K.one()], [-s3, K.zero(), K.one()], []])

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
  return C34CrvDiv(D1.C, [new_f, new_g, new_h])



def add_22_11(D1, D2):
  K = D1.K
  f, g, h = D1.f, D1.g, D1.h
  F, G, H = D2.f, D2.g, D2.h
  new_f, new_g, new_h = [], [], []
  
  # Construct matrix M
  # M = [ a1  a2  a3  a4  a5 ]
  #     [  0   1   0  a6  a7 ]
  a1 = F[0] - f[0]
  a7 = G[0] - g[2]
  # Subtotal : 0I 0M
  # Most of the entries of M aren't needed, since M_rref is so simple.
  
  # Reduce matrix M
  # M_rref = [ 1  0  -f0    0  -w0 ]
  #          [ 0  1    0  -f0  -w1 ]
  if a1 == 0 : # The three points are colinear
    b1 = -g[0] - G[0]*a7
    if b1 == 0 : # D1 and D2 are non-disjoint
      if (g[0] == G[0]*G[0]) and (g[2] == G[0] + G[0]) : # if g = G^2
        # D1 = 2*D2
        return triple(D2)
      else :
        # D1 + D2 = (P + Q) + (P)
        # where P = (F[0], G[0])
        #   and Q = (F[0], g[2] - G[0])
        # Using commutativity/associativity, compute instead (Q) + (2P)
        return add(C34CrvDiv(D1.C, [[f[0], K.one()], [g[2] - G[0], K.zero(), K.one()], []]), double(D2))
    else : # D1 and D2 are disjoint. D1 + D2 = <f> is principal
      return C34CrvDiv(D1.C, [[f[0], K.one()], [], []])
  alpha = 1/a1
  w0 = alpha*(g[0] + G[0]*a7)
  w1 = -a7
  # Subtotal : 1I 2M

  # Find polynomials forming an ideal generating set for D1 + D2 = <f'',g'',h''>
  # f'' = x^2 + ...
  # g'' = xy + ...
  # h'' = y^2 + ...
  new_f = [ F[0]*f[0], f[0] + F[0], K.zero(), K.one() ]
  new_g = [ G[0]*f[0], G[0], f[0], K.zero(), K.one() ]
  new_h = [ F[0]*w0 + G[0]*w1, w0, w1 + G[0], K.zero(), K.zero(), K.one() ]
  # Subtotal : 0I 4M
  # Total : 1I 6M
  # Notes : Can we save a multiplication with Karatsuba's technique?
  return C34CrvDiv(D1.C, [new_f, new_g, new_h])



def add_22_21(D1, D2):
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



def add_31_22(D1, D2):
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

  def gcd_31_31(D1, D2) :
    a1 = D1.f[0] - D2.f[0]
    a2 = D1.g[0] - D2.g[0]
    a3 = D1.h[0] - D2.h[0]
    a4 = D1.f[1] - D2.f[1]
    a5 = D1.g[1] - D2.g[1]
    a6 = D1.h[1] - D2.h[1]
    a7 = D1.f[2] - D2.f[2]
    a8 = D1.g[2] - D2.g[2]
    a9 = D1.h[2] - D2.h[2]
    
    # Given the matrix
    #
    # [ a1  a2  a3 ] 
    # [ a4  a5  a6 ]
    # [ a7  a8  a9 ]
    #
    # Perform column operations to reduce it to one of the following forms
    #
    # [ 1  0  0 ]  [ 0  *  * ]  [ 0  0  * ]     [ 0  0  * ]
    # [ 0  1  0 ], [ 0  1  0 ], [ 0  0  * ], or [ 0  0  1 ]
    # [ 0  0  1 ]  [ 0  0  1 ]  [ 0  0  1 ]     [ 0  0  0 ]
    #
    # In the first case, gcd(D1, D2) is 0
    # In the second, third, and fourth, gcd(D1, D2) is type 11, 21, and 22, respectively.
    
    if (a9 != 0) or (a8 != 0) or (a7 != 0) :
      if (a9 == 0) :
        if (a8 != 0) :
          a2, a5, a8, a3, a6, a9 = a3, a6, a9, a2, a5, a8
        else :
          a1, a4, a7, a3, a6, a9 = a3, a6, a9, a1, a4, a7
      # Compute
      # [ b1  b2  a3 ] 
      # [ b3  b4  a6 ]
      # [ 0   0   a9 ]
      b1 = a1*a9 - a3*a7
      b2 = a2*a9 - a3*a8
      b3 = a4*a9 - a6*a7
      b4 = a5*a9 - a6*a8

      if (b3 != 0) or (b4 != 0) :
        if (b4 == 0) :
          b1, b3, b2, b4 = b2, b4, b1, b3
        # Compute
        # [ c1  b2  a3 ] 
        # [ 0   b4  a6 ]
        # [ 0   0   a9 ]
        c1 = b1*b4 - b2*b3

        if (c1 != 0) :
          # gcd(D1, D2) is zero
          return D1.C.zero_divisor()
        else :
          # gcd(D1, D2) is type 11
          gamma = 1/(b4*a9)
          alpha = gamma*b4
          beta = gamma*a9
          u0 = beta*b2
          v0 = alpha*(a3 - a6*u0)
          return C34CrvDiv(D1.C, [[u0, 1], [v0, 0, 1], []])
      else :
        # gcd(D1, D2) is type 21
        alpha = 1/a9
        u1 = alpha*a6
        u0 = alpha*a3
        v1 = D1.f[1] - D1.f[2]*u1
        v0 = D1.f[0] - D1.f[2]*u0
        return C34CrvDiv(D1.C, [[u0, u1, 1], [v0, v1, 0, 1], []])        
    else :
      # gcd(D1, D2) is type 22
      if (a6 == 0) :
        if (a5 != 0) :
          a2, a5, a3, a6 = a3, a6, a2, a5
        else :
          a1, a4, a3, a6 = a3, a6, a1, a4
      alpha = 1/a6
      u0 = alpha*a3
      v0 = D1.h[0] - D1.h[1]*u0
      v2 = D1.h[2]
      return C34CrvDiv(D1.C, [[u0, 1], [v0, 0, v2, 0, 0, 1], []])        

  # Compute M
  #
  #     [ a1   a2   a3   a4   a5   a6   a7  ]
  # M = [ a8   a9   a10  a11  a12  a13  a14 ]
  #     [ a15  a16  a17  a18  a19  a20  a21 ]
  #
  # Columns 4, 5, and 6 are computed by
  #
  # [ a4   a5   a6  ]   [ 0  -G[0]  -G[0] ] [ a1   a2   a3  ] 
  # [ a11  a12  a13 ] = [ 1  -G[1]  -G[1] ]*[ a8   a9   a10 ]
  # [ a18  a19  a20 ]   [ 0  -G[2]  -G[2] ] [ a15  a16  a17 ]
  #
  # The last column similarly by
  #
  # [ a7  ]   [ 0  -G[0]  -G[0] ] [ a4  ] 
  # [ a14 ] = [ 1  -G[1]  -G[1] ]*[ a11 ]
  # [ a21 ]   [ 0  -G[2]  -G[2] ] [ a18 ]
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

  #print "M = "
  #print Matrix([
  #  [a1,  a2,  a3,  a4,  a5,  a6],
  #  [a8,  a9,  a10, a11, a12, a13],
  #  [a15, a16, a17, a18, a19, a20]])
  #print
  
  if (a1 != 0 or a8 != 0 or a15 != 0) :
    # Swap rows, if necessary, to ensure a1 != 0
    aswap = 0
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
        G = gcd_31_31(D1, D2)
        return add(flip(flip(L)), G)
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
        # If c1, c2, c3 are zero, then c4 is zero and D1, D2 are non-disjoint.
        assert c4 == 0
        raise NotImplementedError("Divisors are non-disjoint.")

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
      u3 =      r0
      v0 = f[0]*s0 + h[0]*s1
      v1 = f[1]*s0 + h[1]*s1 + f[0]*s2
      v2 = f[2]*s0 + h[2]*s1
      v3 =      s0           + f[1]*s2 + f[0]
      v4 =                     f[2]*s2
      v5 =                s1
      v6 =                          s2 + f[1]
      v7 =                               f[2]
      
      # g'' is x^4 + v7*x^2*y + v6*x^3 + v5*y^2 + v4*x*y + v3*x^2 + v2*y + v1*x + v0
      # Kill off monomials divisible by xy by subtracting appropriate multiples of f'' = xy + u3*x^2 + u2*y + u1*x + u0.
      z = v4 - v7*u2
      v0 = v0         - z*u0
      v1 = v1 - v7*u0 - z*u1
      v2 = v2         - z*u2
      v3 = v3 - v7*u1 - z*u3
      #v4 = 0
      #v5 = v5
      v6 = v6 - v7*u3
      #v7 = 0
      new_f = [u0, u1, u2, u3, 1]
      new_g = [v0, v1, v2, v3, 0, v5, v6, 0, 0, 1]
    else :
      # If b1, b2, b7, b8 are all zero, then D1, D2 are non-disjoint.
      raise NotImplementedError("Divisors are non-disjoint.")
  elif (a2 != 0 or a9 != 0 or a16 != 0) :
    if (a2 == 0) :
      if (a9 != 0) :
        a2, a3, a5, a6, a9, a10, a12, a13 = a9, a10, a12, a13, a2, a3, a5, a6
      else :
        a2, a3, a5, a6, a16, a17, a19, a20 = a16, a17, a19, a20, a2, a3, a5, a6

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
      
      assert c1 == 0
      if (c2 != 0) :
        new_f = [f[0], f[1], f[2], 1]
      else :
        raise NotImplementedError("Divisors are non-disjoint.")
    else :
      raise NotImplementedError("Divisors are non-disjoint.")
  else :
    # If first two columns are zero, then D1, D2 are non-disjoint.
    raise NotImplementedError("Divisors are non-disjoint.")

  
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

