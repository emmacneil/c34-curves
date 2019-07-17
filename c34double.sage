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
    
    Input: A typical C34CrvDiv D.
    Output: The C34CrvDiv E equivalent to D + D.
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



def double_0(D):
  return C34CrvDiv(D.C, [[D.K.one()], [], []])



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
  
  return C34CrvDiv(D.C, [new_f, new_g, new_h])



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
    return C34CrvDiv(D.C, [[f[0], f[1], K.one()], [], []])
    #return C34CrvDiv(D.C, [[K.one()], [], []])
  
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
    # print Matrix(K, 2, 5, [-1, u0, 0, v0, w0, 0, 0, -1, v1, w1])
    # print
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
    # print Matrix(K, 2, 5, [-1, u0, v0, 0, w0, 0, 0, 0, -1, w1])
    # print
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
  return C34CrvDiv(D.C, [new_f, new_g, new_h])



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
  return C34CrvDiv(D.C, [new_f, new_g, new_h])



def double_31(D):
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
      return C34CrvDiv(D.C, [[f[0], f[1], K.one()], [], []])
      #return C34CrvDiv(D.C, [[K.one()], [], []])
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

      return C34CrvDiv(D.C, [[u0, u1, u2, u3, 1], [v0, v1, v2, v3, 0, v5, v6, 0, 0, 1], []])
    
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
    A = C34CrvDiv(D.C, [fA, gA, []])
    B = C34CrvDiv(D.C, [fB, gB, []])
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

  return C34CrvDiv(D.C, [new_f, new_g, new_h])

