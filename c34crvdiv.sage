"""
  Divisor of a C34 curve.

  Divisors are represented by up to 3 monic polynomials, f, g, h.
  Internally, these polynomials are represented by lists.
  The i'th element of a list represents the coefficient of the i'th monomial in
  
    [ 1, x, y, x*x, x*y, y*y, x*x*x, x*x*y, x*y*y, x*x*x*x ]
  
  The last element of each list is 1, and the length of the list determines the leading term of the polynomial.
  For example, the list

    [ 4, 3, 2, 1]
  
  represents the polynomial
  
    x^2 + 2y + 3x + 4

  Fields
    C       : The curve on which this divisor is defined.
    K       : The field over which this divisor (and its curve) is defined.
    R       : The polynomial ring containing this divisor's curve's defining polynomial.
              This is K[x,y].
    type    : A string, either "typical" or "semi-typical"
    f,g,h   : Monic polynomials whose vanishing set is the support of this divisor, counting
              multiplicities. Should be of one of the two forms above, for typical or semi-typical
              divisors.
    degree  : The (effective) degree of the divisor.
    reduced : True if the divisor is a reduced divisor. Otherwise False.
"""
class C34CurveDivisor :
  def __init__(self, C, lst, degree = -1, typ = -1, reduced = 2, typical = 2, inv = 0) :
    """
      If divisor is supported by three non-colinear points
      
        P(x1, y1), Q(x2, y2), R(x3, y3)
      
      it is represented by three polynomials
      
        F : x*x + a*y * b*x + c
        G : x*y + d*y + e*x + f
        H : y*y + r*y + s*x + t
      
      where
      
        a = alpha*(x1*x2*(x1 - x2) + x1*x3*(x3 - x1) + x2*x3*(x2 - x3))
        b = alpha*((x2^2 - x3^2)*y1 + (x3^2 - x1^2)*y2 + (x1^2 - x2^2)*y3)
        c = alpha*(x2*x3*(x3 - x2)*y1 + x1*x3*(x1 - x3)*y2 + x1*x2*(x2 - x1)*y3)

        d = alpha*(x1*(x2 - x3)*y1 + x2*(x3 - x1)*y2 + x3*(x1 - x2)*y3)
        e = alpha*((x2 - x1)*y1*y2 + (x1 - x3)*y1*y3 + (x3 - x2)*y2*y3)
        f = alpha*(x1*(x2 - x3)*y2*y3 + x2*(x3 - x1)*y1*y3 + x3*(x1 - x2)*y1*y2)

        r = alpha*((x2 - x3)*y1^2 + (x3 - x1)*y2^2 + (x1 - x2)*y3^2)
        s = alpha*((y2 - y1)*y1*y2 + (y1 - y3)*y1*y3 + (y3 - y2)*y2*y3)
        t = alpha*(x1*(y2 - y3)*y2*y3 + x2*(y3 - y1)*y1*y3 + x3*(y1 - y2)*y1*y2)
      
      and
      
        alpha = ((x3 - x2)*y1 + (x1 - x3)*y2 + (x2 - x1)*y3)^(-1)
    """
    self.C = C
    self.R = C.R
    self.K = C.K
    self.f = []
    self.g = []
    self.h = []
    self.inv = inv
    self.degree = degree
    self.reduced = reduced
    self.typical = typical
    self.type = typ

    # Make sure 'lst' is a list
    if not isinstance(lst, list) :
      raise TypeError("'lst' must be of type 'list'.")
    
    # Makes sure elements in 'lst' are of same type
    if not [type(x) for x in lst] == [type(x) for x in [lst[0]]*len(lst)] :
      raise TypeError("Elements in 'lst' must be of same type.")

    # Check type of elements in 'lst'.
    # Should be one of :
    #   * C34CurvePoints
    #   * Polynomials in self.R
    #   * List of field elements in self.K (representing polynomial coefficients
    if lst[0] in self.R : # If 'lst' elements are polynomials
      polys = copy(lst)
      polys.sort()
      if len(polys) > 3 :
        raise ValueError("Divisor must be given by 3 or fewer polynomials.")
      
      f = polys[0] if len(polys) > 0 else 0
      g = polys[1] if len(polys) > 1 else 0
      h = polys[2] if len(polys) > 2 else 0
      x, y = self.R.gens()
      mmap = {1 : 0, x : 1, y : 2, x*x : 3, x*y : 4, y*y : 5, x*x*x : 6, x*x*y : 7, x*y*y : 8, x*x*x*x : 9}
      if f != 0 :
        LMf = f.monomials()[0]       # Find leading monomial of f
        self.f = [0]*(mmap[LMf] + 1) # Determine length of list representing f
        for m in f.monomials() :     # Construct list representing f
          self.f[mmap[m]] = f.monomial_coefficient(m)
      if g != 0 :
        LMg = g.monomials()[0]
        self.g = [0]*(mmap[LMg] + 1)
        for m in g.monomials() :
          self.g[mmap[m]] = g.monomial_coefficient(m)
      if h != 0 :
        LMh = h.monomials()[0]
        self.h = [0]*(mmap[LMh] + 1)
        for m in h.monomials() :
          self.h[mmap[m]] = h.monomial_coefficient(m)
      
    elif isinstance(lst[0], C34CurvePoint) : # If 'lst' elements are points
      points = lst
      if len(points) == 1 :
        # XXX : Assumes point is not at infinity
        P = points[0]
        if (P[0] not in self.K) or (P[1] not in self.K) :
          raise TypeError("Point must be in curve's base field.")
        x1, y1 = self.K(P[0]), self.K(P[1])
        self.f = [-x1, self.K(1)]
        self.g = [-y1, self.K(0), self.K(1)]
        self.h = []
      else :
        D = C34CurveDivisor(C, [points[0]]) + C34CurveDivisor(C, points[1:])
        self.f, self.g, self.h = D.f, D.g, D.h

    elif isinstance(lst, list) :
      if len(lst) > 0 :
        self.f = copy(lst[0])
        # self.f = [self.K(t) for t in lst[0]]
      if len(lst) > 1 :
        self.g = copy(lst[1])
        # self.g = [self.K(t) for t in lst[1]]
      if len(lst) > 2 :
        self.h = copy(lst[2])
        # self.h = [self.K(t) for t in lst[2]]

    if (self.degree < 0) or (self.type < 0) or (self.typical == 2) or (self.reduced == 2) :
      self.classify()

  
  
  def classify(self) :      
    """
      Determines the type and degree of this divisor, whether it is typical and/or reduced,
      and sets its fields accordingly.

      This determination is done solely on the form of the polynomials generating its ideal.
      If the polynomials have the 'form' of a type 31 divisor, this method does not verify that
      they intersect the parent curve in precisely three places.

      Raises a RuntimeError if the the form of the divisor's generators do not correspond to
      a reduced Groebner basis of any divisor of degree 6 or less.
    """
    L = (len(self.f), len(self.g), len(self.h))
    if   ( L == (1, 0, 0) ) :
      self.type = 0
      self.degree = 0
      self.reduced = True
      self.typical = False
    elif ( L == (2, 3, 0) ) :
      self.type = 11
      self.degree = 1
      self.reduced = True
      self.typical = False
    elif ( L == (3, 4, 0) ) :
      self.type = 21
      self.degree = 2
      self.reduced = True
      self.typical = False
    elif ( L == (2, 6, 0) ) :
      self.type = 22
      self.degree = 2
      self.reduced = True
      self.typical = False
    elif ( L == (4, 5, 6) ) :
      self.type = 31
      self.degree = 3
      self.reduced = True
      if self.f[2] != 0 :
        self.typical = True
      else :
        self.typical = False
    elif ( L == (3, 7, 0) ) :
      self.type = 32
      self.degree = 3
      self.reduced = False
      self.typical = False
    elif ( L == (2, 0, 0) ) :
      self.type = 33
      self.degree = 3
      self.reduced = False
      self.typical = False
    elif ( L == (5, 6, 7) ) :
      self.type = 41
      self.degree = 4
      self.reduced = False
      if self.f[3]^2 + self.g[3] != 0 :
        self.typical = True
      else :
        self.typical = False
    elif ( L == (4, 5, 0) ) :
      self.type = 42
      self.degree = 4
      self.reduced = False
      self.typical = False
    elif ( L == (4, 6, 0) ) :
      self.type = 43
      self.degree = 4
      self.reduced = False
      self.typical = False
    elif ( L == (3, 0, 0) ) :
      self.type = 44
      self.degree = 4
      self.reduced = False
      self.typical = False
    elif ( L == (6, 7, 8) ) :
      self.type = 51
      self.degree = 5
      self.reduced = False
      if self.f[4]*(self.C.c[8] - self.f[4]) + self.f[3] - self.C.c[7] + self.g[4] != 0 :
        self.typical = True
      else :
        self.typical = False
    elif ( L == (5, 6, 0) ) :
      self.type = 52
      self.degree = 5
      self.reduced = False
      self.typical = False
    elif ( L == (5, 7, 0) ) :
      self.type = 53
      self.degree = 5
      self.reduced = False
      self.typical = False
    elif ( L == (4, 9, 0) ) :
      self.type = 54
      self.degree = 5
      self.reduced = False
      self.typical = False
    elif ( L == (7, 8, 9) ) :
      self.type = 61
      self.degree = 6
      self.reduced = False
      if self.f[4] - self.g[5] + self.f[5]*(self.f[5] - self.C.c[8]) != 0 :
        self.typical = True
      else :
        self.typical = False
    elif ( L == (6, 7, 0) ) :
      self.type = 62
      self.degree = 6
      self.reduced = False
      self.typical = False
    elif ( L == (6, 8, 0) ) :
      self.type = 63
      self.degree = 6
      self.reduced = False
      self.typical = False
    elif ( L == (5, 10, 0) ) :
      self.type = 64
      self.degree = 6
      self.reduced = False
      self.typical = False
    elif ( L == (4, 0, 0) ) :
      self.type = 65
      self.degree = 6
      self.reduced = False
      self.typical = False
    else : 
      raise RuntimeError("Could not classify divisor.\nD = {}".format(self))
  
  
  
  def extension(self) :
    """
      Return the degree of the smallest extension L of K such that every point in D is K-rational.
    """
    ret = 1
    for p in self.points() :
      d = parent(p[0]).degree()
      ret = d if d > ret else ret
    return ret



  """
    Returns the formal sum reprentation of point of the divisor as a list of pairs.
    
    If D is the divisor (1 : 1 : 1) + 2*(2 : 2 : 1) + (3 : 3 : 1), then this method returns the
    list
    
      [ ((1 : 1 : 1), 1), ((2 : 2 : 1), 2), ((3 : 3: 1), 1)]
    
    The coordinates of the point may come from an extension of the divisor's base field.
    
    The algorithm proceeds as follows. Let I be the K[x,y]-ideal representing D. Compute the
    primary decomposition of I. This gives I as an intersection of a family of primary ideals Q_i.
    
      I = Q_1 cap Q_2 cap ... cap Q_n
    
    For every Q_i that is a power of a prime K[x,y]-ideal P_i of the form <x - a_i, y - b_i>,
    compute r_i such that (P_i)^(r_i) = Q_i and add r_i*(a_i : b_i : 1) to the formal sum. Then
    divide I out by these Q_i and perform the primary decomposition of the remainder, but viewed
    as a L[x,y]-ideal for an extension L of K. Doing this over higher and higher extensions
    eventually gives I as a product of powers of prime ideals

      I = (P_1)^(r_1) * ... * (P_n)^(r_n).
    
    This is not particularly fast.
  """
  def formal_sum(self) :
    F = self.C.defining_polynomial()
    R = self.R
    K = self.K
    I = self.ideal()
    n = 1
    ret = []
    max_ext = self.degree
    while not I.is_one() :
      assert n <= max_ext
      L = K.extension(n)
      S = PolynomialRing(L, 2, R.variable_names(), order = R.term_order())
      x, y = S.gens()
      J = S.ideal(1)
      for Q, P in S.ideal(I).complete_primary_decomposition() :
        # If P represents a type 11 divisor
        if (P.ngens() == 2) :
          gens = list(P.gens())
          gens.sort()
          f, g = gens
          if (f.lm() == x) and (g.lm() == y) :
            # P represents a type 11 divisor.
            point = self.C.point(-f.constant_coefficient(), -g.constant_coefficient())
            r = 1
            # Find r such that P^r = Q
            while (P^r + S(F)) != (Q + S(F)) :
              r = r + 1
            # Add r*P to the sum
            ret = ret + [(point, r)]
            J = S.ideal(J.intersection(Q).groebner_basis())
      J = R.ideal(J)
      I = I.quotient(J)
      n = n + 1
    ret.sort()
    return ret


  
  def ideal(self) :
    I = self.R.ideal(self.polys() + [self.C.defining_polynomial()])
    G = list(I.groebner_basis())
    G.sort()
    return self.R.ideal(G)

  """
    Returns true if self is equivalent to other in the divisor class group.
    
    Two divisors are equivalent if and only if their flips are equal.
  """
  def is_equivalent_to(self, other) :
    return (- self) == (- other)
  
  

  """
    Returns true if the ideal of D is squarefree.
    
    Equivalently, returns true if all points in the formal sum of D are distinct.
    Otherwise, returns false.
  """
  def is_squarefree(self) :
    I = self.ideal()
    return I.radical() == I

  def matrix(self, polys) :
    """
      Returns the matrix of the given polynomials modulo I_D.
      The divisor must be reduced, otherwise a ValueError is raised.

      Input : A list of polynomials [p_1, p_2, ..., p_n]
      Output : The matrix [ C_1  C_2  ... C_n ] were the i-th column, C_i, is the reduction of the
               i-th polynomial, p_i, modulo I_D.
    """
    # Input : D, a type 31 divisor
    #         polys, a list of polynomials in D.R
    #
    # Return the matrix M = [ C1  C2 ... Cn ]
    # where column i is the reduction of polys[i] modulo the ideal of D
    if (not self.reduced) :
      raise ValueError("Divisor is not reduced.")

    Q = self.R.quotient(self.ideal())
    x, y = self.R.gens()
    ret = Matrix(self.K, 3, len(polys))
    for i in range(len(polys)) : 
      f  = Q(polys[i]).lift()
      ret[0,i] = f.subs(x=0, y=0)
      ret[1,i] = f.subs(x=1, y=0) - ret[0,i]
      ret[2,i] = f.subs(x=0, y=1) - ret[0,i]
    return ret



  def order_at_point(self, point) :
    # Assumes point is affine
    ret = -1
    K = point.base_field()
    AA = self.C.ambient_space().base_extend(K)
    X = self.C.scheme().base_extend(K)
    P = AA(point[0], point[1])
    for f in self.polys() :
      Y = AA.subscheme(f)
      n = X.intersection_multiplicity(Y, P)
      if (ret == -1) or (n < ret) :
        ret = n
    """
    ret = -1
    AA.<x,y> = AffineSpace(2, self.K)
    P = AA(point[0:2])
    X = self.C.scheme()
    for f in self.polys() :
      Y = AA.subscheme(f)
      n = X.intersection_multiplicity(Y, P)
      if (ret == -1) or (n < ret) :
        ret = n
    """
    return ret



  # Return the points in the divisor
  # TODO : parse and return C34CurvePoints
  def points(self) :
    x, y = self.R.gens()
    V = self.variety()
    
    # Convert variety to list of points
    pts = [(t[x].as_finite_field_element(minimal=True)[1], t[y].as_finite_field_element(minimal=True)[1]) for t in V]
    pts.sort()
    pts = [self.C.point(p[0], p[1]) for p in pts]
    #pts = [(t[x], t[y]) for t in V]
    #pts.sort()
    #pts = [self.C.point(p[0], p[1]) for p in pts]
    return pts



  # Return f, g, and h in a list, in the form of polynomials rather than lists.
  def polys(self) :
    x, y = self.R.gens()
    m = [self.R.one(), x, y, x*x, x*y, y*y, x*x*x, x*x*y, x*y*y, x*x*x*x, y*y*y]
    m.sort()
    f = [ self.f[i] * m[i] for i in range(len(self.f))] #+ m[len(self.f)]
    g = [ self.g[i] * m[i] for i in range(len(self.g))] #+ m[len(self.g)]
    h = [ self.h[i] * m[i] for i in range(len(self.h))] #+ m[len(self.h)]
    f = sum(f)
    g = sum(g)
    h = sum(h)
    ret = [f, g, h]
    ret = filter(lambda l : l != 0, ret)
    ret.sort()
    return ret
  
  
  
  """
    Returns the sum of self and other.
    
    More precisely, it returns a reduced divisor that is equivalent in the divisor class group to
    the composition of self with other.

    This method is "slow" in the sense that it uses Sage's provided ideal arithmetic rather than
    explicit formulas.
  """
  def slow_add(self, other) :
    return self.slow_compose(other).slow_reduce()
  
  
  """
    Returns the composition of self with other.
    
    The composition of two divisors D1 and D2 is the divisor whose formal sum of points is simply
    the sum of D1 and D2's formal sums of points. E.g., if
    
      D1 = P + Q + R
      D2 = S + T + U
    
    then D3 = D1.slow_compose(D2) is the divisor
    
      D3 = P + Q + R + S + T + U.
    
    If D1 and D2 are non-disjoint, say
    
      D1 = P + Q
      D2 = Q + R,
    
    then D3 = D1.slow_compose(D2) is
    
      D3 = P + 2*Q + R.
    
    This method is "slow" in the sense that it uses Sage's provided ideal arithmetic rather than
    explicit formulas.
  """
  def slow_compose(self, other) :
    C = self.C
    R = self.R
    F = C.defining_polynomial()
    x, y = R.gens()
    
    I = self.ideal() * other.ideal() + R.ideal(F)
    J = R.ideal(I.groebner_basis())
    polys = J.gens()[0:]
    f = polys[0]
    if f.lm() == y^3 :
      polys = polys[1:]
      
    return C34CurveDivisor(C, polys)

  
  
  """
    Returns the flip of this divisor.
    
    The flip of a divisor D is a reduced divisor that is equivalent in the divisor class group 
    to -D.

    This method is "slow" in the sense that it uses Sage's provided ideal arithmetic rather than
    explicit formulas.
  """
  def slow_flip(self) :
    C = self.C
    R = self.R
    F = C.defining_polynomial()
    x, y = R.gens()
    
    I = self.ideal()
    f = self.polys()[0]
    J = R.ideal(f, F)
    Q = R.ideal(J.quotient(I).groebner_basis())
    polys = Q.gens()[0:]
    f = polys[0]
    if f.lm() == y^3 :
      polys = polys[1:]
    
    return C34CurveDivisor(C, polys)



  """
    Returns the greatest common divisor of self and other.
    
    The gcd of two divisors D1 and D2 is the sum of points that are common to both D1 and D2.
    If the order of D1 at a point P is m and the order of D2 at P is n, then the order of
    gcd(D1, D2) at P is min(m, n). If D1 and D2 are disjoint, then this results in the zero
    divisor. Some examples of gcd's of non-disjoint divisors:
    
      D1 = P + Q + R
      D2 = Q + R + S
      gcd(D1, D2) = Q + R
      
      D1 = P + Q + 2*R
      D2 = Q + R + S
      gcd(D1, D2) = Q + R
      
      D1 = P + Q + 2*R
      D2 = Q + 2*R + S
      gcd(D1, D2) = Q + 2*R

      D1 = 3*P
      D2 = 2*P
      gcd(D1, D2) = 2*P
    
    At the level of ideals, the gcd of two divisors is analogous to the sum of two ideals.

    This method is "slow" in the sense that it uses Sage's provided ideal arithmetic rather than
    explicit formulas.
  """
  def slow_gcd(self, other) :
    C = self.C
    R = self.R
    F = C.defining_polynomial()
    x, y = R.gens()
    
    I = self.ideal() + other.ideal() + R.ideal(F)
    J = R.ideal(I.groebner_basis())
    polys = J.gens()[0:]
    f = polys[0]
    if f.lm() == y^3 :
      polys = polys[1:]
      
    return C34CurveDivisor(C, polys)



  """
    Returns the least common divisor of self and other.
    
    The lcm of two divisors D1 and D2 is the sum of points found in one of either D1 or D2,
    counting multiplicities. If the order of D1 at a point P is m and the order of D2 at P is n,
    then the order of lcm(D1, D2) at P is max(m, n). If D1 and D2 are disjoint, then this the same
    as the composition of the two divisors.
        
      D1 = P + Q + R
      D2 = Q + R + S
      gcd(D1, D2) = P + Q + R + S
      
      D1 = P + Q + 2*R
      D2 = Q + R + S
      gcd(D1, D2) = P + Q + 2*R + S
      
      D1 = P + Q + 2*R
      D2 = Q + 2*R + S
      gcd(D1, D2) = P + Q + 2*R + S

      D1 = 3*P
      D2 = 2*P
      gcd(D1, D2) = 3*P
    
    At the level of ideals, the gcd of two divisors is analogous to the intersection of two ideals.

    This method is "slow" in the sense that it uses Sage's provided ideal arithmetic rather than
    explicit formulas.
  """
  def slow_lcm(self, other) :
    C = self.C
    R = self.R
    F = C.defining_polynomial()
    x, y = R.gens()
    
    I = self.ideal().intersection(other.ideal()) + R.ideal(F)
    J = R.ideal(I.groebner_basis())
    polys = J.gens()[0:]
    f = polys[0]
    if f.lm() == y^3 :
      polys = polys[1:]
      
    return C34CurveDivisor(C, polys)



  """
    Given an integer n (may be positive, zero, or negative), returns the reduced divisor equivalent
    in the divisor class group to the scalar multiple n*D.
    
    If D = P + Q + R, then this returns the unique reduced divisor equivalent to n*P + n*Q + n*R.

    This method is "slow" in the sense that it uses Sage's provided ideal arithmetic rather than
    explicit formulas. It *does* use the additive analogue of fast modular exponentiation.
  """
  def slow_scale(self, n) :
    if (n < 0) :
      return self.slow_flip().slow_scale(-n)

    C = self.C
    R = C.R
    F = C.defining_polynomial()
    x, y = R.gens()

    ret = C.zero_divisor()
    base = self
    while (n > 0) :
      if (n & 1 == 1) :
        ret = ret.slow_add(base)
      n = n >> 1
      base = base.slow_add(base)
    return ret



  """
    Returns the unique reduced divisor equivalent in the divisor class group to this divisor.
    
    If this divisor is already reduced, then a copy of it is returned.

    This method is "slow" in the sense that it uses Sage's provided ideal arithmetic rather than
    explicit formulas.
  """
  def slow_reduce(self) :
    if self.reduced :
      return copy(self)
    return self.slow_flip().slow_flip()
  
  
  
  def variety(self) : 
    return self.ideal().variety(self.K.base().algebraic_closure())
  
  
  
  def __add__(self, other) :
    """
      Input : Two typical C34CurveDivisors, D1 and D2.
      Output : The C34CurveDivisor D3 equivalent to D1 + D2. May be typical or semi-typical (or neither?) 
    """
    # TODO: Make sure divisors come from same curve
    # TODO: Still need to reduce divisors by flipping twice.
    #       E.g. return flip(flip(add(D1, D2)))
    D1, D2 = self, other
    if D1 == D2 :
      return double(D1)
    else :
      return add(D1, D2)
  
  def __eq__(self, other) :
    return (self.f == other.f) and (self.g == other.g) and (self.h == other.h)
  
  def __mul__(self, rhs) :
    # Assumes rhs is of an integer type.
    if (rhs < 0) :
      return self.__mul__(-rhs).__neg__()
    ret = self.C.zero_divisor()
    base = self
    while (rhs > 0) :
      if (rhs & 1 == 1) :
        ret = ret + base
      rhs = rhs >> 1
      base = base + base
    return ret
    
  
  def __ne__(self, other) :
    return not self.__eq__(other)
  
  def __neg__(self) :
    return flip(self)
  
  def __repr__(self) :
    s = ""
    f = self.polys()
    if len(f) == 0 :
      s = "<0>"
    elif len(f) == 1 :
      s = "<{}>".format(str(f[0]))
    elif len(f) == 2 :
      s = "<{}, {}>".format(str(f[0]), str(f[1]))
    elif len(f) == 3 :
      s = "<{}, {}, {}>".format(str(f[0]), str(f[1]), str(f[2]))
    return s
    
  def __rmul__(self, lhs) :
    return self.__mul__(lhs)

  def __sub__(self, rhs) :
    return self.__add__(rhs.__neg__())
