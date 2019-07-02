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
class C34CrvDiv :
  def __init__(self, C, lst) :
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
    self.degree = -1
    self.reduced = False
    self.typical = False
    self.type = 0
    
    def classify_divisor() :      
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
        self.typical = False # TODO : Find condition under which it is typical
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

    # Make sure 'lst' is a list
    if not isinstance(lst, list) :
      raise TypeError("'lst' must be of type 'list'.")
    
    # Makes sure elements in 'lst' are of same type
    if not [type(x) for x in lst] == [type(x) for x in [lst[0]]*len(lst)] :
      raise TypeError("Elements in 'lst' must be of same type.")

    # Check type of elements in 'lst'.
    # Should be one of :
    #   * C34CrvPts
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
        LMf = f.monomials()[0]
        self.f = [0]*(mmap[LMf] + 1)
        for m in f.monomials() :
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
      
    elif isinstance(lst[0], C34CrvPt) : # If 'lst' elements are points
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
        D = C34CrvDiv(C, [points[0]]) + C34CrvDiv(C, points[1:])
        self.f, self.g, self.h = D.f, D.g, D.h

    elif isinstance(lst, list) :
      if len(lst) > 0 :
        self.f = [self.K(t) for t in lst[0]]
      if len(lst) > 1 :
        self.g = [self.K(t) for t in lst[1]]
      if len(lst) > 2 :
        self.h = [self.K(t) for t in lst[2]]

    classify_divisor()

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
    F = self.C.poly()
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
    I = self.R.ideal(self.polys() + [self.C.poly()])
    G = list(I.groebner_basis())
    G.sort()
    return self.R.ideal(G)

  def is_squarefree(self) :
    """
      Returns true if the ideal of D is squarefree.
      
      Equivalently, returns true if all points in the formal sum of D are distinct.
      Otherwise, returns false.
    """
    I = self.ideal()
    return I.radical() == I

  # Return a matrix whose columns are a basis for this divisor's corresponding vector space W_D^10
  def matrix(self) :
    return self.WD10().basis_matrix().transpose()



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
  # TODO : parse and return C34CrvPts
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
  
  def variety(self) : 
    return self.ideal().variety(self.K.base().algebraic_closure())
  
  # Return the vector space W_D^10
  # Assumes D is semi-typical
  def WD10(self) :
    e1 = self.f + [0]*4
    e2 = self.g + [0]*3
    e3 = self.h + [0]*2
    e4 = [0, self.f[0], 0, self.f[1], self.f[2], 0, 1, 0]
    e5 = [0, self.g[0], 0, self.g[1], self.g[2], 0, 0, 1]
    ret = VectorSpace(K, 8).subspace_with_basis([e1, e2, e3, e4, e5])
    return ret
  


  def __add__(self, other) :
    """
      Input : Two typical C34CrvDivs, D1 and D2.
      Output : The C34CrvDiv D3 equivalent to D1 + D2. May be typical or semi-typical (or neither?) 
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
