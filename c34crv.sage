"""
  General Form :
  xxxx + a_3*xxx + a_2*xxy + a_1*xyy + yyy + a_6*xx + a_5*xy + a_4*yy + a_9*x + a_8*y + a_12 

  Reduced Forms :
  In char K = 3
  xxxx + a_2*xxy + a_1*xyy + yyy + a_6*xx + a_5*xy + a_4*yy + a_9*x + a_8*y + a_12 

  In char K = 2
  yyy + (xx + x + 1)y + (xxxx + xxx + xx + x + 1)
  xxxx + a_3*xxx + a_2*xxy + yyy + a_6*xx + a_5*xy + a_9*x + a_8*y + a_12 

  In char K = 0 or char K > 3
  yyy + (xx + x + 1)y + (xxxx + xx + x + 1)
  xxxx + a_2*xxy + yyy + a_6*xx + a_5*xy + a_9*x + a_8*y + a_12 
"""

class C34Crv :
  def __init__(self, fld, coeffs, coords = "x,y") :
    if len(coeffs) != 9 :
      raise TypeError("coeffs must be a list of length 9 of field elements")
    
    self.K = fld
    self.c = copy(coeffs)
    self.c = [fld(t) for t in self.c]
    self.R = PolynomialRing(self.K, 2, coords, order = "m(3,4,0,1)")
    x, y = self.R.gens()
    m = [self.K.one(), x, y, x*x, x*y, y*y, x*x*x, x*x*y, x*y*y]
    self._poly = sum([m[i]*self.c[i] for i in range(9)]) + x^4 + y^3
    self._ambient_space = AffineSpace(2, self.K, coords)
    self._scheme = self._ambient_space.subscheme(self._poly)
    
    # Ensure curve is non-singular
    fx = self.partial_x(x,y)
    fy = self.partial_y(x,y)
    dim = self._scheme.intersection(self._ambient_space.subscheme(fx)).intersection(self._ambient_space.subscheme(fy)).dimension()
    if dim >= 0 :
      raise ValueError("Curve is singular. C = {}".format(self))

  def ambient_space(self) :
    return self._ambient_space

  def coefficients(self) :
    return copy(self.c)
  
  def divisor(self, points) :
    return C34CrvDiv(self, points)
  
  def hessian(self, x0, y0) :
    """
      Returns the determinant of the Hessian matrix of C at the point (x0, y0).
      The Hessian matrix is the matrix
        [ Fxx  Fxy ]
        [ Fyx  Fyy ]
      where, e.g., Fxy is the second partial derivative of the curve equation, obtained by
      taking the partial derivative with respect to x, then with respect to y.
      Returns (Fxx*Fyy - Fxy*Fyx)(x0, y0)
    """
    return self.partial_xx(x0, y0) * self.partial_yy(x0, y0) - self.partial_xy(x0, y0)^2
  
  def partial_x(self, x0, y0) :
    # Returns the value of the partial derivative of the curve equation with respect to its first variable at the point (x0, y0)
    c = self.c
    return 4*x0^3 + c[8]*y0^2 + 2*c[7]*x0*y0 + 3*c[6]*x0^2 + c[4]*y0 + 2*c[3]*x0 + c[1]
    
  def partial_y(self, x0, y0) :
    # Returns the value of the partial derivative of the curve equation with respect to its second variable at the point (x0, y0)
    c = self.c
    return 3*y0^2 + 2*c[8]*x0*y0 + c[7]*x0^2 + 2*c[5]*y0 + c[4]*x0 + c[2]
  
  def partial_xx(self, x0, y0) :
    c = self.c
    return 6*x0*(2*x0 + c[6]) + 2*(c[7]*y0 + c[3])
  
  def partial_xy(self, x0, y0) :
    c = self.c
    return 2*(c[8]*y0 + c[7]*x0) + c[4]
  
  def partial_yy(self, x0, y0) :
    c = self.c
    return 2*(3*y0 + c[8]*x0 + c[5])

  def point(self, *args) :
    return C34CrvPt(self, *args)
  
  def point_at_infinity(self) :
    return C34CrvPt(self, self.K.zero(), self.K.one(), self.K.zero())
  
  # Return all points in the curve's base field
  # Assumes C's base field K is finite
  def points(self) :
    ret = []
    x, y = self.R.gens()
    q = self.K.order()
    f, g = x^q - x, y^q - y
    V = self.R.ideal(f, g, self.poly()).variety(self.K)
    for v in V :
      P = self.point(v[x], v[y])
      ret = ret + [P]
    ret.sort()
    """
    if (not L.is_finite()) :
      print "Error: field must be finite."
      return
    ret = [self.point_at_infinity()]
    f = self.poly()
    for u in L :
      for v in L :
        if (f(u,v) == self.K(0)) :
          ret = ret + [C34CrvPt(self, u, v)]
    """
    return ret
  
  def poly(self) :
    return copy(self._poly)
  
  def poly_to_vec(self, poly) :
    """
      Converts a polynomial into a coefficient vector.
      
      The i'th element of the returned vector represents the coefficient of the i'th polynomial in the list
      [ 1, x, y, x^2, xy, y^2, x^3, x^2y, xy^2, x^4, ... ]
      
      This list does not have any powers of y larger than y^2.
      If the given polynomial has any powers of y larger than y^2, it is reduced modulo the curve equation
      until those terms are gone.
    """
    p = copy(poly)
    x, y = self.R.gens()
    F = self.poly()
    ret = []
    
    # First, reduce p modulo the curve equation to get rid of powers y^3 or larger
    d = p.degree(y)
    while (d > 2) :
      q = p.coefficient(y^d)*(y^(d-3))
      p = p - q*F
      d = p.degree(y)
    
    place = 0
    while (p != 0) :
      coeff = self.K.zero()
      if place == 0 :
        xdeg, ydeg = 0, 0
      elif place == 1 :
        xdeg, ydeg = 1, 0
      elif place == 2 :
        xdeg, ydeg = 0, 1
      else :
        xdeg = 1 + (place // 3) - (place % 3)
        ydeg = place % 3
      coeff = p.coefficient({x:xdeg,y:ydeg})
      p = p - coeff*x^xdeg*y^ydeg
      place = place + 1
      ret = ret + [coeff]
    
    return ret
    
  
  def random_divisor(self) :
    """
      Returns a random reduced divisor of this curve.
      
      The divisor is chosen by taking a large multiple of a randomly selected rational point on the
      curve. The choice of initial point is not made uniformly at random. The "large multiple" is
      a random integer in the range [0, q^3] where q is the order of the base field. This is based
      off a probably-incorrect assumption that the divisor class group has order approximately q^3.
    """
    min_multiple = 0
    max_multiple = self.K.order()^3
    D = self.divisor([self.random_point()])
    n = randint(min_multiple, max_multiple)
    return n*D
  
  def random_point(self) :
    """
      Returns a random affine rational point on the curve.
      
      This does not choose a point uniformly at random from all points on the curve.
      The random point is chosen by selecting a random x-coordinate in the base field of the curve.
      The three points (counting multiplicity) on the curve with this x-coordinate are determined.
      This requires computing cube roots and can be slow.
      If any of these point are rational, one is chosen at random.
      Otherwise, we choose a new x-coordinate and try again.
    """
    x, y = self.R.gens()
    a = self.K.random_element()
    f = self.poly().subs(x=a)
    arr = f.univariate_polynomial().roots()
    roots = [arr[i][0] for i in range(len(arr))]
    shuffle(roots)
    for r in roots :
      if r in self.K :
        return self.point(a, r)
    return self.random_point()
  
  def scheme(self) :
    return self._scheme       
  
  def tangent_line(self, x0, y0) :
    """
      Returns the equation of the tangent line to C at the point (x0, y0)
    """
    x, y = self.R.gens()
    return self.partial_x(x0, y0)*(x - x0) + self.partial_y(x0, y0)*(y - y0)
  
  def zero_divisor(self) :
    """
      Returns the identity element in the curve's divisor class group.
      This is the divisor corresponding to the principal ideal <1>.
    """
    return C34CrvDiv(self, [[self.K.one()], [], []])
  
  def __repr__(self) :
    ret = "C34 curve defined by {} over {}".format(str(self.poly()), str(self.K))
    return ret
 

  # If base field is finite, return all points 
  #def __call__(self, *args) :
  #  if len(args) == 1 :
  #    return self.points(args[0])
  #  elif (len(args) == 2) or (len(args) == 3) :
  #    return self.point(*args)
