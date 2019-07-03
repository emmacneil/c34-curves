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
  # Assumes K is a finite field
  def points(self) :
    ret = []
    x, y = self.R.gens()

    for k in self.K :
      p = self.poly().subs(x=k).univariate_polynomial()
      rts = p.roots()
      for r in rts : 
        P = self.point(k, r[0])
        ret = ret + [P]

    ret.sort()
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
    
  
  
  @staticmethod
  def random_curve(K) :
    """
      Returns a random non-singular C34 curve over the field K.
    """
    while True :
      c = [K.random_element() for i in range(9)]
      try :
        C = C34Crv(K, c)
        return C
      except ValueError :
        continue
    return C
    
  
  
  def random_divisor(self, T = -1, typical = 2) :
    """
      Returns a random divisor on this curve.
      
      If T is specified, then the divisor will be of type T.
      
      If T is not specified, then the random divisor is generated by computing r*P where P is a
      randomly chosen affine rational point on the curve, where r is a randomly chosen integer
      between 0 and q^3, and where q is the order of the base field over which C is defined.
      The divisor generated will be reduced, therefore of type 11, 21, 22, or 31. Divisors of type
      other than 31 are rarely generated over large finite fields.
      
      If typical is 1 or True, then the divisor will be typical.
      If typical is 0 or False, then the divisor will not be typical.
      If typical is unspecified, then the divisor may or may not be typical.
      This is only relevant if T is specified and is either 31, 41, 51, or 61.
    """
    K = self.K
    x, y = self.R.gens()
    F = self.poly()

    ret = self.zero_divisor()
    
    if T == 0 :
      ret = self.zero_divisor()
    elif T == 11 :
      ret = self.divisor([self.random_point()])
    elif T == 21 :
      # Generate a random non-vertical line, f.
      f1 = K.random_element()
      f0 = K.random_element()
      f = y + f1*x + f0
      # Reduce curve equation F modulo f.
      Ff = F.mod(f)
      # Factor F mod f
      factors = []
      for t in Ff.factor() :
        factors = factors + [t[0]]*t[1]
      factors.sort()
      # If there is a cubic or larger factor, then restart
      if factors[-1].degree(x) >= 3 :
        return self.random_divisor(T)
      shuffle(factors)
      g = self.R.one()
      while g.degree(x) < 2 :
        if g.degree(x) + factors[0].degree(x) <= 2 :
          g = g * factors[0]
        factors = factors[1:]
      ret = C34CrvDiv(self, [f, g])
    elif T == 22 :
      ret = - self.random_divisor(11)
    elif T == 31 :
      # Generate a random parabola, y = f2*x^2 + f1*x + f0.
      f2 = K.random_element()
      if (typical == False) :
        f2 = K.zero()
      elif (typical == True) : 
        while f2 == 0 :
          f2 = K.random_element()
      f1 = K.random_element()
      f0 = K.random_element()
      if f2 == 0 :
        ret = sage_compose(self.random_divisor(11), self.random_divisor(22)) # Can sometimes be type 33
        if ret.type == 33 :
          return self.random_divisor(T, typical)
      else :
        f = f2*x*x + f1*x + f0
        # Reduce curve equation F modulo f.
        Ff = F.subs(y = f)
        Ff = Ff / Ff.coefficient(x^6)
        # Factor F mod f
        factors = []
        for t in Ff.factor() :
          factors = factors + [t[0]]*t[1]
        factors.sort()
        # If there is a quartic or larger factor, then restart
        if factors[-1].degree(x) >= 4 :
          return self.random_divisor(T, typical)
        # If there are only quadratic factors, restart
        if factors[0].degree(x) == 2 :
          return self.random_divisor(T, typical)
        # If there are exactly 2 linear factors, then we *must* take one linear and one quadratic factor
        funny_case = False
        if len(factors) >= 3 :
          if (factors[0].degree(x) == 1) and (factors[1].degree(x) == 1) and (factors[2].degree(x) > 1) :
            funny_case = True
        shuffle(factors)
        g = self.R.one()
        while g.degree(x) < 3 :
          if funny_case == True :
            # Throw out the first linear term we encounter
            if factors[0].degree(x) == 1 :
              factors = factors[1:]
              funny_case = False
          if g.degree(x) + factors[0].degree(x) <= 3 :
            g = g * factors[0]
          factors = factors[1:]
        f = y - f
        f = f / f.coefficient(x^2)
        g = g.mod(f)
        g = g / g.coefficient(x*y)
        a = y + g.monomial_coefficient(x)
        b = x + f.monomial_coefficient(x) - g.monomial_coefficient(y)
        f2 = f.monomial_coefficient(y)
        h = (a*f - b*g)/f2
        ret = C34CrvDiv(self, [f, g, h])
    elif T == 32 :
      # Pick a random type 11 divisor D = <f, g>.
      D = self.random_divisor(11)
      f, g = D.polys()
      # Generate a random monic line h = y + ... in I_D
      a0 = self.K.random_element()
      h = a0*f + g
      # Compute the colon ideal h : I_D
      I = self.R.ideal(h, F).quotient(self.R.ideal(f, g, F))
      ret = C34CrvDiv(self, I.gens()[0:2])
    elif T == 33 :
      ret = C34CrvDiv(self, [[self.K.random_element(), self.K.one()], [], []])
    elif T == 41 :
      # Pick a random type 31 divisor D = <f, g, h>.
      D = self.random_divisor(31, typical)
      f, g, h = D.polys()
      # Generate a random monic hyperbola p = xy + ... in I_D
      a0 = self.K.random_element()
      p = a0*f + g
      # Compute the colon ideal h : I_D
      I = self.R.ideal(p, F).quotient(self.R.ideal(f, g, h, F))
      ret = C34CrvDiv(self, I.gens()[0:3])
    elif T == 42 :
      # Pick a random type 22 divisor D = <f, g>.
      D = self.random_divisor(22)
      f, g = D.polys()
      # Generate a random monic parabola h = x^2 + ... in I_D
      a0 = self.K.random_element()
      a = x + a0
      h = a*f
      # Compute the colon ideal h : I_D
      I = self.R.ideal(h, F).quotient(self.R.ideal(f, g, F))
      ret = C34CrvDiv(self, I.gens()[0:2])
    elif T == 43 :
      # Pick a random type 21 divisor D = <f, g>.
      D = self.random_divisor(21)
      f, g = D.polys()
      # Generate a random monic parabola h = x^2 + ... in I_D
      a0 = self.K.random_element()
      h = a0*f + g
      # Compute the colon ideal h : I_D
      I = self.R.ideal(h, F).quotient(self.R.ideal(f, g, F))
      ret = C34CrvDiv(self, I.gens()[0:2])
    elif T == 44 :
      ret = C34CrvDiv(self, [[self.K.random_element(), self.K.random_element(), self.K.one()], [], []])
    elif T == 51 :
      # Pick a random type 31 divisor D = <f, g, h>.
      D = self.random_divisor(31, typical)
      f, g, h = D.polys()
      # Generate a random monic quadratic p = y^2 + ... in I_D
      a0 = self.K.random_element()
      b0 = self.K.random_element()
      p = a0*f + b0*g + h
      # Compute the colon ideal h : I_D
      I = self.R.ideal(p, F).quotient(self.R.ideal(f, g, h, F))
      ret = C34CrvDiv(self, I.gens()[0:3])
    elif T == 52 :
      # Pick a random type 22 divisor D = <f, g>.
      D = self.random_divisor(22)
      f, g = D.polys()
      # Generate a random monic hyperbola h = xy + ... in I_D
      a1 = self.K.random_element()
      a0 = self.K.random_element()
      a = y + a1*x + a0
      h = a*f
      # Compute the colon ideal h : I_D
      I = self.R.ideal(h, F).quotient(self.R.ideal(f, g, F))
      ret = C34CrvDiv(self, I.gens()[0:2])
    elif T == 53 :
      # Pick a random type 21 divisor D = <f, g>.
      D = self.random_divisor(21)
      f, g = D.polys()
      # Generate a random monic hyperbola h = xy + ... in I_D
      a0 = self.K.random_element()
      b0 = self.K.random_element()
      a = x + a0
      h = a*f + b0*g
      # Compute the colon ideal h : I_D
      I = self.R.ideal(h, F).quotient(self.R.ideal(f, g, F))
      ret = C34CrvDiv(self, I.gens()[0:2])
    elif T == 54 :
      # Pick a random type 11 divisor D = <f, g>.
      D = self.random_divisor(11)
      f, g = D.polys()
      # Generate a random monic parabola h = x^2 + ... in I_D
      a0 = self.K.random_element()
      b0 = self.K.random_element()
      a = x + a0
      h = a*f + b0*g
      # Compute the colon ideal h : I_D
      I = self.R.ideal(h, F).quotient(self.R.ideal(f, g, F))
      ret = C34CrvDiv(self, I.gens()[0:2])
    elif T == 61 :
      # Pick a random type 31 divisor D = <f, g, h>.
      D = self.random_divisor(31, typical)
      f, g, h = D.polys()
      # Generate a random monic cubic p = x^3 + ... in I_D
      a0 = self.K.random_element()
      b0 = self.K.random_element()
      c0 = self.K.random_element()
      a = x + a0
      p = a*f + b0*g + c0*h
      # Compute the colon ideal h : I_D
      I = self.R.ideal(p, F).quotient(self.R.ideal(f, g, h, F))
      ret = C34CrvDiv(self, I.gens()[0:3])
    elif T == 62 :
      # Pick a random type 22 divisor D = <f, g>.
      D = self.random_divisor(22)
      f, g = D.polys()
      # Generate a random monic quadratic h = y^2 + ... in I_D
      a1 = self.K.random_element()
      a0 = self.K.random_element()
      a = y + a1*x + a0
      h = a*f + g
      # Compute the colon ideal h : I_D
      I = self.R.ideal(h, F).quotient(self.R.ideal(f, g, F))
      ret = C34CrvDiv(self, I.gens()[0:2])
    elif T == 63 :
      # Pick a random type 21 divisor D = <f, g>.
      D = self.random_divisor(21)
      f, g = D.polys()
      # Generate a random monic quadratic h = y^2 + ... in I_D
      a1 = self.K.random_element()
      a0 = self.K.random_element()
      b0 = self.K.random_element()
      a = y + a1*x + a0
      h = a*f + b0*g
      # Compute the colon ideal h : I_D
      I = self.R.ideal(h, F).quotient(self.R.ideal(f, g, F))
      ret = C34CrvDiv(self, I.gens()[0:2])
    elif T == 64 :
      # Pick a random type 11 divisor D = <f, g>.
      D = self.random_divisor(11)
      f, g = D.polys()
      # Generate a random monic hyperbola h = xy + ... in I_D
      a1 = self.K.random_element()
      a0 = self.K.random_element()
      b0 = self.K.random_element()
      a = y + a1*x + a0
      h = a*f + b0*g
      # Compute the colon ideal h : I_D
      I = self.R.ideal(h, F).quotient(self.R.ideal(f, g, F))
      ret = C34CrvDiv(self, I.gens()[0:2])
    elif T == 65 :
      ret = C34CrvDiv(self, [[self.K.random_element(), self.K.random_element(), self.K.random_element(), self.K.one()], [], []])
    else :
      r = randint(0, self.K.order()^3)
      return r*self.random_divisor(11)
    
    assert ret.type == T, "{} is not of type {}".format(ret, T)
    return ret
  
  
  
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

