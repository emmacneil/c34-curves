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

class C34Curve :
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
  
  

  def differential(self, f) :
    """
      Returns fx*Fy - fy*Fx where F is the curve equation, f is the given polynomial, and, e.g.,
      Fx and fx are the partial formal derivatives of these with respect to x.

      The module of Kahler differentials modulo F is generated by a generator w such that

        dx = Fy*w,   dy = -Fx*w.
       
      This method sends a polynomial f to its differential df expressed in terms of this
      generator. We have

        df = fx*dx + fy*dy = fx*Fy - fy*Fx.
    """
    F = self.defining_polynomial()
    ret = f.derivative(x)*F.derivative(y) - f.derivative(y)*F.derivative(x)
    return ret.mod(F)



  def divisor(self, points) :
    return C34CurveDivisor(self, points)
  
  
  
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
    return C34CurvePoint(self, *args)
  
  
  
  def point_at_infinity(self) :
    """
      Returns the curve's point at infinity. This is always the point (0 : 1 : 0)
    """
    return C34CurvePoint(self, self.K.zero(), self.K.one(), self.K.zero())
  
  
  
  def rational_points(self) :
    ret = []
    x, y = self.R.gens()

    for k in self.K :
      p = self.defining_polynomial().subs(x=k).univariate_polynomial()
      rts = p.roots()
      for r in rts : 
        P = self.point(k, r[0])
        ret = ret + [P]

    ret.sort()
    return ret
  
  
  
  def defining_polynomial(self) :
    return copy(self._poly)
  
  
  
  @staticmethod
  def random_curve(K) :
    """
      Returns a random non-singular C34 curve over the field K.
    """
    while True :
      c = [K.random_element() for i in range(9)]
      try :
        C = C34Curve(K, c)
        return C
      except ValueError :
        continue
    return C
    
  
  
  def random_divisor(self) :
    """
      Returns a random reduced divisor on this curve.

      Any element of the divisor class group may be returned, but they are not chosen uniformly
      at random.

      A polynomial f is chosen uniformly at random. Then a random divisor whose reduced Groebner
      basis begins with f is chosen. Consequently, the more divisors there are beginning with f,
      the less likely each individual one is to be chosen.

      Constructing these divisors calls on Sage's polynomial factoring routines. This may be slow
      if one needs to construct many random divisors.
    """
    # TODO : Not returning any atypical type 31 divisors.
    K = self.K
    x, y = self.R.gens()
    F = self.defining_polynomial()
    
    f3, f2, f1, f0 = [self.K.random_element() for t in range(4)]

    if (f3 != 0) and (f2 != 0) :
      # Set f := y - f3*x^2 - f1*x - f0
      # Compute Fbar = F modulo f = F(x, f3*x^2 + f1*x + f0) = a6*x^6 + a5*x^5 + ... + a1*x + a0
      Fbar = F.subs(y = f3*x*x + f1*x + f0).univariate_polynomial().monic()

      # Factor Fbar into irreducible components. Store the factors in a list
      factors = []
      for fac, pwr in Fbar.factor() :
        factors = factors + [fac]*pwr
      factors.sort()

      # Determine if Fbar can be factored into a product of cubics.
      # If Fbar has a quartic or larger factor, restart.
      # If Fbar's smallest factor is quadratic, restart.
      if (factors[-1].degree() > 3) or (factors[0].degree() == 2) :
        return self.random_divisor()

      # Find a degree 3 polynomial g that divides Fbar.
      # If this is impossible, then restart.
      shuffle(factors)
      g = factors[0]
      for i in range(1, len(factors)) :
        if (g.degree() + factors[i].degree() < 4) :
          g = g*factors[i]
      z = g.parent().gen()
      g0 = g.constant_coefficient()
      g1 = g.monomial_coefficient(z)
      g3 = g.monomial_coefficient(z^2)

      # We have
      #   f = y - f3*x^2 - f1*x - f0
      #   g = x^3 + g3*x^2 + g1*x + g0
      # Rewrite f in the form
      #   f = x^2 + f2*y + f1*x + f0
      # and reduce g modulo f to get g of the form
      #   g = xy + g2*y + g1*x + g0
      f2 = -1/f3
      f1 = -f1*f2
      f0 = -f0*f2
      g2 = f3*(f2*(f1 - g3))
      g1 = f3*(f1*(f1 - g3) + g1 - f0)
      g0 = f3*(f0*(f1 - g3) + g0)
      
      # Compute f = ((y + g1)f - (x + f1 - g2)g)/f2
      h2 = g1 - f3*(g2*(g2 - f1) + f0)
      h1 = -f3*(g1*g2 - g0)
      h0 = -f3*(f0*g1 + g0*(g2 - f1))

      return C34CurveDivisor(self, [[f0, f1, f2, 1], [g0, g1, g2, 0, 1], [h0, h1, h2, 0, 0, 1]])

    elif (f3 != 0) and (f2 == 0) :
      # Set f := x^2 + f1*x + f0
      # Check if f can be factored into
      #   f = pr = (x + p0)(x + r0)
      # for some p0, r0 in the base field.
      f = x^2 + f1*x + f0
      factors = []
      for fac, pwr in f.factor() :
        if fac.degree(x) > 1 :
          return self.random_divisor()
        factors = factors + [fac]*pwr
      shuffle(factors)
      p0 = factors[0].constant_coefficient()
      r0 = factors[1].constant_coefficient()
      p = x + p0
      r = x + r0

      # Check if there are 2 rational points with x-coordinates -p0 and -r0
      if (F.subs(x = -p0).univariate_polynomial().is_irreducible()) :
        return self.random_divisor()
      if (F.subs(x = -r0).univariate_polynomial().is_irreducible()) :
        return self.random_divisor()

      # Generated a type 11 divisor D1 = (x + p0, y + q0)
      # and a type 22 divisor D2 = (x + r0, y^2 + s2*y + s0).
      
      # Generating D1...
      # Compute t = F (mod x + p0) and factor t
      t = F.subs(x = -p0).univariate_polynomial()
      factors = []
      for fac, pwr in t.factor() :
        if fac.degree() > 1 :
          continue
        factors = factors + [fac]*pwr
      shuffle(factors)
      q0 = factors[0].constant_coefficient()
      D1 = C34CurveDivisor(self, [[p0, 1], [q0, 0, 1], []])

      # Generating D2
      # Compute t = F (mod x + r0) and factor t
      t = F.subs(x = -r0).univariate_polynomial()
      factors = []
      for fac, pwr in t.factor() :
        if fac.degree() > 2 :
          continue
        factors = factors + [fac]*pwr
      # t is either a product of three linear terms, or of a quadratic and a linear term.
      if (len(factors) == 3) :
        shuffle(factors)
        u = factors[0]*factors[1]
      else :
        factors.sort() # Sorting puts the linear term first, the quadratic last.
        u = factors[1]
      s2 = u.monomial_coefficient(t.parent().gen()^2)
      s0 = u.constant_coefficient()
      D2 = C34CurveDivisor(self, [[r0, 1], [s0, 0, s2, 0, 0, 1], []])
      D = D1 + D2
      if (D.type == 31) :
        return D if randint(0,1) == 0 else -D
      else :
        return self.random_divisor()
      
    elif (f2 != 0) :
      # Set f := y + f1*x + f0
      # Reduce F modulo f to get a monic quartic t in x
      t = F.subs(y = -f1*x - f0).univariate_polynomial()

      # Factor t.
      # If this factorization has a cubic or quartic term, restart.
      factors = []
      for fac, pwr in t.factor() :
        if fac.degree() > 2 :
          return self.random_divisor()
        factors = factors + [fac]*pwr

      # t is either four linear terms, two linears and a quadratic, or two quadratics.
      if (len(factors) == 2) :
        g = factors[0] if (randint(0,1) == 0) else factors[1]
      elif (len(factors) == 3) :
        factors.sort() # Sorts the quadratic to the end of the list
        g = factors[0]*factors[1] if (randint(0,1) == 0) else factors[2]
      else :
        shuffle(factors)
        g = factors[0]*factors[1]
      g1 = g.monomial_coefficient(g.parent().gen())
      g0 = g.constant_coefficient()
      return C34CurveDivisor(self, [[f0, f1, 1], [g0, g1, 0, 1], []])
    
    elif (f1 != 0) :
      # Set f = x + f0 and compute t = F (mod f)
      t = F.subs(x = -f0).univariate_polynomial()
      
      # Factor t, but throw away any non-linear terms.
      factors = []
      for fac, pwr in t.factor() :
        if fac.degree() > 1 :
          continue
        factors = factors + [fac]*pwr

      # If t is irreducible over K, then retart
      if (len(factors) == 0) :
        return self.random_divisor()

      # Choose a linear term at random
      shuffle(factors)
      g0 = factors[0].constant_coefficient()
      D = C34CurveDivisor(self, [[f0, 1], [g0, 0, 1], []])

      # Randomly return either D or -D
      return D if randint(0,1) == 0 else -D

    elif (f0 != 0) :
      return self.zero_divisor()
    
    else :
      return self.random_divisor()
  
  

  def random_divisor_of_type(self, T, typical = 2) :
    """
      Returns a random divisor on this curve of the specified type T.
      
      If typical is 1 or True, then the divisor will be typical.
      If typical is 0 or False, then the divisor will non-typical.
      If typical is unspecified, then the divisor may or may not be typical.
      This parameter is ignored if T is not one of 31, 41, 51, or 61.
      
      TODO : This method can lead to infinite recursion if the curve does not have a divisor of
             the desired type. E.g., the curve
             
               C = y^3 + x^4 + x*y + y over GF(2)
             
             has no divisors of types 31, 41, 51, 61. The curve 
             
               C = y^3 + x^4 + y + x + 1 over GF(2)
             
             has only one divisor class -- every divisor is principal, equivalent to 0. The curve
             
               C = y^3 + x^4 - 2*y + 2 over GF(5)
             
             has no divisors of types 11, 22, 31, 32, 41, 42, 51, 52, 54, 62, 61, 64.
    """
    K = self.K
    x, y = self.R.gens()
    F = self.defining_polynomial()

    ret = self.zero_divisor()
    
    if T == 0 :
      ret = self.zero_divisor()
    elif T == 11 :
      ret = self.divisor([self.random_rational_point()])
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
        return self.random_divisor_of_type(T)
      shuffle(factors)
      g = self.R.one()
      while g.degree(x) < 2 :
        if g.degree(x) + factors[0].degree(x) <= 2 :
          g = g * factors[0]
        factors = factors[1:]
      ret = C34CurveDivisor(self, [f, g])
    elif T == 22 :
      ret = - self.random_divisor_of_type(11)
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
        ret = self.random_divisor_of_type(11).slow_compose(self.random_divisor_of_type(22)) # Can sometimes be type 33
        if ret.type == 33 :
          return self.random_divisor_of_type(T, typical)
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
          return self.random_divisor_of_type(T, typical)
        # If there are only quadratic factors, restart
        if factors[0].degree(x) == 2 :
          return self.random_divisor_of_type(T, typical)
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
        ret = C34CurveDivisor(self, [f, g, h])
    elif T == 32 :
      # Pick a random type 11 divisor D = <f, g>.
      D = self.random_divisor_of_type(11)
      f, g = D.polys()
      # Generate a random monic line h = y + ... in I_D
      a0 = self.K.random_element()
      h = a0*f + g
      # Compute the colon ideal h : I_D
      I = self.R.ideal(h, F).quotient(self.R.ideal(f, g, F))
      ret = C34CurveDivisor(self, I.gens()[0:2])
    elif T == 33 :
      ret = C34CurveDivisor(self, [[self.K.random_element(), self.K.one()], [], []])
    elif T == 41 :
      # Pick a random type 31 divisor D = <f, g, h>.
      D = self.random_divisor_of_type(31, typical)
      f, g, h = D.polys()
      # Generate a random monic hyperbola p = xy + ... in I_D
      a0 = self.K.random_element()
      p = a0*f + g
      # Compute the colon ideal h : I_D
      I = self.R.ideal(p, F).quotient(self.R.ideal(f, g, h, F))
      ret = C34CurveDivisor(self, I.gens()[0:3])
    elif T == 42 :
      # Pick a random type 22 divisor D = <f, g>.
      D = self.random_divisor_of_type(22)
      f, g = D.polys()
      # Generate a random monic parabola h = x^2 + ... in I_D
      a0 = self.K.random_element()
      a = x + a0
      h = a*f
      # Compute the colon ideal h : I_D
      I = self.R.ideal(h, F).quotient(self.R.ideal(f, g, F))
      ret = C34CurveDivisor(self, I.gens()[0:2])
    elif T == 43 :
      # Pick a random type 21 divisor D = <f, g>.
      D = self.random_divisor_of_type(21)
      f, g = D.polys()
      # Generate a random monic parabola h = x^2 + ... in I_D
      a0 = self.K.random_element()
      h = a0*f + g
      # Compute the colon ideal h : I_D
      I = self.R.ideal(h, F).quotient(self.R.ideal(f, g, F))
      ret = C34CurveDivisor(self, I.gens()[0:2])
    elif T == 44 :
      ret = C34CurveDivisor(self, [[self.K.random_element(), self.K.random_element(), self.K.one()], [], []])
    elif T == 51 :
      # Pick a random type 31 divisor D = <f, g, h>.
      D = self.random_divisor_of_type(31, typical)
      f, g, h = D.polys()
      # Generate a random monic quadratic p = y^2 + ... in I_D
      a0 = self.K.random_element()
      b0 = self.K.random_element()
      p = a0*f + b0*g + h
      # Compute the colon ideal h : I_D
      I = self.R.ideal(p, F).quotient(self.R.ideal(f, g, h, F))
      ret = C34CurveDivisor(self, I.gens()[0:3])
    elif T == 52 :
      # Pick a random type 22 divisor D = <f, g>.
      D = self.random_divisor_of_type(22)
      f, g = D.polys()
      # Generate a random monic hyperbola h = xy + ... in I_D
      a1 = self.K.random_element()
      a0 = self.K.random_element()
      a = y + a1*x + a0
      h = a*f
      # Compute the colon ideal h : I_D
      I = self.R.ideal(h, F).quotient(self.R.ideal(f, g, F))
      ret = C34CurveDivisor(self, I.gens()[0:2])
    elif T == 53 :
      # Pick a random type 21 divisor D = <f, g>.
      D = self.random_divisor_of_type(21)
      f, g = D.polys()
      # Generate a random monic hyperbola h = xy + ... in I_D
      a0 = self.K.random_element()
      b0 = self.K.random_element()
      a = x + a0
      h = a*f + b0*g
      # Compute the colon ideal h : I_D
      I = self.R.ideal(h, F).quotient(self.R.ideal(f, g, F))
      ret = C34CurveDivisor(self, I.gens()[0:2])
    elif T == 54 :
      # Pick a random type 11 divisor D = <f, g>.
      D = self.random_divisor_of_type(11)
      f, g = D.polys()
      # Generate a random monic parabola h = x^2 + ... in I_D
      a0 = self.K.random_element()
      b0 = self.K.random_element()
      a = x + a0
      h = a*f + b0*g
      # Compute the colon ideal h : I_D
      I = self.R.ideal(h, F).quotient(self.R.ideal(f, g, F))
      ret = C34CurveDivisor(self, I.gens()[0:2])
    elif T == 61 :
      # Pick a random type 31 divisor D = <f, g, h>.
      D = self.random_divisor_of_type(31, typical)
      f, g, h = D.polys()
      # Generate a random monic cubic p = x^3 + ... in I_D
      a0 = self.K.random_element()
      b0 = self.K.random_element()
      c0 = self.K.random_element()
      a = x + a0
      p = a*f + b0*g + c0*h
      # Compute the colon ideal h : I_D
      I = self.R.ideal(p, F).quotient(self.R.ideal(f, g, h, F))
      ret = C34CurveDivisor(self, I.gens()[0:3])
    elif T == 62 :
      # Pick a random type 22 divisor D = <f, g>.
      D = self.random_divisor_of_type(22)
      f, g = D.polys()
      # Generate a random monic quadratic h = y^2 + ... in I_D
      a1 = self.K.random_element()
      a0 = self.K.random_element()
      a = y + a1*x + a0
      h = a*f + g
      # Compute the colon ideal h : I_D
      I = self.R.ideal(h, F).quotient(self.R.ideal(f, g, F))
      ret = C34CurveDivisor(self, I.gens()[0:2])
    elif T == 63 :
      # Pick a random type 21 divisor D = <f, g>.
      D = self.random_divisor_of_type(21)
      f, g = D.polys()
      # Generate a random monic quadratic h = y^2 + ... in I_D
      a1 = self.K.random_element()
      a0 = self.K.random_element()
      b0 = self.K.random_element()
      a = y + a1*x + a0
      h = a*f + b0*g
      # Compute the colon ideal h : I_D
      I = self.R.ideal(h, F).quotient(self.R.ideal(f, g, F))
      ret = C34CurveDivisor(self, I.gens()[0:2])
    elif T == 64 :
      # Pick a random type 11 divisor D = <f, g>.
      D = self.random_divisor_of_type(11)
      f, g = D.polys()
      # Generate a random monic hyperbola h = xy + ... in I_D
      a1 = self.K.random_element()
      a0 = self.K.random_element()
      b0 = self.K.random_element()
      a = y + a1*x + a0
      h = a*f + b0*g
      # Compute the colon ideal h : I_D
      I = self.R.ideal(h, F).quotient(self.R.ideal(f, g, F))
      ret = C34CurveDivisor(self, I.gens()[0:2])
    elif T == 65 :
      ret = C34CurveDivisor(self, [[self.K.random_element(), self.K.random_element(), self.K.random_element(), self.K.one()], [], []])
    else :
      raise ValueError("\"Type {}\" is not a valid type".format(T))
    
    assert ret.type == T, "{} is not of type {}".format(ret, T)
    return ret
  
  
  
  def random_rational_point(self) :
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
    f = self.defining_polynomial().subs(x=a)
    arr = f.univariate_polynomial().roots()
    roots = [arr[i][0] for i in range(len(arr))]
    shuffle(roots)
    for r in roots :
      if r in self.K :
        return self.point(a, r)
    return self.random_rational_point()
  
  
  
  def scheme(self) :
    return self._scheme       



  def short_form(self) :
    """
      Performs an invertible change of coordinates and returns a C34 curve isomorphic to C such
      that some curve coefficients of the new curve equation are zero.

      Let K be the field over which this curve is defined.
      If K has characteristic 2, C2 has the form

        C2 : y^3 + x^4 + d7*x^2*y + d6*x^3 + d4*x*y + d3*x^2 + d2*y + d1*x + d0.
      
      If K has characteristic 3, C2 has the form

        C2 : y^3 + x^4 + d8*x*y^2 + d7*x^2*y + d5*y^2 + d4*x*y + d3*x^2 + d2*y + d1*x + d0.

      Otherwise, C2 has the form

        C2 : y^3 + x^4 + d7*x^2*y + d4*x*y + d3*x^2 + d2*y + d1*x + d0.
    """
    K = self.K
    R = self.R
    x, y = R.gens()
    char = K.characteristic()
    F = self.defining_polynomial()
    c = self.coefficients()

    if (char == 2) :
      F2 = F.subs(y = y - (c[8]*x - c[5])/K(3))
    elif (char == 3) :
      F2 = F.subs(x = x - c[6]/K(4))
    else :
      a = (27*c[6] - 9*c[7]*c[8] + 2*c[8]^3)/K(27)
      F2 = F.subs(x = x - a/K(4), y = y - c[8]/K(3)*x + (a*c[8] - 4*c[5])/K(12))
    c2 = [F2.constant_coefficient(),
          F2.monomial_coefficient(x),
          F2.monomial_coefficient(y),
          F2.monomial_coefficient(x*x),
          F2.monomial_coefficient(x*y),
          F2.monomial_coefficient(y*y),
          F2.monomial_coefficient(x*x*x),
          F2.monomial_coefficient(x*x*y),
          F2.monomial_coefficient(x*y*y)]
    return C34Curve(K, c2)
  
  
  
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
    return C34CurveDivisor(self, [[self.K.one()], [], []], degree = 0, typ = 0, reduced = True, typical = False)
  
  
  
  def __eq__(self, other) :
    return (self.K == other.K) and (self.c == other.c)
  
  
  
  def __neq__(self, other) :
    return not self.__eq__(other)
  
  
  
  def __repr__(self) :
    ret = "C34 curve defined by {} over {}".format(str(self.defining_polynomial()), str(self.K))
    return ret

