"""
  A point on a C34 Curve
"""
class C34CurvePoint :
  def __init__(self, curve, *coords) :
    # TODO : This should raise an exception if
    #         * coords are 0, 0, 0
    #         * coords have length other than 2 or 3
    #         * coords do not satisfy curve equation 
    if (len(coords) < 2) or (len(coords) > 3) :
      raise TypeError("coords = {} must be of length 2 or 3.".format(coords))
    
    self.crv = curve
    c = curve.coefficients()

    # If any coordinate comes from an extension field, take the largest extension.    
    self.K = self.crv.K
    for P in coords :
      if is_field(parent(P)) :
        L = parent(P)
        if (L.degree() > self.K.degree()) :
          self.K = L

    self.x = [self.K.zero(), self.K.zero(), self.K.zero()]

    # If two coordinates are given, construct an affine point
    if (len(coords) == 2) :
      self.x[0], self.x[1], self.x[2] = self.K(coords[0]), self.K(coords[1]), self.K.one()
    
    # If three coordinates are given, construct a projective point.
    # Normalize its last non-zero coordinate to be 1.
    elif (len(coords) == 3) :
      if (coords[2] != self.K.zero()) :
        self.x[0], self.x[1], self.x[2] = self.K(coords[0])/self.K(coords[2]), self.K(coords[1])/coords[2], self.K.one()
      elif (coords[1] != self.K.zero()) :
        self.x[0], self.x[1], self.x[2] = self.K(coords[0])/self.K(coords[1]), self.K.one(), self.K.zero()
      elif (coords[0] != self.K.zero()) :
        self.x[0], self.x[1], self.x[2] = self.K.one(), self.K.zero(), self.K.zero()
      else :
        raise ValueError("Error: (0 : 0 : 0) is not a point in the projective plane.")
    
    # Verify this point is actually a point on the curve
    x0, x1, x2 = self.x
    if x2 == 0 :
      if (x0 != 0) or (x1 != 1) :
        raise ValueError("({} : {} : {}) is not a point on the curve.".format(x0, x1, x2))
    else :
      val = x0^4 + x1^3 + c[8]*x0*x1^2 + c[7]*x0^2*x1 + c[6]*x0^3 + c[5]*x1^2 + c[4]*x0*x1 + c[3]*x0^2 + c[2]*x1 + c[1]*x0  + c[0]
      if val != 0 :
        raise ValueError("({} : {} : {}) is not a point on the curve.".format(x0, x1, x2))

  def base_field(self) :
    return self.K

  def __eq__(self, other) :
    return self.x == other.x

  def __getitem__(self, key) :
    return self.x[key]
  
  def __lt__(self, other) :
    if self.x[0] != other.x[0] :
      return self.x[0] < other.x[0]
    else :
      return self.x[1] < other.x[1]
  
  def __repr__(self) :
    ret = "({0} : {1} : {2})".format(str(self.x[0]), str(self.x[1]), str(self.x[2]))
    return ret

