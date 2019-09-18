"""
var("x,y")
var("f0,f1,f2,f3,f4,f5")
var("g0,g1,g2,g3,g4,g5")
var("h0,h1,h2,h3,h4,h5")
var("F0,F1,F2,F3,F4,F5")
var("G0,G1,G2,G3,G4,G5")
var("H0,H1,H2,H3,H4,H5")
var("c0,c1,c2,c3,c4,c5,c6,c7,c8")

c = x*x*x*x + y*y*y + c8*x*y*y + c7*x*x*y + c6*x*x*x + c5*y*y + c4*x*y + c3*x*x + c2*y + c1*x + c0

f = coeff_vec([f0, f1, 1])
g = coeff_vec([g0, g1, 0, 1])
#h = coeff_vec([h0, h1, h2, h3, h4, h5, 0, 0, 1])
F = coeff_vec([F0, F1, 1])
G = coeff_vec([G0, G1, 0, 1])
#H = coeff_vec([H0, H1, H2, H3, H4, H5, 0, 0, 1])

M1 = Matrix([vecpow(f, 0, 0), vecpow(g, 0, 0), vecpow(f, 1, 0), vecpow(f, 0, 1), vecpow(g, 1, 0)]).transpose()
M2 = Matrix([vecpow(F, 0, 0), vecpow(G, 0, 0), vecpow(F, 1, 0), vecpow(F, 0, 1), vecpow(G, 1, 0)]).transpose()
M = M1.augment(M2).submatrix(0, 0, 7, 10)
print(M)

COEFF_VEC_SIZE = 19
def coeff_vec(v) :
  if len(v) > COEFF_VEC_SIZE :
    return v[0:COEFF_VEC_SIZE]
  while len(v) < COEFF_VEC_SIZE :
    v = v + [0]
  return v

# Assumes largest non-zero monomial in v is x^6
# I.e. v represents a polynomial of order <= 18
def times_x(v) :
  return [0, v[0], 0, v[1], v[2], 0] + v[3:COEFF_VEC_SIZE-3]

# Assumes largest non-zero monomial in v is x^3y^2
# I.e. v represents a polynomial of order <= 17
def times_y(v) :
  return [      - c0*v[5],
                - c1*v[5] - c0*v[8],
           v[0] - c2*v[5],
                - c3*v[5] - c1*v[8] - c0*v[11],
           v[1] - c4*v[5] - c2*v[8],
           v[2] - c5*v[5],
                - c6*v[5] - c3*v[8] - c1*v[11] - c0*v[14],
           v[3] - c7*v[5] - c4*v[8] - c2*v[11],
           v[4] - c8*v[5] - c5*v[8],
         - v[5]           - c6*v[8] - c3*v[11] - c1*v[14],
           v[6]           - c7*v[8] - c4*v[11] - c2*v[14],
           v[7]           - c8*v[8] - c5*v[11],
         - v[8]                     - c6*v[11] - c3*v[14],
           v[9]                     - c7*v[11] - c4*v[14],
           v[10]                    - c8*v[11] - c5*v[14],
         - v[11]                               - c6*v[14],
                                               - c7*v[14],
                                               - c8*v[14],
         - v[12]]

# Multiplies the coefficient vector v by the monomial (x^a)*(y^b)
def vecpow(v, a, b) :
  if (a == 0) and (b == 0) :
    return v
  elif (a == 0) :
    return vecpow(times_y(v), a, b - 1)
  else :
    return vecpow(times_x(v), a - 1, b)
"""


# If D is a degree 3 divisor given by <f, g> where
#   f = x^2 + f[2]*y + f[1]*x + f[0]
#   g = x*y + g[2]*y + g[1]*x + g[0]
# computes a third polynomial
#   h = y^2 + h[2]*y + h[1]*x + h[0]
# such that D = <f, g, h>
def compute_h(D) :
  if (D.degree != 3) or (len(D.f) != 4) or (len(D.g) != 5) or (len(D.h) != 0) or (D.f[2] == 0):
    raise ValueError("D is not a divisor of the form <x^2, xy>. D = {}".format(D))

  f, g = D.f, D.g

  # This gives h of the form y^2 + ry + sx + t in 1I 7M
  a = 1/f[2]
  h = [a*(g[1]*f[0] - (f[1] - g[2])*g[0]),
       a*(- g[0] + g[1]*g[2]),
       a*(f[0] - (f[1] - g[2])*g[2]) + g[1],
       K.zero(), K.zero(), K.one()]

  # This gives h of the form y^2 + axy + bx^2 + cy + dx + e in 1I 4M
  # h = [K.zero(), -g[0]*a, f[0]*a, -g[1]*a , (f[1]-g[2])*a , K.one()]
  
  # This gives h of the form y^2 + ay + bx + c in 1I 9M
  # h = [a*(       g[1]*f[0] - (f[1] - g[2])*g[0]),
  #      a*(-g[0] + g[1]*f[1] - (f[1] - g[2])*g[1]),
  #      a*(f[0] + g[1]*f[2] - (f[1] - g[2])*g[2]),
  #      K.zero(), K.zero(), K.one()]

  # This gives h of the form y^2 + ay + bx + c in 1I 10M
  # H = [K.zero(), -g[0]*a, f[0]*a, -g[1]*a , (f[1]-g[2])*a , K.one()]
  # h = [H[0] - H[3]*f[0] - H[4]*g[0],
  #      H[1] - H[3]*f[1] - H[4]*g[1],
  #      H[2] - H[3]*f[2] - H[4]*g[2],
  #      K.zero(), K.zero(), K.one()]
  return h

def val(C, P, f) :
   """
     Returns the valuation of a polynomial f at a point P (on the curve C).
     P should be an affine point.
   """
   A = C.ambient_space()
   X = C.scheme()
   Y = A.subscheme(f)
   PP = A(P[0], P[1])
   return X.intersection_multiplicity(Y, PP)
