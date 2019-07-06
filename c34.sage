import random, timeit

load("c34add.sage")
load("c34crv.sage")
load("c34crvdiv.sage")
load("c34crvpt.sage")
load("c34double.sage")
load("c34flip.sage")
load("c34slow.sage")
load("c34test.sage")
load("c34triple.sage")
load("c34util.sage")

#suite = unittest.TestLoader().loadTestsFromTestCase(TestAdd)
#unittest.TextTestRunner(verbosity=2).run(suite)
#suite = unittest.TestLoader().loadTestsFromTestCase(TestDouble)
#unittest.TextTestRunner(verbosity=2).run(suite)
suite = unittest.TestLoader().loadTestsFromTestCase(TestFlip)
unittest.TextTestRunner(verbosity=2).run(suite)

C = C_1009
K = C.K
R = C.R
x, y = R.gens()
F = C.poly()
c = C.coefficients()
print C

D = C34CrvDiv(C_1009, [[400, 243, 940, 734, 368, 690, 1],
                       [205, 314, 698,  13, 475, 495, 0, 1],
                       [100, 639, 724, 960, 184, 964, 0, 0, 1]])
f, g, h = D.polys()
DL = C34CrvDiv(C, [f,g])
DR = C34CrvDiv(C, [f,h])

f0, f1, f2, f3, f4, f5 = D.f[0:6]
g0, g1, g2, g3, g4, g5 = D.g[0:6]
h0, h1, h2, h3, h4, h5 = D.h[0:6]
c0, c1, c2, c3, c4, c5, c6, c7, c8 = c
r0 = g3 + f5*(c6 - f3)
r1 = f5
s0 = f3 - g4 + f5*(f4 - c7)
t0 = f4 - g5 + f5*(f5 - c8)
r = y + r1*x + r0
s = x + s0
t = t0

tt1 = c8 - f5
rr2 = c7 - f4
kk0 = f5*rr2 + h5
rr1 = kk0 + c6 - f3
tt0 = rr2*(c8*f5 - f4) + f5*(h5 - rr1) + c5 - h4
ss0 = f3*rr2 + f4*rr1 + f2 + h4*tt1 + h3 - c7*kk0 - c4
rr0 = c6*kk0 + c3 - h3*tt1 - f1 - f3*rr1
rr = x*x + rr2*y + rr1*x + rr0
ss = ss0
tt = y + tt1*x + tt0
kk = x + kk0



def flip_72(D) :
  c0, c1, c2, c3, c4, c5, c6, c7, c8 = D.C.coefficients()
  f0, f1, f2, f3, f4, f5 = D.f[0:6]
  g0, g1, g2, g3, g4, g5 = D.g[0:6]
  h0, h1, h2, h3, h4, h5 = D.h[0:6]
  r0 = g3 + f5*(c6 - f3)
  r1 = f5
  s0 = f3 - g4 + f5*(f4 - c7)
  t0 = f4 - g5 + f5*(f5 - c8)

  tt1 = c8 - f5
  rr2 = c7 - f4
  kk0 = f5*rr2 + h5
  rr1 = kk0 + c6 - f3
  tt0 = rr2*(c8*f5 - f4) + f5*(h5 - rr1) + c5 - h4
  ss0 = f3*rr2 + f4*rr1 + f2 + h4*tt1 + h3 - c7*kk0 - c4
  rr0 = c6*kk0 + c3 - h3*tt1 - f1 - f3*rr1
  
  z0 = g4 - h5 + g5*(tt1 - c8)
  z1 = - h3 - tt1*(h4 - g3) + g5*(g4 + s0 - c6)
  z2 = g2 + tt1*(h3 + g5*s0) - g5*(c5 + g3) - z0*h5
  
  #a1 = g3 - g5*(c7 - f4)
  #a2 = tt1*c7*h3 - (tt1*g5 - c8*g5 + g4 - h5)*(tt1*h4 + h3) + (g3*s0 + g1)*tt1 - (tt1*g5*s0 - (c5 + g3)*g5 + tt1*h3 - (tt1*g5 - c8*g5 + g4 - h5)*h5 + g2)*(c7 - f4) - (tt1*g3 - (c6 - g4 - s0)*g5 - tt1*h4 - h3)*(g4 + s0) - (c7*g3 - g4*s0 + c3 - g2)*g5 - tt1*h2 - tt0*h4 - h1
  #a5 = -(tt1*g5 - c8*g5 + g4 - h5)*tt1 + tt1*(g4 + s0) - (c7 - g5)*g5 - tt1*h5 - tt0 + g3 - h4
  a1 = g3 + g5*(f4 - c7)
  a2 = - h1 - tt0*h4 + tt1*(c7*h3 - h2 + g1 + g3*s0) + g5*(g2 + g4*s0 - c3 - c7*g3) - z0*(h3 + tt1*h4) - z1*(g4 + s0) + z2*(f4 - c7)
  a5 = g3 - h4 - tt0 + tt1*(g4 + s0 - h5 - z0) + g5*(g5 - c7)

  v0 = a1*a5 - a2
  v2 = - a5
  
  return C34CrvDiv(D.C, [[s0, 1], [v0, 0, v2, 0, 0, 1], []])

def flip_73(D) :
  c0, c1, c2, c3, c4, c5, c6, c7, c8 = D.C.coefficients()
  f0, f1, f2, f3, f4, f5 = D.f[0:6]
  g0, g1, g2, g3, g4, g5 = D.g[0:6]
  h0, h1, h2, h3, h4, h5 = D.h[0:6]
  r0 = g3 + f5*(c6 - f3)
  r1 = f5
  s0 = f3 - g4 + f5*(f4 - c7)
  t0 = f4 - g5 + f5*(f5 - c8)

  tt1 = c8 - f5
  rr2 = c7 - f4
  kk0 = f5*rr2 + h5
  rr1 = kk0 + c6 - f3
  tt0 = rr2*(c8*f5 - f4) + f5*(h5 - rr1) + c5 - h4
  ss0 = f3*rr2 + f4*rr1 + f2 + h4*tt1 + h3 - c7*kk0 - c4
  rr0 = c6*kk0 + c3 - h3*tt1 - f1 - f3*rr1
  
  z0 = h4 - g3 + g5*(c7 - f4)
  z1 = h3 + g5*(c6 - f3 + s0)
  
  a2 = - g2 - g4*s0 + g5*(c5 + c8*s0 - z0) - z1*f5
  a3 = h5 - g4 - s0 + g5*(c8 - f5)
  
  v0 = h5*a3 - a2
  v1 = - a3
  
  return C34CrvDiv(D.C, [[tt0, tt1, 1], [v0, v1, 0, 1], []])

def find_super_rare_case(C) :
  max_time = 30 # seconds
  t0 = timeit.default_timer()
  while timeit.default_timer() - t0 < max_time : 
    # Get a random atypical degree 3 divisor
    D = C.random_divisor(61, False)
    if (D.f[5] != 0) :
      continue
    # Check if <f, h> = <f, g, h>
    # If so, restart
    # If not, return
    f, g, h = D.polys()
    F = C.poly()
    R = C.R
    if not (g in R.ideal(f, h, F)) :
      if not (f in R.ideal(g, h, F)) :
        return D
  return C.zero_divisor()

"""
def f(T, typical = 2) :
  for C in [C_2, C_2_4, C_3, C_3_3, C_31, C_31_2, C_1009] :
    K = C.K
    D = C.random_divisor(T, typical)
    A = sage_flip(D)
    crv_str = "C_{}".format(K.prime_subfield().order())
    if K.degree() > 1 :
      crv_str = crv_str + "_{}".format(K.degree())
    print("    D = C34CrvDiv({}, {})".format(crv_str, [D.f, D.g, D.h]))
    print("    A = C34CrvDiv({}, {})".format(crv_str, [A.f, A.g, A.h]))
    print("    self.assertEqual(flip_{}(D), A)".format(T))
    print

def g(T) :
  for C in [C_2, C_2_4, C_3, C_3_3, C_31, C_31_2, C_1009] :
    K = C.K
    D = find_rare_case(C, T)
    A = sage_flip(D)
    crv_str = "C_{}".format(K.prime_subfield().order())
    if K.degree() > 1 :
      crv_str = crv_str + "_{}".format(K.degree())
    print("    D = C34CrvDiv({}, {})".format(crv_str, [D.f, D.g, D.h]))
    print("    A = C34CrvDiv({}, {})".format(crv_str, [A.f, A.g, A.h]))
    print("    self.assertEqual(flip_{}(D), A)".format(T))
    print

def find_rare_case(C, T) :
  # XXX : Requires T in [31, 41, 51, 61]
  D = C.zero_divisor()
  max_time = 60 # seconds
  t0 = timeit.default_timer()
  while timeit.default_timer() - t0 < max_time : 
    # Get a random atypical degree 3 divisor
    D = C.random_divisor(T, False)
    # Check if <f, h> = <f, g, h>
    # If so, restart
    # If not, return
    f, g, h = D.polys()
    F = C.poly()
    R = C.R
    if not (g in R.ideal(f, h, F)) :
      break
  return D
"""
