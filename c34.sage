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

C = C_2
K = C.K
R = C.R
x, y = R.gens()
F = C.poly()
print C

D = C34CrvDiv(C_2, [[1, 0, 1, 1, 0, 1], [1, 0, 0, 0, 0, 0, 1], [1, 1, 1, 1, 1, 0, 0, 1]])
f, g, h = D.polys()
print D

def intersect_22_21(D1, D2) :
  if D1.type != 22 :
    raise ValueError("D1 not of expected type. D1 = {}".format(D1))
  if D2.type != 21 :
    raise ValueError("D2 not of expected type. D2 = {}".format(D2))
  p0 = D1.f[0]
  q0, q2 = D1.g[0], D1.g[2]
  r0, r1 = D2.f[0], D2.f[1]
  s0, s1 = D2.g[0], D2.g[1]
  t0 = r0 + r1*(p0 - s1)
  t1 = r1*(r1*s1 + q2 - 2*r0)
  new_f = [s0, s1, K.zero(), K.one()]
  new_g = [t0*p0, t0, p0, K.zero(), K.one()]
  new_h = [q0 + t1*p0, t1, q2, K.zero(), K.zero(), K.one()]
  
  return C34CrvDiv(D1.C, [new_f, new_g, new_h])

def test_intersect(D1, D2) :
  I1 = D1.ideal()
  I2 = D2.ideal()
  J = I1.intersection(I2)
  J = D1.R.ideal(J.groebner_basis())
  E = C34CrvDiv(D1.C, list(J.gens()))
  G = intersect_22_21(D1, D2)
  if E == G :
    print "Pass"
  else :
    print "Fail"

"""
def find_super_rare_case(C) :
  max_time = 30 # seconds
  t0 = timeit.default_timer()
  while timeit.default_timer() - t0 < max_time : 
    # Get a random atypical degree 3 divisor
    D = C.random_divisor(31, False)
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
