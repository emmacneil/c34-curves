import random, timeit

load("c34add.sage")
load("c34crv.sage")
load("c34crvdiv.sage")
load("c34crvpt.sage")
load("c34double.sage")
load("c34flip.sage")
load("c34test.sage")
load("c34triple.sage")
load("c34util.sage")

#suite = unittest.TestLoader().loadTestsFromTestCase(TestAdd)
#unittest.TextTestRunner(verbosity=2).run(suite)
#suite = unittest.TestLoader().loadTestsFromTestCase(TestDouble)
#unittest.TextTestRunner(verbosity=2).run(suite)
#suite = unittest.TestLoader().loadTestsFromTestCase(TestFlip)
#unittest.TextTestRunner(verbosity=2).run(suite)

C = C_3
K = C.K
R = C.R
x, y = R.gens()
F = C.poly()
c = C.coefficients()
print C


def test_add(C, T1, T2) :
  n_trials = 1000
  t0 = timeit.default_timer()
  for i in range(n_trials) :
    if (i > 0) and (i % 100 == 0) :
      print("{} trials passed.".format(i))
    # Get two random disjoint divisors.
    D1 = C.random_divisor(T1)
    D2 = C.random_divisor(T2)
    #while D1.slow_gcd(D2) != C.zero_divisor() :
    #  D1 = C.random_divisor(T1)
    #  D2 = C.random_divisor(T2)
    try :
      D3 = D1 + D2
    except:
      print("{} trials passed.".format(i))
      print("Done.")
      return D1, D2
    if (D1.slow_add(D2) != D3) :
      print("{} trials passed.".format(i))
      print("Done.")
      return D1, D2
  t1 = timeit.default_timer()
  print("{} trials passed in {} seconds.".format(n_trials, t1 - t0))
  print("Done.")
  return C.zero_divisor(), C.zero_divisor()

"""
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
