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

C = C_31
K = C.K
R = C.R
x, y = R.gens()
F = C.poly()
c = C.coefficients()
print C



def fib(C) :
  D1 = C.random_divisor(11)
  D2 = C.random_divisor(11)
  print D1
  print D2

  t0 = timeit.default_timer()
  t1 = timeit.default_timer()
  a = 0
  
  while (t1 - t0) < 60 :
    D3 = D1 + D2
    a = a + 1
    D1 = D2
    D2 = D3
    #print D3
    t1 = timeit.default_timer()
  print("Performed {} additions in {} seconds.".format(a, t1 - t0))


def get_lin_comb(D) :
  """
    Returns polynomials r, s, t such that f*r + g*s + h*t = C.
  """
  C = D.C
  c = C.coefficients()
  x, y = C.R.gens()

  t1 = c[8]
  r2 = c[7] - D.f[2]
  r1 = c[6] - D.f[1]
  t0 = c[5] - D.h[2] - D.f[2]*r2
  s0 = c[4] - D.h[2]*t1 - D.h[1] - D.f[2]*r1 - D.f[1]*r2
  r0 = c[3] - D.h[1]*t1 - D.f[1]*r1 - D.f[0]
  r = x^2 + r2*y + r1*x + r0
  s = s0
  t = y + t1*x + t0
  
  return r, s, t



def test_add(C, T1, T2, disjoint = False, n_trials = 1000) :
  t0 = timeit.default_timer()
  for i in range(n_trials) :
    if (i > 0) and (i % 100 == 0) :
      print("{} trials passed.".format(i))
    # Get two random disjoint divisors.
    D1 = C.random_divisor(T1)
    D2 = C.random_divisor(T2)
    if (disjoint == True) :
      while D1.slow_gcd(D2) != C.zero_divisor() :
        D1 = C.random_divisor(T1)
        D2 = C.random_divisor(T2)
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



def test_double(C, T) :
  n_trials = 1000
  t0 = timeit.default_timer()
  for i in range(n_trials) :
    if (i > 0) and (i % 100 == 0) :
      print("{} trials passed.".format(i))
    # Get a random divisors.
    D = C.random_divisor(T)
    try :
      DD = double(D)
    except:
      print("{} trials passed.".format(i))
      print("Done.")
      return D1, D2
    if (D.slow_add(D) != DD) :
      print("{} trials passed.".format(i))
      print("Done.")
      return D
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
