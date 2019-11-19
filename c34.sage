import random, timeit

load("c34add.sage")
load("c34crv.sage")
load("c34crvdiv.sage")
load("c34crvpt.sage")
load("c34double.sage")
load("c34flip.sage")
load("c34reduce.sage")
load("c34test.sage")
load("c34triple.sage")
load("c34util.sage")

#suite = unittest.TestLoader().loadTestsFromTestCase(TestAdd)
#unittest.TextTestRunner(verbosity=2).run(suite)
#suite = unittest.TestLoader().loadTestsFromTestCase(TestDouble)
#unittest.TextTestRunner(verbosity=2).run(suite)
#suite = unittest.TestLoader().loadTestsFromTestCase(TestFlip)
#unittest.TextTestRunner(verbosity=2).run(suite)

C = C34Curve.random_curve(GF(previous_prime(2^28)))
K = C.K
R = C.R
x, y = R.gens()
F = C.defining_polynomial()
print(C)



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



def test_fast_add(C, n_trials = 1000) :
  t0 = timeit.default_timer()
  i = 0
  restarts = 0
  while (i < n_trials) :
    print(i)
    D1 = C.random_divisor_of_type(31, True)
    D2 = C.random_divisor_of_type(31, True)
    D3 = D1 + D2
    L = D1.slow_lcm(D2)
    if (D3.typical) and (L.degree == 6) :
      try :
        got = fast_add(D1, D2)
      except : 
        print("{} trials passed.".format(i))
        print("Done.")
        return D1, D2
      if (got != D3) :
        print("{} trials passed.".format(i))
        print("Done.")
        return D1, D2
    else :
      restarts = restarts + 1
      continue
    i = i + 1
  t1 = timeit.default_timer()
  print("{} trials passed in {} seconds.".format(n_trials, t1 - t0))
  print("{} restarts.".format(restarts))
  print("Done.")
  return C.zero_divisor(), C.zero_divisor()


def compare_adds(C, t) :
  # Time old addition
  D1 = C.random_divisor()
  D2 = C.random_divisor()
  t0 = timeit.default_timer()
  ctr = 0
  while (timeit.default_timer() - t0 < t) :
    D3 = D1 + D2
    D1 = D2
    D2 = D3
    ctr = ctr + 1
  print("Performed {} old additions in {} seconds.".format(ctr, timeit.default_timer() - t0))
  
  t0 = timeit.default_timer()
  ctr = 0
  while (timeit.default_timer() - t0 < t) :
    if (not D1.typical) or (not D2.typical) :
      D3 = D1 + D2
    else :
      try :
        D3 = fast_add(D1, D2)
      except : 
        D3 = D1 + D2
    D1 = D2
    D2 = D3
    ctr = ctr + 1
  print("Performed {} NEW additions in {} seconds.".format(ctr, timeit.default_timer() - t0))



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



def time_random_divisor_generation(C, t) :
  """
    Prints the number of random divisors on the curve C that can be generated in t seconds.
  """
  t0 = timeit.default_timer()
  ctr = 0
  while (timeit.default_timer() - t0 < t) :
    D = C.random_divisor()
    ctr = ctr + 1
  print("Generated {} divisors in {} seconds.".format(ctr, timeit.default_timer() - t0))



def time_divisor_addition(C, t = 10, algo = fast_add_31_31_high_char, initial_seed = 0) :
  """
    Prints the number of additions that can be performed in Div(C) in t seconds.
  """
  set_random_seed(initial_seed)
  D1 = C.random_divisor()
  D2 = C.random_divisor()
  t0 = timeit.default_timer()
  ctr = 0
  restarts = 0
  while (timeit.default_timer() - t0 < t) :
    try :
      D3 = algo(D1, D2)
    except :
      D1 = C.random_divisor()
      D2 = C.random_divisor()
      restarts = restarts + 1
      continue
    D1 = D2
    D2 = D3
    ctr = ctr + 1
  print("Performed {} additions in {} seconds.".format(ctr, timeit.default_timer() - t0))
  print("{} exception(s) raised.".format(restarts))


def time_divisor_doubling(C, t = 10, algo = fast_double_31_high_char, initial_seed = 0) :
  """
    Prints the number of additions that can be performed in Div(C) in t seconds.
  """
  set_random_seed(initial_seed)
  D = C.random_divisor()
  t0 = timeit.default_timer()
  ctr = 0
  restarts = 0
  while (timeit.default_timer() - t0 < t) :
    try :
      D2 = algo(D)
    except :
      D = C.random_divisor()
      restarts = restarts + 1
      continue
    D = D2
    ctr = ctr + 1
  print("Performed {} doublings in {} seconds.".format(ctr, timeit.default_timer() - t0))
  print("{} exception(s) raised.".format(restarts))


