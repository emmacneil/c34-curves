import random, timeit

load("c34crv.sage")


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



def time_divisor_addition(C, t = 10, algo = fast_add_31_31_high_char, initial_seed = 0, verbose = False) :
  """
    Prints the number of additions that can be performed in Div(C) in t seconds.
  """
  set_random_seed(initial_seed)
  D1 = C.random_divisor()
  D2 = C.random_divisor()
  ret = [D1.type] + D1.f[0:3] + D1.g[0:3] + D1.h[0:3]
  ret = ret + [D2.type] + D2.f[0:3] + D2.g[0:3] + D2.h[0:3]
  if verbose :
    print D1
    print D2
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
  ret = ret + [ctr, restarts]
  return ret



def time_divisor_doubling(C, t = 10, algo = fast_double_31_high_char, initial_seed = 0, verbose = False) :
  """
    Prints the number of additions that can be performed in Div(C) in t seconds.
  """
  set_random_seed(initial_seed)
  D = C.random_divisor()
  ret = [D.type] + D.f[0:3] + D.g[0:3] + D.h[0:3]
  if verbose :
    print D
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
  ret = ret + [ctr, restarts]
  return ret



def timing_script(t = 10, n_primes = 20) :
  """
    Compares addition and doubling algorithms across several curves
  """
  PRIMES = [2^28 + 81, 2^28 + 105, 2^28 + 357]
  i = 0
  table = []
  total = [0, 0, 0, 0]
  for p in PRIMES :
    i = i + 1
    C = C34Curve.random_curve(GF(p))
    print("")
    print(C)
    res = C.c[0:5] + [C.c[7]]

    tmp = time_divisor_addition(C, t, algo = fast_add_31_31_high_char, verbose = True)
    if (tmp[0] != 31) or (tmp[10] != 31) or (tmp[-1] > 0) :
      continue
    res = res + tmp[1:10] + tmp[11:21]
    my_add = tmp[-2]

    tmp = time_divisor_addition(C, t, algo = km_add_31_31, verbose = True)
    if tmp[-1] > 0 :
      continue
    res = res + [tmp[-2]]
    res = res + [RR(res[-2]/res[-1])]
    km_add = tmp[-2]

    tmp = time_divisor_doubling(C, t, algo = fast_double_31_high_char, verbose = True)
    if tmp[-1] > 0 :
      continue
    res = res + [tmp[-2]]
    my_dbl = tmp[-2]

    tmp = time_divisor_doubling(C, t, algo = km_double_31, verbose = True)
    if tmp[-1] > 0 :
      continue
    res = res + [tmp[-2]]
    res = res + [RR(res[-2]/res[-1])]
    km_dbl = tmp[-2]
    total = [total[0] + my_add, total[1] + km_add, total[2] + my_dbl, total[3] + km_dbl]

    add_adv = (float(my_add/km_add) - 1.0)*100.0
    dbl_adv = (float(my_dbl/km_dbl) - 1.0)*100.0
    print("{} & $2^{{28}} - {}$ & {} & {} & {:.1f}\% & {} & {} & {:.1f}\% \\\\".format(
           i, p + 2^28, my_add, km_add, float(add_adv), my_dbl, km_dbl, float(dbl_adv)))
    table = table + [[i, p - 2^28, my_add, km_add, float(add_adv), my_dbl, km_dbl, float(dbl_adv)]]
  
  print("")
  for j in range(i) :
    t = table[j]
    print("{} & $2^{{28}} + {}$ & {} & {} & {:.1f}\% & {} & {} & {:.1f}\% \\\\".format(
         t[0], t[1], t[2], t[3], t[4], t[5], t[6], t[7]))
  ttl_add_adv = (float(total[0]/total[1]) - 1.0)*100.0
  ttl_dbl_adv = (float(total[2]/total[3]) - 1.0)*100.0
  print("&& {} & {} & {:.1f}\% & {} & {} & {:.1f}\% \\\\".format(
       total[0], total[1], float(ttl_add_adv), total[2], total[3], float(ttl_dbl_adv)))


