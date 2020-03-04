import random, timeit

def time_divisor_addition(C, t = 10, algo = fast_add_31_31_high_char, initial_seed = 0, verbose = True) :
  """
    Prints the number of additions that can be performed in Div(C) in t seconds.

    Does so by generating random divisors D_1, D_2, then computing the sequence D_{i+2} = D_i + D_{i+1}.
  """
  set_random_seed(initial_seed)

  D1 = C.random_divisor()
  D2 = C.random_divisor()
  ret = [D1.type] + D1.f[0:3] + D1.g[0:3] + D1.h[0:3] + [D2.type] + D2.f[0:3] + D2.g[0:3] + D2.h[0:3]
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
  if verbose :
    print("Performed {} additions in {} seconds.".format(ctr, timeit.default_timer() - t0))
    print("{} exception(s) raised.".format(restarts))
  ret = ret + [ctr, restarts]
  return ret



def time_divisor_doubling(C, t = 10, algo = faster_double_31_high_char, initial_seed = 0, verbose = True) :
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
  if verbose :
    print("Performed {} doublings in {} seconds.".format(ctr, timeit.default_timer() - t0))
    print("{} exception(s) raised.".format(restarts))
  ret = ret + [ctr, restarts]
  return ret



def timing_script(first_p = 2^28, t = 10, n_primes = 3, filename = "out.csv") :
  """
    Compares addition and doubling algorithms across several curves
  """
  fout = open(filename, 'w')
  header = "p, c7, c4, c3, c2, c1, c0, T1, f2, f1, f0, g2, g1, g0, h2, h1, h0, T2, F2, F1, F0, G2, G1, G0, H2, H1, H0, EMM adds, KM adds, EMM doubles, KM doubles\n"
  fout.write(header)
  p = first_p
  v = False
  for i in range(n_primes) :
    p = next_prime(p)
    print("p = {}".format(p))

    set_random_seed() # Picks a random seed, for curve generation. The time_divisor_XXX subroutines reset the seed to 0 so that they always add the same sequences.
    C = C34Curve.random_curve(GF(p))

    line = "{}, {}, {}, {}, {}, {}, {}, ".format(p, C.c[7], C.c[4], C.c[3], C.c[2], C.c[1], C.c[0])

    # Khuri-Makdisi's algorithms are more likely to encounter exceptions, so run those first
    # KM-Addition
    tmp = time_divisor_addition(C, t, algo = km_add_31_31, verbose = v)
    T1, f0, f1, f2, g0, g1, g2, h0, h1, h2, T2, F0, F1, F2, G0, G1, G2, H0, H1, H2, km_adds, r = tmp
    line = line + "{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, ".format(T1, f2, f1, f0, g2, g1, g0, h2, h1, h0)
    line = line + "{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, ".format(T2, F2, F1, F0, G2, G1, G0, H2, H1, H0)
    if r > 0 :
      p = p - 1
      i = i - 1
      continue

    # KM-Doubling
    tmp = time_divisor_doubling(C, t, algo = km_double_31, verbose = v)
    km_dubs, r = tmp[-2], tmp[-1]
    if r > 0 :
      p = p - 1
      i = i - 1
      continue

    # My addition
    tmp = time_divisor_addition(C, t, algo = fast_add_31_31_high_char, verbose = v)
    emm_adds, r = tmp[-2], tmp[-1]
    if r > 0 :
      p = p - 1
      i = i - 1
      continue

    # My doubling
    tmp = time_divisor_doubling(C, t, algo = faster_double_31_high_char, verbose = v)
    emm_dubs, r = tmp[-2], tmp[-1]
    if r > 0 :
      p = p - 1
      i = i - 1
      continue

    line = line + "{}, {}, {}, {}\n".format(emm_adds, km_adds, emm_dubs, km_dubs)

    fout.write(line)
  fout.close()
