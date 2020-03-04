import random
import unittest
import timeit
import datetime

load("c34testadd.sage")
load("c34testdouble.sage")
load("c34testflip.sage")
load("c34testreduce.sage")

suite = unittest.TestLoader().loadTestsFromTestCase(TestAdd)
unittest.TextTestRunner(verbosity=2).run(suite)
suite = unittest.TestLoader().loadTestsFromTestCase(TestDouble)
unittest.TextTestRunner(verbosity=2).run(suite)
suite = unittest.TestLoader().loadTestsFromTestCase(TestFlip)
unittest.TextTestRunner(verbosity=2).run(suite)
suite = unittest.TestLoader().loadTestsFromTestCase(TestReduce)
unittest.TextTestRunner(verbosity=2).run(suite)

"""
z2 = GF(31^2).gen()
z3 = GF(3^3).gen()
z4 = GF(2^4).gen()

# Curves over small primes
C_2 = C34Curve(GF(2), [0, 0, 1, 0, 1, 1, 1, 1, 1])
C_3 = C34Curve(GF(3), [0, 2, 1, 1, 0, 2, 2, 0, 2])
C_11 = C34Curve(GF(11), [6, 5, 1, 8, 6, 3, 7, 10, 4])
C_31 = C34Curve(GF(31), [2, 29, 13, 7, 28, 25, 17, 13, 17])
C_41 = C34Curve(GF(41), [26, 22, 16, 37, 30, 22, 7, 22, 29])
C_1009 = C34Curve(GF(1009), [715, 286, 748, 7, 92, 446, 122, 935, 314])

# Curves over small prime powers
C_2_4 = C34Curve(GF(2^4), [z4^3 + 1, z4, z4^3 + z4^2 + 1, z4^2, z4^3 + z4, z4^3 + z4 + 1, z4^3 + z4, z4^3 + z4^2 + z4, z4^3 + 1])
C_3_3 = C34Curve(GF(3^3), [z3 + 2, z3 + 1, z3^2 + z3 + 2, 2*z3^2 + z3 + 1, 2*z3^2 + z3, z3^2 + 2, 2*z3^2 + z3, 2*z3^2 + 1, 2])
C_31_2 = C34Curve(GF(31^2), [9*z2 + 11, 11*z2 + 28, 17*z2 + 2, 24*z2 + 22, 6*z2 + 29, 21*z2 + 8, 18*z2 + 5, 18*z2 + 15, 13*z2 + 10])
"""


def trial_c34curve_random_divisor() :
  """
    Tests random divisor generation by performing and verifying several trials.

    Specifically, tests the method C34Curve.random_divisor
    (not C34Curve.random_divisor_of_type).
    For several randomly chosen curves over small finite fields, generates
    random divisors via C34Curve.random_divisor and verifies the result.
    Verification is done calling C34CurveDivisor.formal_sum and counting the
    number of points in the divisor (counting multiplicities), and comparing
    it to the expected degree of the divisor.

    Takes about 15 minutes to run.
  """
  MAXQ = 31 # Test over all finite fields of order q <= MAXQ
  CURVES_PER_FIELD = 10 # For each finite field, try this many random curves
  DIVISORS_PER_CURVE = 100 # For each curve, and each combination of types, try adding this many divisors
  passes = 0
  fails = 0
  ret = [None]*2
  PP = prime_powers(MAXQ + 1)
  t0 = timeit.default_timer()
  for q in PP :
    print("Testing over fields of order {}.".format(q))
    K = GF(q)
    for i in range(CURVES_PER_FIELD) :
      C = C34Curve.random_curve(K)
      print("  C = {}".format(C))
      #print("    Adding divisors of types {} and {}.".format(T1, T2))
      for i in range(DIVISORS_PER_CURVE) :
        D = C.random_divisor()
        d = 0
        for (_, n) in D.formal_sum() :
          d = d + n
        if (d == D.degree) :
          passes = passes + 1
        else :
          fails = fails + 1
          ret = (C, D)
          return ret
  t1 = timeit.default_timer()
  print("{} trials. {} passes. {} fails.".format(passes + fails, passes, fails))
  print("Time taken (h:mm:ss) -- {}".format(str(datetime.timedelta(seconds = t1 - t0))))
  return ret



def trial_add(max_q = 31, curves_per_field = 10, divisors_per_curve = 100, force_short_form = False) :
  """
    Tests divisor addition by performing and verifying many random trials.

    Divisors come from many curves over many small finite fields.

    Verification is done by comparing the result of the divisor addition to
    the result of multiplying and reducing their ideals using Sage's ideal
    arithmetic.

    The optional arguments control the number of trials to perform. MAX_Q
    controls the number of fields over which to test. Random curves will be
    chosen from amongst all finite fields of order q <= MAX_Q. For each
    q, generate CURVES_PER_FIELD random curves; and for each curve, perform
    DIVISORS_PER_CURVE trials.

    Using the default values of MAX_Q = 10, CURVES_PER_FIELD = 10,
    DIVISORS_PER_CURVE = 100 --- 17000 trials in all --- takes about 3 minutes
    to run.

    The defining polynomials of the randomly chosen curves will not
    necessarily be in short form. Set force_short_form = True to make curves
    be generated in short form.
  """
  def icg_op(C, I1, I2) :
    # Ideal class group operation
    # Multiplies I1 and I2 in the ideal class group.
    # Returns a reduced representative of the ideal class [I1*I2]
    J = I1*I2 + C.defining_polynomial()
    G = list(J.groebner_basis())
    G.sort()
    Q = C.R.ideal(G[0], C.defining_polynomial())
    JJ = Q.quotient(J)
    G = list(JJ.groebner_basis())
    G.sort()
    Q = C.R.ideal(G[0], C.defining_polynomial())
    I3 = Q.quotient(JJ)
    return I3

  passes = 0
  fails = 0
  doubles = 0
  ret = [None]*3
  
  PP = prime_powers(max_q + 1)
  t0 = timeit.default_timer()
  for q in PP :
    print("Testing over fields of order {}.".format(q))
    K = GF(q)
    for i in range(curves_per_field) :
      C = C34Curve.random_curve(K)
      if (force_short_form) :
        C = C.short_form()
      print("  C = {}".format(C))
      #print("    Adding divisors of types {} and {}.".format(T1, T2))
      for i in range(divisors_per_curve) :
        D1 = C.random_divisor()
        D2 = C.random_divisor()
        if (D1 == D2) :
          doubles = doubles + 1
        D3 = D1 + D2
        I1 = D1.ideal()
        I2 = D2.ideal()
        I3 = icg_op(C, I1, I2)
        if D3.ideal() == I3 :
          passes = passes + 1
        else :
          fails = fails + 1
          ret = (C, D1, D2)
          
  t1 = timeit.default_timer()
  print("{} trials. {} passes. {} fails.".format(passes + fails, passes, fails))
  print("Tested doubling {} times".format(doubles))
  print("Time taken (h:mm:ss) -- {}".format(str(datetime.timedelta(seconds = t1 - t0))))
  return ret
  


def trial_double(max_q = 31, curves_per_field = 10, divisors_per_curve = 100, force_short_form = False) :
  """
    Tests divisor doubling by performing and verifying many random trials.

    Divisors come from many curves over many small finite fields.

    Verification is done by comparing the result of the divisor doubling to
    the result of multiplying and reducing their ideals using Sage's ideal
    arithmetic.

    The optional arguments control the number of trials to perform. MAX_Q
    controls the number of fields over which to test. Random curves will be
    chosen from amongst all finite fields of order q <= MAX_Q. For each
    q, generate CURVES_PER_FIELD random curves; and for each curve, perform
    DIVISORS_PER_CURVE trials.

    Using the default values of MAX_Q = 10, CURVES_PER_FIELD = 10,
    DIVISORS_PER_CURVE = 100 --- 17000 trials in all --- takes about 3 minutes
    to run.

    The defining polynomials of the randomly chosen curves will not
    necessarily be in short form. Set force_short_form = True to make curves
    be generated in short form.
  """
  def icg_op(C, I1, I2) :
    # Ideal class group operation
    # Multiplies I1 and I2 in the ideal class group.
    # Returns a reduced representative of the ideal class [I1*I2]
    J = I1*I2 + C.defining_polynomial()
    G = list(J.groebner_basis())
    G.sort()
    Q = C.R.ideal(G[0], C.defining_polynomial())
    JJ = Q.quotient(J)
    G = list(JJ.groebner_basis())
    G.sort()
    Q = C.R.ideal(G[0], C.defining_polynomial())
    I3 = Q.quotient(JJ)
    return I3

  passes = 0
  fails = 0
  ret = [None]*2
  
  PP = prime_powers(max_q + 1)
  t0 = timeit.default_timer()
  for q in PP :
    print("Testing over fields of order {}.".format(q))
    K = GF(q)
    for i in range(curves_per_field) :
      C = C34Curve.random_curve(K)
      if (force_short_form) :
        C = C.short_form()
      print("  C = {}".format(C))
      for i in range(divisors_per_curve) :
        D1 = C.random_divisor()
        try :
          D2 = 2*D1
        except :
          return (C, D1)
        I1 = D1.ideal()
        I2 = icg_op(C, I1, I1)
        if D2.ideal() == I2 :
          passes = passes + 1
        else :
          fails = fails + 1
          ret = (C, D1)
          return ret
 
  t1 = timeit.default_timer()
  print("{} trials. {} passes. {} fails.".format(passes + fails, passes, fails))
  print("Time taken (h:mm:ss) -- {}".format(str(datetime.timedelta(seconds = t1 - t0))))
  return ret



def gen_add_test_case(C, type1, type2, type3, cname = "C") :
  """
    Generate an addition unit test case.

    Given a curve C and integers type1, type2, type3, finds divisors D1 and
    D2 with type(D1) = type1, type(D2) = type2 and type(lcm(D1, D2)) = type3.
    Then prints out a unit test case.
  """
  TIME_LIMIT = 60 # Number of seconds after which to give up finding a case
  D1 = C.random_divisor_of_type(type1)
  D2 = C.random_divisor_of_type(type2)
  L = D1.slow_lcm(D2)
  t0 = timeit.default_timer()
  while (L.type != type3) :
    if (timeit.default_timer() - t0 > TIME_LIMIT) :
      print("Case not found within {} seconds.".format(TIME_LIMIT))
      return C.zero_divisor(), C.zero_divisor(), C.zero_divisor()
    D1 = C.random_divisor_of_type(type1)
    D2 = C.random_divisor_of_type(type2)
    L = D1.slow_lcm(D2)
  D3 = D1.slow_add(D2)
  print("    D1 = C34CurveDivisor({}, {})".format(cname, [D1.f, D1.g, D1.h]))
  print("    D2 = C34CurveDivisor({}, {})".format(cname, [D2.f, D2.g, D2.h]))
  print("    D3 = C34CurveDivisor({}, {})".format(cname, [D3.f, D3.g, D3.h]))
  print("    self.assertEqual(D1 + D2, D3)")
  return D1, D2, D3



def gen_double_test_case(C, type1, type2, cname = "C") :
  """
    Generate an doubling unit test case.

    Given a curve C and integers type1, type2, finds divisors D1 and D2 with
    type(D1) = type1, type(2*D2) = type2 and D3 is the reduction of 2*D2.
    Then prints out a unit test case.
  """
  MAX_TIME = 60
  D1 = C.random_divisor_of_type(type1)
  L = D1.slow_compose(D1)
  t0 = timeit.default_timer()
  while (L.type != type2) :
    if (timeit.default_timer() - t0 > MAX_TIME) :
      print("No case found after {} second(s).".format(MAX_TIME))
      return C.zero_divisor(), C.zero_divisor(), C.zero_divisor()
    D1 = C.random_divisor_of_type(type1)
    L = D1.slow_compose(D1)
  D2 = D1.slow_add(D1)
  print("    D1 = C34CurveDivisor({}, {})".format(cname, [D1.f, D1.g, D1.h]))
  print("    D2 = C34CurveDivisor({}, {})".format(cname, [D2.f, D2.g, D2.h]))
  print("    self.assertEqual(2*D1, D2)")
  return D1, D2



def gen_reduce_test_cases(typical = True) :
  """
    Generate an entire suite of reduction unit tests.
  """
  #types = [11, 21, 22, 31, 32, 33, 41, 42, 43, 44, 51, 52, 53, 54, 61, 62, 63, 64, 65]
  types = [41, 51, 61]
  curves = [C_2, C_2_4, C_3, C_3_3, C_31, C_31_2, C_1009]
  strings = ["C_2", "C_2_4", "C_3", "C_3_3", "C_31", "C_31_2", "C_1009"]
  for T in types :
    print("  def test_reduce_{}(self)".format(T))
    print("    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009")
    print("    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2")
    print("    z2, z3, z4 = self.z2, self.z3, self.z4")
    print("")
    for i in range(len(curves)) :
      C = curves[i]
      D = C.random_divisor_of_type(T, typical)
      A = D.slow_reduce()
      print("    D = C34CurveDivisor({}, {})".format(strings[i], [D.f, D.g, D.h]))
      print("    A = C34CurveDivisor({}, {})".format(strings[i], [A.f, A.g, A.h]))
      print("    self.assertEqual(reduce(D), A)")
      print("")
    print("")
    print("")
