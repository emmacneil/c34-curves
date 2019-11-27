import random
import unittest
import timeit

load("c34testadd.sage")
load("c34testdouble.sage")
load("c34testflip.sage")
load("c34testreduce.sage")

suite = unittest.TestLoader().loadTestsFromTestCase(TestAdd)
unittest.TextTestRunner(verbosity=2).run(suite)
#suite = unittest.TestLoader().loadTestsFromTestCase(TestDouble)
#unittest.TextTestRunner(verbosity=2).run(suite)
#suite = unittest.TestLoader().loadTestsFromTestCase(TestFlip)
#unittest.TextTestRunner(verbosity=2).run(suite)
suite = unittest.TestLoader().loadTestsFromTestCase(TestReduce)
unittest.TextTestRunner(verbosity=2).run(suite)

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

def test_script() :
  MAXQ = 31
  CURVES_PER_FIELD = 10
  DIVISORS_PER_CURVE = 10
  TYPES = [11, 21, 22, 31]
  
  passes = 0
  fails = 0
  
  PP = prime_powers(MAXQ + 1)
  for q in PP :
    print("Testing over fields of order {}.".format(q))
    K = GF(q)
    for i in range(CURVES_PER_FIELD) :
      C = C34Curve.random_curve(K)
      print("  C = {}".format(C))
      for T1 in TYPES :
        for T2 in TYPES :
          print("    Adding divisors of types {} and {}.".format(T1, T2))
          for i in range(DIVISORS_PER_CURVE) :
            D1 = C.random_divisor_of_type(T1)
            D2 = C.random_divisor_of_type(T2)
            E = D1.slow_add(D2)
            G = D1 + D2
            if E == G :
              passes = passes + 1
            else :
              fails = fails + 1
  print("{} trials. {} passes. {} fails.".format(passes + fails, passes, fails))
  


def gen_add_test_case(C, type1, type2, type3) :
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
  print("    D1 = C34CurveDivisor(C, {})".format([D1.f, D1.g, D1.h]))
  print("    D2 = C34CurveDivisor(C, {})".format([D2.f, D2.g, D2.h]))
  print("    D3 = C34CurveDivisor(C, {})".format([D3.f, D3.g, D3.h]))
  print("    self.assertEqual(D1 + D2, D3)")
  return D1, D2, D3


def gen_reduce_test_cases() :
  """
    Generate a reduction unit test case.

    Given a curve C and integer type1, finds a random divisor of type type1 and
    prints out a test case.
  """
  types = [11, 21, 22, 31, 32, 33, 41, 42, 43, 44, 51, 52, 53, 54, 61, 62, 63, 64, 65]
  curves = [C_2, C_2_4, C_3, C_3_3, C_31, C_31_2, C_1009]
  strings = ["C_2", "C_2_4", "C_3", "C_3_3", "C_31", "C_31_2", "C_1009"]
  for T in types :
    print("  def test_reduce_{}(self)".format(T))
    print("    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009")
    print("    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2")
    print("    z2, z3, z4 = self.z2, self.z3, self.z4")
    print()
    for i in range(len(curves)) :
      C = curves[i]
      D = C.random_divisor_of_type(T, True)
      A = D.slow_reduce()
      print("    D = C34CurveDivisor({}, {})".format(strings[i], [D.f, D.g, D.h]))
      print("    A = C34CurveDivisor({}, {})".format(strings[i], [A.f, A.g, A.h]))
      print("    self.assertEqual(reduce(D), A)")
      print()
    print()
    print()
