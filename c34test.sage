import random
import unittest

load("c34testadd.sage")
load("c34testdouble.sage")
load("c34testflip.sage")

z2 = GF(31^2).gen()
z3 = GF(3^3).gen()
z4 = GF(2^4).gen()

# Curves over small primes
C_2 = C34Crv(GF(2), [0, 0, 1, 0, 1, 1, 1, 1, 1])
C_3 = C34Crv(GF(3), [0, 2, 1, 1, 0, 2, 2, 0, 2])
C_11 = C34Crv(GF(11), [6, 5, 1, 8, 6, 3, 7, 10, 4])
C_31 = C34Crv(GF(31), [2, 29, 13, 7, 28, 25, 17, 13, 17])
C_41 = C34Crv(GF(41), [26, 22, 16, 37, 30, 22, 7, 22, 29])
C_1009 = C34Crv(GF(1009), [715, 286, 748, 7, 92, 446, 122, 935, 314])

# Curves over small prime powers
C_2_4 = C34Crv(GF(2^4), [z4^3 + 1, z4, z4^3 + z4^2 + 1, z4^2, z4^3 + z4, z4^3 + z4 + 1, z4^3 + z4, z4^3 + z4^2 + z4, z4^3 + 1])
C_3_3 = C34Crv(GF(3^3), [z3 + 2, z3 + 1, z3^2 + z3 + 2, 2*z3^2 + z3 + 1, 2*z3^2 + z3, z3^2 + 2, 2*z3^2 + z3, 2*z3^2 + 1, 2])
C_31_2 = C34Crv(GF(31^2), [9*z2 + 11, 11*z2 + 28, 17*z2 + 2, 24*z2 + 22, 6*z2 + 29, 21*z2 + 8, 18*z2 + 5, 18*z2 + 15, 13*z2 + 10])

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
      C = C34Crv.random_curve(K)
      print ("  C = {}".format(C))
      for T1 in TYPES :
        for T2 in TYPES :
          print("    Adding divisors of types {} and {}.".format(T1, T2))
          for i in range(DIVISORS_PER_CURVE) :
            D1 = C.random_divisor(T1)
            D2 = C.random_divisor(T2)
            E = D1.slow_add(D2)
            G = D1 + D2
            if E == G :
              passes = passes + 1
            else :
              fails = fails + 1
  print("{} trials. {} passes. {} fails.".format(passes + fails, passes, fails))
  


def gen_add_test_case(C, type1, type2) :
  """
    Generate a generic addition test case.
  """
  D1 = find_reduced_divisor(C, type1)
  D2 = find_reduced_divisor(C, type2)
  D3 = sage_compose(D1, D2)
  print("D1 = C34CrvDiv(C, {})".format([D1.f, D1.g, D1.h]))
  print("D2 = C34CrvDiv(C, {})".format([D2.f, D2.g, D2.h]))
  print("D3 = C34CrvDiv(C, {})".format([D3.f, D3.g, D3.h]))
  print("self.assertEqual(add_{}_{}(D1, D2), D3)".format(type1, type2))
  return D1, D2, D3


