
class TestAdd(unittest.TestCase) :
  def setUp(self) :
    self.z2 = GF(313^2).gen()
    self.z3 = GF(3^3).gen()
    self.z4 = GF(2^4).gen()
    
    self.C_2 = C34Curve(GF(2), [0, 0, 1, 0, 1, 1, 1, 1, 1])
    self.C_3 = C34Curve(GF(3), [0, 2, 1, 1, 0, 2, 2, 0, 2])
    self.C_11 = C34Curve(GF(11), [6, 5, 1, 8, 6, 3, 7, 10, 4]) # The vertical tangent line at P = (9 : 6 : 1) intersects P with multiplicity 2
    self.C_31 = C34Curve(GF(31), [2, 29, 13, 7, 28, 25, 17, 13, 17])
    self.C_41 = C34Curve(GF(41), [26, 22, 16, 37, 30, 22, 7, 22, 29]) # Vertical tangents at (31 : 33 : 1) and (28 : 22 : 1) intersect with m = 2
    self.C_2_4 = C34Curve(GF(2^4), [z4^3 + 1, z4, z4^3 + z4^2 + 1, z4^2, z4^3 + z4, z4^3 + z4 + 1, z4^3 + z4, z4^3 + z4^2 + z4, z4^3 + 1])
    self.C_3_3 = C34Curve(GF(3^3), [z3 + 2, z3 + 1, z3^2 + z3 + 2, 2*z3^2 + z3 + 1, 2*z3^2 + z3, z3^2 + 2, 2*z3^2 + z3, 2*z3^2 + 1, 2])
    self.C_31_2 = C34Curve(GF(31^2), [9*z2 + 11, 11*z2 + 28, 17*z2 + 2, 24*z2 + 22, 6*z2 + 29, 21*z2 + 8, 18*z2 + 5, 18*z2 + 15, 13*z2 + 10])
  
  def test_add_11_11(self) :
    C_2, C_3, C_31 = self.C_2, self.C_3, self.C_31
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4
    
    # Test most common case, when sum has type 21
    D1 = C34CurveDivisor(C_2, [[0, 1], [0, 0, 1], []])
    D2 = C34CurveDivisor(C_2, [[1, 1], [0, 0, 1], []])
    D3 = C34CurveDivisor(C_2, [[0, 0, 1], [0, 1, 0, 1], []])
    self.assertEqual(add_11_11(D1, D2), D3)

    D1 = C34CurveDivisor(C_2_4, [[z4, 1], [z4 + 1, 0, 1], []])
    D2 = C34CurveDivisor(C_2_4, [[z4^2 + 1, 1], [z4^2 + 1, 0, 1], []])
    D3 = C34CurveDivisor(C_2_4, [[z4^3 + z4^2 + 1, z4^2 + z4 + 1, 1], [z4^3 + z4, z4^2 + z4 + 1, 0, 1], []])
    self.assertEqual(add_11_11(D1, D2), D3)

    D1 = C34CurveDivisor(C_3, [[1, 1], [2, 0, 1], []])
    D2 = C34CurveDivisor(C_3, [[0, 1], [1, 0, 1], []])
    D3 = C34CurveDivisor(C_3, [[1, 2, 1], [0, 1, 0, 1], []])
    self.assertEqual(add_11_11(D1, D2), D3)

    D1 = C34CurveDivisor(C_3_3, [[2*z3, 1], [2*z3 + 1, 0, 1], []])
    D2 = C34CurveDivisor(C_3_3, [[1, 1], [2*z3 + 2, 0, 1], []])
    D3 = C34CurveDivisor(C_3_3, [[z3^2 + z3 + 2, z3^2 + 2*z3, 1], [2*z3, 2*z3 + 1, 0, 1], []])
    self.assertEqual(add_11_11(D1, D2), D3)

    D1 = C34CurveDivisor(C_31, [[17, 1], [24, 0, 1], []])
    D2 = C34CurveDivisor(C_31, [[8, 1], [23, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[29, 24, 1], [12, 25, 0, 1], []])
    self.assertEqual(add_11_11(D1, D2), D3)

    D1 = C34CurveDivisor(C_31_2, [[29*z2, 1], [27*z2 + 4, 0, 1], []])
    D2 = C34CurveDivisor(C_31_2, [[20, 1], [7*z2 + 28, 0, 1], []])
    D3 = C34CurveDivisor(C_31_2, [[30*z2 + 14, 12*z2 + 21, 1], [22*z2, 29*z2 + 20, 0, 1], []])
    self.assertEqual(add_11_11(D1, D2), D3)

    # Test rarer case, when sum has type 22
    D1 = C34CurveDivisor(C_2, [[1, 1], [0, 0, 1], []])
    D2 = C34CurveDivisor(C_2, [[1, 1], [1, 0, 1], []])
    D3 = C34CurveDivisor(C_2, [[1, 1], [0, 0, 1, 0, 0, 1], []])
    self.assertEqual(add_11_11(D1, D2), D3)

    D1 = C34CurveDivisor(C_2_4, [[z4^2 + 1, 1], [1, 0, 1], []])
    D2 = C34CurveDivisor(C_2_4, [[z4^2 + 1, 1], [z4^2, 0, 1], []])
    D3 = C34CurveDivisor(C_2_4, [[z4^2 + 1, 1], [z4^2, 0, z4^2 + 1, 0, 0, 1], []])
    self.assertEqual(add_11_11(D1, D2), D3)

    D1 = C34CurveDivisor(C_3, [[2, 1], [0, 0, 1], []])
    D2 = C34CurveDivisor(C_3, [[2, 1], [2, 0, 1], []])
    D3 = C34CurveDivisor(C_3, [[2, 1], [0, 0, 2, 0, 0, 1], []])
    self.assertEqual(add_11_11(D1, D2), D3)

    D1 = C34CurveDivisor(C_3_3, [[z3^2 + 1, 1], [2*z3 + 1, 0, 1], []])
    D2 = C34CurveDivisor(C_3_3, [[z3^2 + 1, 1], [z3^2 + 2*z3, 0, 1], []])
    D3 = C34CurveDivisor(C_3_3, [[z3^2 + 1, 1], [2*z3^2 + z3 + 1, 0, z3^2 + z3 + 1, 0, 0, 1], []])
    self.assertEqual(add_11_11(D1, D2), D3)

    D1 = C34CurveDivisor(C_31, [[19, 1], [22, 0, 1], []])
    D2 = C34CurveDivisor(C_31, [[19, 1], [18, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[19, 1], [24, 0, 9, 0, 0, 1], []])
    self.assertEqual(add_11_11(D1, D2), D3)



  def test_add_21_11(self) :
    C_2, C_3, C_31 = self.C_2, self.C_3, self.C_31
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4
    
    # Test case where D1 + D2 = (P + Q) + (R) under various fields
    D1 = C34CurveDivisor(C_2, [[0, 0, 1], [0, 1, 0, 1], []])
    D2 = C34CurveDivisor(C_2, [[1, 1], [1, 0, 1], []])
    D3 = C34CurveDivisor(C_2, [[0, 1, 0, 1], [0, 0, 1, 0, 1], [0, 0, 1, 0, 0, 1]])
    self.assertEqual(add_21_11(D1, D2), D3)
    
    D1 = C34CurveDivisor(C_2_4, [[z4^3 + 1, z4^2 + 1, 1], [z4^2, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_2_4, [[z4^3 + z4^2 + z4 + 1, 1], [z4^3 + z4 + 1, 0, 1], []])
    D3 = C34CurveDivisor(C_2_4, [[1, z4^2, z4^3 + z4, 1], [z4^3 + z4 + 1, z4^3, z4^3 + z4 + 1, 0, 1], [z4^3 + z4^2 + 1, z4^3 + z4^2 + z4 + 1, z4 + 1, 0, 0, 1]])
    self.assertEqual(add_21_11(D1, D2), D3)

    D1 = C34CurveDivisor(C_3, [[1, 1, 1], [0, 2, 0, 1], []])
    D2 = C34CurveDivisor(C_3, [[0, 1], [0, 0, 1], []])
    D3 = C34CurveDivisor(C_3, [[0, 2, 0, 1], [0, 2, 0, 0, 1], [0, 1, 1, 0, 0, 1]])
    self.assertEqual(add_21_11(D1, D2), D3)

    D1 = C34CurveDivisor(C_3_3, [[z3, z3, 1], [2*z3^2 + 2*z3, 2*z3^2 + 2*z3 + 1, 0, 1], []])
    D2 = C34CurveDivisor(C_3_3, [[z3^2 + z3 + 2, 1], [2*z3 + 1, 0, 1], []])
    D3 = C34CurveDivisor(C_3_3, [[2*z3^2 + z3 + 1, 2*z3^2 + z3 + 2, 2*z3^2, 1], [1, 1, z3^2 + 2*z3 + 1, 0, 1], [2*z3^2, 2*z3^2, z3^2 + z3 + 2, 0, 0, 1]])
    self.assertEqual(add_21_11(D1, D2), D3)

    D1 = C34CurveDivisor(C_31, [[22, 9, 1], [20, 25, 0, 1], []])
    D2 = C34CurveDivisor(C_31, [[28, 1], [17, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[26, 19, 20, 1], [10, 10, 3, 0, 1], [5, 1, 12, 0, 0, 1]])
    self.assertEqual(add_21_11(D1, D2), D3)
    
    D1 = C34CurveDivisor(C_31_2, [[7*z2 + 26, 13*z2 + 23, 1], [19*z2 + 26, 15*z2 + 12, 0, 1], []])
    D2 = C34CurveDivisor(C_31_2, [[3, 1], [6*z2 + 22, 0, 1], []])
    D3 = C34CurveDivisor(C_31_2, [[30*z2 + 18, 9*z2 + 11, 7*z2 + 28, 1], [22*z2 + 28, 20*z2 + 7, 6*z2 + 4, 0, 1], [6*z2 + 9, 5*z2 + 23, 8*z2 + 4, 0, 0, 1]])
    self.assertEqual(add_21_11(D1, D2), D3)
    
    # Test other cases where points in D3 are non-distinct
    # D1 + D2 = (P + P) + (Q)
    D1 = C34CurveDivisor(C_31, [[30, 28, 1], [2, 16, 0, 1], []])
    D2 = C34CurveDivisor(C_31, [[30, 1], [25, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[27, 29, 6, 1], [20, 27, 17, 0, 1], [4, 6, 13, 0, 0, 1]])
    self.assertEqual(add_21_11(D1, D2), D3)

    # D1 + D2 = (P + Q) + (P)
    D1 = C34CurveDivisor(C_31, [[20, 30, 1], [6, 26, 0, 1], []])
    D2 = C34CurveDivisor(C_31, [[28, 1], [17, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[20, 16, 10, 1], [22, 8, 7, 0, 1], [21, 22, 13, 0, 0, 1]])
    self.assertEqual(add_21_11(D1, D2), D3)

    # D1 + D2 = (P + P) + (P)
    D1 = C34CurveDivisor(C_31, [[4, 24, 1], [8, 30, 0, 1], []])
    D2 = C34CurveDivisor(C_31, [[15, 1], [16, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[2, 25, 14, 1], [12, 12, 20, 0, 1], [24, 3, 5, 0, 0, 1]])
    self.assertEqual(add_21_11(D1, D2), D3)
    
    # Test case where D3 is of type 32
    D1 = C34CurveDivisor(C_31, [[30, 4, 1], [13, 14, 0, 1], []])
    D2 = C34CurveDivisor(C_31, [[17, 1], [24, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[30, 4, 1], [4, 3, 0, 0, 0, 0, 1], []])
    self.assertEqual(add_21_11(D1, D2), D3)

  
  
  def test_add_21_21(self) :
    C_2, C_3, C_31 = self.C_2, self.C_3, self.C_31
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4
    
    # Test case where D1 + D2 = (P + Q) + (R + S) under various fields
    # where P, Q, R, S are distinct
    D1 = C34CurveDivisor(C_2, [[0, 0, 1], [0, 1, 0, 1], []])
    D2 = C34CurveDivisor(C_2, [[0, 1, 1], [1, 1, 0, 1], []]) # D2 defined over F_2, but points defined over F_4
    D3 = C34CurveDivisor(C_2, [[0, 1, 1, 1, 1], [0, 1, 1, 1, 0, 1], [0, 0, 1, 1, 0, 0, 1]])
    self.assertEqual(add_21_21(D1, D2), D3)
    
    D1 = C34CurveDivisor(C_2_4, [[z4^2 + z4 + 1, z4^2 + z4 + 1, 1], [z4^3, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_2_4, [[z4^2 + z4 + 1, z4^3 + z4^2, 1], [z4^3 + z4 + 1, z4^3, 0, 1], []])
    D3 = C34CurveDivisor(C_2_4, [[z4, z4^3, z4^2, 0, 1], [z4 + 1, z4 + 1, z4^3 + z4, z4, 0, 1], [z4^2, z4^2 + z4, z4, z4^3 + z4^2, 0, 0, 1]])
    self.assertEqual(add_21_21(D1, D2), D3)
    
    D1 = C34CurveDivisor(C_3, [[0, 1, 1], [0, 1, 0, 1], []])
    D2 = C34CurveDivisor(C_3, [[1, 1, 1], [0, 2, 0, 1], []])
    D3 = C34CurveDivisor(C_3, [[0, 2, 0, 0, 1], [0, 0, 1, 1, 0, 1], [0, 2, 0, 0, 0, 0, 1]])
    self.assertEqual(add_21_21(D1, D2), D3)

    D1 = C34CurveDivisor(C_3_3, [[z3 + 2, 1, 1], [2*z3^2 + z3, z3^2 + 1, 0, 1], []])
    D2 = C34CurveDivisor(C_3_3, [[2*z3^2 + 2*z3, 2*z3^2 + 2*z3 + 1, 1], [2*z3^2 + z3 + 1, z3^2 + z3 + 2, 0, 1], []])
    D3 = C34CurveDivisor(C_3_3, [[2*z3^2 + z3, z3 + 1, 2*z3^2 + z3 + 1, 2*z3^2 + 2*z3 + 1, 1],
                           [2*z3^2 + 2*z3 + 1, z3^2, 2*z3^2, z3 + 1, 0, 1],
                           [2*z3^2 + 2*z3 + 1, 2*z3^2 + 2*z3, z3^2 + 2, 2, 0, 0, 1]])
    self.assertEqual(add_21_21(D1, D2), D3)

    D1 = C34CurveDivisor(C_31, [[21, 23, 1], [8, 9, 0, 1], []])
    D2 = C34CurveDivisor(C_31, [[11, 14, 1], [9, 21, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[11, 25, 1, 14, 1], [10, 25, 26, 21, 0, 1], [0, 29, 16, 29, 0, 0, 1]])
    self.assertEqual(add_21_21(D1, D2), D3)

    D1 = C34CurveDivisor(C_31_2, [[22*z2 + 12, 23*z2 + 1, 1], [20*z2 + 13, z2 + 12, 0, 1], []])
    D2 = C34CurveDivisor(C_31_2, [[16*z2 + 15, 9*z2 + 5, 1], [13*z2 + 9, 7*z2 + 4, 0, 1], []])
    D3 = C34CurveDivisor(C_31_2, [[15*z2 + 18, 5*z2 + 11, 7*z2 + 22, 14*z2 + 19, 1],
                            [3*z2 + 22, 17*z2, 22*z2 + 9, 4*z2 + 25, 0, 1],
                            [11*z2 + 2, 3*z2 + 1, 8*z2 + 12, z2 + 25, 0, 0, 1]])
    self.assertEqual(add_21_21(D1, D2), D3)    

    # Test other cases
    # D1 + D2 = (P + P) + (Q + S)
    D1 = C34CurveDivisor(C_31, [[0, 9, 1], [4, 27, 0, 1], []])
    D2 = C34CurveDivisor(C_31, [[7, 9, 1], [21, 22, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[18, 1, 9, 29, 1], [17, 14, 0, 24, 0, 1], [21, 20, 11, 9, 0, 0, 1]])
    self.assertEqual(add_21_21(D1, D2), D3)

    # D1 + D2 = (P + Q) + (Q + S)
    # XXX : This test is failing
    #       add_21_21 is producing a type 31 divisor (the reduction of D3).
    D1 = C34CurveDivisor(C_31, [[11, 14, 1], [9, 21, 0, 1], []])
    D2 = C34CurveDivisor(C_31, [[10, 21, 1], [16, 3, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[8, 24, 12, 14, 1], [16, 27, 4, 21, 0, 1], [7, 25, 8, 12, 0, 0, 1]])
    self.assertEqual(add_21_21(D1, D2), D3)

    # D1 + D2 = (P + P) + (P + S)
    D1 = C34CurveDivisor(C_31, [[10, 21, 1], [16, 23, 0, 1], []])
    D2 = C34CurveDivisor(C_31, [[15, 12, 1], [6, 10, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[5, 3, 2, 1], [6, 20, 17, 0, 2, 1], []])
    self.assertEqual(add_21_21(D1, D2), D3)

    # D1 + D2 = (P + P) + (S + S)
    D1 = C34CurveDivisor(C_31, [[4, 6, 1], [9, 6, 0, 1], []])
    D2 = C34CurveDivisor(C_31, [[17, 16, 1], [5, 19, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[19, 14, 24, 25, 1], [22, 13, 20, 11, 0, 1], [7, 26, 13, 1, 0, 0, 1]])
    self.assertEqual(add_21_21(D1, D2), D3)
    
    # Test cases where D1 + D2 is type 41 (atypical), 43, 44.
    # (Two type 21 divisors cannot sum to a type 42 divisor.)
    
    # Case where D3 is type 41 but not typical
    D1 = C34CurveDivisor(C_31, [[2, 13, 1], [4, 29, 0, 1], []])
    D2 = C34CurveDivisor(C_31, [[28, 28, 1], [13, 17, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[5, 2, 19, 28, 1], [6, 28, 26, 22, 0, 1], [10, 1, 9, 27, 0, 0, 1]])
    self.assertEqual(add_21_21(D1, D2), D3)
    
    # Case where D3 is type 43
    D1 = C34CurveDivisor(C_31, [[19, 11, 1], [19, 30, 0, 1], []])
    D2 = C34CurveDivisor(C_31, [[4, 17, 1], [18, 18, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[26, 21, 2, 1], [19, 5, 21, 0, 28, 1], []])
    self.assertEqual(add_21_21(D1, D2), D3)

    # Case where D3 is type 44
    D1 = C34CurveDivisor(C_31, [[10, 12, 1], [3, 9, 0, 1], []])
    D2 = C34CurveDivisor(C_31, [[10, 12, 1], [13, 14, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[10, 12, 1], [], []])
    self.assertEqual(add_21_21(D1, D2), D3)



  def test_add_22_11(self) :
    C_2, C_3, C_11, C_31 = self.C_2, self.C_3, self.C_11, self.C_31
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    # Test case where D1 + D2 = (P + Q) + (R) under various fields
    # where P, Q, R are distinct
    D1 = C34CurveDivisor(C_2, [[1, 1], [0, 0, 1, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_2, [[0, 1], [0, 0, 1], []])
    D3 = C34CurveDivisor(C_2, [[0, 1, 0, 1], [0, 0, 1, 0, 1], [0, 0, 1, 0, 0, 1]])
    self.assertEqual(add_22_11(D1, D2), D3)

    D1 = C34CurveDivisor(C_2_4, [[z4^3 + z4^2 + 1, 1], [z4^3, 0, 1, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_2_4, [[z4^3 + z4^2 + z4, 1], [z4^3 + z4^2 + 1, 0, 1], []])
    D3 = C34CurveDivisor(C_2_4, [[z4^3 + z4, z4 + 1, 0, 1],
                           [z4^3 + z4^2 + z4, z4^3 + z4^2 + 1, z4^3 + z4^2 + 1, 0, 1],
                           [z4^3 + z4, z4^3, 1, 0, 0, 1]])
    self.assertEqual(add_22_11(D1, D2), D3)

    D1 = C34CurveDivisor(C_3, [[0, 1], [0, 0, 1, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_3, [[1, 1], [2, 0, 1], []])
    D3 = C34CurveDivisor(C_3, [[0, 1, 0, 1], [0, 2, 0, 0, 1], [0, 2, 1, 0, 0, 1]])
    self.assertEqual(add_22_11(D1, D2), D3)
    
    D1 = C34CurveDivisor(C_3_3, [[2*z3, 1], [z3^2 + 1, 0, z3^2 + 1, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_3_3, [[2*z3^2 + z3 + 2, 1], [2*z3^2 + z3 + 1, 0, 1], []])
    D3 = C34CurveDivisor(C_3_3, [[2*z3^2 + 2*z3 + 2, 2*z3^2 + 2, 0, 1],
                       [2*z3^2 + 2, 2*z3^2 + z3 + 1, 2*z3, 0, 1],
                       [2*z3, 2*z3^2 + z3 + 2, z3^2 + 1, 0, 0, 1]])
    self.assertEqual(add_22_11(D1, D2), D3)

    D1 = C34CurveDivisor(C_31, [[1, 1], [10, 0, 24, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_31, [[12, 1], [6, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[12, 13, 0, 1], [6, 6, 1, 0, 1], [18, 8, 24, 0, 0, 1]])
    self.assertEqual(add_22_11(D1, D2), D3)
    
    D1 = C34CurveDivisor(C_31_2, [[4*z2 + 25, 1], [18*z2, 0, 3*z2 + 23, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_31_2, [[4*z2 + 8, 1], [10*z2 + 19, 0, 1], []])
    D3 = C34CurveDivisor(C_31_2, [[9*z2 + 28, 8*z2 + 2, 0, 1],
                            [3*z2 + 14, 10*z2 + 19, 4*z2 + 25, 0, 1],
                            [12*z2 + 27, 29*z2 + 15, 3*z2 + 23, 0, 0, 1]])
    self.assertEqual(add_22_11(D1, D2), D3)

    # Test other cases where points in D3 are non-distinct
    # D1 + D2 = (P + P) + (Q)
    D1 = C34CurveDivisor(C_11, [[2, 1], [3, 0, 10, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_11, [[3, 1], [6, 0, 1], []])
    D3 = C34CurveDivisor(C_11, [[6, 5, 0, 1], [1, 6, 2, 0, 1], [5, 1, 10, 0, 0, 1]])
    self.assertEqual(add_22_11(D1, D2), D3)

    # D1 + D2 = (P + Q) + (P)
    D1 = C34CurveDivisor(C_31, [[19, 1], [4, 0, 25, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_31, [[19, 1], [3, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[20, 7, 0, 1], [26, 3, 19, 0, 1], [10, 15, 25, 0, 0, 1]])
    self.assertEqual(add_22_11(D1, D2), D3)

    # D1 + D2 = (P + P) + (P)
    D1 = C34CurveDivisor(C_11, [[2, 1], [3, 0, 10, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_11, [[2, 1], [5, 0, 1], []])
    D3 = C34CurveDivisor(C_11, [[4, 4, 0, 1], [10, 5, 2, 0, 1], [1, 10, 10, 0, 0, 1]])
    self.assertEqual(add_22_11(D1, D2), D3)
    
    # Test case where D3 is of type 33
    # (TODO: Verify type 22 and 11 divisors cannot sum to a type 32?)
    D1 = C34CurveDivisor(C_31, [[16, 1], [12, 0, 30, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_31, [[16, 1], [2, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[16, 1], [], []])
    self.assertEqual(add_22_11(D1, D2), D3)



  def test_add_22_21(self) :
    C_2, C_3, C_11, C_31 = self.C_2, self.C_3, self.C_11, self.C_31
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    # Test case where D1 + D2 = (P + Q) + (R + S) under various fields
    # where P, Q, R, S are distinct
    D1 = C34CurveDivisor(C_2, [[0, 1], [1, 0, 1, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_2, [[0, 1, 1], [1, 1, 0, 1], []])
    D3 = C34CurveDivisor(C_2, [[0, 0, 0, 1, 1], [1, 0, 1, 0, 0, 1], [0, 1, 0, 1, 0, 0, 1]])
    self.assertEqual(add_22_21(D1, D2), D3)

    D1 = C34CurveDivisor(C_2_4, [[0, 1], [z4^2, 0, z4^2, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_2_4, [[z4^2 + z4 + 1, z4^3 + z4, 1], [z4^2, z4^3 + 1, 0, 1], []])
    D3 = C34CurveDivisor(C_2_4, [[0, z4^2 + z4 + 1, 0, z4^3 + z4, 1],
                           [z4^2, z4^3 + 1, z4^2, z4^2 + z4, 0, 1],
                           [0, z4^2, 0, z4^3 + 1, 0, 0, 1]])
    self.assertEqual(add_22_21(D1, D2), D3)

    D1 = C34CurveDivisor(C_3, [[1, 1], [2, 0, 1, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_3, [[1, 2, 1], [0, 2, 0, 1], []])
    D3 = C34CurveDivisor(C_3, [[1, 0, 1, 2, 1], [0, 2, 1, 1, 0, 1], [0, 2, 0, 0, 0, 0, 1]])
    self.assertEqual(add_22_21(D1, D2), D3)

    D1 = C34CurveDivisor(C_3_3, [[z3^2 + 1, 1], [2*z3^2 + z3 + 1, 0, z3^2 + z3 + 1, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_3_3, [[z3, z3^2 + 1, 1], [z3^2 + 2*z3 + 1, 2*z3^2 + 2*z3, 0, 1], []])
    D3 = C34CurveDivisor(C_3_3, [[2*z3 + 2, 1, z3^2 + 1, z3^2 + 1, 1],
                           [z3^2 + 2, 2*z3^2 + 2, z3^2 + z3 + 1, z3^2, 0, 1],
                           [2, 2*z3^2 + z3 + 2, 0, 2*z3 + 1, 0, 0, 1]])
    self.assertEqual(add_22_21(D1, D2), D3)
    
    D1 = C34CurveDivisor(C_31, [[1, 1], [1, 0, 13, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_31, [[0, 15, 1], [20, 20, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[0, 15, 1, 15, 1], [19, 28, 13, 10, 0, 1], [20, 9, 0, 21, 0, 0, 1]])
    self.assertEqual(add_22_21(D1, D2), D3)
    
    D1 = C34CurveDivisor(C_31_2, [[10*z2 + 26, 1], [4*z2 + 1, 0, 26*z2 + 1, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_31_2, [[30*z2 + 2, 27*z2 + 1, 1], [28*z2 + 11, 2*z2 + 23, 0, 1], []])
    D3 = C34CurveDivisor(C_31_2, [[5*z2 + 20, 11*z2 + 24, 10*z2 + 26, 27*z2 + 1, 1],
                            [28*z2 + 25, 26*z2 + 11, 26*z2 + 1, 6*z2 + 26, 0, 1],
                            [3*z2 + 4, 9*z2 + 22, 0, 12*z2 + 18, 0, 0, 1]])
    self.assertEqual(add_22_21(D1, D2), D3)
    
    # Test other cases where points in D3 are non-distinct
    # D1 + D2 = (P + P) + (Q + R)
    D1 = C34CurveDivisor(C_11, [[2, 1], [3, 0, 10, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_11, [[5, 10, 1], [7, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_11, [[7, 0, 0, 1], [6, 3, 2, 0, 1], []])
    self.assertEqual(add_22_21(D1, D2), D3)

    # D1 + D2 = (P + Q) + (Q + R)
    D1 = C34CurveDivisor(C_31, [[23, 1], [7, 0, 2, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_31, [[20, 28, 1], [6, 30, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[26, 13, 23, 28, 1], [18, 20, 2, 23, 0, 1], [14, 14, 0, 22, 0, 0, 1]]) # XXX: Fail
    self.assertEqual(add_22_21(D1, D2), D3)

    # D1 + D2 = (P + Q) + (R + R)
    D1 = C34CurveDivisor(C_31, [[13, 1], [2, 0, 12, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_31, [[16, 9, 1], [1, 29, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[22, 9, 13, 9, 1], [8, 17, 12, 22, 0, 1], [13, 6, 0, 11, 0, 0, 1]])
    self.assertEqual(add_22_21(D1, D2), D3)

    # D1 + D2 = (P + P) + (P + Q)
    D1 = C34CurveDivisor(C_11, [[2, 1], [3, 0, 10, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_11, [[4, 5, 1], [10, 7, 0, 1], []])
    D3 = C34CurveDivisor(C_11, [[8, 3, 2, 5, 1], [6, 4, 10, 4, 0, 1], [9, 2, 0, 9, 0, 0, 1]])
    self.assertEqual(add_22_21(D1, D2), D3)

    # D1 + D2 = (P + P) + (Q + Q)
    D1 = C34CurveDivisor(C_11, [[2, 1], [3, 0, 10, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_11, [[10, 7, 1], [0, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_11, [[9, 2, 2, 7, 1], [0, 7, 10, 7, 0, 1], [0, 0, 0, 2, 0, 0, 1]])
    self.assertEqual(add_22_21(D1, D2), D3)

    # D1 + D2 = (P + Q) + (Q + Q)
    D1 = C34CurveDivisor(C_31, [[11, 1], [26, 0, 27, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_31, [[2, 17, 1], [28, 22, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[22, 3, 11, 17, 1], [12, 29, 27, 13, 0, 1], [29, 22, 0, 2, 0, 0, 1]])
    self.assertEqual(add_22_21(D1, D2), D3)

    # Test case where D3 is of type 41 (atypical) or 42
    # TODO : Verify D3 cannot be of type 43 or 44.
    
    # Test case where D3 is of type 41 (atypical)
    D1 = C34CurveDivisor(C_31, [[19, 1], [4, 0, 25, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_31, [[17, 3, 1], [7, 4, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[13, 12, 19, 3, 1], [12, 4, 25, 22, 0, 1], [9, 21, 0, 23, 0, 0, 1]])
    self.assertEqual(add_22_21(D1, D2), D3)
    
    # Test case where D3 is of type 42
    D1 = C34CurveDivisor(C_31, [[26, 1], [11, 0, 28, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_31, [[20, 0, 1], [3, 13, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[3, 13, 0, 1], [24, 20, 26, 0, 1], []])
    self.assertEqual(add_22_21(D1, D2), D3)



  def test_add_22_22(self) :
    C_2, C_3, C_31, C_41 = self.C_2, self.C_3, self.C_31, self.C_41
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4
    
    # Test case where D1 + D2 = (P + Q) + (R + S) under various fields
    # where P, Q, R, S are distinct
    D1 = C34CurveDivisor(C_2, [[0, 1], [1, 0, 1, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_2, [[1, 1], [0, 0, 1, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_2, [[0, 1, 0, 1], [1, 1, 1, 0, 0, 1], []])
    self.assertEqual(add_22_22(D1, D2), D3)

    D1 = C34CurveDivisor(C_2_4, [[z4^3 + z4 + 1, 1], [z4^3 + z4^2 + z4, 0, 1, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_2_4, [[z4^3, 1], [z4^3, 0, z4^3 + z4^2 + z4, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_2_4, [[z4^2 + z4 + 1, z4 + 1, 0, 1], [z4^3 + z4 + 1, z4, 0, 0, z4^2 + 1, 1], []])
    self.assertEqual(add_22_22(D1, D2), D3)
    
    D1 = C34CurveDivisor(C_3, [[2, 1], [0, 0, 2, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_3, [[1, 1], [2, 0, 1, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_3, [[2, 0, 0, 1], [1, 2, 0, 0, 2, 1], []])
    self.assertEqual(add_22_22(D1, D2), D3)
    
    D1 = C34CurveDivisor(C_3_3, [[2*z3^2 + z3 + 1, 1], [1, 0, 2*z3^2 + z3 + 2, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_3_3, [[2*z3^2, 1], [2*z3^2 + 2*z3, 0, z3^2 + z3 + 1, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_3_3, [[z3 + 1, z3^2 + z3 + 1, 0, 1], [z3^2 + 2*z3 + 1, z3^2 + z3, z3^2 + 1, 0, 2*z3^2 + 1, 1], []])
    self.assertEqual(add_22_22(D1, D2), D3)
    
    D1 = C34CurveDivisor(C_31, [[25, 1], [28, 0, 14, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_31, [[27, 1], [24, 0, 30, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[24, 21, 0, 1], [16, 2, 0, 0, 23, 1], []])
    self.assertEqual(add_22_22(D1, D2), D3)

    D1 = C34CurveDivisor(C_31_2, [[20*z2 + 1, 1], [28*z2 + 7, 0, 17*z2 + 21, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_31_2, [[6*z2 + 24, 1], [6*z2 + 2, 0, 19*z2 + 27, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_31_2, [[13*z2 + 5, 26*z2 + 25, 0, 1],
                            [3*z2 + 13, 10*z2 + 17, 26, 0, 7*z2 + 22, 1], []])
    self.assertEqual(add_22_22(D1, D2), D3)

    # Test other cases where points in D3 are non-distinct
    # D1 + D2 = (P + P) + (Q + R)
    D1 = C34CurveDivisor(C_41, [[13, 1], [33, 0, 38, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_41, [[37, 1], [27, 0, 22, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_41, [[30, 9, 0, 1], [26, 31, 33, 0, 28, 1], []])

    # D1 + D2 = (P + Q) + (Q + R)
    D1 = C34CurveDivisor(C_41, [[32, 1], [5, 0, 12, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_41, [[32, 1], [35, 0, 10, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_41, [[40, 23, 0, 1], [12, 26, 32, 0, 1], []])

    # D1 + D2 = (P + P) + (P + Q)
    D1 = C34CurveDivisor(C_41, [[13, 1], [33, 0, 38, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_41, [[13, 1], [36, 0, 36, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_41, [[5, 26, 0, 1], [1, 19, 13, 0, 1], []])

    # D1 + D2 = (P + P) + (Q + Q)
    D1 = C34CurveDivisor(C_41, [[13, 1], [33, 0, 38, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_41, [[10, 1], [23, 0, 16, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_41, [[7, 23, 0, 1], [17, 24, 11, 0, 20, 1], []])

    # Test case where D3 is of type 42 or 43.
    # These cases are hit in the above tests, but test explicitly for this below.
    # TODO : Verify D3 cannot be of type 41
    # D3 cannot be of type 44.
    
    # Test case where D3 is of type 42
    D1 = C34CurveDivisor(C_31, [[1, 1], [1, 0, 13, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_31, [[1, 1], [10, 0, 24, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[1, 2, 0, 1], [29, 29, 1, 0, 1], []])
    self.assertEqual(add_22_22(D1, D2), D3)
    
    # Test case where D3 is of type 43
    D1 = C34CurveDivisor(C_31, [[15, 1], [6, 0, 2, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_31, [[23, 1], [7, 0, 2, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[4, 7, 0, 1], [8, 27, 2, 0, 0, 1], []])
    self.assertEqual(add_22_22(D1, D2), D3)



  def test_add_31_11(self) :
    C_2, C_3, C_31 = self.C_2, self.C_3, self.C_31
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    # Test case where D1 and D2 are disjoint under various fields
    D1 = C34CurveDivisor(C_2, [[1, 0, 0, 1], [0, 0, 1, 0, 1], [1, 1, 1, 0, 0, 1]])
    D2 = C34CurveDivisor(C_2, [[0, 1], [0, 0, 1], []])
    D3 = C34CurveDivisor(C_2, [[0, 0, 1, 0, 1], [0, 1, 1, 1, 0, 1], [0, 1, 0, 0, 0, 0, 1]])
    self.assertEqual(add_31_11(D1, D2), D3)
    
    D1 = C34CurveDivisor(C_2_4, [[z4^3 + z4, z4^3 + z4^2 + 1, z4^3 + z4^2 + z4 + 1, 1],
                           [z4^3 + z4^2 + z4 + 1, z4^3 + 1, 1, 0, 1],
                           [z4, z4^2 + 1, z4^3 + z4^2, 0, 0, 1]])
    D2 = C34CurveDivisor(C_2_4, [[z4, 1], [z4 + 1, 0, 1], []])
    D3 = C34CurveDivisor(C_2_4, [[z4^3 + z4, z4^2 + z4, z4^3 + z4^2 + z4 + 1, z4^3 + 1, 1],
                           [z4^2 + 1, z4^3 + z4^2, 1, z4, 0, 1],
                           [z4^3 + z4 + 1, z4^2 + z4 + 1, z4^2 + z4 + 1, 1, 0, 0, 1]])
    self.assertEqual(add_31_11(D1, D2), D3)
    
    D1 = C34CurveDivisor(C_3, [[1, 1, 0, 1], [0, 0, 2, 0, 1], [2, 1, 2, 0, 0, 1]])
    D2 = C34CurveDivisor(C_3, [[1, 1], [2, 0, 1], []])
    D3 = C34CurveDivisor(C_3, [[2, 2, 2, 2, 1], [1, 0, 2, 2, 0, 1], [1, 2, 0, 2, 0, 0, 1]])
    self.assertEqual(add_31_11(D1, D2), D3)
    
    D1 = C34CurveDivisor(C_3_3, [[z3^2 + z3 + 2, z3^2 + 1, z3^2 + 2*z3 + 2, 1],
                           [2*z3^2 + z3 + 1, z3, 1, 0, 1],
                           [z3 + 2, 2*z3^2 + z3 + 2, 2, 0, 0, 1]])
    D2 = C34CurveDivisor(C_3_3, [[2*z3^2 + z3 + 1, 1], [z3^2 + 1, 0, 1], []])
    D3 = C34CurveDivisor(C_3_3, [[z3^2 + 2*z3 + 1, 2*z3, 0, z3^2 + z3 + 1, 1],
                           [2*z3 + 2, 2*z3 + 1, 2*z3^2 + 2*z3 + 2, 2*z3 + 1, 0, 1],
                           [z3 + 1, 2*z3^2 + z3 + 1, z3, z3, 0, 0, 1]])
    self.assertEqual(add_31_11(D1, D2), D3)

    D1 = C34CurveDivisor(C_31, [[20, 4, 1, 1], [5, 19, 26, 0, 1], [25, 24, 22, 0, 0, 1]])
    D2 = C34CurveDivisor(C_31, [[29, 1], [18, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[17, 9, 8, 13, 1], [28, 6, 2, 11, 0, 1], [5, 3, 21, 20, 0, 0, 1]])
    self.assertEqual(add_31_11(D1, D2), D3)

    D1 = C34CurveDivisor(C_31_2, [[10*z2 + 22, 7*z2 + 8, 27*z2 + 4, 1],
                            [13*z2 + 1, 14*z2 + 1, 12*z2 + 7, 0, 1],
                            [8*z2 + 29, 19*z2 + 20, 28*z2 + 3, 0, 0, 1]])
    D2 = C34CurveDivisor(C_31_2, [[20*z2 + 6, 1], [6*z2 + 14, 0, 1], []])
    D3 = C34CurveDivisor(C_31_2, [[6*z2 + 3, 10*z2 + 10, 30*z2 + 14, 7*z2 + 4, 1],
                            [11*z2 + 14, 28*z2, 9*z2 + 7, 2*z2 + 26, 0, 1],
                            [6*z2 + 6, 14*z2 + 17, 10*z2 + 3, 9*z2 + 7, 0, 0, 1]])
    self.assertEqual(add_31_11(D1, D2), D3)

    # Test case where D1 is not typical
    D1 = C34CurveDivisor(C_31, [[16, 21, 0, 1], [11, 18, 23, 0, 1], [27, 13, 2, 0, 0, 1]])
    D2 = C34CurveDivisor(C_31, [[19, 1], [3, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[18, 2, 23, 14, 1], [27, 13, 2, 0, 0, 1], [25, 12, 0, 9, 0, 0, 1]])
    self.assertEqual(add_31_11(D1, D2), D3)    
    
    # Test other cases where points in D3 are non-distinct
    # D1 + D2 = (P + P + Q) + (R)
    D1 = C34CurveDivisor(C_31, [[17, 23, 8, 1], [26, 17, 19, 0, 1], [27, 10, 29, 0, 0, 1]])
    D2 = C34CurveDivisor(C_31, [[29, 1], [18, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[19, 13, 23, 16, 1], [15, 12, 27, 23, 0, 1], [0, 22, 17, 17, 0, 0, 1]])
    self.assertEqual(add_31_11(D1, D2), D3)    

    # D1 + D2 = (P + Q + R) + (R)
    D1 = C34CurveDivisor(C_31, [[10, 7, 23, 1], [16, 5, 18, 0, 1], [26, 14, 10, 0, 0, 1]])
    D2 = C34CurveDivisor(C_31, [[17, 1], [24, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[23, 13, 0, 10, 1], [1, 12, 30, 13, 0, 1], [13, 16, 19, 11, 0, 0, 1]])
    self.assertEqual(add_31_11(D1, D2), D3)

    # D1 + D2 = (P + P + P) + (Q)
    D1 = C34CurveDivisor(C_31, [[26, 11, 27, 1], [6, 25, 25, 0, 1], [18, 8, 24, 0, 0, 1]])
    D2 = C34CurveDivisor(C_31, [[12, 1], [6, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[5, 21, 18, 25, 1], [24, 1, 4, 5, 0, 1], [22, 25, 24, 30, 0, 0, 1]])
    self.assertEqual(add_31_11(D1, D2), D3)

    # D1 + D2 = (P + P + Q) + (P)
    D1 = C34CurveDivisor(C_31, [[6, 14, 20, 1], [16, 21, 2, 0, 1], [6, 23, 17, 0, 0, 1]])
    D2 = C34CurveDivisor(C_31, [[23, 1], [6, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[15, 29, 9, 5, 1], [18, 20, 26, 2, 0, 1], [24, 27, 1, 30, 0, 0, 1]])
    self.assertEqual(add_31_11(D1, D2), D3)

    # D1 + D2 = (P + P + Q) + (Q)
    D1 = C34CurveDivisor(C_31, [[30, 22, 13, 1], [13, 22, 9, 0, 1], [2, 19, 1, 0, 0, 1]])
    D2 = C34CurveDivisor(C_31, [[25, 1], [20, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[19, 14, 24, 25, 1], [22, 13, 20, 11, 0, 1], [7, 26, 13, 1, 0, 0, 1]])
    self.assertEqual(add_31_11(D1, D2), D3)

    # D1 + D2 = (P + P + P) + (P)
    D1 = C34CurveDivisor(C_31, [[22, 24, 13, 1], [7, 22, 8, 0, 1], [0, 13, 21, 0, 0, 1]])
    D2 = C34CurveDivisor(C_31, [[21, 1], [4, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[1, 7, 27, 11, 1], [3, 5, 27, 10, 0, 1], [15, 1, 15, 26, 0, 0, 1]])
    self.assertEqual(add_31_11(D1, D2), D3)

    # Test case where D3 is of type 41 (typical/atypical), 42, 43.
    # D3 cannot be type 44
    
    # D3 is type 41 but not typical
    D1 = C34CurveDivisor(C_31, [[15, 2, 21, 1], [26, 10, 20, 0, 1], [25, 26, 19, 0, 0, 1]])
    D2 = C34CurveDivisor(C_31, [[8, 1], [23, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[12, 4, 19, 28, 1], [14, 8, 16, 22, 0, 1], [23, 9, 17, 11, 0, 0, 1]])
    self.assertEqual(add_31_11(D1, D2), D3)
    
    # D3 is type 42
    D1 = C34CurveDivisor(C_31, [[29, 14, 0, 1], [8, 20, 19, 0, 1], [9, 27, 21, 0, 0, 1]])
    D2 = C34CurveDivisor(C_31, [[19, 1], [22, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[29, 14, 0, 1], [8, 20, 19, 0, 1], []])
    self.assertEqual(add_31_11(D1, D2), D3)

    # D3 is type 42
    D1 = C34CurveDivisor(C_31, [[20, 16, 10, 1], [17, 29, 25, 0, 1], [2, 15, 7, 0, 0, 1]])
    D2 = C34CurveDivisor(C_31, [[15, 1], [16, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[20, 16, 10, 1], [17, 26, 9, 0, 10, 1], []])
    self.assertEqual(add_31_11(D1, D2), D3)


    
  def test_add_31_21(self) :
    C_2, C_3, C_31 = self.C_2, self.C_3, self.C_31
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    # Test case where D1 and D2 are disjoint under various fields
    D1 = C34CurveDivisor(C_2, [[0, 0, 1, 1], [1, 1, 0, 0, 1], [0, 1, 1, 0, 0, 1]])
    D2 = C34CurveDivisor(C_2, [[0, 0, 1], [0, 1, 0, 1], []])
    D3 = C34CurveDivisor(C_2, [[0, 1, 0, 1, 0, 1], [0, 0, 1, 1, 1, 0, 1], [0, 1, 0, 1, 0, 0, 0, 1]])
    self.assertEqual(add_31_11(D1, D2), D3)

    D1 = C34CurveDivisor(C_2_4, [[z4^3 + 1, z4^3 + z4, z4^3 + z4^2 + z4 + 1, 1],
                           [z4 + 1, 1, z4 + 1, 0, 1],
                           [z4^3, 0, z4^3 + 1, 0, 0, 1]])
    D2 = C34CurveDivisor(C_2_4, [[z4^3 + z4^2 + z4 + 1, z4, 1], [0, z4^2 + 1, 0, 1], []])
    D3 = C34CurveDivisor(C_2_4, [[z4^3 + z4 + 1, z4^3 + z4 + 1, z4^3, z4^3 + z4^2 + z4, z4^3 + z4^2 + 1, 1],
                           [z4^3, z4^2 + z4 + 1, z4^3 + z4^2, z4^2 + 1, z4^3 + z4^2 + 1, 0, 1],
                           [z4^3 + z4^2 + z4, 0, z4^3 + 1, z4^2 + z4 + 1, z4^3 + 1, 0, 0, 1]])
    self.assertEqual(add_31_21(D1, D2), D3)

    D1 = C34CurveDivisor(C_3, [[0, 2, 0, 1], [0, 0, 2, 0, 1], [0, 0, 2, 0, 0, 1]])
    D2 = C34CurveDivisor(C_3, [[1, 2, 1], [0, 1, 0, 1], []])
    D3 = C34CurveDivisor(C_3, [[0, 2, 1, 1, 1, 1], [0, 2, 0, 0, 0, 0, 1], [0, 1, 0, 2, 2, 0, 0, 1]])
    self.assertEqual(add_31_21(D1, D2), D3)

    D1 = C34CurveDivisor(C_3_3, [[2*z3^2 + z3 + 1, 2*z3^2 + 2*z3, z3^2 + z3, 1],
                           [2*z3, 2*z3^2 + z3 + 1, z3^2 + z3 + 2, 0, 1],
                           [z3, z3 + 1, 2*z3^2 + 2, 0, 0, 1]])
    D2 = C34CurveDivisor(C_3_3, [[2*z3^2 + 2, 2*z3^2 + 1, 1], [2*z3^2 + z3 + 1, 2*z3^2 + 2*z3 + 2, 0, 1], []])
    D3 = C34CurveDivisor(C_3_3, [[2, z3^2 + z3 + 2, z3^2 + z3, z3^2 + 2*z3 + 2, 2*z3^2, 1],
                           [2*z3^2 + 2*z3, 2*z3^2 + z3 + 2, 2, 2*z3^2 + 2, z3^2 + z3, 0, 1],
                           [2*z3^2 + 2*z3 + 1, 2*z3^2 + z3, 2*z3 + 2, 2*z3^2 + 2*z3 + 2, 2*z3 + 2, 0, 0, 1]])
    self.assertEqual(add_31_21(D1, D2), D3)

    D1 = C34CurveDivisor(C_31, [[13, 21, 20, 1], [18, 18, 24, 0, 1], [2, 30, 30, 0, 0, 1]])
    D2 = C34CurveDivisor(C_31, [[14, 2, 1], [11, 8, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[23, 14, 25, 7, 3, 1], [4, 0, 10, 15, 28, 0, 1], [24, 9, 21, 10, 7, 0, 0, 1]])
    self.assertEqual(add_31_21(D1, D2), D3)

    D1 = C34CurveDivisor(C_31_2, [[21*z2 + 23, 29*z2 + 29, 23, 1],
                            [26*z2, 25*z2 + 20, 11*z2 + 2, 0, 1],
                            [15*z2 + 22, 17*z2 + 9, 5*z2, 0, 0, 1]])
    D2 = C34CurveDivisor(C_31_2, [[22*z2 + 12, 18*z2 + 23, 1], [10*z2 + 14, 23*z2 + 22, 0, 1], []])
    D3 = C34CurveDivisor(C_31_2, [[15*z2 + 11, 23*z2 + 21, 18*z2 + 14, 2*z2 + 6, 15*z2 + 15, 1],
                            [z2 + 16, 27*z2 + 10, 30*z2 + 27, 20*z2 + 21, 21*z2 + 10, 0, 1],
                            [29*z2 + 28, z2, 27*z2 + 15, 6, 12*z2 + 1, 0, 0, 1]])
    self.assertEqual(add_31_21(D1, D2), D3)

    # Test case where D1 is not typical
    D1 = C34CurveDivisor(C_31, [[28, 29, 0, 1], [17, 17, 1, 0, 1], [7, 6, 13, 0, 0, 1]])
    D2 = C34CurveDivisor(C_31, [[4, 17, 1], [4, 29, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[10, 6, 9, 28, 27, 1], [0, 12, 30, 13, 30, 0, 1], [19, 16, 14, 28, 15, 0, 0, 1]])
    self.assertEqual(add_31_21(D1, D2), D3)

    # Test other cases where points in D3 are non-distinct
    # D1 + D2 = (P + Q + R) + (P + S)
    D1 = C34CurveDivisor(C_31, [[18, 3, 24, 1], [13, 16, 17, 0, 1], [17, 25, 6, 0, 0, 1]])
    D2 = C34CurveDivisor(C_31, [[1, 17, 1], [9, 10, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[8, 4, 18, 5, 21, 1], [18, 21, 24, 4, 24, 0, 1], [16, 11, 8, 19, 26, 0, 0, 1]])
    self.assertEqual(add_31_21(D1, D2), D3)

    # D1 + D2 = (P + Q + R) + (S + S)
    D1 = C34CurveDivisor(C_31, [[29, 28, 25, 1], [29, 23, 22, 0, 1], [16, 29, 4, 0, 0, 1]])
    D2 = C34CurveDivisor(C_31, [[11, 14, 1], [19, 13, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[20, 15, 4, 25, 4, 1], [14, 26, 9, 10, 5, 0, 1], [26, 6, 16, 10, 22, 0, 0, 1]])
    self.assertEqual(add_31_21(D1, D2), D3)

    # D1 + D2 = (P + P + Q) + (R + S)
    D1 = C34CurveDivisor(C_31, [[23, 22, 24, 1], [4, 3, 16, 0, 1], [29, 7, 9, 0, 0, 1]])
    D2 = C34CurveDivisor(C_31, [[10, 2, 1], [29, 14, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[22, 7, 7, 6, 18, 1], [24, 28, 8, 16, 18, 0, 1], [26, 11, 17, 2, 5, 0, 0, 1]])
    self.assertEqual(add_31_21(D1, D2), D3)

    # D1 + D2 = (P + P + P) + (Q + R)
    D1 = C34CurveDivisor(C_31, [[4, 18, 30, 1], [26, 6, 7, 0, 1], [14, 15, 17, 0, 0, 1]])
    D2 = C34CurveDivisor(C_31, [[21, 5, 1], [27, 5, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[6, 0, 11, 3, 4, 1], [18, 3, 26, 17, 7, 0, 1], [0, 11, 19, 29, 13, 0, 0, 1]])
    self.assertEqual(add_31_21(D1, D2), D3)

    # D1 + D2 = (P + P + Q) + (P + R)
    D1 = C34CurveDivisor(C_31, [[3, 30, 29, 1], [1, 2, 16, 0, 1], [4, 0, 4, 0, 0, 1]])
    D2 = C34CurveDivisor(C_31, [[13, 22, 1], [9, 3, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[1, 13, 12, 15, 14, 1], [29, 16, 23, 1, 21, 0, 1], [23, 6, 22, 30, 17, 0, 0, 1]])
    self.assertEqual(add_31_21(D1, D2), D3)

    # D1 + D2 = (P + P + Q) + (Q + R)
    D1 = C34CurveDivisor(C_31, [[27, 5, 24, 1], [9, 16, 1, 0, 1], [1, 30, 26, 0, 0, 1]])
    D2 = C34CurveDivisor(C_31, [[13, 22, 1], [14, 13, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[24, 1, 26, 13, 29, 1], [4, 14, 21, 15, 22, 0, 1], [28, 17, 28, 8, 4, 0, 0, 1]])
    self.assertEqual(add_31_21(D1, D2), D3)

    # D1 + D2 = (P + P + Q) + (R + R)
    D1 = C34CurveDivisor(C_31, [[15, 20, 9, 1], [4, 19, 23, 0, 1], [2, 24, 18, 0, 0, 1]])
    D2 = C34CurveDivisor(C_31, [[9, 11, 1], [19, 18, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[16, 16, 29, 1, 23, 1], [8, 6, 24, 25, 0, 0, 1], [27, 11, 5, 7, 5, 0, 0, 1]])
    self.assertEqual(add_31_21(D1, D2), D3)

    # D1 + D2 = (P + Q + R) + (P + P)
    D1 = C34CurveDivisor(C_31, [[24, 30, 27, 1], [20, 24, 2, 0, 1], [27, 24, 1, 0, 0, 1]])
    D2 = C34CurveDivisor(C_31, [[16, 9, 1], [1, 29, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[23, 4, 30, 8, 15, 1], [2, 0, 17, 8, 7, 0, 1], [24, 30, 20, 29, 22, 0, 0, 1]])
    self.assertEqual(add_31_21(D1, D2), D3)

    # D1 + D2 = (P + Q + R) + (P + Q)
    D1 = C34CurveDivisor(C_31, [[1, 23, 1, 1], [16, 21, 29, 0, 1], [24, 4, 10, 0, 0, 1]])
    D2 = C34CurveDivisor(C_31, [[26, 30, 1], [6, 24, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[6, 27, 1, 21, 15, 1], [20, 6, 19, 18, 20, 0, 1], [9, 18, 14, 0, 27, 0, 0, 1]])
    self.assertEqual(add_31_21(D1, D2), D3)

    # D1 + D2 = (P + P + P) + (P + Q)
    D1 = C34CurveDivisor(C_31, [[2, 25, 14, 1], [12, 12, 20, 0, 1], [24, 3, 5, 0, 0, 1]])
    D2 = C34CurveDivisor(C_31, [[26, 11, 1], [15, 16, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[5, 27, 9, 14, 9, 1], [0, 20, 7, 15, 26, 0, 1], [13, 12, 8, 2, 15, 0, 0, 1]])
    self.assertEqual(add_31_21(D1, D2), D3)

    # D1 + D2 = (P + P + P) + (Q + Q)
    D1 = C34CurveDivisor(C_31, [[2, 25, 14, 1], [12, 12, 20, 0, 1], [24, 3, 5, 0, 0, 1]])
    D2 = C34CurveDivisor(C_31, [[3, 17, 1], [10, 3, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[28, 7, 22, 0, 21, 1], [17, 9, 21, 3, 1, 0, 1], [7, 13, 30, 5, 14, 0, 0, 1]])
    self.assertEqual(add_31_21(D1, D2), D3)

    # D1 + D2 = (P + P + Q) + (P + P)
    D1 = C34CurveDivisor(C_31, [[4, 6, 10, 1], [4, 0, 22, 0, 1], [25, 12, 17, 0, 0, 1]])
    D2 = C34CurveDivisor(C_31, [[17, 3, 1], [20, 7, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[8, 0, 25, 29, 21, 1], [23, 19, 19, 24, 21, 0, 1], [12, 26, 22, 14, 11, 0, 0, 1]])
    self.assertEqual(add_31_21(D1, D2), D3)

    # D1 + D2 = (P + P + Q) + (P + Q)
    D1 = C34CurveDivisor(C_31, [[7, 9, 13, 1], [16, 16, 20, 0, 1], [15, 21, 12, 0, 0, 1]])
    D2 = C34CurveDivisor(C_31, [[14, 24, 1], [11, 7, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[2, 14, 24, 3, 25, 1], [24, 12, 9, 15, 8, 0, 1], [23, 10, 2, 17, 21, 0, 0, 1]])
    self.assertEqual(add_31_21(D1, D2), D3)

    # D1 + D2 = (P + P + Q) + (Q + Q)
    D1 = C34CurveDivisor(C_31, [[27, 27, 16, 1], [28, 20, 4, 0, 1], [9, 11, 14, 0, 0, 1]])
    D2 = C34CurveDivisor(C_31, [[29, 27, 1], [9, 25, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[22, 5, 5, 25, 14, 1], [3, 6, 26, 24, 19, 0, 1], [9, 30, 19, 20, 1, 0, 0, 1]])
    self.assertEqual(add_31_21(D1, D2), D3)

    # D1 + D2 = (P + P + P) + (P + P)
    D1 = C34CurveDivisor(C_31, [[0, 7, 7, 1], [22, 26, 6, 0, 1], [19, 28, 3, 0, 0, 1]])
    D2 = C34CurveDivisor(C_31, [[17, 16, 1], [5, 19, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[22, 11, 22, 20, 1], [22, 10, 0, 6, 0, 6, 1], []])
    self.assertEqual(add_31_21(D1, D2), D3)
    
    # TODO: Test case where D3 is of type 51 (typical/atypical), 52, 53, 54.

    self.fail("Test not implemented.")

  def test_add_31_22(self) :
    # Test case where D1 + D2 = (P + Q + R) + (S + T) under various fields
    # where P, Q, R, S, T are distinct
    
    # Test case where D1 is not typical
    
    # Test other cases where points in D3 are non-distinct
    # D1 + D2 = (P + Q + R) + (P + S)
    # D1 + D2 = (P + Q + R) + (S + S)
    # D1 + D2 = (P + P + Q) + (R + S)
    # D1 + D2 = (P + P + P) + (Q + R)
    # D1 + D2 = (P + P + Q) + (P + R)
    # D1 + D2 = (P + P + Q) + (Q + R)
    # D1 + D2 = (P + P + Q) + (R + R)
    # D1 + D2 = (P + Q + R) + (P + P)
    # D1 + D2 = (P + Q + R) + (P + Q)
    # D1 + D2 = (P + P + P) + (P + Q)
    # D1 + D2 = (P + P + P) + (Q + Q)
    # D1 + D2 = (P + P + Q) + (P + P)
    # D1 + D2 = (P + P + Q) + (P + Q)
    # D1 + D2 = (P + P + Q) + (Q + Q)
    # D1 + D2 = (P + P + P) + (P + P)
    
    # Test case where D3 is of type 51 (typical/atypical), 52, 53, 54.

    self.fail("Test not implemented.")

  def test_add_31_31(self) :
    # Test case where D1 + D2 = (P + Q + R) + (S + T + U) under various fields
    # where P, Q, R, S, T, U are distinct
    
    # Test case where D1 or D2 is not typical
    
    # Test other cases where points in D3 are non-distinct
    # D1 + D2 = (P + P + Q) + (R + S + T)
    # D1 + D2 = (P + Q + R) + (R + S + T)
    # D1 + D2 = (P + P + P) + (Q + R + S)
    # D1 + D2 = (P + P + Q) + (R + R + S)
    # D1 + D2 = (P + P + Q) + (Q + R + S)
    # D1 + D2 = (P + P + Q) + (P + R + S)
    # D1 + D2 = (P + P + Q) + (P + Q + S)
    # D1 + D2 = (P + Q + R) + (P + Q + S)
    # D1 + D2 = (P + Q + R) + (P + P + S)
    # D1 + D2 = (P + P + P) + (Q + Q + R)
    # D1 + D2 = (P + P + P) + (P + Q + R)
    # D1 + D2 = (P + P + Q) + (P + R + R)
    # D1 + D2 = (P + P + Q) + (Q + R + R)
    # D1 + D2 = (P + P + Q) + (P + P + R)
    # D1 + D2 = (P + P + Q) + (P + Q + R)
    # D1 + D2 = (P + P + Q) + (Q + Q + R)
    # D1 + D2 = (P + P + P) + (P + P + Q)
    # D1 + D2 = (P + P + P) + (P + Q + Q)
    # D1 + D2 = (P + P + P) + (Q + Q + Q)
    # D1 + D2 = (P + P + Q) + (P + Q + Q)
    
    # Test case where D3 is of type 51 (typical/atypical), 52, 53, 54.

    self.fail("Test not implemented.")

