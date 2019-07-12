class TestDouble(unittest.TestCase) :
  def setUp(self) :
    self.z2 = GF(313^2).gen()
    self.z3 = GF(3^3).gen()
    self.z4 = GF(2^4).gen()
    
    self.C_2 = C34Crv(GF(2), [0, 0, 1, 0, 1, 1, 1, 1, 1])
    self.C_3 = C34Crv(GF(3), [0, 2, 1, 1, 0, 2, 2, 0, 2])
    self.C_11 = C34Crv(GF(11), [6, 5, 1, 8, 6, 3, 7, 10, 4]) # The vertical tangent line at P = (9 : 6 : 1) intersects P with multiplicity 2
    self.C_31 = C34Crv(GF(31), [2, 29, 13, 7, 28, 25, 17, 13, 17])
    self.C_41 = C34Crv(GF(41), [26, 22, 16, 37, 30, 22, 7, 22, 29]) # Vertical tangents at (31 : 33 : 1) and (28 : 22 : 1) intersect with m = 2
    self.C_1009 = C34Crv(GF(1009), [715, 286, 748, 7, 92, 446, 122, 935, 314])
    self.C_2_4 = C34Crv(GF(2^4), [z4^3 + 1, z4, z4^3 + z4^2 + 1, z4^2, z4^3 + z4, z4^3 + z4 + 1, z4^3 + z4, z4^3 + z4^2 + z4, z4^3 + 1])
    self.C_3_3 = C34Crv(GF(3^3), [z3 + 2, z3 + 1, z3^2 + z3 + 2, 2*z3^2 + z3 + 1, 2*z3^2 + z3, z3^2 + 2, 2*z3^2 + z3, 2*z3^2 + 1, 2])
    self.C_31_2 = C34Crv(GF(31^2), [9*z2 + 11, 11*z2 + 28, 17*z2 + 2, 24*z2 + 22, 6*z2 + 29, 21*z2 + 8, 18*z2 + 5, 18*z2 + 15, 13*z2 + 10])
  


  def test_double_0(self) :
    C_31 = self.C_31
    D1 = C34CrvDiv(C_31, [[1], [], []])
    D2 = C34CrvDiv(C_31, [[1], [], []])
    self.assertEqual(double_0(D1), D2)



  def test_double_11(self) :
    C_2, C_3, C_31, C_41, C_1009 = self.C_2, self.C_3, self.C_31, self.C_41, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4
    
    # Test cases where 2*D1 is type 21
    D1 = C34CrvDiv(C_2, [[1, 1], [0, 0, 1], []])
    D2 = C34CrvDiv(C_2, [[1, 1, 1], [1, 0, 0, 1], []])
    self.assertEqual(double_11(D1), D2)

    D1 = C34CrvDiv(C_2_4, [[z4, 1], [z4 + 1, 0, 1], []])
    D2 = C34CrvDiv(C_2_4, [[z4^3 + 1, z4^2 + 1, 1], [z4^2, 0, 0, 1], []])
    self.assertEqual(double_11(D1), D2)
    
    D1 = C34CrvDiv(C_3, [[1, 1], [2, 0, 1], []])
    D2 = C34CrvDiv(C_3, [[0, 1, 1], [1, 2, 0, 1], []])
    self.assertEqual(double_11(D1), D2)

    D1 = C34CrvDiv(C_3_3, [[z3^2 + z3 + 2, 1], [2*z3 + 1, 0, 1], []])
    D2 = C34CrvDiv(C_3_3, [[0, z3^2 + 2*z3 + 1, 1], [2*z3 + 2, 2*z3^2 + 2*z3 + 1, 0, 1], []])
    self.assertEqual(double_11(D1), D2)
    
    D1 = C34CrvDiv(C_31, [[18, 1], [20, 0, 1], []])
    D2 = C34CrvDiv(C_31, [[1, 11, 1], [14, 5, 0, 1], []])
    self.assertEqual(double_11(D1), D2)

    D1 = C34CrvDiv(C_1009, [[306, 1], [347, 0, 1], []])
    D2 = C34CrvDiv(C_1009, [[945, 701, 1], [808, 612, 0, 1], []])
    self.assertEqual(double_11(D1), D2)

    # Test cases where 2*D1 is type 22
    # (No type 11 divisor on C_2_4, C_3_3, C_31 doubles to a type 22)
    D1 = C34CrvDiv(C_2, [[1, 1], [1, 0, 1], []])
    D2 = C34CrvDiv(C_2, [[1, 1], [1, 0, 0, 0, 0, 1], []])
    self.assertEqual(double_11(D1), D2)

    D1 = C34CrvDiv(C_3, [[2, 1], [2, 0, 1], []])
    D2 = C34CrvDiv(C_3, [[2, 1], [1, 0, 1, 0, 0, 1], []])
    self.assertEqual(double_11(D1), D2)

    D1 = C34CrvDiv(C_41, [[13, 1], [19, 0, 1], []])
    D2 = C34CrvDiv(C_41, [[13, 1], [33, 0, 38, 0, 0, 1], []])
    self.assertEqual(double_11(D1), D2)

    D1 = C34CrvDiv(C_1009, [[647, 1], [26, 0, 1], []])
    D2 = C34CrvDiv(C_1009, [[647, 1], [676, 0, 52, 0, 0, 1], []])
    self.assertEqual(double_11(D1), D2)



  def test_double_21(self) :
    C_2, C_3, C_31, C_41, C_1009 = self.C_2, self.C_3, self.C_31, self.C_41, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    # Test typical cases where 2*D1 is type 41
    D1 = C34CrvDiv(C_2, [[0, 0, 1], [0, 1, 0, 1], []])
    D2 = C34CrvDiv(C_2, [[0, 0, 1, 0, 1], [0, 0, 0, 0, 0, 1], [0, 0, 1, 1, 0, 0, 1]]) # Atypical divisor
    self.assertEqual(double_21(D1), D2)

    D1 = C34CrvDiv(C_2_4, [[z4^3 + z4^2 + z4 + 1, z4^2 + z4, 1], [0, z4, 0, 1], []])
    D2 = C34CrvDiv(C_2_4, [[z4^2 + 1, z4^3 + z4, z4^3 + z4^2 + z4, z4^3 + z4^2, 1],
                           [z4^3 + z4, 0, 0, z4^2 + z4 + 1, 0, 1],
                           [z4^2 + z4, z4^2 + z4, z4^2 + 1, z4^3 + z4^2 + z4, 0, 0, 1]])
    self.assertEqual(double_21(D1), D2)
    
    D1 = C34CrvDiv(C_3, [[0, 1, 1], [0, 1, 0, 1], []])
    D2 = C34CrvDiv(C_3, [[0, 0, 0, 1, 1], [0, 0, 0, 2, 0, 1], [0, 2, 1, 2, 0, 0, 1]])
    self.assertEqual(double_21(D1), D2)

    D1 = C34CrvDiv(C_3_3, [[z3, z3^2, 1], [2*z3^2 + 2*z3 + 2, 2*z3 + 2, 0, 1], []])
    D2 = C34CrvDiv(C_3_3, [[2*z3 + 2, 2*z3^2, z3^2 + 1, z3^2, 1],
                           [2*z3 + 1, 2*z3^2 + 1, 2*z3^2 + z3, 2*z3^2 + z3, 0, 1],
                           [z3^2, 2*z3 + 1, 2*z3 + 2, 2*z3^2 + z3, 0, 0, 1]])
    self.assertEqual(double_21(D1), D2)
    
    D1 = C34CrvDiv(C_31, [[15, 17, 1], [29, 30, 0, 1], []])
    D2 = C34CrvDiv(C_31, [[7, 30, 11, 3, 1], [18, 17, 28, 1, 0, 1], [26, 4, 27, 18, 0, 0, 1]])
    self.assertEqual(double_21(D1), D2)

    D1 = C34CrvDiv(C_31_2, [[16*z2 + 14, 20*z2 + 28, 1], [24*z2 + 27, 21*z2 + 12, 0, 1], []])
    D2 = C34CrvDiv(C_31_2, [[7*z2 + 16, 12*z2 + 4, 14*z2 + 13, 9*z2 + 22, 1],
                            [19*z2 + 23, 14*z2 + 18, 26*z2 + 19, 2*z2 + 21, 0, 1],
                            [14*z2 + 1, 9*z2 + 6, 10*z2 + 28, 28*z2 + 11, 0, 0, 1]])
    self.assertEqual(double_21(D1), D2)

    D1 = C34CrvDiv(C_1009, [[487, 900, 1], [388, 506, 0, 1], []])
    D2 = C34CrvDiv(C_1009, [[620, 784, 870, 751, 1], [8, 170, 942, 33, 0, 1], [907, 762, 584, 142, 0, 0, 1]])
    self.assertEqual(double_21(D1), D2)

    # Test cases where 2*D1 is type 42
    # TODO : Verify this case even exists
    
    # Test cases where 2*D1 is type 43
    D1 = C34CrvDiv(C_41, [[26, 10, 1], [7, 23, 0, 1], []])
    D2 = C34CrvDiv(C_41, [[7, 23, 0, 1], [17, 24, 11, 0, 20, 1], []])
    self.assertEqual(double_21(D1), D2)

    # Test cases where 2*D1 is type 44
    C = C34Crv(GF(31), [4, 26, 5, 4, 24, 30, 3, 1, 26])
    D1 = C34CrvDiv(C, [[0, 30, 1], [29, 0, 0, 1], []]) # D1 = P + Q
    D2 = C34CrvDiv(C, [[0, 30, 1], [], []])            # D2 = 2P + 2Q
    self.assertEqual(double_21(D1), D2)

    C = C34Crv(GF(1009), [105, 392, 370, 162, 517, 810, 26, 260, 148])
    D1 = C34CrvDiv(C, [[666, 460, 1], [671, 588, 0, 1], []]) # D1 = 2P
    D2 = C34CrvDiv(C, [[666, 460, 1], [], []])               # D2 = 4P
    self.assertEqual(double_21(D1), D2)

    # Test cases where D1 = 2*P
    D1 = C34CrvDiv(C_1009, [[227, 543, 1], [40, 361, 0, 1], []])
    D2 = C34CrvDiv(C_1009, [[323, 529, 285, 649, 1], [424, 962, 707, 698, 0, 1], [184, 242, 223, 437, 0, 0, 1]])
    self.assertEqual(double_21(D1), D2)



  def test_double_22(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    # Test most common case where 2*D1 is type 43
    D1 = C34CrvDiv(C_2, [[0, 1], [1, 0, 1, 0, 0, 1], []])
    D2 = C34CrvDiv(C_2, [[0, 0, 0, 1], [1, 1, 1, 0, 1, 1], []])
    self.assertEqual(double_22(D1), D2)

    D1 = C34CrvDiv(C_2_4, [[z4^3 + z4^2 + z4, 1], [z4^3 + z4^2 + z4, 0, 1, 0, 0, 1], []])
    D2 = C34CrvDiv(C_2_4, [[z4^3 + z4 + 1, 0, 0, 1], [z4^3 + z4, z4^3 + z4^2, z4^2 + 1, 0, z4^3 + z4^2, 1], []])
    self.assertEqual(double_22(D1), D2)

    D1 = C34CrvDiv(C_3, [[2, 1], [0, 0, 2, 0, 0, 1], []])
    D2 = C34CrvDiv(C_3, [[1, 1, 0, 1], [0, 0, 2, 0, 1], []])
    self.assertEqual(double_22(D1), D2)
    
    D1 = C34CrvDiv(C_3_3, [[2*z3^2 + z3 + 2, 1], [z3^2 + z3 + 2, 0, z3^2, 0, 0, 1], []])
    D2 = C34CrvDiv(C_3_3, [[z3^2 + z3, z3^2 + 2*z3 + 1, 0, 1], [2*z3^2, 2, z3 + 1, 0, 2*z3^2, 1], []])
    self.assertEqual(double_22(D1), D2)

    D1 = C34CrvDiv(C_31, [[8, 1], [14, 0, 21, 0, 0, 1], []])
    D2 = C34CrvDiv(C_31, [[2, 16, 0, 1], [17, 12, 26, 0, 20, 1], []])
    self.assertEqual(double_22(D1), D2)

    D1 = C34CrvDiv(C_31_2, [[7*z2 + 6, 1], [11*z2 + 15, 0, 25*z2 + 22, 0, 0, 1], []])
    D2 = C34CrvDiv(C_31_2, [[27*z2 + 13, 14*z2 + 12, 0, 1],
                            [6*z2 + 29, 21*z2 + 19, 5*z2 + 29, 0, 27*z2 + 13, 1], []])
    self.assertEqual(double_22(D1), D2)

    D1 = C34CrvDiv(C_1009, [[172, 1], [426, 0, 779, 0, 0, 1], []])
    D2 = C34CrvDiv(C_1009, [[323, 344, 0, 1], [478, 540, 241, 0, 79, 1], []])
    self.assertEqual(double_22(D1), D2)

    # Test a case where 2*D1 is type 42
    # Here, P = (42 : 42 : 1) is an inflection point and x = 42 is the tangent line at P.
    # D1 = 2P and D2 = 4P
    C = C34Crv(GF(43), [38, 15, 36, 13, 14, 4, 36, 24, 1])
    D1 = C34CrvDiv(C, [[1, 1], [1, 0, 2, 0, 0, 1], []])
    D2 = C34CrvDiv(C, [[1, 2, 0, 1], [1, 1, 1, 0, 1], []])
    self.assertEqual(double_22(D1), D2)
    
    # Test another case where D1 = 2P.
    # D2 = 4P has type 43
    D1 = C34CrvDiv(C_1009, [[647, 1], [676, 0, 52, 0, 0, 1], []])
    D2 = C34CrvDiv(C_1009, [[883, 285, 0, 1], [927, 593, 877, 0, 112, 1], []])
    self.assertEqual(double_22(D1), D2)



  def test_double_31(self) :
    C_2, C_3, C_11, C_31, C_1009 = self.C_2, self.C_3, self.C_11, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    # Test typical case over various fields
    D1 = C34CrvDiv(C_2, [[0, 0, 1, 1], [0, 0, 1, 0, 1], [0, 0, 1, 0, 0, 1]])
    D2 = C34CrvDiv(C_2, [[0, 0, 1, 0, 1, 1, 1], [0, 0, 0, 0, 1, 0, 0, 1], [0, 0, 0, 0, 0, 1, 0, 0, 1]])
    self.assertEqual(double_31(D1), D2)

    D1 = C34CrvDiv(C_2, [[1, 0, 1, 1], [1, 1, 0, 0, 1], [1, 1, 0, 0, 0, 1]])
    D2 = C34CrvDiv(C_2, [[1, 0, 1, 1, 1, 1], [0, 0, 1, 1, 0, 0, 1], []]) # Type 62
    self.assertEqual(double_31(D1), D2)
    
    D1 = C34CrvDiv(C_2_4, [[z4^2 + 1, z4^3 + z4^2, z4^2 + z4 + 1, 1],
                           [z4, z4^3 + z4^2 + z4 + 1, z4^3 + z4^2, 0, 1],
                           [z4^2 + z4 + 1, z4^3 + 1, z4, 0, 0, 1]])
    D2 = C34CrvDiv(C_2_4, [[z4^3 + z4 + 1, z4^3 + z4^2, z4^3 + z4^2, z4^3 + z4^2 + z4, z4^2 + 1, z4^3 + z4^2 + 1, 1],
                           [z4^3 + z4^2, z4^3 + z4^2 + z4, z4 + 1, z4^3 + 1, z4^3 + z4^2 + z4, z4^3 + z4^2 + 1, 0, 1],
                           [z4^3 + z4^2, z4^2 + 1, z4^3 + z4, z4^3 + z4^2, z4^2 + 1, 0, 0, 0, 1]])
    self.assertEqual(double_31(D1), D2)

    D1 = C34CrvDiv(C_3, [[1, 0, 1, 1], [0, 2, 0, 0, 1], [2, 0, 0, 0, 0, 1]])
    D2 = C34CrvDiv(C_3, [[1, 1, 2, 2, 2, 1], [0, 0, 0, 2, 2, 0, 1, 1], []]) # Type 63
    self.assertEqual(double_31(D1), D2)

    D1 = C34CrvDiv(C_3_3, [[z3, z3^2, 0, 1],
                           [2*z3^2 + z3 + 1, z3^2 + 2*z3, 2*z3 + 1, 0, 1],
                           [z3 + 1, z3^2 + z3 + 2, 2*z3^2 + z3 + 1, 0, 0, 1]])
    D2 = C34CrvDiv(C_3_3, [[2*z3^2, 2*z3^2 + 2*z3, z3^2 + z3, 2*z3 + 1, 2*z3^2 + 2, 2*z3^2 + 1, 1],
                           [2*z3^2 + 2, z3^2 + 2*z3, z3^2 + 2*z3, 2*z3 + 2, 2*z3 + 2, z3^2 + 2*z3 + 1, 0, 1],
                           [2*z3^2 + z3 + 2, z3^2 + 2*z3 + 1, 2*z3 + 2, 2*z3^2 + 2*z3 + 2, 2*z3^2 + 2*z3, 2*z3^2 + 2, 0, 0, 1]])
    self.assertEqual(double_31(D1), D2) # XXX : Fail
    
    D1 = C34CrvDiv(C_11, [[4, 2, 10, 1], [10, 10, 2, 0, 1], [4, 1, 6, 0, 0, 1]])
    D2 = C34CrvDiv(C_11, [[4, 7, 9, 4, 1], [4, 3, 9, 3, 0, 1, 1, 0, 0, 1], []]) # Type 64
    self.assertEqual(double_31(D1), D2) # XXX : Fail

    D1 = C34CrvDiv(C_31, [[27, 18, 13, 1], [1, 26, 10, 0, 1], [20, 8, 10, 0, 0, 1]])
    D2 = C34CrvDiv(C_31, [[21, 26, 22, 14, 14, 25, 1], [23, 18, 3, 11, 23, 13, 0, 1], [14, 22, 20, 6, 8, 19, 0, 0, 1]])
    self.assertEqual(double_31(D1), D2)

    D1 = C34CrvDiv(C_31_2, [[28*z2 + 16, 23*z2 + 2, 22*z2 + 16, 1],
                            [8*z2 + 30, 15*z2 + 28, 21*z2 + 3, 0, 1],
                            [2*z2 + 5, 10*z2 + 3, 5*z2 + 21, 0, 0, 1]])
    D2 = C34CrvDiv(C_31_2, [[z2 + 14, 28*z2 + 1, 29*z2 + 23, 7*z2 + 5, 13*z2 + 30, 12*z2 + 27, 1],
                            [22*z2 + 9, 17*z2 + 12, 6*z2, 4*z2, 9*z2 + 23, 10*z2 + 20, 0, 1],
                            [12*z2 + 24, 11*z2 + 1, 5*z2 + 13, 21*z2 + 1, 24*z2 + 16, 10*z2 + 13, 0, 0, 1]])
    self.assertEqual(double_31(D1), D2)

    D1 = C34CrvDiv(C_1009, [[788, 432, 202, 1], [900, 690, 797, 0, 1], [202, 670, 141, 0, 0, 1]])
    D2 = C34CrvDiv(C_1009, [[696, 629, 726, 743, 4, 92, 1], [53, 743, 310, 872, 317, 315, 0, 1], [933, 339, 970, 273, 964, 389, 0, 0, 1]])
    self.assertEqual(double_31(D1), D2)

    # Test case where D1 is atypical
    D1 = C34CrvDiv(C_2, [[0, 1, 0, 1], [0, 0, 1, 0, 1], [0, 0, 1, 0, 0, 1]])
    D2 = C34CrvDiv(C_2, [[0, 0, 1, 0, 1], [0, 0, 0, 1, 0, 0, 0, 0, 0, 1], []]) # Type 64
    self.assertEqual(double_31(D1), D2) # XXX : Fail

    D1 = C34CrvDiv(C_3, [[2, 0, 0, 1], [1, 2, 2, 0, 1], [0, 0, 2, 0, 0, 1]])
    D2 = C34CrvDiv(C_3, [[2, 0, 1, 0, 2, 0, 1], [1, 1, 0, 1, 2, 0, 0, 1], [2, 2, 0, 2, 0, 2, 0, 0, 1]])
    self.assertEqual(double_31(D1), D2) # XXX : Fail
    
    D1 = C34CrvDiv(C_31, [[19, 20, 0, 1], [29, 26, 19, 0, 1], [24, 19, 25, 0, 0, 1]])
    D2 = C34CrvDiv(C_31, [[25, 27, 15, 24, 24, 27, 1], [11, 27, 12, 3, 19, 29, 0, 1], [2, 13, 12, 17, 22, 4, 0, 0, 1]])
    self.assertEqual(double_31(D1), D2) # XXX : Fail

    D1 = C34CrvDiv(C_1009, [[406, 311, 0, 1], [676, 908, 303, 0, 1], [200, 59, 601, 0, 0, 1]])
    D2 = C34CrvDiv(C_1009, [[194, 925, 436, 195, 333, 130, 1], [754, 287, 54, 180, 479, 917, 0, 1], [232, 331, 440, 422, 103, 554, 0, 0, 1]])
    self.assertEqual(double_31(D1), D2) # XXX : Fail

    # TODO : Find cases where D2 is type 62, 63, 64, but not characteristic 2 or 3.
    
    # Cases where D2 is type 65 (principal)
    C = C34Crv(GF(97), [80, 10, 22, 25, 68, 53, 87, 66, 27])
    D1 = C34CrvDiv(C, [[4, 6, 33, 1], [18, 88, 91, 0, 1], [10, 54, 8, 0, 0, 1]])
    D2 = C34CrvDiv(C, [[4, 6, 33, 1], [], []])

    C = C34Crv(GF(97), [1, 33, 36, 7, 5, 37, 46, 95, 19])
    D1 = C34CrvDiv(C, [[3, 19, 57, 1], [13, 21, 50, 0, 1], [32, 25, 4, 0, 0, 1]]) # D1 = 2P + Q
    D2 = C34CrvDiv(C, [[3, 19, 57, 1], [], []])                                   # D2 = 4P + 2Q
