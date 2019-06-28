class TestFlip(unittest.TestCase) :
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
  


  def test_flip_0(self) :
    D = C34CrvDiv(self.C_31, [[1], [], []])
    self.assertEqual(flip_0(D), D)



  def test_flip_11(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    D = C34CrvDiv(C_2, [[0, 1], [0, 0, 1], []])
    A = C34CrvDiv(C_2, [[0, 1], [1, 0, 1, 0, 0, 1], []])
    self.assertEqual(flip_11(D), A)

    D = C34CrvDiv(C_2_4, [[z4, 1], [z4 + 1, 0, 1], []])
    A = C34CrvDiv(C_2_4, [[z4, 1], [z4^3 + z4^2 + z4 + 1, 0, z4^3 + 1, 0, 0, 1], []])
    self.assertEqual(flip_11(D), A)

    D = C34CrvDiv(C_3, [[2, 1], [2, 0, 1], []])
    A = C34CrvDiv(C_3, [[2, 1], [0, 0, 2, 0, 0, 1], []])
    self.assertEqual(flip_11(D), A)

    D = C34CrvDiv(C_3_3, [[z3 + 2, 1], [z3^2 + 2, 0, 1], []])
    A = C34CrvDiv(C_3_3, [[z3 + 2, 1], [2*z3^2 + z3, 0, z3 + 2, 0, 0, 1], []])
    self.assertEqual(flip_11(D), A)

    D = C34CrvDiv(C_31, [[13, 1], [9, 0, 1], []])
    A = C34CrvDiv(C_31, [[13, 1], [2, 0, 12, 0, 0, 1], []])
    self.assertEqual(flip_11(D), A)
    
    D = C34CrvDiv(C_31_2, [[6*z2 + 6, 1], [28*z2 + 22, 0, 1], []])
    A = C34CrvDiv(C_31_2, [[6*z2 + 6, 1], [12*z2 + 3, 0, 9*z2 + 5, 0, 0, 1], []])
    self.assertEqual(flip_11(D), A)
    
    D = C34CrvDiv(C_1009, [[731, 1], [830, 0, 1], []])
    A = C34CrvDiv(C_1009, [[731, 1], [865, 0, 134, 0, 0, 1], []])
    self.assertEqual(flip_11(D), A)



  def test_flip_21(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    D = C34CrvDiv(C_2, [[0, 1, 1], [0, 1, 0, 1], []])
    A = C34CrvDiv(C_2, [[0, 1, 1], [1, 1, 0, 1], []])
    self.assertEqual(flip_21(D), A)

    D = C34CrvDiv(C_2_4, [[z4^3 + 1, z4^2 + z4, 1], [z4^3 + z4^2 + z4 + 1, z4^3 + z4^2 + 1, 0, 1], []])
    A = C34CrvDiv(C_2_4, [[z4^3 + 1, z4^2 + z4, 1], [1, z4^3 + z4^2 + z4, 0, 1], []])
    self.assertEqual(flip_21(D), A)

    D = C34CrvDiv(C_3, [[0, 1, 1], [0, 1, 0, 1], []])
    A = C34CrvDiv(C_3, [[0, 1, 1], [1, 2, 0, 1], []])
    self.assertEqual(flip_21(D), A)

    D = C34CrvDiv(C_3_3, [[z3^2 + 2*z3, z3^2, 1], [1, 2*z3 + 2, 0, 1], []])
    A = C34CrvDiv(C_3_3, [[z3^2 + 2*z3, z3^2, 1], [2*z3^2 + 2*z3 + 1, z3, 0, 1], []])
    self.assertEqual(flip_21(D), A)

    D = C34CrvDiv(C_31, [[30, 1, 1], [10, 11, 0, 1], []])
    A = C34CrvDiv(C_31, [[30, 1, 1], [1, 9, 0, 1], []])
    self.assertEqual(flip_21(D), A)

    D = C34CrvDiv(C_31_2, [[12*z2 + 12, 5*z2 + 2, 1], [19*z2 + 27, 12*z2 + 18, 0, 1], []])
    A = C34CrvDiv(C_31_2, [[12*z2 + 12, 5*z2 + 2, 1], [21*z2 + 20, 21*z2 + 25, 0, 1], []])
    self.assertEqual(flip_21(D), A)

    D = C34CrvDiv(C_1009, [[412, 379, 1], [78, 558, 0, 1], []])
    A = C34CrvDiv(C_1009, [[412, 379, 1], [260, 988, 0, 1], []])
    self.assertEqual(flip_21(D), A)



  def test_flip_22(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    D = C34CrvDiv(C_2, [[1, 1], [0, 0, 1, 0, 0, 1], []])
    A = C34CrvDiv(C_2, [[1, 1], [1, 0, 1], []])
    self.assertEqual(flip_22(D), A)

    D = C34CrvDiv(C_2_4, [[z4^3 + z4, 1], [z4, 0, z4^3 + z4, 0, 0, 1], []])
    A = C34CrvDiv(C_2_4, [[z4^3 + z4, 1], [z4^2, 0, 1], []])
    self.assertEqual(flip_22(D), A)

    D = C34CrvDiv(C_3, [[0, 1], [1, 0, 2, 0, 0, 1], []])
    A = C34CrvDiv(C_3, [[0, 1], [0, 0, 1], []])
    self.assertEqual(flip_22(D), A)

    D = C34CrvDiv(C_3_3, [[z3 + 2, 1], [2*z3^2 + z3, 0, z3 + 2, 0, 0, 1], []])
    A = C34CrvDiv(C_3_3, [[z3 + 2, 1], [z3^2 + 2, 0, 1], []])
    self.assertEqual(flip_22(D), A)

    D = C34CrvDiv(C_31, [[28, 1], [17, 0, 28, 0, 0, 1], []])
    A = C34CrvDiv(C_31, [[28, 1], [17, 0, 1], []])
    self.assertEqual(flip_22(D), A)

    D = C34CrvDiv(C_31_2, [[21*z2 + 19, 1], [28*z2 + 26, 0, 8*z2 + 10, 0, 0, 1], []])
    A = C34CrvDiv(C_31_2, [[21*z2 + 19, 1], [2*z2 + 7, 0, 1], []])
    self.assertEqual(flip_22(D), A)

    D = C34CrvDiv(C_1009, [[306, 1], [928, 0, 434, 0, 0, 1], []])
    A = C34CrvDiv(C_1009, [[306, 1], [792, 0, 1], []])
    self.assertEqual(flip_22(D), A)



  def test_flip_31(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    D = C34CrvDiv(C_2, [[1, 0, 0, 1], [0, 0, 1, 0, 1], [1, 1, 1, 0, 0, 1]])
    A = C34CrvDiv(C_2, [[1, 0, 0, 1], [1, 1, 1, 0, 1], [0, 1, 0, 0, 0, 1]])
    self.assertEqual(flip_31(D), A)

    D = C34CrvDiv(C_2_4, [[z4^3 + z4, z4^2 + z4, z4^2 + z4, 1], [z4 + 1, z4^3 + z4^2 + 1, z4^3, 0, 1], [z4 + 1, z4^2 + z4 + 1, z4^2, 0, 0, 1]])
    A = C34CrvDiv(C_2_4, [[z4^3 + z4, z4^2 + z4, z4^2 + z4, 1], [z4^3, z4^3 + z4^2 + z4, z4^3 + z4^2 + 1, 0, 1], [z4^2 + z4 + 1, z4^3 + z4^2 + z4, z4^3 + z4^2, 0, 0, 1]])
    self.assertEqual(flip_31(D), A)

    D = C34CrvDiv(C_3, [[1, 1, 2, 1], [2, 1, 0, 0, 1], [1, 2, 0, 0, 0, 1]])
    A = C34CrvDiv(C_3, [[1, 1, 2, 1], [2, 0, 0, 0, 1], [2, 2, 2, 0, 0, 1]])
    self.assertEqual(flip_31(D), A)

    D = C34CrvDiv(C_3_3, [[2*z3^2 + 2*z3 + 1, z3^2 + 2*z3 + 1, z3, 1], [z3^2 + 2*z3 + 1, z3 + 2, z3^2 + z3, 0, 1], [0, 2*z3^2 + 2*z3 + 2, z3^2 + z3 + 1, 0, 0, 1]])
    A = C34CrvDiv(C_3_3, [[2*z3^2 + 2*z3 + 1, z3^2 + 2*z3 + 1, z3, 1], [z3^2 + z3 + 2, 2*z3^2, 2*z3 + 1, 0, 1], [0, z3, 2*z3^2 + z3, 0, 0, 1]])
    self.assertEqual(flip_31(D), A)

    D = C34CrvDiv(C_31, [[1, 18, 14, 1], [4, 10, 3, 0, 1], [23, 24, 29, 0, 0, 1]])
    A = C34CrvDiv(C_31, [[1, 18, 14, 1], [21, 21, 25, 0, 1], [12, 5, 7, 0, 0, 1]])
    self.assertEqual(flip_31(D), A)

    D = C34CrvDiv(C_31_2, [[13*z2 + 10, 4*z2 + 21, 20, 1], [23*z2 + 27, 4*z2 + 10, 12*z2 + 9, 0, 1], [26*z2 + 3, 13*z2 + 13, 6*z2, 0, 0, 1]])
    A = C34CrvDiv(C_31_2, [[13*z2 + 10, 4*z2 + 21, 20, 1], [23*z2, 5*z2 + 21, 11*z2 + 29, 0, 1], [3*z2 + 20, 3*z2 + 16, 20, 0, 0, 1]])
    self.assertEqual(flip_31(D), A)

    D = C34CrvDiv(C_1009, [[50, 418, 205, 1], [102, 1004, 698, 0, 1], [901, 391, 924, 0, 0, 1]])
    A = C34CrvDiv(C_1009, [[50, 418, 205, 1], [735, 1003, 935, 0, 1], [917, 855, 162, 0, 0, 1]])
    self.assertEqual(flip_31(D), A)
    
    # Test case where <f, g, h> = <f, h>, but <f, g, h> =/= <f, g>
    D = C34CrvDiv(C_31, [[17, 21, 0, 1], [18, 27, 11, 0, 1], [26, 1, 23, 0, 0, 1]])
    A = C34CrvDiv(C_31, [[17, 21, 0, 1], [10, 1, 10, 0, 1], [7, 22, 14, 0, 0, 1]])
    self.assertEqual(flip_31(D), A)
    
    # TODO : Test case where <f, g> =/= <f, g, h> =/= <f, h>
    self.fail("Test not implemented.")


  def test_flip_32(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    D = C34CrvDiv(C_2, [[1, 1, 1], [1, 1, 0, 1, 0, 0, 1], []])
    A = C34CrvDiv(C_2, [[1, 1], [0, 0, 1], []])
    self.assertEqual(flip_32(D), A)

    D = C34CrvDiv(C_2_4, [[z4^3 + z4^2, z4^2 + z4, 1], [z4^3 + z4 + 1, z4, 0, z4^2 + z4, 0, 0, 1], []])
    A = C34CrvDiv(C_2_4, [[z4^2 + 1, 1], [1, 0, 1], []])
    self.assertEqual(flip_32(D), A)

    D = C34CrvDiv(C_3, [[1, 2, 1], [1, 2, 0, 2, 0, 0, 1], []])
    A = C34CrvDiv(C_3, [[0, 1], [1, 0, 1], []])
    self.assertEqual(flip_32(D), A)

    D = C34CrvDiv(C_3_3, [[2*z3^2 + z3 + 1, 2*z3^2 + z3 + 1, 1], [2*z3^2 + 2, z3^2, 0, z3^2 + 1, 0, 0, 1], []])
    A = C34CrvDiv(C_3_3, [[1, 1], [0, 0, 1], []])
    self.assertEqual(flip_32(D), A)

    D = C34CrvDiv(C_31, [[29, 12, 1], [9, 6, 0, 1, 0, 0, 1], []])
    A = C34CrvDiv(C_31, [[22, 1], [13, 0, 1], []])
    self.assertEqual(flip_32(D), A)

    D = C34CrvDiv(C_31_2, [[14*z2 + 24, 13*z2 + 11, 1], [26*z2 + 13, 15*z2 + 29, 0, 21*z2 + 26, 0, 0, 1], []])
    A = C34CrvDiv(C_31_2, [[3, 1], [6*z2 + 22, 0, 1], []])
    self.assertEqual(flip_32(D), A)

    D = C34CrvDiv(C_1009, [[964, 689, 1], [889, 163, 0, 587, 0, 0, 1], []])
    A = C34CrvDiv(C_1009, [[593, 1], [23, 0, 1], []])
    self.assertEqual(flip_32(D), A)



  def test_flip_33(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    D = C34CrvDiv(C_2, [[1, 1], [], []])
    A = C34CrvDiv(C_2, [[1], [], []])
    self.assertEqual(flip_33(D), A)

    D = C34CrvDiv(C_2_4, [[z4^2 + 1, 1], [], []])
    A = C34CrvDiv(C_2_4, [[1], [], []])
    self.assertEqual(flip_33(D), A)

    D = C34CrvDiv(C_3, [[0, 1], [], []])
    A = C34CrvDiv(C_3, [[1], [], []])
    self.assertEqual(flip_33(D), A)

    D = C34CrvDiv(C_3_3, [[2*z3^2 + 2*z3 + 2, 1], [], []])
    A = C34CrvDiv(C_3_3, [[1], [], []])
    self.assertEqual(flip_33(D), A)

    D = C34CrvDiv(C_31, [[14, 1], [], []])
    A = C34CrvDiv(C_31, [[1], [], []])
    self.assertEqual(flip_33(D), A)

    D = C34CrvDiv(C_31_2, [[14*z2 + 2, 1], [], []])
    A = C34CrvDiv(C_31_2, [[1], [], []])
    self.assertEqual(flip_33(D), A)

    D = C34CrvDiv(C_1009, [[24, 1], [], []])
    A = C34CrvDiv(C_1009, [[1], [], []])
    self.assertEqual(flip_33(D), A)



  def test_flip_41(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    # Test case where <f, g, h> = <f, g>
    D = C34CrvDiv(C_2, [[0, 0, 1, 0, 1], [0, 0, 0, 1, 0, 1], [0, 0, 0, 1, 0, 0, 1]])
    A = C34CrvDiv(C_2, [[0, 1, 1, 1], [0, 0, 1, 0, 1], [0, 0, 0, 0, 0, 1]])
    self.assertEqual(flip_41(D), A)

    D = C34CrvDiv(C_2_4, [[z4^3 + z4^2 + 1, z4^2 + 1, z4 + 1, 1, 1], [z4^3 + z4^2, 0, z4^3 + z4^2, z4^3 + z4, 0, 1], [z4^3 + 1, z4, z4^3, z4^3 + z4, 0, 0, 1]])
    A = C34CrvDiv(C_2_4, [[z4^2 + z4 + 1, z4^3 + z4^2 + z4, z4^3 + z4 + 1, 1], [z4^3 + z4, z4^3 + z4 + 1, z4^3, 0, 1], [z4^3 + z4^2, z4^3 + z4^2, 1, 0, 0, 1]])
    self.assertEqual(flip_41(D), A)

    D = C34CrvDiv(C_3, [[0, 2, 2, 1, 1], [0, 0, 2, 0, 0, 1], [0, 2, 0, 0, 0, 0, 1]])
    A = C34CrvDiv(C_3, [[1, 0, 1, 1], [2, 2, 1, 0, 1], [1, 0, 1, 0, 0, 1]])
    self.assertEqual(flip_41(D), A)

    D = C34CrvDiv(C_3_3, [[z3 + 2, 2*z3^2 + z3 + 2, 2*z3^2 + z3, z3 + 2, 1], [z3^2 + z3 + 2, 2*z3^2 + 2, 2*z3^2, z3^2 + 2*z3 + 2, 0, 1], [z3, z3^2 + 2*z3 + 1, z3^2 + z3, z3^2, 0, 0, 1]])
    A = C34CrvDiv(C_3_3, [[2*z3^2 + z3 + 2, 2*z3^2 + 2*z3 + 2, 2*z3^2, 1], [z3^2 + z3, 2*z3^2 + 2*z3, z3^2 + 2*z3 + 2, 0, 1], [z3^2 + z3 + 2, z3^2 + 2, 0, 0, 0, 1]])
    self.assertEqual(flip_41(D), A)

    D = C34CrvDiv(C_31, [[21, 22, 10, 8, 1], [18, 27, 17, 20, 0, 1], [0, 19, 22, 15, 0, 0, 1]])
    A = C34CrvDiv(C_31, [[22, 9, 22, 1], [0, 12, 20, 0, 1], [12, 25, 23, 0, 0, 1]])
    self.assertEqual(flip_41(D), A)

    D = C34CrvDiv(C_31_2, [[15*z2 + 18, 20*z2 + 9, 7*z2 + 24, 30*z2 + 24, 1], [5*z2 + 16, 30*z2 + 1, 26*z2 + 17, 29*z2 + 4, 0, 1], [17*z2 + 28, 28*z2 + 21, 30*z2 + 12, 16*z2 + 9, 0, 0, 1]])
    A = C34CrvDiv(C_31_2, [[23*z2 + 21, 10*z2 + 28, 14*z2 + 19, 1], [26*z2 + 3, 14*z2 + 20, 28*z2 + 22, 0, 1], [z2 + 21, 3*z2 + 2, 29*z2 + 2, 0, 0, 1]])
    self.assertEqual(flip_41(D), A)

    D = C34CrvDiv(C_1009, [[380, 975, 197, 929, 1], [448, 240, 92, 210, 0, 1], [478, 585, 588, 164, 0, 0, 1]])
    A = C34CrvDiv(C_1009, [[646, 314, 556, 1], [601, 870, 281, 0, 1], [193, 1003, 918, 0, 0, 1]])
    self.assertEqual(flip_41(D), A)

    # TODO : Test case where <f, g, h> = <f, h> =/= <f, g>
    # TODO : Test case where <f, h> =/= <f, g, h> =/= <f, g>
    self.fail("Test not implemented.")



  def test_flip_42(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    D = C34CrvDiv(C_2, [[0, 0, 0, 1], [0, 0, 0, 0, 1], []])
    A = C34CrvDiv(C_2, [[0, 1], [1, 0, 1, 0, 0, 1], []])
    self.assertEqual(flip_42(D), A)

    D = C34CrvDiv(C_2_4, [[1, z4^2 + z4 + 1, 0, 1], [z4^3 + z4^2 + z4 + 1, 1, z4^3 + z4^2 + z4 + 1, 0, 1], []])
    A = C34CrvDiv(C_2_4, [[z4^3, 1], [z4^3, 0, z4^3 + z4^2 + z4, 0, 0, 1], []])
    self.assertEqual(flip_42(D), A)

    D = C34CrvDiv(C_3, [[0, 1, 0, 1], [0, 2, 0, 0, 1], []])
    A = C34CrvDiv(C_3, [[1, 1], [2, 0, 1, 0, 0, 1], []])
    self.assertEqual(flip_42(D), A)

    D = C34CrvDiv(C_3_3, [[2, 2*z3^2 + 2, 0, 1], [2*z3 + 1, 2*z3^2 + z3 + 1, 2*z3 + 2, 0, 1], []])
    A = C34CrvDiv(C_3_3, [[2*z3^2 + z3, 1], [z3^2 + 2, 0, z3^2 + 1, 0, 0, 1], []])
    self.assertEqual(flip_42(D), A)

    D = C34CrvDiv(C_31, [[10, 22, 0, 1], [18, 20, 4, 0, 1], []])
    A = C34CrvDiv(C_31, [[18, 1], [7, 0, 9, 0, 0, 1], []])
    self.assertEqual(flip_42(D), A)

    D = C34CrvDiv(C_31_2, [[10, 15*z2 + 30, 0, 1], [27*z2 + 22, 13*z2 + 9, 7*z2 + 19, 0, 1], []])
    A = C34CrvDiv(C_31_2, [[8*z2 + 11, 1], [z2, 0, 11*z2 + 15, 0, 0, 1], []])
    self.assertEqual(flip_42(D), A)

    D = C34CrvDiv(C_1009, [[119, 884, 0, 1], [11, 565, 734, 0, 1], []])
    A = C34CrvDiv(C_1009, [[150, 1], [690, 0, 204, 0, 0, 1], []])
    self.assertEqual(flip_42(D), A)



  def test_flip_43(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    D = C34CrvDiv(C_2, [[1, 0, 0, 1], [0, 0, 1, 0, 1], []])
    A = C34CrvDiv(C_2, [[1, 1], [1, 0, 0, 0, 0, 1], []])
    self.assertEqual(flip_43(D), A)

    D = C34CrvDiv(C_2_4, [[z4^3 + z4^2, 0, 0, 1], [z4 + 1, z4 + 1, z4^3 + z4, 0, z4^3 + 1, 1], []])
    A = C34CrvDiv(C_2_4, [[1, 0, 1], [z4^3 + z4^2, 0, 0, 1], []])
    self.assertEqual(flip_43(D), A)

    D = C34CrvDiv(C_3, [[1, 1, 0, 1], [1, 2, 2, 0, 1], []])
    A = C34CrvDiv(C_3, [[2, 1], [0, 0, 2, 0, 0, 1], []])
    self.assertEqual(flip_43(D), A)

    D = C34CrvDiv(C_3_3, [[z3^2 + 2*z3 + 2, z3^2 + 2*z3 + 1, 0, 1], [z3^2 + z3, 2*z3^2 + 2, 2*z3^2 + z3 + 2, 0, 2*z3^2 + z3, 1], []])
    A = C34CrvDiv(C_3_3, [[2*z3^2 + 2*z3, z3^2 + 2*z3 + 2, 1], [z3^2 + 2*z3 + 2, z3^2 + 2*z3 + 1, 0, 1], []])
    self.assertEqual(flip_43(D), A)

    D = C34CrvDiv(C_31, [[17, 17, 0, 1], [5, 9, 25, 0, 24, 1], []])
    A = C34CrvDiv(C_31, [[0, 24, 1], [17, 17, 0, 1], []])
    self.assertEqual(flip_43(D), A)

    D = C34CrvDiv(C_31_2, [[21*z2 + 16, 5*z2 + 4, 0, 1], [10*z2 + 21, 22*z2 + 12, 27*z2 + 14, 0, 30, 1], []])
    A = C34CrvDiv(C_31_2, [[25*z2 + 25, 13*z2 + 11, 1], [21*z2 + 16, 5*z2 + 4, 0, 1], []])
    self.assertEqual(flip_43(D), A)

    D = C34CrvDiv(C_1009, [[762, 264, 0, 1], [689, 37, 546, 0, 225, 1], []])
    A = C34CrvDiv(C_1009, [[909, 89, 1], [762, 264, 0, 1], []])
    self.assertEqual(flip_43(D), A)



  def test_flip_44(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    D = C34CrvDiv(C_2, [[1, 0, 1], [], []])
    A = C34CrvDiv(C_2, [[1], [], []])
    self.assertEqual(flip_44(D), A)

    D = C34CrvDiv(C_2_4, [[z4^2, z4^3 + z4^2 + z4 + 1, 1], [], []])
    A = C34CrvDiv(C_2_4, [[1], [], []])
    self.assertEqual(flip_44(D), A)

    D = C34CrvDiv(C_3, [[1, 0, 1], [], []])
    A = C34CrvDiv(C_3, [[1], [], []])
    self.assertEqual(flip_44(D), A)

    D = C34CrvDiv(C_3_3, [[z3^2, 2*z3^2 + z3, 1], [], []])
    A = C34CrvDiv(C_3_3, [[1], [], []])
    self.assertEqual(flip_44(D), A)

    D = C34CrvDiv(C_31, [[29, 9, 1], [], []])
    A = C34CrvDiv(C_31, [[1], [], []])
    self.assertEqual(flip_44(D), A)

    D = C34CrvDiv(C_31_2, [[19*z2 + 23, 27*z2 + 19, 1], [], []])
    A = C34CrvDiv(C_31_2, [[1], [], []])
    self.assertEqual(flip_44(D), A)

    D = C34CrvDiv(C_1009, [[728, 532, 1], [], []])
    A = C34CrvDiv(C_1009, [[1], [], []])
    self.assertEqual(flip_44(D), A)



  def test_flip_51(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    # Test case where <f, g, h> = <f, g>
    # Test case where <f, g, h> = <f, h> =/= <f, g>
    # Test case where <f, h> =/= <f, g, h> =/= <f, g>
    self.fail("Test not implemented.")



  def test_flip_52(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    D = C34CrvDiv(C_2, [[1, 0, 1, 1, 1], [1, 0, 0, 1, 0, 1], [0, 1, 1, 0, 0, 0, 0, 0, 0, 1]])
    A = C34CrvDiv(C_2, [[1, 1], [1, 0, 0, 0, 0, 1], []])
    self.assertEqual(flip_52(D), A)

    D = C34CrvDiv(C_2_4, [[0, z4^2 + z4 + 1, z4^3 + z4^2 + z4, z4^3 + 1, 1], [0, z4^2 + 1, z4^3 + z4, z4^3 + z4^2 + 1, 0, 1], [z4^3 + 1, z4 + 1, z4^3 + z4^2 + z4 + 1, z4^2 + z4 + 1, 0, 0, z4^3 + z4^2 + 1, 0, 0, 1]])
    A = C34CrvDiv(C_2_4, [[z4^3 + z4^2 + z4, 1], [z4^3 + z4^2 + z4, 0, 1, 0, 0, 1], []])
    self.assertEqual(flip_52(D), A)

    D = C34CrvDiv(C_3, [[0, 2, 2, 1, 1], [0, 1, 1, 2, 0, 1], [0, 2, 1, 0, 0, 0, 0, 0, 0, 1]])
    A = C34CrvDiv(C_3, [[2, 1], [1, 0, 1, 0, 0, 1], []])
    self.assertEqual(flip_52(D), A)

    D = C34CrvDiv(C_3_3, [[2*z3^2, z3^2 + 2*z3 + 1, 2*z3^2 + z3 + 2, z3^2 + 2*z3, 1], [2*z3^2 + 1, 2*z3^2 + 2, z3, z3^2 + 1, 0, 1], [z3 + 1, z3^2 + 1, z3, 0, 0, 0, 2*z3^2, 0, 0, 1]])
    A = C34CrvDiv(C_3_3, [[2*z3^2 + z3 + 2, 1], [z3^2 + z3 + 2, 0, z3^2, 0, 0, 1], []])
    self.assertEqual(flip_52(D), A)

    D = C34CrvDiv(C_31, [[30, 17, 14, 27, 1], [18, 29, 7, 15, 0, 1], [20, 19, 10, 26, 0, 0, 2, 0, 0, 1]])
    A = C34CrvDiv(C_31, [[14, 1], [26, 0, 2, 0, 0, 1], []])
    self.assertEqual(flip_52(D), A)

    D = C34CrvDiv(C_31_2, [[2*z2 + 28, 2*z2 + 9, 9*z2 + 12, 25*z2 + 24, 1], [26*z2 + 10, 27*z2 + 26, 12*z2 + 2, 30*z2 + 28, 0, 1], [2*z2 + 2, 10*z2 + 23, 19*z2 + 16, 18*z2 + 16, 0, 0, 4*z2 + 28, 0, 0, 1]])
    A = C34CrvDiv(C_31_2, [[9*z2 + 12, 1], [1, 0, 17*z2 + 28, 0, 0, 1], []])
    self.assertEqual(flip_52(D), A)

    D = C34CrvDiv(C_1009, [[692, 560, 313, 179, 1], [692, 565, 341, 247, 0, 1], [243, 397, 991, 589, 0, 0, 203, 0, 0, 1]])
    A = C34CrvDiv(C_1009, [[313, 1], [839, 0, 201, 0, 0, 1], []])
    self.assertEqual(flip_52(D), A)



  def test_flip_53(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    D = C34CrvDiv(C_2, [[0, 0, 1, 0, 1], [0, 0, 0, 0, 0, 1, 1], []])
    A = C34CrvDiv(C_2, [[1, 1, 1], [1, 0, 0, 1], []])
    self.assertEqual(flip_53(D), A)

    D = C34CrvDiv(C_2_4, [[z4^2, z4 + 1, z4^3 + z4^2 + 1, 0, 1], [z4^2 + 1, z4^2 + z4, z4^2 + 1, z4^3, 0, z4^2 + 1, 1], []])
    A = C34CrvDiv(C_2_4, [[z4^3 + 1, z4^2 + 1, 1], [z4^3 + 1, z4^3 + z4^2 + z4 + 1, 0, 1], []])
    self.assertEqual(flip_53(D), A)

    D = C34CrvDiv(C_3, [[2, 0, 1, 1, 1], [0, 0, 1, 0, 0, 1, 1], []])
    A = C34CrvDiv(C_3, [[1, 2, 1], [2, 0, 0, 1], []])
    self.assertEqual(flip_53(D), A)

    D = C34CrvDiv(C_3_3, [[2*z3^2 + z3, 2*z3^2 + z3, 2*z3^2 + z3, 2*z3, 1], [z3 + 2, z3^2 + 2*z3 + 2, 2*z3^2 + 2*z3, 2*z3^2 + z3 + 1, 0, z3^2, 1], []])
    A = C34CrvDiv(C_3_3, [[z3 + 2, z3^2 + 2*z3, 1], [2*z3^2 + 2*z3 + 1, 2, 0, 1], []])
    self.assertEqual(flip_53(D), A)

    D = C34CrvDiv(C_31, [[4, 27, 8, 13, 1], [0, 25, 2, 26, 0, 6, 1], []])
    A = C34CrvDiv(C_31, [[20, 19, 1], [26, 19, 0, 1], []])
    self.assertEqual(flip_53(D), A)

    D = C34CrvDiv(C_31_2, [[24*z2 + 14, 7*z2 + 10, 23*z2 + 24, 11*z2 + 16, 1], [14*z2 + 23, 14*z2 + 4, 24*z2 + 3, 14*z2 + 7, 0, 6*z2, 1], []])
    A = C34CrvDiv(C_31_2, [[6, 17*z2 + 16, 1], [10*z2 + 30, z2 + 24, 0, 1], []])
    self.assertEqual(flip_53(D), A)

    D = C34CrvDiv(C_1009, [[482, 936, 246, 432, 1], [325, 327, 189, 634, 0, 3, 1], []])
    A = C34CrvDiv(C_1009, [[888, 435, 1], [680, 339, 0, 1], []])
    self.assertEqual(flip_53(D), A)



  def test_flip_54(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    D = C34CrvDiv(C_2, [[1, 1, 1, 1], [1, 0, 1, 0, 1, 1, 0, 0, 1], []])
    A = C34CrvDiv(C_2, [[1, 1], [1, 0, 1], []])
    self.assertEqual(flip_54(D), A)

    D = C34CrvDiv(C_2_4, [[z4^3, z4^2 + z4, z4^3 + z4^2, 1], [z4^2 + 1, 0, z4^3 + z4^2 + 1, 0, z4^3 + z4, z4^3 + z4 + 1, 0, 0, 1], []])
    A = C34CrvDiv(C_2_4, [[z4^3 + z4 + 1, 1], [z4^2 + z4, 0, 1], []])
    self.assertEqual(flip_54(D), A)

    D = C34CrvDiv(C_3, [[2, 0, 2, 1], [1, 0, 1, 0, 2, 0, 0, 0, 1], []])
    A = C34CrvDiv(C_3, [[2, 1], [0, 0, 1], []])
    self.assertEqual(flip_54(D), A)

    D = C34CrvDiv(C_3_3, [[z3^2, z3 + 1, 2*z3 + 1, 1], [z3^2 + 1, z3^2, 2*z3, 0, 2, 2*z3^2 + 2*z3, 0, 0, 1], []])
    A = C34CrvDiv(C_3_3, [[z3^2 + z3 + 2, 1], [z3^2 + 2*z3, 0, 1], []])
    self.assertEqual(flip_54(D), A)

    D = C34CrvDiv(C_31, [[6, 25, 19, 1], [7, 27, 28, 0, 23, 29, 0, 0, 1], []])
    A = C34CrvDiv(C_31, [[14, 1], [2, 0, 1], []])
    self.assertEqual(flip_54(D), A)

    D = C34CrvDiv(C_31_2, [[12*z2 + 5, z2 + 5, 18*z2 + 28, 1], [8*z2 + 26, 25*z2 + 22, 15*z2 + 2, 0, 2*z2 + 24, 20*z2 + 19, 0, 0, 1], []])
    A = C34CrvDiv(C_31_2, [[23*z2 + 5, 1], [11*z2 + 17, 0, 1], []])
    self.assertEqual(flip_54(D), A)

    D = C34CrvDiv(C_1009, [[924, 198, 558, 1], [950, 186, 218, 0, 216, 208, 0, 0, 1], []])
    A = C34CrvDiv(C_1009, [[344, 1], [692, 0, 1], []])
    self.assertEqual(flip_54(D), A)



  def test_flip_61(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    # Test case where <f, g, h> = <f, g>
    # Test case where <f, g, h> = <f, h> =/= <f, g>
    # Test case where <f, h> =/= <f, g, h> =/= <f, g>
    self.fail("Test not implemented.")



  def test_flip_62(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    D = C34CrvDiv(C_2, [[1, 1, 1, 1, 0, 1], [1, 1, 1, 0, 0, 0, 1], []])
    A = C34CrvDiv(C_2, [[0, 1], [1, 0, 1, 0, 0, 1], []])
    self.assertEqual(flip_62(D), A)

    D = C34CrvDiv(C_2_4, [[z4^3 + z4^2 + z4 + 1, z4^2 + 1, z4 + 1, 0, z4^3 + z4, 1], [z4^3 + z4^2 + z4, 0, z4^3 + z4 + 1, 1, z4 + 1, 0, 1], []])
    A = C34CrvDiv(C_2_4, [[z4^3 + z4 + 1, 1], [z4^3 + z4^2 + z4, 0, 1, 0, 0, 1], []])
    self.assertEqual(flip_62(D), A)

    D = C34CrvDiv(C_3, [[1, 2, 2, 1, 0, 1], [0, 0, 0, 0, 2, 0, 1], []])
    A = C34CrvDiv(C_3, [[0, 1], [1, 0, 2, 0, 0, 1], []])
    self.assertEqual(flip_62(D), A)

    D = C34CrvDiv(C_3_3, [[2*z3^2, 2*z3^2 + 1, 2*z3^2 + z3 + 1, 1, 2*z3^2 + 2*z3 + 2, 1], [z3^2 + 2*z3 + 2, z3, 2*z3^2 + 2*z3 + 1, 2*z3^2 + 2*z3 + 2, 2*z3^2 + 2*z3 + 1, 0, 1], []])
    A = C34CrvDiv(C_3_3, [[2*z3^2 + z3 + 1, 1], [1, 0, 2*z3^2 + z3 + 2, 0, 0, 1], []])
    self.assertEqual(flip_62(D), A)

    D = C34CrvDiv(C_31, [[18, 1, 20, 21, 22, 1], [14, 1, 16, 22, 9, 0, 1], []])
    A = C34CrvDiv(C_31, [[7, 1], [17, 0, 21, 0, 0, 1], []])
    self.assertEqual(flip_62(D), A)

    D = C34CrvDiv(C_31_2, [[28*z2 + 18, 17*z2 + 1, 12*z2 + 4, 7*z2 + 8, 6*z2 + 8, 1], [25*z2, z2 + 9, 6*z2 + 26, 6*z2 + 11, 14*z2 + 24, 0, 1], []])
    A = C34CrvDiv(C_31_2, [[30*z2 + 1, 1], [z2 + 12, 0, 26*z2 + 9, 0, 0, 1], []])
    self.assertEqual(flip_62(D), A)

    D = C34CrvDiv(C_1009, [[304, 682, 916, 35, 610, 1], [76, 208, 255, 348, 849, 0, 1], []])
    A = C34CrvDiv(C_1009, [[44, 1], [723, 0, 310, 0, 0, 1], []])
    self.assertEqual(flip_62(D), A)



  def test_flip_63(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    self.fail("Test not implemented.")



  def test_flip_64(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    D = C34CrvDiv(C_2, [[0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 0, 1, 0, 0, 1], []])
    A = C34CrvDiv(C_2, [[0, 1], [0, 0, 1], []])
    self.assertEqual(flip_64(D), A)

    D = C34CrvDiv(C_2_4, [[z4^3 + 1, 1, z4^3 + z4^2 + z4 + 1, z4^3 + z4^2 + z4, 1], [z4^3 + z4, z4^3 + z4^2 + z4 + 1, 0, z4^3 + z4^2, 0, 0, z4^2 + 1, 0, 0, 1], []])
    A = C34CrvDiv(C_2_4, [[z4^3 + z4^2 + z4 + 1, 1], [1, 0, 1], []])
    self.assertEqual(flip_64(D), A)

    D = C34CrvDiv(C_3, [[0, 2, 2, 1, 1], [2, 1, 0, 0, 0, 0, 2, 0, 0, 1], []])
    A = C34CrvDiv(C_3, [[0, 1], [0, 0, 1], []])
    self.assertEqual(flip_64(D), A)

    D = C34CrvDiv(C_3_3, [[2*z3^2 + z3 + 1, z3 + 2, z3^2 + z3 + 2, 2*z3^2, 1], [z3^2 + 2*z3 + 2, 2*z3^2, 0, 2*z3^2 + z3 + 2, 0, 0, z3^2 + 2, 0, 0, 1], []])
    A = C34CrvDiv(C_3_3, [[2*z3^2 + 2*z3 + 1, 1], [z3 + 2, 0, 1], []])
    self.assertEqual(flip_64(D), A)

    D = C34CrvDiv(C_31, [[5, 27, 24, 3, 1], [20, 28, 0, 2, 0, 0, 22, 0, 0, 1], []])
    A = C34CrvDiv(C_31, [[13, 1], [9, 0, 1], []])
    self.assertEqual(flip_64(D), A)

    D = C34CrvDiv(C_31_2, [[30*z2 + 4, 19*z2 + 8, 8*z2 + 24, 12*z2 + 5, 1], [22*z2 + 24, 30*z2 + 12, 0, 7*z2 + 4, 0, 0, 28*z2 + 4, 0, 0, 1], []])
    A = C34CrvDiv(C_31_2, [[9*z2 + 5, 1], [15*z2 + 10, 0, 1], []])
    self.assertEqual(flip_64(D), A)

    D = C34CrvDiv(C_1009, [[398, 421, 424, 3, 1], [145, 876, 0, 302, 0, 0, 1, 0, 0, 1], []])
    A = C34CrvDiv(C_1009, [[539, 1], [559, 0, 1], []])
    self.assertEqual(flip_64(D), A)



  def test_flip_65(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    D = C34CrvDiv(C_2, [[0, 0, 0, 1], [], []])
    A = C34CrvDiv(C_2, [[1], [], []])
    self.assertEqual(flip_65(D), A)

    D = C34CrvDiv(C_2_4, [[z4, z4^3 + z4, z4^3 + z4^2 + z4, 1], [], []])
    A = C34CrvDiv(C_2_4, [[1], [], []])
    self.assertEqual(flip_65(D), A)

    D = C34CrvDiv(C_3, [[2, 1, 1, 1], [], []])
    A = C34CrvDiv(C_3, [[1], [], []])
    self.assertEqual(flip_65(D), A)

    D = C34CrvDiv(C_3_3, [[z3^2 + 2*z3 + 2, z3 + 2, 2, 1], [], []])
    A = C34CrvDiv(C_3_3, [[1], [], []])
    self.assertEqual(flip_65(D), A)

    D = C34CrvDiv(C_31, [[14, 23, 22, 1], [], []])
    A = C34CrvDiv(C_31, [[1], [], []])
    self.assertEqual(flip_65(D), A)

    D = C34CrvDiv(C_31_2, [[6*z2 + 10, 5*z2 + 23, 3*z2 + 27, 1], [], []])
    A = C34CrvDiv(C_31_2, [[1], [], []])
    self.assertEqual(flip_65(D), A)

    D = C34CrvDiv(C_1009, [[408, 893, 618, 1], [], []])
    A = C34CrvDiv(C_1009, [[1], [], []])
    self.assertEqual(flip_65(D), A)

