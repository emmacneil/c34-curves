class TestReduce(unittest.TestCase) :
  def setUp(self) :
    self.z2 = GF(31^2).gen()
    self.z3 = GF(3^3).gen()
    self.z4 = GF(2^4).gen()
    z2, z3, z4 = self.z2, self.z3, self.z4
    
    self.C_2 = C34Curve(GF(2), [0, 0, 1, 0, 1, 1, 1, 1, 1])
    self.C_3 = C34Curve(GF(3), [0, 2, 1, 1, 0, 2, 2, 0, 2])
    self.C_11 = C34Curve(GF(11), [6, 5, 1, 8, 6, 3, 7, 10, 4]) # The vertical tangent line at P = (9 : 6 : 1) intersects P with multiplicity 2
    self.C_31 = C34Curve(GF(31), [2, 29, 13, 7, 28, 25, 17, 13, 17])
    self.C_2_4 = C34Curve(GF(2^4), [z4^3 + 1, z4, z4^3 + z4^2 + 1, z4^2, z4^3 + z4, z4^3 + z4 + 1, z4^3 + z4, z4^3 + z4^2 + z4, z4^3 + 1])
    self.C_3_3 = C34Curve(GF(3^3), [z3 + 2, z3 + 1, z3^2 + z3 + 2, 2*z3^2 + z3 + 1, 2*z3^2 + z3, z3^2 + 2, 2*z3^2 + z3, 2*z3^2 + 1, 2])
    self.C_31_2 = C34Curve(GF(31^2), [9*z2 + 11, 11*z2 + 28, 17*z2 + 2, 24*z2 + 22, 6*z2 + 29, 21*z2 + 8, 18*z2 + 5, 18*z2 + 15, 13*z2 + 10])
    self.C_1009 = C34Curve(GF(1009), [715, 286, 748, 7, 92, 446, 122, 935, 314])



  def test_reduce_11(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    D = C34CurveDivisor(C_2, [[1, 1], [0, 0, 1], []])
    A = C34CurveDivisor(C_2, [[1, 1], [0, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_2_4, [[z4^3 + z4^2 + z4 + 1, 1], [z4^3 + z4^2 + z4 + 1, 0, 1], []])
    A = C34CurveDivisor(C_2_4, [[z4^3 + z4^2 + z4 + 1, 1], [z4^3 + z4^2 + z4 + 1, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_3, [[2, 1], [0, 0, 1], []])
    A = C34CurveDivisor(C_3, [[2, 1], [0, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_3_3, [[z3 + 1, 1], [2*z3, 0, 1], []])
    A = C34CurveDivisor(C_3_3, [[z3 + 1, 1], [2*z3, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_31, [[12, 1], [6, 0, 1], []])
    A = C34CurveDivisor(C_31, [[12, 1], [6, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_31_2, [[23*z2 + 25, 1], [15*z2 + 30, 0, 1], []])
    A = C34CurveDivisor(C_31_2, [[23*z2 + 25, 1], [15*z2 + 30, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_1009, [[217, 1], [414, 0, 1], []])
    A = C34CurveDivisor(C_1009, [[217, 1], [414, 0, 1], []])
    self.assertEqual(reduce(D), A)



  def test_reduce_21(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    D = C34CurveDivisor(C_2, [[1, 1, 1], [1, 0, 0, 1], []])
    A = C34CurveDivisor(C_2, [[1, 1, 1], [1, 0, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_2_4, [[z4^2 + z4, z4^3 + z4^2 + 1, 1], [z4^2, z4^3 + 1, 0, 1], []])
    A = C34CurveDivisor(C_2_4, [[z4^2 + z4, z4^3 + z4^2 + 1, 1], [z4^2, z4^3 + 1, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_3, [[0, 2, 1], [0, 0, 0, 1], []])
    A = C34CurveDivisor(C_3, [[0, 2, 1], [0, 0, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_3_3, [[z3^2, z3^2 + 2*z3 + 2, 1], [z3^2 + z3, 2*z3^2 + z3 + 2, 0, 1], []])
    A = C34CurveDivisor(C_3_3, [[z3^2, z3^2 + 2*z3 + 2, 1], [z3^2 + z3, 2*z3^2 + z3 + 2, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_31, [[0, 29, 1], [11, 25, 0, 1], []])
    A = C34CurveDivisor(C_31, [[0, 29, 1], [11, 25, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_31_2, [[23*z2 + 7, 2*z2 + 30, 1], [28*z2 + 13, 26*z2 + 25, 0, 1], []])
    A = C34CurveDivisor(C_31_2, [[23*z2 + 7, 2*z2 + 30, 1], [28*z2 + 13, 26*z2 + 25, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_1009, [[43, 352, 1], [842, 751, 0, 1], []])
    A = C34CurveDivisor(C_1009, [[43, 352, 1], [842, 751, 0, 1], []])
    self.assertEqual(reduce(D), A)



  def test_reduce_22(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    D = C34CurveDivisor(C_2, [[0, 1], [1, 0, 1, 0, 0, 1], []])
    A = C34CurveDivisor(C_2, [[0, 1], [1, 0, 1, 0, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_2_4, [[z4^2 + z4 + 1, 1], [z4^3 + z4 + 1, 0, z4^2, 0, 0, 1], []])
    A = C34CurveDivisor(C_2_4, [[z4^2 + z4 + 1, 1], [z4^3 + z4 + 1, 0, z4^2, 0, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_3, [[0, 1], [0, 0, 1, 0, 0, 1], []])
    A = C34CurveDivisor(C_3, [[0, 1], [0, 0, 1, 0, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_3_3, [[z3, 1], [z3^2 + 2, 0, z3^2 + 2*z3 + 1, 0, 0, 1], []])
    A = C34CurveDivisor(C_3_3, [[z3, 1], [z3^2 + 2, 0, z3^2 + 2*z3 + 1, 0, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_31, [[27, 1], [24, 0, 30, 0, 0, 1], []])
    A = C34CurveDivisor(C_31, [[27, 1], [24, 0, 30, 0, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_31_2, [[15*z2 + 6, 1], [7*z2 + 5, 0, 24*z2 + 24, 0, 0, 1], []])
    A = C34CurveDivisor(C_31_2, [[15*z2 + 6, 1], [7*z2 + 5, 0, 24*z2 + 24, 0, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_1009, [[524, 1], [910, 0, 902, 0, 0, 1], []])
    A = C34CurveDivisor(C_1009, [[524, 1], [910, 0, 902, 0, 0, 1], []])
    self.assertEqual(reduce(D), A)



  def test_reduce_31(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    D = C34CurveDivisor(C_2, [[0, 1, 1, 1], [0, 0, 1, 0, 1], [0, 0, 0, 0, 0, 1]])
    A = C34CurveDivisor(C_2, [[0, 1, 1, 1], [0, 0, 1, 0, 1], [0, 0, 0, 0, 0, 1]])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_2_4, [[z4 + 1, z4^2 + 1, z4^2 + z4, 1], [z4^3 + z4^2, z4^2, z4^3 + 1, 0, 1], [z4^3 + 1, z4^3 + z4^2, z4^3 + z4^2, 0, 0, 1]])
    A = C34CurveDivisor(C_2_4, [[z4 + 1, z4^2 + 1, z4^2 + z4, 1], [z4^3 + z4^2, z4^2, z4^3 + 1, 0, 1], [z4^3 + 1, z4^3 + z4^2, z4^3 + z4^2, 0, 0, 1]])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_3, [[0, 2, 1, 1], [0, 0, 1, 0, 1], [0, 0, 2, 0, 0, 1]])
    A = C34CurveDivisor(C_3, [[0, 2, 1, 1], [0, 0, 1, 0, 1], [0, 0, 2, 0, 0, 1]])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_3_3, [[z3^2 + z3 + 1, 2*z3^2 + 1, z3^2 + 1, 1], [2*z3^2 + 2*z3 + 1, 2*z3, 2*z3^2 + 2, 0, 1], [2*z3^2 + z3, 1, z3^2 + 1, 0, 0, 1]])
    A = C34CurveDivisor(C_3_3, [[z3^2 + z3 + 1, 2*z3^2 + 1, z3^2 + 1, 1], [2*z3^2 + 2*z3 + 1, 2*z3, 2*z3^2 + 2, 0, 1], [2*z3^2 + z3, 1, z3^2 + 1, 0, 0, 1]])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_31, [[3, 14, 1, 1], [8, 29, 10, 0, 1], [24, 3, 23, 0, 0, 1]])
    A = C34CurveDivisor(C_31, [[3, 14, 1, 1], [8, 29, 10, 0, 1], [24, 3, 23, 0, 0, 1]])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_31_2, [[6*z2 + 1, 10*z2 + 6, 12*z2 + 12, 1], [3*z2 + 3, 18*z2, 13*z2 + 13, 0, 1], [24*z2 + 29, 4*z2 + 23, 14*z2 + 4, 0, 0, 1]])
    A = C34CurveDivisor(C_31_2, [[6*z2 + 1, 10*z2 + 6, 12*z2 + 12, 1], [3*z2 + 3, 18*z2, 13*z2 + 13, 0, 1], [24*z2 + 29, 4*z2 + 23, 14*z2 + 4, 0, 0, 1]])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_1009, [[983, 603, 748, 1], [914, 88, 887, 0, 1], [808, 899, 476, 0, 0, 1]])
    A = C34CurveDivisor(C_1009, [[983, 603, 748, 1], [914, 88, 887, 0, 1], [808, 899, 476, 0, 0, 1]])
    self.assertEqual(reduce(D), A)



  def test_reduce_32(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    D = C34CurveDivisor(C_2, [[0, 0, 1], [0, 0, 0, 1, 0, 0, 1], []])
    A = C34CurveDivisor(C_2, [[0, 1], [1, 0, 1, 0, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_2_4, [[z4^2 + z4, z4^3 + z4 + 1, 1], [z4^3 + z4^2 + 1, z4^3 + z4 + 1, 0, 0, 0, 0, 1], []])
    A = C34CurveDivisor(C_2_4, [[z4 + 1, 1], [z4^2, 0, z4^3 + z4 + 1, 0, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_3, [[1, 2, 1], [0, 1, 0, 1, 0, 0, 1], []])
    A = C34CurveDivisor(C_3, [[1, 1], [2, 0, 1, 0, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_3_3, [[2*z3^2 + z3, 2, 1], [z3 + 2, z3^2 + z3 + 1, 0, z3^2 + 2*z3, 0, 0, 1], []])
    A = C34CurveDivisor(C_3_3, [[2*z3 + 1, 1], [z3^2 + z3 + 1, 0, 2*z3^2 + 2*z3 + 2, 0, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_31, [[8, 17, 1], [6, 14, 0, 14, 0, 0, 1], []])
    A = C34CurveDivisor(C_31, [[30, 1], [1, 0, 17, 0, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_31_2, [[9, 20*z2 + 8, 1], [2*z2 + 18, 15*z2 + 24, 0, 14*z2 + 30, 0, 0, 1], []])
    A = C34CurveDivisor(C_31_2, [[17*z2 + 18, 1], [21*z2 + 14, 0, 10*z2 + 9, 0, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_1009, [[887, 436, 1], [785, 298, 0, 408, 0, 0, 1], []])
    A = C34CurveDivisor(C_1009, [[852, 1], [806, 0, 585, 0, 0, 1], []])
    self.assertEqual(reduce(D), A)



  def test_reduce_33(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    D = C34CurveDivisor(C_2, [[0, 1], [], []])
    A = C34CurveDivisor(C_2, [[1], [], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_2_4, [[z4^3 + z4^2 + z4, 1], [], []])
    A = C34CurveDivisor(C_2_4, [[1], [], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_3, [[1, 1], [], []])
    A = C34CurveDivisor(C_3, [[1], [], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_3_3, [[2*z3^2 + 2*z3 + 2, 1], [], []])
    A = C34CurveDivisor(C_3_3, [[1], [], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_31, [[0, 1], [], []])
    A = C34CurveDivisor(C_31, [[1], [], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_31_2, [[20*z2 + 24, 1], [], []])
    A = C34CurveDivisor(C_31_2, [[1], [], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_1009, [[867, 1], [], []])
    A = C34CurveDivisor(C_1009, [[1], [], []])
    self.assertEqual(reduce(D), A)



  def test_reduce_41(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    # Test cases where D is typical
    D = C34CurveDivisor(C_2, [[1, 1, 1, 1, 1], [1, 0, 1, 0, 0, 1], [0, 1, 1, 0, 0, 0, 1]])
    A = C34CurveDivisor(C_2, [[0, 0, 1, 1], [0, 0, 1, 0, 1], [0, 0, 1, 0, 0, 1]])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_2_4, [[z4^3, z4^3 + z4^2 + z4 + 1, 0, z4^3, 1], [z4^3 + z4^2 + z4, z4^2, z4^3 + 1, 0, 0, 1], [z4^2, z4^3 + 1, z4^3 + z4^2 + z4 + 1, 0, 0, 0, 1]])
    A = C34CurveDivisor(C_2_4, [[z4^3 + z4 + 1, z4^2 + 1, z4^3 + z4^2, 1], [z4, z4^2 + z4 + 1, z4^3 + 1, 0, 1], [z4^3 + z4^2, z4^3 + z4^2 + z4 + 1, z4^3 + z4^2, 0, 0, 1]])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_3, [[0, 2, 2, 1, 1], [2, 0, 2, 1, 0, 1], [2, 2, 0, 1, 0, 0, 1]])
    A = C34CurveDivisor(C_3, [[0, 0, 2, 1], [1, 1, 2, 0, 1], [1, 2, 0, 0, 0, 1]])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_3_3, [[z3^2 + z3 + 1, z3^2 + 2*z3 + 1, z3^2 + z3 + 2, z3, 1], [2*z3^2 + z3 + 2, z3^2 + z3 + 2, 2*z3^2 + z3, 2*z3^2 + 2*z3 + 2, 0, 1], [z3 + 2, 2*z3^2 + z3, 2, z3, 0, 0, 1]])
    A = C34CurveDivisor(C_3_3, [[z3^2 + 2*z3 + 1, z3^2 + 2*z3 + 2, 2*z3 + 2, 1], [2*z3^2 + z3, 2*z3^2 + 2, 2*z3^2 + 2*z3 + 2, 0, 1], [z3, z3^2 + z3 + 1, z3^2 + 2*z3 + 2, 0, 0, 1]])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_31, [[14, 23, 16, 30, 1], [1, 23, 16, 16, 0, 1], [29, 21, 14, 21, 0, 0, 1]])
    A = C34CurveDivisor(C_31, [[29, 26, 17, 1], [26, 2, 14, 0, 1], [27, 22, 23, 0, 0, 1]])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_31_2, [[15*z2 + 17, 13*z2 + 14, 20*z2 + 19, 27*z2 + 4, 1], [23*z2 + 20, 12*z2 + 4, 21*z2 + 10, 30*z2 + 11, 0, 1], [10*z2 + 15, 5*z2 + 4, 14*z2 + 9, 24*z2 + 29, 0, 0, 1]])
    A = C34CurveDivisor(C_31_2, [[2*z2 + 13, 20, 30*z2 + 10, 1], [26*z2 + 11, 11*z2 + 14, 5*z2 + 14, 0, 1], [z2 + 6, 2, 8*z2 + 22, 0, 0, 1]])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_1009, [[24, 778, 274, 749, 1], [436, 147, 860, 826, 0, 1], [125, 427, 798, 212, 0, 0, 1]])
    A = C34CurveDivisor(C_1009, [[933, 339, 823, 1], [125, 564, 884, 0, 1], [130, 971, 871, 0, 0, 1]])
    self.assertEqual(reduce(D), A)

    # Test cases where D is atypical
    D = C34CurveDivisor(C_2, [[0, 1, 1, 1, 1], [0, 1, 1, 1, 0, 1], [1, 0, 0, 0, 0, 0, 1]])
    A = C34CurveDivisor(C_2, [[0, 1, 0, 1], [0, 0, 0, 0, 1], [1, 1, 1, 0, 0, 1]])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_2_4, [[z4, 0, z4, z4^3 + 1, 1], [z4, z4^3, z4 + 1, z4^3 + z4^2 + 1, 0, 1], [z4^3 + z4^2 + z4, z4, z4^3 + z4, z4^2 + z4 + 1, 0, 0, 1]])
    A = C34CurveDivisor(C_2_4, [[z4^2 + z4 + 1, z4^3, 0, 1], [z4^3 + z4^2 + 1, z4 + 1, z4^3 + z4, 0, 1], [z4 + 1, z4^3 + z4^2, z4^3 + z4, 0, 0, 1]])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_3, [[2, 0, 2, 1, 1], [1, 0, 2, 2, 0, 1], [1, 0, 1, 1, 0, 0, 1]])
    A = C34CurveDivisor(C_3, [[1, 1, 0, 1], [0, 0, 2, 0, 1], [2, 1, 2, 0, 0, 1]])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_3_3, [[2*z3^2 + 2*z3, 0, 2*z3^2 + z3 + 2, 1, 1], [2*z3^2 + z3 + 2, 2*z3^2 + z3 + 1, z3^2 + 2*z3, 2, 0, 1], [z3^2 + 1, 2*z3^2 + 2, 2*z3^2 + z3 + 1, 2*z3^2 + 2, 0, 0, 1]])
    A = C34CurveDivisor(C_3_3, [[z3^2 + z3 + 2, 2*z3, 0, 1], [1, 2*z3^2 + z3 + 1, z3^2 + z3 + 1, 0, 1], [z3^2 + 2, z3, 2*z3^2, 0, 0, 1]])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_31, [[16, 21, 26, 15, 1], [6, 16, 5, 23, 0, 1], [29, 29, 5, 28, 0, 0, 1]])
    A = C34CurveDivisor(C_31, [[20, 22, 0, 1], [13, 20, 27, 0, 1], [2, 21, 30, 0, 0, 1]])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_31_2, [[23*z2 + 29, 13*z2 + 22, 24*z2 + 5, 7*z2 + 6, 1], [15*z2 + 4, 16*z2 + 17, 8*z2 + 16, 4*z2 + 18, 0, 1], [27*z2 + 29, 20*z2 + 11, 27*z2 + 1, 3*z2 + 28, 0, 0, 1]])
    A = C34CurveDivisor(C_31_2, [[23*z2 + 12, 22*z2 + 22, 0, 1], [12*z2 + 7, 2*z2 + 7, 29*z2 + 17, 0, 1], [5*z2 + 1, 17*z2 + 19, 19*z2 + 28, 0, 0, 1]])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_1009, [[827, 256, 944, 1002, 1], [874, 458, 978, 960, 0, 1], [24, 675, 633, 492, 0, 0, 1]])
    A = C34CurveDivisor(C_1009, [[948, 650, 0, 1], [476, 623, 715, 0, 1], [334, 199, 173, 0, 0, 1]])
    self.assertEqual(reduce(D), A)



  def test_reduce_42(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    D = C34CurveDivisor(C_2, [[0, 0, 0, 1], [0, 0, 0, 0, 1], []])
    A = C34CurveDivisor(C_2, [[0, 1], [0, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_2_4, [[z4^2 + z4 + 1, z4, 0, 1], [z4^2 + z4, z4^3 + z4 + 1, z4^3 + z4^2 + 1, 0, 1], []])
    A = C34CurveDivisor(C_2_4, [[z4^3 + z4^2 + z4 + 1, 1], [z4^3 + z4 + 1, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_3, [[2, 0, 0, 1], [1, 2, 2, 0, 1], []])
    A = C34CurveDivisor(C_3, [[1, 1], [2, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_3_3, [[z3^2 + 2*z3 + 2, 2*z3^2 + z3 + 2, 0, 1], [z3 + 1, 2*z3 + 1, 2*z3^2 + 2*z3 + 2, 0, 1], []])
    A = C34CurveDivisor(C_3_3, [[2*z3, 1], [2*z3 + 1, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_31, [[5, 11, 0, 1], [13, 9, 29, 0, 1], []])
    A = C34CurveDivisor(C_31, [[13, 1], [9, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_31_2, [[21, 9*z2 + 18, 0, 1], [29*z2 + 1, 24*z2 + 14, 22*z2 + 18, 0, 1], []])
    A = C34CurveDivisor(C_31_2, [[18*z2, 1], [24*z2 + 14, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_1009, [[361, 114, 0, 1], [535, 538, 528, 0, 1], []])
    A = C34CurveDivisor(C_1009, [[595, 1], [538, 0, 1], []])
    self.assertEqual(reduce(D), A)



  def test_reduce_43(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    D = C34CurveDivisor(C_2, [[0, 0, 1, 1], [0, 1, 1, 0, 0, 1], []])
    A = C34CurveDivisor(C_2, [[0, 1, 1], [1, 1, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_2_4, [[z4^2 + z4 + 1, 1, z4^3 + z4, 1], [z4^3 + z4^2 + 1, z4^3 + z4^2 + z4 + 1, z4^2, 0, z4^3 + z4^2 + z4, 1], []])
    A = C34CurveDivisor(C_2_4, [[0, z4^2 + z4 + 1, 1], [z4 + 1, z4^2 + z4, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_3, [[1, 1, 1, 1], [2, 2, 2, 0, 0, 1], []])
    A = C34CurveDivisor(C_3, [[1, 2, 1], [2, 0, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_3_3, [[2*z3, z3^2, z3^2 + z3 + 2, 1], [z3^2 + 1, 2, 2*z3, 0, z3^2 + z3 + 2, 1], []])
    A = C34CurveDivisor(C_3_3, [[2*z3^2 + z3, 2*z3^2 + 2*z3, 1], [2*z3 + 2, 2*z3^2 + 2, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_31, [[19, 16, 26, 1], [3, 25, 16, 0, 2, 1], []])
    A = C34CurveDivisor(C_31, [[11, 15, 1], [8, 26, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_31_2, [[15*z2 + 12, 22*z2 + 19, 19*z2 + 1, 1], [18*z2 + 25, 26*z2 + 14, 30*z2 + 25, 0, 20*z2 + 27, 1], []])
    A = C34CurveDivisor(C_31_2, [[22*z2 + 13, 24*z2 + 14, 1], [7*z2 + 15, 26*z2 + 8, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_1009, [[332, 725, 919, 1], [29, 293, 554, 0, 744, 1], []])
    A = C34CurveDivisor(C_1009, [[299, 579, 1], [226, 781, 0, 1], []])
    self.assertEqual(reduce(D), A)



  def test_reduce_44(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    D = C34CurveDivisor(C_2, [[1, 0, 1], [], []])
    A = C34CurveDivisor(C_2, [[1], [], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_2_4, [[z4^2 + z4 + 1, z4^3, 1], [], []])
    A = C34CurveDivisor(C_2_4, [[1], [], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_3, [[0, 2, 1], [], []])
    A = C34CurveDivisor(C_3, [[1], [], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_3_3, [[z3^2 + z3 + 1, z3^2 + z3 + 2, 1], [], []])
    A = C34CurveDivisor(C_3_3, [[1], [], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_31, [[3, 30, 1], [], []])
    A = C34CurveDivisor(C_31, [[1], [], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_31_2, [[4*z2 + 27, 15*z2 + 29, 1], [], []])
    A = C34CurveDivisor(C_31_2, [[1], [], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_1009, [[823, 108, 1], [], []])
    A = C34CurveDivisor(C_1009, [[1], [], []])
    self.assertEqual(reduce(D), A)



  def test_reduce_51(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    # Test cases where D is typical
    D = C34CurveDivisor(C_2, [[1, 0, 0, 1, 1, 1], [0, 0, 1, 1, 1, 0, 1], [1, 0, 0, 1, 1, 0, 0, 1]])
    A = C34CurveDivisor(C_2, [[0, 1, 1, 1], [0, 0, 1, 0, 1], [0, 0, 0, 0, 0, 1]])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_2_4, [[z4^2 + z4, z4^2 + z4, z4^2 + z4 + 1, z4^2 + z4 + 1, z4 + 1, 1], [z4^3 + z4^2 + 1, z4^2 + z4 + 1, z4^3 + z4^2 + z4 + 1, z4^3 + z4^2 + z4 + 1, z4^3, 0, 1], [z4^3 + z4^2, z4^2 + 1, z4^3 + z4 + 1, 1, z4^2 + z4 + 1, 0, 0, 1]])
    A = C34CurveDivisor(C_2_4, [[z4^2 + 1, z4^2 + 1, z4^3 + z4^2, 1], [z4^3 + z4 + 1, z4^3 + z4^2 + z4 + 1, z4, 0, 1], [z4^2 + z4 + 1, z4^3 + 1, z4^3 + z4^2 + 1, 0, 0, 1]])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_3, [[1, 0, 2, 2, 1, 1], [0, 1, 0, 2, 1, 0, 1], [0, 0, 0, 2, 0, 0, 0, 1]])
    A = C34CurveDivisor(C_3, [[0, 2, 2, 1], [0, 0, 2, 0, 1], [0, 0, 0, 0, 0, 1]])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_3_3, [[2*z3 + 2, z3, 2*z3^2 + z3 + 2, z3 + 2, z3^2, 1], [z3^2 + 2*z3 + 1, 2*z3^2 + 2*z3, 2*z3^2, z3^2 + 2*z3 + 2, 2*z3^2 + z3, 0, 1], [2*z3, z3^2 + 2, 1, 1, z3^2 + z3 + 2, 0, 0, 1]])
    A = C34CurveDivisor(C_3_3, [[z3^2, z3^2, 2*z3^2 + 2, 1], [2*z3^2, 2*z3^2 + z3, z3^2, 0, 1], [z3^2 + z3 + 2, 2*z3, z3^2 + 2*z3 + 1, 0, 0, 1]])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_31, [[4, 5, 10, 2, 24, 1], [23, 22, 15, 27, 23, 0, 1], [5, 26, 16, 6, 29, 0, 0, 1]])
    A = C34CurveDivisor(C_31, [[6, 8, 1, 1], [25, 20, 11, 0, 1], [9, 9, 28, 0, 0, 1]])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_31_2, [[12*z2 + 15, 25*z2 + 16, 27*z2 + 27, 13*z2 + 17, 27*z2 + 14, 1], [25*z2 + 5, 23*z2 + 17, 26*z2 + 4, 17*z2 + 24, 25*z2 + 5, 0, 1], [4*z2 + 21, 23*z2 + 4, 30*z2 + 2, 27*z2 + 27, 28*z2 + 17, 0, 0, 1]])
    A = C34CurveDivisor(C_31_2, [[12, 18*z2 + 8, 17*z2, 1], [30*z2 + 29, 10*z2 + 30, 8, 0, 1], [25*z2 + 16, 22*z2 + 10, 28*z2 + 22, 0, 0, 1]])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_1009, [[808, 890, 561, 122, 741, 1], [140, 152, 818, 458, 967, 0, 1], [371, 700, 862, 152, 483, 0, 0, 1]])
    A = C34CurveDivisor(C_1009, [[617, 117, 436, 1], [389, 574, 334, 0, 1], [388, 256, 457, 0, 0, 1]])
    self.assertEqual(reduce(D), A)

    # Test cases where D is atypical
    D = C34CurveDivisor(C_2, [[0, 1, 0, 0, 0, 1], [0, 1, 1, 0, 1, 0, 1], [1, 1, 1, 0, 0, 0, 0, 1]])
    A = C34CurveDivisor(C_2, [[0, 1, 0, 1], [0, 0, 0, 0, 1], [1, 1, 1, 0, 0, 1]])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_2_4, [[0, z4^3 + 1, z4^3, z4 + 1, z4^3 + z4, 1], [z4^2, z4^3 + 1, z4^2 + 1, 0, 0, 0, 1], [z4^3 + z4^2 + z4 + 1, z4^3 + z4^2, z4^2 + z4 + 1, z4^3 + z4 + 1, z4^3 + z4^2 + z4 + 1, 0, 0, 1]])
    A = C34CurveDivisor(C_2_4, [[z4^3 + z4, 0, 0, 1], [z4^3 + z4^2 + z4 + 1, 1, z4^3 + z4^2 + z4 + 1, 0, 1], [z4^3 + z4, z4^3, z4^3 + z4, 0, 0, 1]])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_3, [[0, 2, 0, 2, 1, 1], [0, 2, 2, 1, 0, 0, 1], [2, 2, 0, 0, 1, 0, 0, 1]])
    A = C34CurveDivisor(C_3, [[0, 2, 0, 1], [0, 0, 0, 0, 1], [1, 2, 2, 0, 0, 1]])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_3_3, [[2*z3^2 + z3 + 1, 2*z3^2 + 1, 0, 2*z3, 2*z3^2 + 2, 1], [2*z3, z3^2 + 2, 2*z3, z3^2 + z3, z3^2 + 1, 0, 1], [z3^2 + 2*z3 + 1, 2*z3^2 + z3 + 2, 2*z3 + 2, z3^2 + 2*z3 + 1, 2*z3^2 + 2*z3, 0, 0, 1]])
    A = C34CurveDivisor(C_3_3, [[2*z3^2 + z3 + 2, 2*z3, 0, 1], [2*z3 + 2, z3^2, 2*z3^2 + z3 + 1, 0, 1], [2*z3^2, 2*z3, 2*z3^2 + z3 + 2, 0, 0, 1]])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_31, [[30, 23, 15, 29, 11, 1], [14, 30, 17, 6, 11, 0, 1], [26, 15, 19, 11, 25, 0, 0, 1]])
    A = C34CurveDivisor(C_31, [[4, 7, 0, 1], [2, 27, 15, 0, 1], [8, 27, 2, 0, 0, 1]])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_31_2, [[28*z2 + 2, 17*z2 + 10, 4*z2 + 23, 14*z2 + 30, 15*z2 + 18, 1], [27*z2 + 30, 17*z2 + 1, 29*z2 + 26, 26*z2 + 13, 3*z2 + 8, 0, 1], [13*z2 + 12, z2 + 23, 7*z2 + 17, 13*z2 + 24, 22*z2 + 27, 0, 0, 1]])
    A = C34CurveDivisor(C_31_2, [[5*z2 + 21, 22*z2 + 1, 0, 1], [10*z2 + 10, 8*z2 + 4, 19*z2 + 8, 0, 1], [23*z2 + 22, 9*z2 + 1, 2*z2 + 11, 0, 0, 1]])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_1009, [[460, 115, 693, 237, 650, 1], [104, 804, 256, 331, 145, 0, 1], [21, 243, 739, 682, 943, 0, 0, 1]])
    A = C34CurveDivisor(C_1009, [[572, 705, 0, 1], [398, 858, 993, 0, 1], [878, 519, 163, 0, 0, 1]])
    self.assertEqual(reduce(D), A)



  def test_reduce_52(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    D = C34CurveDivisor(C_2, [[0, 1, 0, 1, 1], [0, 1, 1, 1, 0, 1], []])
    A = C34CurveDivisor(C_2, [[0, 1], [0, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_2_4, [[z4^2 + z4, z4, z4^3, z4^2, 1], [z4^3 + z4^2 + z4 + 1, z4^3 + z4^2, z4 + 1, z4 + 1, 0, 1], []])
    A = C34CurveDivisor(C_2_4, [[z4^3, 1], [1, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_3, [[0, 1, 1, 1, 1], [0, 1, 1, 2, 0, 1], []])
    A = C34CurveDivisor(C_3, [[1, 1], [2, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_3_3, [[2*z3^2, z3^2, 2*z3^2 + 2*z3 + 1, 2*z3^2 + 1, 1], [2*z3^2 + 2*z3 + 1, z3^2 + 2*z3, 2, z3^2 + z3 + 2, 0, 1], []])
    A = C34CurveDivisor(C_3_3, [[2*z3^2 + 2*z3 + 1, 1], [z3 + 2, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_31, [[5, 27, 7, 24, 1], [3, 19, 12, 13, 0, 1], []])
    A = C34CurveDivisor(C_31, [[7, 1], [11, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_31_2, [[2*z2 + 12, 24*z2 + 10, 18*z2 + 29, 7*z2 + 25, 1], [2*z2 + 9, 10*z2 + 14, 20*z2 + 10, 17*z2 + 18, 0, 1], []])
    A = C34CurveDivisor(C_31_2, [[18*z2 + 29, 1], [8*z2 + 12, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_1009, [[267, 682, 161, 298, 1], [414, 665, 456, 997, 0, 1], []])
    A = C34CurveDivisor(C_1009, [[161, 1], [884, 0, 1], []])
    self.assertEqual(reduce(D), A)



  def test_reduce_53(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    D = C34CurveDivisor(C_2, [[0, 0, 1, 0, 1], [0, 0, 1, 1, 0, 1, 1], []])
    A = C34CurveDivisor(C_2, [[0, 1, 1], [1, 1, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_2_4, [[z4 + 1, z4^3 + z4^2 + z4 + 1, z4, z4^3 + z4, 1], [z4^3 + z4^2 + z4, z4^3, 0, z4^3 + z4^2, 0, 0, 1], []])
    A = C34CurveDivisor(C_2_4, [[z4^3, z4^3 + z4, 1], [z4^2 + z4 + 1, z4^3 + z4^2 + z4, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_3, [[1, 1, 2, 2, 1], [0, 2, 0, 0, 0, 2, 1], []])
    A = C34CurveDivisor(C_3, [[0, 1, 1], [0, 1, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_3_3, [[z3^2 + z3 + 1, 1, z3^2 + 2*z3 + 1, z3 + 2, 1], [z3^2 + 1, 2*z3^2 + 2*z3, 2*z3 + 2, 2*z3^2 + z3 + 1, 0, 2*z3^2, 1], []])
    A = C34CurveDivisor(C_3_3, [[2*z3 + 1, 2*z3^2 + z3 + 2, 1], [z3^2 + 2*z3 + 1, z3, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_31, [[21, 21, 12, 15, 1], [7, 24, 21, 9, 0, 28, 1], []])
    A = C34CurveDivisor(C_31, [[3, 12, 1], [5, 3, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_31_2, [[2*z2 + 12, 24*z2 + 14, 20*z2 + 4, 9*z2 + 21, 1], [29*z2 + 16, 20*z2 + 20, 22*z2 + 11, 6*z2 + 12, 0, 9*z2 + 13, 1], []])
    A = C34CurveDivisor(C_31_2, [[2*z2 + 22, 18*z2 + 3, 1], [13*z2 + 24, 17*z2, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_1009, [[543, 837, 41, 508, 1], [261, 904, 509, 465, 0, 739, 1], []])
    A = C34CurveDivisor(C_1009, [[542, 238, 1], [905, 243, 0, 1], []])
    self.assertEqual(reduce(D), A)



  def test_reduce_54(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    D = C34CurveDivisor(C_2, [[0, 1, 0, 1], [1, 1, 1, 0, 1, 1, 0, 0, 1], []])
    A = C34CurveDivisor(C_2, [[0, 1], [1, 0, 1, 0, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_2_4, [[z4^3 + z4 + 1, z4^3 + z4 + 1, z4^2 + z4, 1], [z4, z4 + 1, z4^3 + 1, 0, z4^3 + z4^2, z4^3 + z4, 0, 0, 1], []])
    A = C34CurveDivisor(C_2_4, [[z4, 1], [z4^3 + z4^2 + z4 + 1, 0, z4^3 + 1, 0, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_3, [[0, 1, 0, 1], [0, 2, 0, 0, 1, 0, 0, 0, 1], []])
    A = C34CurveDivisor(C_3, [[1, 1], [2, 0, 1, 0, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_3_3, [[2*z3^2 + 2*z3 + 2, 1, 2*z3, 1], [z3^2 + 2*z3, 2*z3^2 + z3, 2*z3^2 + 2*z3 + 2, 0, 2, 2*z3^2 + z3, 0, 0, 1], []])
    A = C34CurveDivisor(C_3_3, [[z3^2 + z3 + 1, 1], [2*z3^2 + 1, 0, 2*z3^2 + 1, 0, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_31, [[13, 12, 8, 1], [4, 15, 0, 0, 19, 5, 0, 0, 1], []])
    A = C34CurveDivisor(C_31, [[26, 1], [11, 0, 28, 0, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_31_2, [[24*z2 + 14, 14*z2 + 4, 18*z2 + 28, 1], [14*z2 + 5, 8*z2 + 7, 22*z2 + 18, 0, 21*z2 + 6, 21*z2 + 7, 0, 0, 1], []])
    A = C34CurveDivisor(C_31_2, [[4*z2 + 16, 1], [6*z2 + 6, 0, 4*z2 + 14, 0, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_1009, [[939, 222, 386, 1], [494, 593, 10, 0, 123, 78, 0, 0, 1], []])
    A = C34CurveDivisor(C_1009, [[20, 1], [955, 0, 147, 0, 0, 1], []])
    self.assertEqual(reduce(D), A)



  def test_reduce_61(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    # Test cases where D is typical
    D = C34CurveDivisor(C_2, [[0, 0, 1, 0, 0, 0, 1], [0, 0, 1, 0, 0, 1, 0, 1], [0, 0, 0, 0, 1, 1, 0, 0, 1]])
    A = C34CurveDivisor(C_2, [[0, 0, 1, 1], [1, 1, 0, 0, 1], [0, 1, 1, 0, 0, 1]])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_2_4, [[z4^3, z4^3, z4^3 + 1, z4 + 1, z4 + 1, z4^2 + z4 + 1, 1], [z4^3 + z4, z4^3 + z4^2, z4 + 1, z4^3 + z4, 1, z4^3 + z4^2 + z4, 0, 1], [z4^3 + 1, z4^3 + z4 + 1, z4^2, z4^3 + z4^2 + 1, z4^2 + 1, 1, 0, 0, 1]])
    A = C34CurveDivisor(C_2_4, [[1, z4^3 + z4^2, 1, 1], [z4^3 + z4^2 + 1, z4^3, z4^3 + z4^2, 0, 1], [z4^3, z4^2 + z4 + 1, z4^3 + 1, 0, 0, 1]])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_3, [[0, 2, 2, 0, 2, 1, 1], [1, 1, 0, 1, 2, 0, 0, 1], [0, 0, 2, 0, 1, 0, 0, 0, 1]])
    A = C34CurveDivisor(C_3, [[0, 2, 1, 1], [0, 0, 1, 0, 1], [0, 0, 2, 0, 0, 1]])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_3_3, [[2*z3 + 2, z3^2, z3^2 + z3 + 2, 2*z3^2 + z3 + 1, z3^2 + 2*z3 + 1, z3^2 + z3 + 2, 1], [z3^2, z3^2 + 2*z3 + 2, 0, z3^2 + z3 + 1, 2*z3^2 + 2*z3 + 2, z3^2 + z3 + 1, 0, 1], [2*z3^2 + 2, z3^2 + z3 + 1, z3^2 + 2*z3 + 1, 0, z3^2 + z3, 2*z3^2 + 2, 0, 0, 1]])
    A = C34CurveDivisor(C_3_3, [[0, 2*z3^2 + 2, z3^2 + z3 + 1, 1], [z3^2 + 2*z3 + 1, 2*z3^2 + 2*z3 + 2, z3^2 + 1, 0, 1], [2*z3^2 + z3 + 1, z3^2, z3^2 + z3 + 2, 0, 0, 1]])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_31, [[9, 11, 29, 12, 27, 5, 1], [24, 10, 0, 5, 2, 28, 0, 1], [14, 26, 8, 1, 15, 15, 0, 0, 1]])
    A = C34CurveDivisor(C_31, [[1, 10, 1, 1], [17, 5, 18, 0, 1], [17, 11, 26, 0, 0, 1]])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_31_2, [[22*z2 + 10, 11*z2 + 11, 22*z2 + 9, 28*z2 + 20, 2*z2 + 5, 13*z2 + 13, 1], [14*z2 + 7, 3*z2 + 30, 3*z2 + 4, 24*z2 + 12, 7*z2 + 25, 29*z2 + 22, 0, 1], [11*z2 + 3, 0, 26*z2 + 15, 11*z2 + 19, 27*z2 + 3, 13*z2 + 28, 0, 0, 1]])
    A = C34CurveDivisor(C_31_2, [[16*z2 + 22, 30*z2 + 25, 12*z2 + 22, 1], [19*z2 + 30, 14*z2 + 16, 9*z2 + 28, 0, 1], [5*z2, 28*z2 + 11, z2 + 7, 0, 0, 1]])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_1009, [[219, 928, 182, 533, 167, 449, 1], [745, 821, 588, 123, 688, 539, 0, 1], [456, 860, 265, 683, 331, 955, 0, 0, 1]])
    A = C34CurveDivisor(C_1009, [[128, 950, 712, 1], [498, 393, 15, 0, 1], [899, 522, 681, 0, 0, 1]])
    self.assertEqual(reduce(D), A)

    # Test cases where D is typical
    D = C34CurveDivisor(C_2, [[0, 0, 0, 0, 0, 1, 1], [0, 0, 1, 0, 0, 0, 0, 1], [0, 0, 0, 1, 0, 0, 0, 0, 1]])
    A = C34CurveDivisor(C_2, [[0, 1, 0, 1], [0, 0, 0, 0, 1], [1, 1, 1, 0, 0, 1]])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_2_4, [[z4^3 + 1, z4^3 + 1, z4, z4^3 + 1, z4^2 + 1, z4 + 1, 1], [z4^3, z4^3 + z4^2 + z4 + 1, z4^3 + z4, 0, z4^2 + 1, z4^3, 0, 1], [z4^3 + z4, 0, z4 + 1, z4^2, z4^3 + z4^2, z4^3 + 1, 0, 0, 1]])
    A = C34CurveDivisor(C_2_4, [[z4^2 + z4, 1, 0, 1], [z4^2 + 1, z4 + 1, z4 + 1, 0, 1], [z4^2 + z4, z4^3 + z4^2 + z4 + 1, z4^3 + z4 + 1, 0, 0, 1]])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_3, [[0, 2, 2, 1, 1, 2, 1], [1, 0, 1, 2, 0, 1, 0, 1], [1, 2, 0, 1, 2, 1, 0, 0, 1]])
    A = C34CurveDivisor(C_3, [[0, 2, 0, 1], [2, 1, 2, 0, 1], [1, 2, 2, 0, 0, 1]])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_3_3, [[z3 + 1, 2*z3^2, z3, z3^2 + z3 + 1, z3 + 1, z3, 1], [1, 2*z3^2 + z3 + 2, z3^2 + 2, 2*z3^2 + 2*z3 + 1, 2*z3^2 + 2, z3^2 + 2*z3 + 1, 0, 1], [2*z3^2 + z3 + 2, 2*z3^2 + 2*z3, 2*z3^2 + 2*z3 + 1, 2*z3^2 + z3, 2*z3^2 + 2*z3 + 2, 2*z3 + 2, 0, 0, 1]])
    A = C34CurveDivisor(C_3_3, [[2, 2*z3^2 + z3 + 1, 0, 1], [2*z3^2 + 2*z3 + 1, z3, 2*z3^2 + 2*z3, 0, 1], [2*z3^2 + z3 + 2, 1, 2*z3^2, 0, 0, 1]])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_31, [[20, 8, 15, 14, 9, 1, 1], [19, 13, 30, 17, 7, 24, 0, 1], [25, 8, 3, 28, 7, 14, 0, 0, 1]])
    A = C34CurveDivisor(C_31, [[5, 15, 0, 1], [18, 17, 12, 0, 1], [22, 5, 1, 0, 0, 1]])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_31_2, [[6*z2 + 25, 3*z2, 23, 15*z2 + 1, 25*z2 + 9, 12*z2 + 5, 1], [8*z2 + 30, 20*z2 + 10, 4*z2 + 29, 20*z2 + 28, 15*z2 + 20, 29*z2 + 20, 0, 1], [5*z2 + 2, 25*z2 + 14, 7*z2 + 13, 25*z2 + 10, 7*z2 + 17, 9*z2 + 16, 0, 0, 1]])
    A = C34CurveDivisor(C_31_2, [[30*z2 + 12, 7*z2, 0, 1], [9*z2 + 10, 30*z2 + 23, 22, 0, 1], [24*z2 + 29, 10*z2 + 10, 2*z2 + 21, 0, 0, 1]])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_1009, [[521, 499, 451, 471, 246, 306, 1], [661, 348, 97, 411, 1008, 825, 0, 1], [250, 678, 935, 406, 152, 123, 0, 0, 1]])
    A = C34CurveDivisor(C_1009, [[375, 858, 0, 1], [117, 170, 339, 0, 1], [593, 336, 863, 0, 0, 1]])
    self.assertEqual(reduce(D), A)



  def test_reduce_62(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    D = C34CurveDivisor(C_2, [[1, 1, 1, 1, 1, 1], [0, 0, 0, 1, 0, 0, 1], []])
    A = C34CurveDivisor(C_2, [[0, 1], [0, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_2_4, [[z4^3 + 1, z4, z4^3, z4^3 + 1, 1, 1], [z4^3 + z4^2 + z4, z4^3 + z4 + 1, z4^2 + 1, z4^3 + z4^2 + 1, z4^3 + z4^2 + z4 + 1, 0, 1], []])
    A = C34CurveDivisor(C_2_4, [[z4 + 1, 1], [z4^3, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_3, [[1, 0, 1, 2, 1, 1], [1, 0, 1, 1, 0, 0, 1], []])
    A = C34CurveDivisor(C_3, [[2, 1], [2, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_3_3, [[z3^2 + z3 + 1, z3^2 + 2*z3, z3^2, z3^2 + 1, 1, 1], [2*z3^2 + 2*z3 + 1, z3^2 + z3 + 1, 2*z3^2, 1, z3^2 + 2, 0, 1], []])
    A = C34CurveDivisor(C_3_3, [[z3^2 + z3 + 1, 1], [2*z3^2 + 2*z3 + 1, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_31, [[1, 27, 27, 28, 1, 1], [5, 0, 5, 27, 0, 0, 1], []])
    A = C34CurveDivisor(C_31, [[7, 1], [10, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_31_2, [[28*z2 + 22, 29*z2 + 10, 10*z2 + 3, 13*z2 + 10, 1, 1], [25*z2 + 23, 4*z2 + 12, 8*z2 + 21, 2*z2 + 19, 23*z2 + 27, 0, 1], []])
    A = C34CurveDivisor(C_31_2, [[20*z2, 1], [24*z2 + 10, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_1009, [[486, 349, 47, 371, 1, 1], [225, 714, 438, 507, 251, 0, 1], []])
    A = C34CurveDivisor(C_1009, [[536, 1], [125, 0, 1], []])
    self.assertEqual(reduce(D), A)



  def test_reduce_63(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    D = C34CurveDivisor(C_2, [[0, 1, 0, 1, 1, 1], [1, 1, 0, 0, 0, 0, 0, 1], []])
    A = C34CurveDivisor(C_2, [[0, 1, 1], [1, 1, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_2_4, [[z4^2, z4, z4, z4^3 + z4^2, z4^2 + z4 + 1, 1], [z4^3 + z4^2, z4^3 + 1, z4^2 + 1, z4^3 + z4 + 1, z4^3 + z4^2 + z4, 0, z4^3, 1], []])
    A = C34CurveDivisor(C_2_4, [[z4^3 + z4, z4^3 + z4^2 + z4 + 1, 1], [1, z4 + 1, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_3, [[0, 2, 1, 1, 0, 1], [0, 2, 0, 2, 1, 0, 0, 1], []])
    A = C34CurveDivisor(C_3, [[0, 0, 1], [1, 0, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_3_3, [[z3^2 + 2*z3, z3^2 + z3, z3^2 + 1, z3, z3^2 + 2*z3, 1], [2*z3^2 + z3 + 1, z3 + 2, 2*z3^2 + z3, 2*z3 + 2, 2, 0, z3^2 + z3 + 1, 1], []])
    A = C34CurveDivisor(C_3_3, [[2*z3^2 + z3 + 1, z3 + 2, 1], [z3^2 + 2*z3, z3^2 + 2, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_31, [[25, 22, 8, 29, 15, 1], [17, 10, 1, 21, 21, 0, 14, 1], []])
    A = C34CurveDivisor(C_31, [[10, 1, 1], [13, 5, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_31_2, [[25*z2 + 2, 16, 12*z2 + 20, 17*z2 + 28, 9*z2 + 2, 1], [21*z2 + 15, 5*z2 + 7, 21*z2 + 23, 30*z2 + 6, 29*z2 + 9, 0, 17*z2 + 23, 1], []])
    A = C34CurveDivisor(C_31_2, [[z2 + 7, 23*z2 + 10, 1], [29*z2 + 1, 26*z2 + 24, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_1009, [[57, 108, 417, 671, 719, 1], [649, 703, 887, 223, 312, 0, 585, 1], []])
    A = C34CurveDivisor(C_1009, [[429, 134, 1], [188, 200, 0, 1], []])
    self.assertEqual(reduce(D), A)



  def test_reduce_64(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    D = C34CurveDivisor(C_2, [[0, 1, 1, 0, 1], [1, 0, 0, 0, 0, 1, 0, 0, 0, 1], []])
    A = C34CurveDivisor(C_2, [[0, 1], [1, 0, 1, 0, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_2_4, [[0, z4^2 + 1, 0, z4^3 + z4^2, 1], [0, 1, 0, z4^2, 0, 0, z4, 0, 0, 1], []])
    A = C34CurveDivisor(C_2_4, [[z4^3 + z4, 1], [z4, 0, z4^3 + z4, 0, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_3, [[2, 0, 2, 1, 1], [0, 0, 1, 0, 0, 1, 0, 0, 0, 1], []])
    A = C34CurveDivisor(C_3, [[2, 1], [1, 0, 1, 0, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_3_3, [[z3^2, z3, z3^2 + z3 + 2, z3^2 + 1, 1], [1, z3^2 + z3 + 1, 2*z3^2 + z3 + 2, z3^2 + 1, 0, z3^2 + 2*z3 + 1, 2, 0, 0, 1], []])
    A = C34CurveDivisor(C_3_3, [[z3 + 2, 1], [2*z3^2 + z3, 0, z3 + 2, 0, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_31, [[8, 11, 13, 14, 1], [19, 9, 23, 4, 0, 4, 12, 0, 0, 1], []])
    A = C34CurveDivisor(C_31, [[21, 1], [23, 0, 5, 0, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_31_2, [[12*z2 + 13, 29*z2, 19*z2 + 6, 22*z2 + 13, 1], [18*z2 + 20, 9*z2 + 17, 15*z2 + 28, 19*z2 + 30, 0, 12*z2 + 25, 23*z2 + 8, 0, 0, 1], []])
    A = C34CurveDivisor(C_31_2, [[30*z2 + 18, 1], [22*z2 + 21, 0, 3*z2 + 4, 0, 0, 1], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_1009, [[139, 967, 271, 711, 1], [648, 637, 239, 844, 0, 26, 402, 0, 0, 1], []])
    A = C34CurveDivisor(C_1009, [[418, 1], [983, 0, 894, 0, 0, 1], []])
    self.assertEqual(reduce(D), A)



  def test_reduce_65(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    D = C34CurveDivisor(C_2, [[1, 1, 0, 1], [], []])
    A = C34CurveDivisor(C_2, [[1], [], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_2_4, [[0, z4^3 + z4^2 + z4 + 1, z4^3 + 1, 1], [], []])
    A = C34CurveDivisor(C_2_4, [[1], [], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_3, [[2, 2, 1, 1], [], []])
    A = C34CurveDivisor(C_3, [[1], [], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_3_3, [[2*z3^2 + 2*z3 + 1, 2*z3^2 + z3, 0, 1], [], []])
    A = C34CurveDivisor(C_3_3, [[1], [], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_31, [[2, 30, 19, 1], [], []])
    A = C34CurveDivisor(C_31, [[1], [], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_31_2, [[10*z2 + 10, 11*z2 + 30, 27*z2 + 5, 1], [], []])
    A = C34CurveDivisor(C_31_2, [[1], [], []])
    self.assertEqual(reduce(D), A)

    D = C34CurveDivisor(C_1009, [[145, 974, 141, 1], [], []])
    A = C34CurveDivisor(C_1009, [[1], [], []])
    self.assertEqual(reduce(D), A)




