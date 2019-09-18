class TestFlip(unittest.TestCase) :
  def setUp(self) :
    self.z2 = GF(313^2).gen()
    self.z3 = GF(3^3).gen()
    self.z4 = GF(2^4).gen()
    
    self.C_2 = C34Curve(GF(2), [0, 0, 1, 0, 1, 1, 1, 1, 1])
    self.C_3 = C34Curve(GF(3), [0, 2, 1, 1, 0, 2, 2, 0, 2])
    self.C_11 = C34Curve(GF(11), [6, 5, 1, 8, 6, 3, 7, 10, 4]) # The vertical tangent line at P = (9 : 6 : 1) intersects P with multiplicity 2
    self.C_31 = C34Curve(GF(31), [2, 29, 13, 7, 28, 25, 17, 13, 17])
    self.C_41 = C34Curve(GF(41), [26, 22, 16, 37, 30, 22, 7, 22, 29]) # Vertical tangents at (31 : 33 : 1) and (28 : 22 : 1) intersect with m = 2
    self.C_1009 = C34Curve(GF(1009), [715, 286, 748, 7, 92, 446, 122, 935, 314])
    self.C_2_4 = C34Curve(GF(2^4), [z4^3 + 1, z4, z4^3 + z4^2 + 1, z4^2, z4^3 + z4, z4^3 + z4 + 1, z4^3 + z4, z4^3 + z4^2 + z4, z4^3 + 1])
    self.C_3_3 = C34Curve(GF(3^3), [z3 + 2, z3 + 1, z3^2 + z3 + 2, 2*z3^2 + z3 + 1, 2*z3^2 + z3, z3^2 + 2, 2*z3^2 + z3, 2*z3^2 + 1, 2])
    self.C_31_2 = C34Curve(GF(31^2), [9*z2 + 11, 11*z2 + 28, 17*z2 + 2, 24*z2 + 22, 6*z2 + 29, 21*z2 + 8, 18*z2 + 5, 18*z2 + 15, 13*z2 + 10])
  


  def test_flip_0(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    D = C34CurveDivisor(C_2, [[1], [], []])
    A = C34CurveDivisor(C_2, [[1], [], []])
    self.assertEqual(flip_0(D), A)

    D = C34CurveDivisor(C_2_4, [[1], [], []])
    A = C34CurveDivisor(C_2_4, [[1], [], []])
    self.assertEqual(flip_0(D), A)

    D = C34CurveDivisor(C_3, [[1], [], []])
    A = C34CurveDivisor(C_3, [[1], [], []])
    self.assertEqual(flip_0(D), A)

    D = C34CurveDivisor(C_3_3, [[1], [], []])
    A = C34CurveDivisor(C_3_3, [[1], [], []])
    self.assertEqual(flip_0(D), A)

    D = C34CurveDivisor(C_31, [[1], [], []])
    A = C34CurveDivisor(C_31, [[1], [], []])
    self.assertEqual(flip_0(D), A)

    D = C34CurveDivisor(C_31_2, [[1], [], []])
    A = C34CurveDivisor(C_31_2, [[1], [], []])
    self.assertEqual(flip_0(D), A)

    D = C34CurveDivisor(C_1009, [[1], [], []])
    A = C34CurveDivisor(C_1009, [[1], [], []])
    self.assertEqual(flip_0(D), A)



  def test_flip_11(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    D = C34CurveDivisor(C_2, [[0, 1], [0, 0, 1], []])
    A = C34CurveDivisor(C_2, [[0, 1], [1, 0, 1, 0, 0, 1], []])
    self.assertEqual(flip_11(D), A)

    D = C34CurveDivisor(C_2_4, [[z4^2 + z4 + 1, 1], [z4^2 + 1, 0, 1], []])
    A = C34CurveDivisor(C_2_4, [[z4^2 + z4 + 1, 1], [z4^3 + z4 + 1, 0, z4^2, 0, 0, 1], []])
    self.assertEqual(flip_11(D), A)

    D = C34CurveDivisor(C_3, [[2, 1], [0, 0, 1], []])
    A = C34CurveDivisor(C_3, [[2, 1], [1, 0, 1, 0, 0, 1], []])
    self.assertEqual(flip_11(D), A)

    D = C34CurveDivisor(C_3_3, [[z3^2 + 2*z3 + 2, 1], [2, 0, 1], []])
    A = C34CurveDivisor(C_3_3, [[z3^2 + 2*z3 + 2, 1], [z3^2 + z3 + 2, 0, 2*z3^2 + 2*z3 + 2, 0, 0, 1], []])
    self.assertEqual(flip_11(D), A)

    D = C34CurveDivisor(C_31, [[30, 1], [25, 0, 1], []])
    A = C34CurveDivisor(C_31, [[30, 1], [1, 0, 17, 0, 0, 1], []])
    self.assertEqual(flip_11(D), A)

    D = C34CurveDivisor(C_31_2, [[29*z2, 1], [27*z2 + 4, 0, 1], []])
    A = C34CurveDivisor(C_31_2, [[29*z2, 1], [11*z2 + 5, 0, 4*z2 + 19, 0, 0, 1], []])
    self.assertEqual(flip_11(D), A)

    D = C34CurveDivisor(C_1009, [[125, 1], [437, 0, 1], []])
    A = C34CurveDivisor(C_1009, [[125, 1], [773, 0, 110, 0, 0, 1], []])
    self.assertEqual(flip_11(D), A)



  def test_flip_21(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    D = C34CurveDivisor(C_2, [[0, 0, 1], [0, 0, 0, 1], []])
    A = C34CurveDivisor(C_2, [[0, 0, 1], [0, 1, 0, 1], []])
    self.assertEqual(flip_21(D), A)

    D = C34CurveDivisor(C_2_4, [[1, z4^3 + 1, 1], [z4^3 + z4^2, z4^2 + 1, 0, 1], []])
    A = C34CurveDivisor(C_2_4, [[1, z4^3 + 1, 1], [z4^2 + z4, z4^3, 0, 1], []])
    self.assertEqual(flip_21(D), A)

    D = C34CurveDivisor(C_3, [[1, 1, 1], [0, 2, 0, 1], []])
    A = C34CurveDivisor(C_3, [[1, 1, 1], [2, 1, 0, 1], []])
    self.assertEqual(flip_21(D), A)

    D = C34CurveDivisor(C_3_3, [[z3^2 + 2*z3, z3 + 1, 1], [z3 + 2, z3 + 1, 0, 1], []])
    A = C34CurveDivisor(C_3_3, [[z3^2 + 2*z3, z3 + 1, 1], [z3^2 + 1, 2*z3^2 + 2, 0, 1], []])
    self.assertEqual(flip_21(D), A)

    D = C34CurveDivisor(C_31, [[14, 10, 1], [12, 29, 0, 1], []])
    A = C34CurveDivisor(C_31, [[14, 10, 1], [20, 0, 0, 1], []])
    self.assertEqual(flip_21(D), A)

    D = C34CurveDivisor(C_31_2, [[10*z2 + 23, 9*z2 + 23, 1], [12*z2 + 12, 13*z2 + 23, 0, 1], []])
    A = C34CurveDivisor(C_31_2, [[10*z2 + 23, 9*z2 + 23, 1], [23*z2 + 12, 2, 0, 1], []])
    self.assertEqual(flip_21(D), A)

    D = C34CurveDivisor(C_1009, [[950, 71, 1], [173, 58, 0, 1], []])
    A = C34CurveDivisor(C_1009, [[950, 71, 1], [523, 310, 0, 1], []])
    self.assertEqual(flip_21(D), A)



  def test_flip_22(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    D = C34CurveDivisor(C_2, [[1, 1], [0, 0, 1, 0, 0, 1], []])
    A = C34CurveDivisor(C_2, [[1, 1], [1, 0, 1], []])
    self.assertEqual(flip_22(D), A)

    D = C34CurveDivisor(C_2_4, [[z4^3, 1], [z4^3, 0, z4^3 + z4^2 + z4, 0, 0, 1], []])
    A = C34CurveDivisor(C_2_4, [[z4^3, 1], [1, 0, 1], []])
    self.assertEqual(flip_22(D), A)

    D = C34CurveDivisor(C_3, [[2, 1], [0, 0, 2, 0, 0, 1], []])
    A = C34CurveDivisor(C_3, [[2, 1], [2, 0, 1], []])
    self.assertEqual(flip_22(D), A)

    D = C34CurveDivisor(C_3_3, [[z3^2 + z3 + 2, 1], [2*z3^2 + z3 + 1, 0, z3^2 + z3 + 1, 0, 0, 1], []])
    A = C34CurveDivisor(C_3_3, [[z3^2 + z3 + 2, 1], [z3^2, 0, 1], []])
    self.assertEqual(flip_22(D), A)

    D = C34CurveDivisor(C_31, [[7, 1], [28, 0, 19, 0, 0, 1], []])
    A = C34CurveDivisor(C_31, [[7, 1], [11, 0, 1], []])
    self.assertEqual(flip_22(D), A)

    D = C34CurveDivisor(C_31_2, [[15*z2 + 2, 1], [16*z2 + 26, 0, 10*z2 + 29, 0, 0, 1], []])
    A = C34CurveDivisor(C_31_2, [[15*z2 + 2, 1], [3*z2 + 17, 0, 1], []])
    self.assertEqual(flip_22(D), A)

    D = C34CurveDivisor(C_1009, [[590, 1], [801, 0, 867, 0, 0, 1], []])
    A = C34CurveDivisor(C_1009, [[590, 1], [984, 0, 1], []])
    self.assertEqual(flip_22(D), A)



  def test_flip_31(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    D = C34CurveDivisor(C_2, [[0, 1, 1, 1], [1, 0, 1, 0, 1], [0, 1, 0, 0, 0, 1]])
    A = C34CurveDivisor(C_2, [[0, 1, 1, 1], [0, 0, 1, 0, 1], [0, 0, 0, 0, 0, 1]])
    self.assertEqual(flip_31(D), A)

    D = C34CurveDivisor(C_2_4, [[z4^3 + z4, z4^2 + z4, z4^3 + z4^2 + z4 + 1, 1], [z4^2, z4^2, z4^2 + z4, 0, 1], [z4^3 + 1, 1, z4^3 + z4 + 1, 0, 0, 1]])
    A = C34CurveDivisor(C_2_4, [[z4^3 + z4, z4^2 + z4, z4^3 + z4^2 + z4 + 1, 1], [z4^3 + z4 + 1, z4^3 + z4^2, z4^3 + z4^2 + z4, 0, 1], [z4^2 + 1, 1, z4^2 + z4 + 1, 0, 0, 1]])
    self.assertEqual(flip_31(D), A)

    D = C34CurveDivisor(C_3, [[0, 2, 2, 1], [2, 1, 1, 0, 1], [2, 1, 2, 0, 0, 1]])
    A = C34CurveDivisor(C_3, [[0, 2, 2, 1], [0, 2, 0, 0, 1], [0, 0, 2, 0, 0, 1]])
    self.assertEqual(flip_31(D), A)

    D = C34CurveDivisor(C_3_3, [[z3^2 + z3 + 2, z3^2 + z3 + 1, z3^2 + z3 + 1, 1], [z3^2 + 2*z3 + 2, 2*z3^2 + z3, z3^2 + 1, 0, 1], [2*z3^2 + 2*z3 + 2, z3^2 + 1, z3^2 + 1, 0, 0, 1]])
    A = C34CurveDivisor(C_3_3, [[z3^2 + z3 + 2, z3^2 + z3 + 1, z3^2 + z3 + 1, 1], [2, 2*z3 + 2, z3^2 + 2*z3 + 1, 0, 1], [2*z3^2 + z3 + 2, z3^2 + z3, 2*z3^2, 0, 0, 1]])
    self.assertEqual(flip_31(D), A)

    D = C34CurveDivisor(C_31, [[6, 9, 6, 1], [18, 28, 17, 0, 1], [21, 4, 0, 0, 0, 1]])
    A = C34CurveDivisor(C_31, [[6, 9, 6, 1], [29, 4, 14, 0, 1], [23, 20, 27, 0, 0, 1]])
    self.assertEqual(flip_31(D), A)

    D = C34CurveDivisor(C_31_2, [[28*z2 + 15, 5*z2 + 18, z2 + 26, 1], [5*z2 + 4, 26*z2 + 28, 15*z2 + 30, 0, 1], [3*z2 + 18, 9*z2 + 18, 21*z2 + 15, 0, 0, 1]])
    A = C34CurveDivisor(C_31_2, [[28*z2 + 15, 5*z2 + 18, z2 + 26, 1], [5, 26, 19*z2 + 15, 0, 1], [21*z2 + 24, 5*z2 + 13, 15*z2 + 9, 0, 0, 1]])
    self.assertEqual(flip_31(D), A)

    D = C34CurveDivisor(C_1009, [[842, 319, 784, 1], [106, 961, 986, 0, 1], [556, 879, 864, 0, 0, 1]])
    A = C34CurveDivisor(C_1009, [[842, 319, 784, 1], [964, 419, 362, 0, 1], [490, 97, 512, 0, 0, 1]])
    self.assertEqual(flip_31(D), A)
    
    # Test case where <f, g, h> =/= <f, g>, but <f, g, h> = <f, h>
    D = C34CurveDivisor(C_2, [[0, 1, 0, 1], [0, 0, 1, 0, 1], [0, 0, 1, 0, 0, 1]])
    A = C34CurveDivisor(C_2, [[0, 1, 0, 1], [0, 1, 0, 0, 1], [1, 1, 1, 0, 0, 1]])
    self.assertEqual(flip_31(D), A)

    D = C34CurveDivisor(C_2_4, [[z4^3 + 1, z4^3 + z4^2 + z4 + 1, 0, 1], [z4^2, z4 + 1, z4^3 + z4^2 + 1, 0, 1], [z4^2 + z4 + 1, z4^3 + 1, 1, 0, 0, 1]])
    A = C34CurveDivisor(C_2_4, [[z4^3 + 1, z4^3 + z4^2 + z4 + 1, 0, 1], [z4^3 + z4, z4^2 + 1, z4, 0, 1], [z4^2 + 1, z4^2 + 1, z4^3 + 1, 0, 0, 1]])
    self.assertEqual(flip_31(D), A)

    D = C34CurveDivisor(C_3, [[0, 1, 0, 1], [0, 2, 0, 0, 1], [0, 2, 1, 0, 0, 1]])
    A = C34CurveDivisor(C_3, [[0, 1, 0, 1], [1, 1, 1, 0, 1], [0, 1, 1, 0, 0, 1]])
    self.assertEqual(flip_31(D), A)

    D = C34CurveDivisor(C_3_3, [[2*z3^2 + z3 + 2, 2*z3, 0, 1], [z3^2 + z3 + 1, z3^2 + 1, z3^2 + z3 + 2, 0, 1], [2, z3^2 + 2*z3 + 2, z3^2 + z3 + 1, 0, 0, 1]])
    A = C34CurveDivisor(C_3_3, [[2*z3^2 + z3 + 2, 2*z3, 0, 1], [2*z3 + 2, z3^2, 2*z3^2 + z3 + 1, 0, 1], [2*z3^2, 2*z3, 2*z3^2 + z3 + 2, 0, 0, 1]])
    self.assertEqual(flip_31(D), A)

    D = C34CurveDivisor(C_31, [[11, 0, 0, 1], [21, 6, 19, 0, 1], [10, 14, 21, 0, 0, 1]])
    A = C34CurveDivisor(C_31, [[11, 0, 0, 1], [16, 22, 12, 0, 1], [29, 3, 1, 0, 0, 1]])
    self.assertEqual(flip_31(D), A)

    D = C34CurveDivisor(C_31_2, [[14*z2 + 15, 15, 0, 1], [2*z2 + 27, 7*z2 + 28, 4*z2 + 25, 0, 1], [5*z2 + 25, 20*z2 + 8, 13*z2, 0, 0, 1]])
    A = C34CurveDivisor(C_31_2, [[14*z2 + 15, 15, 0, 1], [24*z2 + 9, 4*z2 + 7, 27*z2 + 21, 0, 1], [11*z2 + 27, 12*z2 + 30, 9*z2 + 17, 0, 0, 1]])
    self.assertEqual(flip_31(D), A)

    D = C34CurveDivisor(C_1009, [[437, 526, 0, 1], [997, 792, 107, 0, 1], [791, 629, 501, 0, 0, 1]])
    A = C34CurveDivisor(C_1009, [[437, 526, 0, 1], [168, 653, 419, 0, 1], [905, 365, 267, 0, 0, 1]])
    self.assertEqual(flip_31(D), A)
    
    # Test case where <f, g> =/= <f, g, h> =/= <f, h>
    D = C34CurveDivisor(C_2, [[0, 1, 0, 1], [0, 0, 0, 0, 1], [1, 1, 1, 0, 0, 1]])
    A = C34CurveDivisor(C_2, [[0, 1, 0, 1], [0, 0, 1, 0, 1], [0, 1, 0, 0, 0, 1]])
    self.assertEqual(flip_31(D), A)

    D = C34CurveDivisor(C_2_4, [[1, z4^3 + z4^2 + z4, 0, 1], [z4^3 + z4, z4^2, z4^3 + z4 + 1, 0, 1], [z4^3 + 1, z4^3, 1, 0, 0, 1]])
    A = C34CurveDivisor(C_2_4, [[1, z4^3 + z4^2 + z4, 0, 1], [z4^3 + z4^2 + 1, z4^2 + z4, z4^2 + 1, 0, 1], [z4^3 + z4 + 1, z4^3, z4^2, 0, 0, 1]])
    self.assertEqual(flip_31(D), A)

    D = C34CurveDivisor(C_3, [[2, 0, 0, 1], [2, 2, 1, 0, 1], [0, 1, 1, 0, 0, 1]])
    A = C34CurveDivisor(C_3, [[2, 0, 0, 1], [1, 2, 2, 0, 1], [0, 0, 2, 0, 0, 1]])
    self.assertEqual(flip_31(D), A)

    D = C34CurveDivisor(C_3_3, [[z3^2 + z3, 1, 0, 1], [2*z3^2 + z3 + 2, 2*z3 + 1, 2*z3^2, 0, 1], [z3^2 + z3 + 2, z3^2 + 2*z3, z3^2 + z3 + 1, 0, 0, 1]])
    A = C34CurveDivisor(C_3_3, [[z3^2 + z3, 1, 0, 1], [2*z3^2 + 2*z3 + 2, 2*z3^2 + 2*z3 + 1, z3^2 + 1, 0, 1], [z3^2, z3^2 + 1, 2*z3^2 + z3 + 2, 0, 0, 1]])
    self.assertEqual(flip_31(D), A)

    D = C34CurveDivisor(C_31, [[23, 24, 0, 1], [4, 15, 23, 0, 1], [29, 11, 10, 0, 0, 1]])
    A = C34CurveDivisor(C_31, [[23, 24, 0, 1], [27, 27, 1, 0, 1], [24, 14, 24, 0, 0, 1]])
    self.assertEqual(flip_31(D), A)

    D = C34CurveDivisor(C_31_2, [[22*z2 + 18, 14*z2 + 2, 0, 1], [23*z2 + 3, 18*z2 + 16, 11*z2 + 16, 0, 1], [11*z2 + 27, 6*z2 + 10, 11*z2 + 12, 0, 0, 1]])
    A = C34CurveDivisor(C_31_2, [[22*z2 + 18, 14*z2 + 2, 0, 1], [29*z2 + 24, 26*z2 + 17, 3*z2 + 17, 0, 1], [11*z2 + 1, 27*z2 + 23, 15*z2 + 1, 0, 0, 1]])
    self.assertEqual(flip_31(D), A)

    D = C34CurveDivisor(C_1009, [[161, 370, 0, 1], [828, 867, 790, 0, 1], [951, 982, 862, 0, 0, 1]])
    A = C34CurveDivisor(C_1009, [[161, 370, 0, 1], [59, 747, 589, 0, 1], [707, 987, 289, 0, 0, 1]])
    self.assertEqual(flip_31(D), A)



  def test_flip_32(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    D = C34CurveDivisor(C_2, [[0, 1, 1], [1, 0, 0, 0, 0, 0, 1], []])
    A = C34CurveDivisor(C_2, [[0, 1], [0, 0, 1], []])
    self.assertEqual(flip_32(D), A)

    D = C34CurveDivisor(C_2_4, [[z4, z4^3 + z4^2 + z4, 1], [z4^3 + z4 + 1, z4^3 + z4^2 + 1, 0, z4^3 + z4^2 + z4 + 1, 0, 0, 1], []])
    A = C34CurveDivisor(C_2_4, [[z4^3 + z4, 1], [z4^2, 0, 1], []])
    self.assertEqual(flip_32(D), A)

    D = C34CurveDivisor(C_3, [[1, 2, 1], [1, 2, 0, 2, 0, 0, 1], []])
    A = C34CurveDivisor(C_3, [[0, 1], [1, 0, 1], []])
    self.assertEqual(flip_32(D), A)

    D = C34CurveDivisor(C_3_3, [[2*z3^2 + z3 + 1, 0, 1], [2*z3^2 + z3, 2*z3 + 2, 0, 0, 0, 0, 1], []])
    A = C34CurveDivisor(C_3_3, [[2*z3^2 + z3, 1], [2*z3^2 + z3 + 1, 0, 1], []])
    self.assertEqual(flip_32(D), A)

    D = C34CurveDivisor(C_31, [[22, 22, 1], [25, 12, 0, 27, 0, 0, 1], []])
    A = C34CurveDivisor(C_31, [[12, 1], [6, 0, 1], []])
    self.assertEqual(flip_32(D), A)

    D = C34CurveDivisor(C_31_2, [[2*z2 + 29, 4*z2 + 6, 1], [8*z2 + 30, 6*z2 + 22, 0, 26*z2 + 13, 0, 0, 1], []])
    A = C34CurveDivisor(C_31_2, [[30*z2 + 18, 1], [6*z2 + 2, 0, 1], []])
    self.assertEqual(flip_32(D), A)

    D = C34CurveDivisor(C_1009, [[953, 707, 1], [841, 449, 0, 801, 0, 0, 1], []])
    A = C34CurveDivisor(C_1009, [[724, 1], [648, 0, 1], []])
    self.assertEqual(flip_32(D), A)



  def test_flip_33(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    D = C34CurveDivisor(C_2, [[1, 1], [], []])
    A = C34CurveDivisor(C_2, [[1], [], []])
    self.assertEqual(flip_33(D), A)

    D = C34CurveDivisor(C_2_4, [[1, 1], [], []])
    A = C34CurveDivisor(C_2_4, [[1], [], []])
    self.assertEqual(flip_33(D), A)

    D = C34CurveDivisor(C_3, [[0, 1], [], []])
    A = C34CurveDivisor(C_3, [[1], [], []])
    self.assertEqual(flip_33(D), A)

    D = C34CurveDivisor(C_3_3, [[0, 1], [], []])
    A = C34CurveDivisor(C_3_3, [[1], [], []])
    self.assertEqual(flip_33(D), A)

    D = C34CurveDivisor(C_31, [[21, 1], [], []])
    A = C34CurveDivisor(C_31, [[1], [], []])
    self.assertEqual(flip_33(D), A)

    D = C34CurveDivisor(C_31_2, [[29*z2, 1], [], []])
    A = C34CurveDivisor(C_31_2, [[1], [], []])
    self.assertEqual(flip_33(D), A)

    D = C34CurveDivisor(C_1009, [[959, 1], [], []])
    A = C34CurveDivisor(C_1009, [[1], [], []])
    self.assertEqual(flip_33(D), A)



  def test_flip_41(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    # Test case where <f, g, h> = <f, g>
    D = C34CurveDivisor(C_2, [[0, 0, 1, 0, 1], [0, 0, 0, 1, 0, 1], [0, 0, 0, 1, 0, 0, 1]])
    A = C34CurveDivisor(C_2, [[0, 1, 1, 1], [0, 0, 1, 0, 1], [0, 0, 0, 0, 0, 1]])
    self.assertEqual(flip_41(D), A)

    D = C34CurveDivisor(C_2_4, [[1, z4^3 + z4, z4^2 + z4, z4^2 + 1, 1], [z4^3 + z4^2 + z4, z4, z4^2 + z4, z4^2 + 1, 0, 1], [z4^3 + z4, z4^3, 0, z4, 0, 0, 1]])
    A = C34CurveDivisor(C_2_4, [[z4^3 + 1, z4^3 + z4^2 + z4, z4^2 + z4 + 1, 1], [z4^3 + z4, z4^3 + 1, z4^3 + z4^2 + z4, 0, 1], [z4^3, z4^3, z4^3 + z4, 0, 0, 1]])
    self.assertEqual(flip_41(D), A)

    D = C34CurveDivisor(C_3, [[1, 1, 2, 0, 1], [1, 1, 2, 2, 0, 1], [2, 1, 2, 1, 0, 0, 1]])
    A = C34CurveDivisor(C_3, [[0, 0, 2, 1], [1, 1, 2, 0, 1], [1, 2, 0, 0, 0, 1]])
    self.assertEqual(flip_41(D), A)

    D = C34CurveDivisor(C_3_3, [[z3^2 + z3, 2*z3^2 + 2, 0, z3 + 1, 1], [2*z3^2 + 2, 2*z3^2 + 1, 2*z3, z3^2 + z3 + 1, 0, 1], [z3^2 + 2*z3 + 2, 2*z3 + 2, 2*z3^2, 2, 0, 0, 1]])
    A = C34CurveDivisor(C_3_3, [[2*z3^2 + 2*z3 + 2, z3^2 + 2*z3 + 1, 2*z3^2 + 2, 1], [z3, 2*z3^2 + 2*z3 + 2, z3^2 + 2*z3, 0, 1], [2*z3^2 + 1, z3 + 1, z3 + 2, 0, 0, 1]])
    self.assertEqual(flip_41(D), A)

    D = C34CurveDivisor(C_31, [[21, 8, 4, 11, 1], [1, 18, 6, 5, 0, 1], [26, 15, 26, 18, 0, 0, 1]])
    A = C34CurveDivisor(C_31, [[24, 19, 2, 1], [5, 16, 13, 0, 1], [22, 24, 20, 0, 0, 1]])
    self.assertEqual(flip_41(D), A)

    D = C34CurveDivisor(C_31_2, [[z2 + 23, 2*z2 + 13, 2*z2 + 5, 16*z2 + 9, 1], [29*z2 + 16, 18*z2 + 16, 12*z2 + 15, 27*z2 + 23, 0, 1], [30*z2 + 24, 23*z2 + 4, 6*z2 + 3, 26*z2 + 7, 0, 0, 1]])
    A = C34CurveDivisor(C_31_2, [[12*z2 + 11, 24*z2 + 13, 21*z2 + 18, 1], [15*z2 + 4, 19*z2 + 25, 14, 0, 1], [29*z2 + 15, 22*z2 + 17, 15*z2 + 21, 0, 0, 1]])
    self.assertEqual(flip_41(D), A)

    D = C34CurveDivisor(C_1009, [[930, 542, 132, 786, 1], [670, 525, 754, 180, 0, 1], [206, 306, 653, 402, 0, 0, 1]])
    A = C34CurveDivisor(C_1009, [[368, 345, 468, 1], [256, 793, 569, 0, 1], [427, 444, 83, 0, 0, 1]])
    self.assertEqual(flip_41(D), A)

    # Test case where <f, g> =/= <f, g, h> = <f, h>
    D = C34CurveDivisor(C_2, [[0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 1], [0, 0, 1, 0, 0, 0, 1]])
    A = C34CurveDivisor(C_2, [[0, 1, 0, 1], [0, 0, 0, 0, 1], [1, 1, 1, 0, 0, 1]])
    self.assertEqual(flip_41(D), A)

    D = C34CurveDivisor(C_2_4, [[z4^3 + 1, z4^3 + z4^2 + z4, z4^3 + z4 + 1, z4, 1], [z4^3 + z4^2 + z4, z4 + 1, z4^3, z4^2, 0, 1], [z4^3 + z4 + 1, 0, z4^3 + z4^2 + z4, z4^3 + z4, 0, 0, 1]])
    A = C34CurveDivisor(C_2_4, [[1, z4^3 + z4^2 + z4, 0, 1], [z4^3 + z4 + 1, 1, z4^3 + z4 + 1, 0, 1], [z4^2 + 1, 1, 1, 0, 0, 1]])
    self.assertEqual(flip_41(D), A)

    D = C34CurveDivisor(C_3, [[2, 2, 2, 2, 1], [2, 2, 0, 2, 0, 1], [2, 2, 2, 1, 0, 0, 1]])
    A = C34CurveDivisor(C_3, [[2, 0, 0, 1], [1, 2, 2, 0, 1], [1, 0, 1, 0, 0, 1]])
    self.assertEqual(flip_41(D), A)

    D = C34CurveDivisor(C_3_3, [[z3^2 + 2*z3 + 2, z3^2 + 1, z3^2 + 2*z3 + 1, z3^2 + 2*z3 + 2, 1], [2*z3^2, z3^2, 2*z3^2 + 1, z3, 0, 1], [z3, z3^2 + 2, 1, z3, 0, 0, 1]])
    A = C34CurveDivisor(C_3_3, [[z3^2 + 2*z3 + 1, z3^2 + 2*z3 + 2, 0, 1], [2*z3^2 + 2*z3 + 1, z3^2 + z3 + 1, z3^2 + 2*z3 + 1, 0, 1], [2*z3^2 + 1, z3 + 1, 2*z3 + 1, 0, 0, 1]])
    self.assertEqual(flip_41(D), A)

    D = C34CurveDivisor(C_31, [[8, 6, 23, 3, 1], [11, 1, 19, 22, 0, 1], [20, 1, 19, 2, 0, 0, 1]])
    A = C34CurveDivisor(C_31, [[21, 1, 0, 1], [7, 3, 23, 0, 1], [8, 2, 10, 0, 0, 1]])
    self.assertEqual(flip_41(D), A)

    D = C34CurveDivisor(C_31_2, [[22*z2 + 26, z2 + 20, 19*z2 + 25, 15*z2 + 29, 1], [15*z2 + 15, 19*z2 + 10, 19*z2 + 2, 13*z2 + 20, 0, 1], [9*z2 + 10, 13*z2 + 30, 10*z2 + 23, 13*z2 + 16, 0, 0, 1]])
    A = C34CurveDivisor(C_31_2, [[24*z2, 6*z2 + 9, 0, 1], [z2 + 21, 8*z2 + 29, 19*z2 + 25, 0, 1], [22*z2 + 6, 29, z2 + 2, 0, 0, 1]])
    self.assertEqual(flip_41(D), A)

    D = C34CurveDivisor(C_1009, [[86, 142, 131, 801, 1], [111, 499, 539, 123, 0, 1], [511, 818, 366, 771, 0, 0, 1]])
    A = C34CurveDivisor(C_1009, [[414, 604, 0, 1], [433, 658, 131, 0, 1], [527, 740, 294, 0, 0, 1]])
    self.assertEqual(flip_41(D), A)

    # Test case where <f, g> =/= <f, g, h> =/= <f, h>
    D = C34CurveDivisor(C_2, [[0, 0, 1, 0, 1], [0, 0, 1, 0, 0, 1], [0, 0, 0, 1, 0, 0, 1]])
    A = C34CurveDivisor(C_2, [[0, 1, 0, 1], [0, 0, 1, 0, 1], [0, 0, 1, 0, 0, 1]])
    self.assertEqual(flip_41(D), A)

    D = C34CurveDivisor(C_2_4, [[1, z4^2 + 1, z4^2 + 1, z4^3, 1], [z4 + 1, z4^2 + z4, z4^2, z4^3 + z4^2, 0, 1], [z4^3 + z4 + 1, 0, 0, z4^3 + 1, 0, 0, 1]])
    A = C34CurveDivisor(C_2_4, [[z4^3 + z4^2 + 1, z4 + 1, 0, 1], [z4 + 1, z4^3 + z4^2 + z4, z4^2 + 1, 0, 1], [1, z4^3 + z4^2 + z4 + 1, 1, 0, 0, 1]])
    self.assertEqual(flip_41(D), A)

    D = C34CurveDivisor(C_3, [[0, 1, 0, 1, 1], [1, 0, 2, 2, 0, 1], [0, 2, 0, 1, 0, 0, 1]])
    A = C34CurveDivisor(C_3, [[0, 2, 0, 1], [0, 2, 0, 0, 1], [0, 1, 1, 0, 0, 1]])
    self.assertEqual(flip_41(D), A)

    D = C34CurveDivisor(C_3_3, [[z3^2 + z3 + 1, 2, 2*z3 + 1, 2*z3^2 + 2*z3 + 2, 1], [2*z3^2, z3^2 + z3, 2*z3^2 + 2*z3 + 2, 2*z3^2 + 1, 0, 1], [z3^2 + 2*z3, z3^2 + z3 + 2, 0, 2*z3^2 + 1, 0, 0, 1]])
    A = C34CurveDivisor(C_3_3, [[z3^2 + 2*z3 + 2, 2*z3^2 + 1, 0, 1], [z3^2 + z3, 2*z3^2 + z3 + 1, 2*z3 + 1, 0, 1], [2*z3^2 + 2*z3 + 1, z3^2 + 2*z3 + 2, 2*z3^2 + z3 + 1, 0, 0, 1]])
    self.assertEqual(flip_41(D), A)

    D = C34CurveDivisor(C_31, [[29, 24, 11, 25, 1], [28, 7, 29, 26, 0, 1], [4, 2, 0, 5, 0, 0, 1]])
    A = C34CurveDivisor(C_31, [[23, 30, 0, 1], [12, 18, 11, 0, 1], [25, 14, 27, 0, 0, 1]])
    self.assertEqual(flip_41(D), A)

    D = C34CurveDivisor(C_31_2, [[6*z2 + 10, 23*z2 + 2, 16*z2, 19*z2 + 11, 1], [23*z2 + 2, 4*z2 + 20, 13*z2 + 2, 7*z2 + 1, 0, 1], [16*z2 + 25, 6*z2 + 12, 0, 24*z2 + 19, 0, 0, 1]])
    A = C34CurveDivisor(C_31_2, [[29*z2 + 2, 25*z2 + 9, 0, 1], [4*z2 + 29, 22*z2 + 26, 16*z2, 0, 1], [17*z2 + 14, 19*z2 + 5, 26*z2 + 7, 0, 0, 1]])
    self.assertEqual(flip_41(D), A)

    D = C34CurveDivisor(C_1009, [[5, 153, 631, 933, 1], [615, 703, 606, 278, 0, 1], [826, 654, 0, 313, 0, 0, 1]])
    A = C34CurveDivisor(C_1009, [[675, 413, 0, 1], [855, 262, 631, 0, 1], [64, 661, 689, 0, 0, 1]])
    self.assertEqual(flip_41(D), A)



  def test_flip_42(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    D = C34CurveDivisor(C_2, [[1, 0, 0, 1], [0, 0, 1, 0, 1], []])
    A = C34CurveDivisor(C_2, [[1, 1], [1, 0, 0, 0, 0, 1], []])
    self.assertEqual(flip_42(D), A)

    D = C34CurveDivisor(C_2_4, [[0, z4^3 + z4^2, 0, 1], [z4^3, z4^3 + z4^2 + z4 + 1, z4^3 + z4^2, 0, 1], []])
    A = C34CurveDivisor(C_2_4, [[0, 1], [z4^2, 0, z4^2, 0, 0, 1], []])
    self.assertEqual(flip_42(D), A)

    D = C34CurveDivisor(C_3, [[2, 0, 0, 1], [0, 0, 1, 0, 1], []])
    A = C34CurveDivisor(C_3, [[2, 1], [1, 0, 1, 0, 0, 1], []])
    self.assertEqual(flip_42(D), A)

    D = C34CurveDivisor(C_3_3, [[z3^2 + 2*z3, 2*z3^2 + z3 + 2, 0, 1], [z3^2 + 2*z3 + 2, 2*z3^2 + z3 + 1, 2, 0, 1], []])
    A = C34CurveDivisor(C_3_3, [[2*z3^2 + z3, 1], [z3^2 + 2, 0, z3^2 + 1, 0, 0, 1], []])
    self.assertEqual(flip_42(D), A)

    D = C34CurveDivisor(C_31, [[6, 30, 0, 1], [3, 27, 7, 0, 1], []])
    A = C34CurveDivisor(C_31, [[23, 1], [24, 0, 10, 0, 0, 1], []])
    self.assertEqual(flip_42(D), A)

    D = C34CurveDivisor(C_31_2, [[20*z2 + 29, z2 + 9, 0, 1], [14*z2 + 14, 21, 11*z2 + 11, 0, 1], []])
    A = C34CurveDivisor(C_31_2, [[21*z2 + 29, 1], [13*z2 + 9, 0, 4*z2 + 20, 0, 0, 1], []])
    self.assertEqual(flip_42(D), A)

    D = C34CurveDivisor(C_1009, [[75, 272, 0, 1], [566, 991, 417, 0, 1], []])
    A = C34CurveDivisor(C_1009, [[864, 1], [502, 0, 589, 0, 0, 1], []])
    self.assertEqual(flip_42(D), A)



  def test_flip_43(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    D = C34CurveDivisor(C_2, [[1, 0, 1, 1], [1, 1, 1, 0, 0, 1], []])
    A = C34CurveDivisor(C_2, [[0, 1, 1], [1, 1, 0, 1], []])
    self.assertEqual(flip_43(D), A)

    D = C34CurveDivisor(C_2_4, [[z4^3 + z4^2 + z4, 0, z4^3 + 1, 1], [z4^2, z4^2 + z4 + 1, z4^3 + 1, 0, 1, 1], []])
    A = C34CurveDivisor(C_2_4, [[z4^3 + z4^2, z4^3, 1], [z4^3, z4^2, 0, 1], []])
    self.assertEqual(flip_43(D), A)

    D = C34CurveDivisor(C_3, [[0, 1, 0, 1], [1, 2, 2, 0, 1, 1], []])
    A = C34CurveDivisor(C_3, [[0, 1, 1], [0, 1, 0, 1], []])
    self.assertEqual(flip_43(D), A)

    D = C34CurveDivisor(C_3_3, [[z3 + 2, 2*z3^2 + 2*z3 + 2, z3^2 + 2*z3, 1], [2*z3^2 + 2*z3 + 2, 0, 2*z3, 0, 2*z3^2, 1], []])
    A = C34CurveDivisor(C_3_3, [[2*z3^2 + z3 + 1, z3^2 + 2, 1], [z3^2 + 2*z3 + 1, 2*z3^2 + 1, 0, 1], []])
    self.assertEqual(flip_43(D), A)

    D = C34CurveDivisor(C_31, [[18, 15, 21, 1], [13, 22, 19, 0, 19, 1], []])
    A = C34CurveDivisor(C_31, [[27, 29, 1], [9, 26, 0, 1], []])
    self.assertEqual(flip_43(D), A)

    D = C34CurveDivisor(C_31_2, [[26*z2 + 9, 30*z2 + 28, 8*z2 + 14, 1], [14*z2 + 16, 23*z2 + 24, z2 + 22, 0, 19*z2 + 2, 1], []])
    A = C34CurveDivisor(C_31_2, [[29*z2 + 20, 25*z2 + 8, 1], [19*z2 + 22, 22*z2 + 20, 0, 1], []])
    self.assertEqual(flip_43(D), A)

    D = C34CurveDivisor(C_1009, [[711, 254, 811, 1], [323, 304, 787, 0, 517, 1], []])
    A = C34CurveDivisor(C_1009, [[947, 806, 1], [543, 420, 0, 1], []])
    self.assertEqual(flip_43(D), A)



  def test_flip_44(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    D = C34CurveDivisor(C_2, [[1, 1, 1], [], []])
    A = C34CurveDivisor(C_2, [[1], [], []])
    self.assertEqual(flip_44(D), A)

    D = C34CurveDivisor(C_2_4, [[z4 + 1, z4^2 + z4, 1], [], []])
    A = C34CurveDivisor(C_2_4, [[1], [], []])
    self.assertEqual(flip_44(D), A)

    D = C34CurveDivisor(C_3, [[0, 0, 1], [], []])
    A = C34CurveDivisor(C_3, [[1], [], []])
    self.assertEqual(flip_44(D), A)

    D = C34CurveDivisor(C_3_3, [[z3^2, 2*z3^2 + z3 + 2, 1], [], []])
    A = C34CurveDivisor(C_3_3, [[1], [], []])
    self.assertEqual(flip_44(D), A)

    D = C34CurveDivisor(C_31, [[24, 28, 1], [], []])
    A = C34CurveDivisor(C_31, [[1], [], []])
    self.assertEqual(flip_44(D), A)

    D = C34CurveDivisor(C_31_2, [[6*z2 + 11, 5*z2 + 3, 1], [], []])
    A = C34CurveDivisor(C_31_2, [[1], [], []])
    self.assertEqual(flip_44(D), A)

    D = C34CurveDivisor(C_1009, [[990, 764, 1], [], []])
    A = C34CurveDivisor(C_1009, [[1], [], []])
    self.assertEqual(flip_44(D), A)



  def test_flip_51(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    # Test case where <f, g, h> = <f, g>
    D = C34CurveDivisor(C_2, [[0, 1, 0, 1, 0, 1], [0, 0, 1, 1, 1, 0, 1], [0, 1, 0, 1, 0, 0, 0, 1]])
    A = C34CurveDivisor(C_2, [[1, 0, 1, 1], [0, 1, 1, 0, 1], [1, 1, 1, 0, 0, 1]])
    self.assertEqual(flip_51(D), A)

    D = C34CurveDivisor(C_2_4, [[z4^3 + z4 + 1, z4^3 + z4^2 + z4 + 1, z4^3, z4^2 + z4, z4^3 + z4^2, 1], [z4^3 + z4^2 + z4, 0, 0, z4 + 1, z4^2 + z4 + 1, 0, 1], [z4 + 1, z4^3 + z4 + 1, z4^2 + 1, z4^3 + z4^2 + z4 + 1, z4^3 + 1, 0, 0, 1]])
    A = C34CurveDivisor(C_2_4, [[z4^2, z4^3 + z4^2, z4^2 + z4, 1], [z4^3 + z4^2, z4^3 + z4 + 1, z4^3 + z4, 0, 1], [z4^3 + z4^2 + z4 + 1, z4^3 + z4^2, z4^3 + z4^2 + z4, 0, 0, 1]])
    self.assertEqual(flip_51(D), A)

    D = C34CurveDivisor(C_3, [[0, 2, 1, 1, 1, 1], [0, 0, 0, 2, 0, 0, 1], [0, 0, 0, 0, 2, 0, 0, 1]])
    A = C34CurveDivisor(C_3, [[0, 2, 1, 1], [0, 0, 1, 0, 1], [0, 0, 2, 0, 0, 1]])
    self.assertEqual(flip_51(D), A)

    D = C34CurveDivisor(C_3_3, [[2*z3^2 + z3 + 1, z3, 2*z3^2 + 2, z3^2 + 2, z3^2 + 2*z3 + 2, 1], [z3^2 + 2*z3 + 2, z3 + 2, 2*z3^2 + 2*z3 + 1, z3^2, 2*z3^2 + 1, 0, 1], [2*z3^2 + 2*z3 + 1, 2*z3^2 + 2*z3 + 1, z3^2 + z3 + 1, 2*z3^2 + 2*z3 + 2, 2*z3^2 + 2*z3 + 2, 0, 0, 1]])
    A = C34CurveDivisor(C_3_3, [[2*z3^2, 2*z3^2 + 1, z3, 1], [2*z3^2 + 2*z3 + 1, 2*z3^2 + 2*z3 + 1, 2*z3^2 + 2, 0, 1], [2*z3 + 2, 2*z3 + 2, 0, 0, 0, 1]])
    self.assertEqual(flip_51(D), A)

    D = C34CurveDivisor(C_31, [[15, 4, 19, 24, 23, 1], [8, 12, 15, 22, 8, 0, 1], [23, 0, 0, 29, 27, 0, 0, 1]])
    A = C34CurveDivisor(C_31, [[28, 12, 26, 1], [28, 18, 6, 0, 1], [1, 15, 1, 0, 0, 1]])
    self.assertEqual(flip_51(D), A)

    D = C34CurveDivisor(C_31_2, [[14*z2, 3*z2 + 17, 29*z2 + 30, 18*z2 + 19, 9*z2 + 28, 1], [19*z2 + 29, 28*z2 + 3, 18*z2 + 9, 15*z2 + 7, 29*z2 + 17, 0, 1], [21*z2 + 23, 16*z2 + 20, 15*z2 + 13, 12*z2 + 27, 20*z2 + 5, 0, 0, 1]])
    A = C34CurveDivisor(C_31_2, [[23*z2 + 19, 10*z2 + 28, 11*z2 + 2, 1], [19*z2 + 30, 24*z2 + 17, 17*z2 + 22, 0, 1], [22*z2 + 27, 17*z2 + 11, 20*z2 + 26, 0, 0, 1]])
    self.assertEqual(flip_51(D), A)

    D = C34CurveDivisor(C_1009, [[759, 73, 399, 895, 471, 1], [82, 869, 481, 850, 817, 0, 1], [609, 376, 438, 69, 138, 0, 0, 1]])
    A = C34CurveDivisor(C_1009, [[733, 879, 522, 1], [150, 707, 879, 0, 1], [554, 361, 57, 0, 0, 1]])
    self.assertEqual(flip_51(D), A)

    # Test case where <f, g> =/= <f, g, h> = <f, h>
    D = C34CurveDivisor(C_2, [[1, 0, 1, 1, 0, 1], [1, 0, 0, 0, 0, 0, 1], [1, 1, 1, 1, 1, 0, 0, 1]])
    A = C34CurveDivisor(C_2, [[0, 1, 0, 1], [0, 0, 0, 0, 1], [1, 1, 1, 0, 0, 1]])
    self.assertEqual(flip_51(D), A)

    D = C34CurveDivisor(C_2_4, [[z4^2, z4^2 + z4, z4^2, z4, 1, 1], [0, z4^3 + z4^2, z4^2 + 1, z4^3 + 1, z4^2, 0, 1], [z4, z4 + 1, z4^2 + z4 + 1, z4, z4^2 + z4 + 1, 0, 0, 1]])
    A = C34CurveDivisor(C_2_4, [[0, z4^3 + z4, 0, 1], [0, z4^2, 0, 0, 1], [z4^2, z4^2 + 1, z4^2, 0, 0, 1]])
    self.assertEqual(flip_51(D), A)

    D = C34CurveDivisor(C_3, [[1, 2, 2, 1, 2, 1], [2, 1, 1, 2, 2, 0, 1], [0, 2, 1, 0, 2, 0, 0, 1]])
    A = C34CurveDivisor(C_3, [[0, 2, 0, 1], [0, 2, 0, 0, 1], [1, 2, 2, 0, 0, 1]])
    self.assertEqual(flip_51(D), A)

    D = C34CurveDivisor(C_3_3, [[2*z3^2 + z3 + 1, z3^2 + 2*z3, 2*z3^2 + 2*z3 + 1, z3^2 + 2, z3^2, 1], [2*z3^2 + 2*z3 + 2, z3^2 + 2*z3 + 1, z3^2 + 2*z3 + 2, 2*z3, 2*z3 + 2, 0, 1], [z3^2 + 2*z3 + 2, 2*z3 + 1, 2*z3 + 1, 2*z3^2 + 2, 2*z3^2 + 2*z3 + 1, 0, 0, 1]])
    A = C34CurveDivisor(C_3_3, [[2*z3^2 + 2*z3, 2*z3^2 + 2*z3, 0, 1], [2*z3^2 + 2*z3 + 1, z3^2 + 1, z3 + 2, 0, 1], [2*z3^2 + 2, 2*z3^2 + 2*z3 + 2, z3 + 2, 0, 0, 1]])
    self.assertEqual(flip_51(D), A)

    D = C34CurveDivisor(C_31, [[11, 4, 26, 25, 13, 1], [19, 11, 29, 24, 29, 0, 1], [28, 12, 10, 3, 18, 0, 0, 1]])
    A = C34CurveDivisor(C_31, [[17, 18, 0, 1], [28, 29, 17, 0, 1], [28, 14, 22, 0, 0, 1]])
    self.assertEqual(flip_51(D), A)

    D = C34CurveDivisor(C_31_2, [[9*z2 + 2, 8*z2, 16*z2 + 24, 2*z2 + 25, 22*z2 + 20, 1], [24*z2 + 11, 3*z2 + 2, 22*z2 + 10, 20*z2 + 27, 6*z2 + 30, 0, 1], [6*z2 + 29, 8*z2 + 10, 14*z2 + 15, 24*z2 + 12, 15*z2 + 3, 0, 0, 1]])
    A = C34CurveDivisor(C_31_2, [[6*z2 + 21, 15*z2 + 25, 0, 1], [12*z2 + 30, 30*z2 + 13, 19, 0, 1], [8*z2 + 15, 14*z2 + 7, z2 + 16, 0, 0, 1]])
    self.assertEqual(flip_51(D), A)

    D = C34CurveDivisor(C_1009, [[410, 466, 162, 217, 757, 1], [654, 485, 157, 811, 72, 0, 1], [914, 992, 708, 9, 538, 0, 0, 1]])
    A = C34CurveDivisor(C_1009, [[630, 851, 0, 1], [882, 357, 596, 0, 1], [199, 609, 13, 0, 0, 1]])
    self.assertEqual(flip_51(D), A)

    # Test case where <f, g> =/= <f, g, h> =/= <f, h>
    D = C34CurveDivisor(C_2, [[0, 0, 1, 1, 1, 1], [0, 0, 1, 0, 0, 0, 1], [0, 1, 1, 1, 0, 0, 0, 1]])
    A = C34CurveDivisor(C_2, [[0, 1, 0, 1], [0, 0, 1, 0, 1], [0, 1, 0, 0, 0, 1]])
    self.assertEqual(flip_51(D), A)

    D = C34CurveDivisor(C_2_4, [[z4^3 + 1, z4^2 + z4, z4^3 + z4, z4^3 + 1, z4^3, 1], [0, z4^2 + z4, z4^3, z4^2 + z4, z4^3 + z4^2 + z4 + 1, 0, 1], [z4^2 + 1, z4^3 + z4, z4, z4^3 + z4^2 + z4 + 1, z4^3 + z4^2 + 1, 0, 0, 1]])
    A = C34CurveDivisor(C_2_4, [[z4^2, z4^3 + z4^2 + z4 + 1, 0, 1], [z4^2 + z4 + 1, z4^2, z4^2 + 1, 0, 1], [z4^2 + z4, z4^3 + z4^2 + z4, z4^2, 0, 0, 1]])
    self.assertEqual(flip_51(D), A)

    D = C34CurveDivisor(C_3, [[0, 1, 1, 1, 2, 1], [1, 0, 0, 2, 2, 0, 1], [2, 2, 2, 2, 2, 0, 0, 1]])
    A = C34CurveDivisor(C_3, [[0, 2, 0, 1], [0, 2, 0, 0, 1], [0, 1, 1, 0, 0, 1]])
    self.assertEqual(flip_51(D), A)

    D = C34CurveDivisor(C_3_3, [[z3^2 + 2*z3, 2*z3 + 2, z3 + 1, 2*z3^2 + 2, z3^2 + 2*z3, 1], [2*z3^2 + 2, 1, 2*z3^2 + 2*z3, z3^2, 2*z3 + 1, 0, 1], [z3^2 + 1, 2*z3, 2*z3^2 + 2*z3 + 2, 0, 2*z3^2 + 2*z3 + 1, 0, 0, 1]])
    A = C34CurveDivisor(C_3_3, [[z3^2 + z3 + 1, z3^2 + z3 + 2, 0, 1], [0, 0, z3^2 + z3 + 1, 0, 1], [z3^2, z3^2, 2*z3^2 + 1, 0, 0, 1]])
    self.assertEqual(flip_51(D), A)

    D = C34CurveDivisor(C_31, [[1, 3, 10, 22, 30, 1], [21, 3, 25, 1, 9, 0, 1], [10, 29, 2, 25, 3, 0, 0, 1]])
    A = C34CurveDivisor(C_31, [[6, 30, 0, 1], [5, 11, 23, 0, 1], [29, 5, 2, 0, 0, 1]])
    self.assertEqual(flip_51(D), A)

    D = C34CurveDivisor(C_31_2, [[19*z2 + 25, 19*z2 + 14, 14*z2 + 15, 14*z2 + 7, 5*z2 + 19, 1], [8*z2 + 22, 21*z2 + 16, 19, 20*z2 + 22, 3*z2 + 20, 0, 1], [9*z2 + 23, 11*z2 + 7, 20*z2 + 8, 12*z2 + 20, 20*z2 + 7, 0, 0, 1]])
    A = C34CurveDivisor(C_31_2, [[2*z2 + 21, 29*z2, 0, 1], [16*z2 + 28, 6*z2 + 1, 30*z2 + 10, 0, 1], [12*z2 + 11, 3*z2 + 1, 24*z2 + 27, 0, 0, 1]])
    self.assertEqual(flip_51(D), A)

    D = C34CurveDivisor(C_1009, [[308, 153, 744, 443, 758, 1], [109, 259, 722, 710, 38, 0, 1], [798, 590, 956, 237, 635, 0, 0, 1]])
    A = C34CurveDivisor(C_1009, [[764, 890, 0, 1], [125, 161, 358, 0, 1], [976, 453, 801, 0, 0, 1]])
    self.assertEqual(flip_51(D), A)



  def test_flip_52(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    D = C34CurveDivisor(C_2, [[0, 1, 0, 1, 1], [0, 1, 1, 1, 0, 1], []])
    A = C34CurveDivisor(C_2, [[0, 1], [1, 0, 1, 0, 0, 1], []])
    self.assertEqual(flip_52(D), A)

    D = C34CurveDivisor(C_2_4, [[z4^3 + z4^2 + z4, 1, z4^3, z4^3 + 1, 1], [z4, 0, 0, z4^3 + z4^2 + 1, 0, 1], []])
    A = C34CurveDivisor(C_2_4, [[z4^3, 1], [z4^3, 0, z4^3 + z4^2 + z4, 0, 0, 1], []])
    self.assertEqual(flip_52(D), A)

    D = C34CurveDivisor(C_3, [[2, 2, 2, 2, 1], [1, 0, 2, 2, 0, 1], []])
    A = C34CurveDivisor(C_3, [[2, 1], [0, 0, 2, 0, 0, 1], []])
    self.assertEqual(flip_52(D), A)

    D = C34CurveDivisor(C_3_3, [[z3^2 + 2*z3 + 2, 2*z3^2, z3^2 + z3 + 2, z3 + 1, 1], [z3, z3^2 + z3 + 1, 2*z3^2 + 1, 2*z3^2 + z3 + 2, 0, 1], []])
    A = C34CurveDivisor(C_3_3, [[z3^2 + z3 + 2, 1], [z3^2 + 2*z3 + 1, 0, z3^2 + 2*z3 + 1, 0, 0, 1], []])
    self.assertEqual(flip_52(D), A)

    D = C34CurveDivisor(C_31, [[20, 8, 29, 22, 1], [0, 3, 21, 12, 0, 1], []])
    A = C34CurveDivisor(C_31, [[29, 1], [3, 0, 10, 0, 0, 1], []])
    self.assertEqual(flip_52(D), A)

    D = C34CurveDivisor(C_31_2, [[9*z2 + 4, 24*z2, 5*z2 + 4, 30*z2 + 16, 1], [23*z2 + 26, 21*z2 + 13, 16*z2 + 12, 30*z2 + 26, 0, 1], []])
    A = C34CurveDivisor(C_31_2, [[5*z2 + 4, 1], [2*z2 + 9, 0, 6*z2 + 24, 0, 0, 1], []])
    self.assertEqual(flip_52(D), A)

    D = C34CurveDivisor(C_1009, [[917, 221, 542, 611, 1], [216, 507, 581, 9, 0, 1], []])
    A = C34CurveDivisor(C_1009, [[542, 1], [465, 0, 1008, 0, 0, 1], []])
    self.assertEqual(flip_52(D), A)



  def test_flip_53(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    D = C34CurveDivisor(C_2, [[0, 0, 1, 0, 1], [0, 0, 0, 0, 0, 1, 1], []])
    A = C34CurveDivisor(C_2, [[1, 1, 1], [1, 0, 0, 1], []])
    self.assertEqual(flip_53(D), A)

    D = C34CurveDivisor(C_2_4, [[z4^3 + z4 + 1, z4^3 + z4^2 + 1, z4^3 + z4^2 + z4 + 1, 0, 1], [z4^3 + z4, z4^3, z4^3 + z4^2 + z4, z4, 0, z4^3 + z4^2 + 1, 1], []])
    A = C34CurveDivisor(C_2_4, [[z4^3 + z4^2 + z4 + 1, z4^3 + z4^2 + 1, 1], [z4^2, z4^2 + z4 + 1, 0, 1], []])
    self.assertEqual(flip_53(D), A)

    D = C34CurveDivisor(C_3, [[1, 1, 1, 1, 1], [0, 2, 2, 0, 0, 2, 1], []])
    A = C34CurveDivisor(C_3, [[2, 0, 1], [2, 2, 0, 1], []])
    self.assertEqual(flip_53(D), A)

    D = C34CurveDivisor(C_3_3, [[2*z3^2, 2*z3^2 + 2, z3 + 2, z3^2 + 1, 1], [2, 2*z3 + 1, 2, 2*z3 + 2, 0, 2*z3^2 + 2*z3 + 1, 1], []])
    A = C34CurveDivisor(C_3_3, [[z3^2 + z3 + 2, 2*z3 + 2, 1], [z3^2 + z3 + 1, 2*z3^2 + 2, 0, 1], []])
    self.assertEqual(flip_53(D), A)

    D = C34CurveDivisor(C_31, [[14, 23, 24, 10, 1], [30, 17, 7, 9, 0, 21, 1], []])
    A = C34CurveDivisor(C_31, [[30, 0, 1], [10, 21, 0, 1], []])
    self.assertEqual(flip_53(D), A)

    D = C34CurveDivisor(C_31_2, [[27, 23*z2 + 10, 21*z2 + 18, 13*z2 + 4, 1], [3*z2 + 2, 16*z2 + 1, 21*z2 + 30, 4*z2 + 4, 0, 3*z2 + 23, 1], []])
    A = C34CurveDivisor(C_31_2, [[6*z2 + 29, 16*z2 + 27, 1], [30*z2 + 2, 28, 0, 1], []])
    self.assertEqual(flip_53(D), A)

    D = C34CurveDivisor(C_1009, [[700, 889, 918, 822, 1], [710, 992, 431, 300, 0, 51, 1], []])
    A = C34CurveDivisor(C_1009, [[711, 873, 1], [518, 714, 0, 1], []])
    self.assertEqual(flip_53(D), A)



  def test_flip_54(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    D = C34CurveDivisor(C_2, [[1, 0, 0, 1], [0, 0, 1, 0, 1, 1, 0, 0, 1], []])
    A = C34CurveDivisor(C_2, [[1, 1], [1, 0, 1], []])
    self.assertEqual(flip_54(D), A)

    D = C34CurveDivisor(C_2_4, [[z4^3 + z4^2 + z4 + 1, z4^3 + z4 + 1, 1, 1], [z4^2 + z4 + 1, z4^2, 1, 0, z4^3 + z4 + 1, z4, 0, 0, 1], []])
    A = C34CurveDivisor(C_2_4, [[0, 1], [z4^3 + z4^2 + z4 + 1, 0, 1], []])
    self.assertEqual(flip_54(D), A)

    D = C34CurveDivisor(C_3, [[1, 2, 0, 1], [2, 2, 1, 0, 1, 1, 0, 0, 1], []])
    A = C34CurveDivisor(C_3, [[1, 1], [2, 0, 1], []])
    self.assertEqual(flip_54(D), A)

    D = C34CurveDivisor(C_3_3, [[2*z3^2 + 1, 0, 2*z3^2 + z3, 1], [2*z3 + 2, z3 + 1, 1, 0, 2*z3^2 + 2*z3, 2*z3^2 + 1, 0, 0, 1], []])
    A = C34CurveDivisor(C_3_3, [[z3 + 2, 1], [z3^2 + 2, 0, 1], []])
    self.assertEqual(flip_54(D), A)

    D = C34CurveDivisor(C_31, [[10, 14, 8, 1], [26, 12, 12, 0, 14, 22, 0, 0, 1], []])
    A = C34CurveDivisor(C_31, [[11, 1], [1, 0, 1], []])
    self.assertEqual(flip_54(D), A)

    D = C34CurveDivisor(C_31_2, [[23*z2 + 15, 13*z2 + 26, 29*z2 + 20, 1], [5*z2 + 15, 30*z2 + 22, 13*z2 + 22, 0, 21*z2 + 2, 19*z2 + 11, 0, 0, 1], []])
    A = C34CurveDivisor(C_31_2, [[23*z2 + 16, 1], [29*z2 + 10, 0, 1], []])
    self.assertEqual(flip_54(D), A)

    D = C34CurveDivisor(C_1009, [[32, 174, 906, 1], [873, 62, 180, 0, 19, 194, 0, 0, 1], []])
    A = C34CurveDivisor(C_1009, [[34, 1], [810, 0, 1], []])
    self.assertEqual(flip_54(D), A)



  def test_flip_61(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    # Test case where <f, g, h> = <f, g>
    D = C34CurveDivisor(C_2, [[1, 1, 1, 1, 0, 0, 1], [0, 1, 0, 1, 0, 1, 0, 1], [1, 0, 0, 1, 1, 0, 0, 0, 1]])
    A = C34CurveDivisor(C_2, [[0, 0, 1, 1], [1, 1, 0, 0, 1], [0, 1, 1, 0, 0, 1]])
    self.assertEqual(flip_61(D), A)

    D = C34CurveDivisor(C_2_4, [[z4 + 1, z4^3 + z4^2, 1, z4^2 + z4, z4^2 + 1, z4^3 + 1, 1], [z4^2 + z4, z4^3 + z4^2 + z4, 1, z4^3 + z4^2, z4^3 + z4^2 + z4 + 1, z4 + 1, 0, 1], [z4^2 + 1, z4^2 + 1, z4^3 + 1, z4 + 1, z4^3 + 1, z4^2 + z4, 0, 0, 1]])
    A = C34CurveDivisor(C_2_4, [[z4^2, z4^3 + z4^2, z4^2 + z4, 1], [z4^3 + z4^2, z4^2 + 1, z4^2 + 1, 0, 1], [z4^2 + z4 + 1, z4^3 + z4^2, z4^3 + z4^2 + z4, 0, 0, 1]])
    self.assertEqual(flip_61(D), A)

    D = C34CurveDivisor(C_3, [[1, 2, 1, 2, 2, 1, 1], [1, 2, 0, 0, 0, 0, 0, 1], [2, 2, 0, 2, 0, 2, 0, 0, 1]])
    A = C34CurveDivisor(C_3, [[1, 2, 1, 1], [2, 0, 0, 0, 1], [2, 1, 1, 0, 0, 1]])
    self.assertEqual(flip_61(D), A)

    D = C34CurveDivisor(C_3_3, [[z3^2 + 2*z3, 2*z3^2 + z3, z3^2 + z3, z3^2 + 2*z3 + 2, 2*z3, 2*z3^2 + 1, 1], [2, 2*z3^2 + 2*z3 + 2, 2*z3^2 + z3 + 1, 2*z3^2 + 2*z3 + 2, 2*z3 + 1, z3^2 + z3 + 1, 0, 1], [z3^2, 2*z3^2 + 2, z3^2 + 1, 2*z3^2 + z3 + 1, 2*z3^2 + 2*z3 + 2, 1, 0, 0, 1]])
    A = C34CurveDivisor(C_3_3, [[z3^2 + z3 + 2, 2*z3^2 + 2*z3 + 2, 1, 1], [z3^2 + 2*z3 + 1, z3^2 + 2, z3^2 + z3 + 1, 0, 1], [2, 0, z3^2 + z3 + 2, 0, 0, 1]])
    self.assertEqual(flip_61(D), A)

    D = C34CurveDivisor(C_31, [[16, 26, 0, 23, 12, 8, 1], [19, 6, 25, 24, 17, 5, 0, 1], [7, 28, 23, 16, 4, 24, 0, 0, 1]])
    A = C34CurveDivisor(C_31, [[25, 8, 28, 1], [29, 8, 25, 0, 1], [17, 5, 13, 0, 0, 1]])
    self.assertEqual(flip_61(D), A)

    D = C34CurveDivisor(C_31_2, [[26*z2 + 21, 25*z2 + 28, 29*z2 + 7, 5*z2 + 15, 7*z2 + 13, 18*z2 + 30, 1], [29*z2 + 23, 2*z2 + 27, 19*z2 + 29, 11*z2 + 11, 10*z2 + 27, 21*z2 + 3, 0, 1], [23, 24*z2 + 22, 20*z2 + 9, 22*z2 + 6, 9*z2 + 26, 14*z2 + 17, 0, 0, 1]])
    A = C34CurveDivisor(C_31_2, [[29*z2 + 18, 29*z2 + 29, 25*z2 + 30, 1], [21*z2 + 19, 11, 9*z2 + 3, 0, 1], [18*z2 + 11, 29*z2 + 12, 20*z2 + 15, 0, 0, 1]])
    self.assertEqual(flip_61(D), A)

    D = C34CurveDivisor(C_1009, [[925, 136, 99, 615, 118, 925, 1], [772, 252, 607, 139, 83, 998, 0, 1], [845, 404, 934, 468, 314, 906, 0, 0, 1]])
    A = C34CurveDivisor(C_1009, [[147, 67, 264, 1], [85, 970, 412, 0, 1], [548, 493, 993, 0, 0, 1]])
    self.assertEqual(flip_61(D), A)

    # Test case where <f, g> =/= <f, g, h> = <f, h>
    D = C34CurveDivisor(C_2, [[1, 1, 1, 0, 1, 1, 1], [1, 1, 1, 0, 1, 1, 0, 1], [0, 1, 0, 0, 1, 0, 0, 0, 1]])
    A = C34CurveDivisor(C_2, [[1, 0, 0, 1], [1, 1, 1, 0, 1], [0, 1, 0, 0, 0, 1]])
    self.assertEqual(flip_61(D), A)

    D = C34CurveDivisor(C_2_4, [[z4^3 + z4^2, z4^3 + z4 + 1, z4^3 + z4^2, z4, z4^3 + z4, z4^2 + 1, 1], [z4^2 + z4 + 1, 0, z4^3 + z4^2 + z4 + 1, z4^3 + 1, z4^2 + z4, z4 + 1, 0, 1], [z4 + 1, z4^3 + z4^2 + 1, z4^3, 1, z4^2 + z4, 1, 0, 0, 1]])
    A = C34CurveDivisor(C_2_4, [[0, z4 + 1, 0, 1], [z4, z4^3 + z4^2 + z4 + 1, z4 + 1, 0, 1], [z4^3 + 1, z4^3 + z4, z4^3 + z4 + 1, 0, 0, 1]])
    self.assertEqual(flip_61(D), A)

    D = C34CurveDivisor(C_3, [[2, 0, 2, 1, 2, 1, 1], [2, 2, 1, 1, 1, 1, 0, 1], [0, 0, 1, 2, 0, 2, 0, 0, 1]])
    A = C34CurveDivisor(C_3, [[2, 0, 0, 1], [1, 2, 2, 0, 1], [1, 0, 1, 0, 0, 1]])
    self.assertEqual(flip_61(D), A)

    D = C34CurveDivisor(C_3_3, [[z3^2 + z3 + 1, 2*z3 + 2, z3^2 + 1, 2*z3^2 + 2*z3, z3^2 + 2, z3^2 + 2*z3 + 2, 1], [2*z3 + 1, z3^2 + z3 + 2, z3^2 + z3 + 1, 2*z3, z3, 2*z3^2 + z3 + 1, 0, 1], [z3 + 2, z3^2 + 2*z3, 2, z3^2 + 2*z3, 2*z3, 2*z3 + 2, 0, 0, 1]])
    A = C34CurveDivisor(C_3_3, [[z3^2 + 1, 2*z3^2, 0, 1], [z3^2 + z3, 2*z3^2 + z3 + 1, 2*z3 + 1, 0, 1], [2*z3^2 + 1, 2*z3, 2*z3^2 + 2*z3 + 2, 0, 0, 1]])
    self.assertEqual(flip_61(D), A)

    D = C34CurveDivisor(C_31, [[6, 17, 30, 10, 15, 9, 1], [20, 29, 0, 12, 25, 5, 0, 1], [18, 3, 24, 18, 13, 2, 0, 0, 1]])
    A = C34CurveDivisor(C_31, [[11, 17, 0, 1], [6, 2, 3, 0, 1], [13, 14, 19, 0, 0, 1]])
    self.assertEqual(flip_61(D), A)

    D = C34CurveDivisor(C_31_2, [[7*z2 + 8, 8*z2 + 30, 7*z2 + 18, 10*z2, 21*z2 + 11, 29*z2 + 12, 1], [12*z2 + 16, 8*z2 + 26, 21*z2 + 27, 12*z2 + 24, 14*z2 + 9, 21*z2 + 7, 0, 1], [13*z2 + 7, 30*z2 + 10, 15*z2 + 24, 11*z2 + 1, 5*z2 + 26, 10*z2 + 14, 0, 0, 1]])
    A = C34CurveDivisor(C_31_2, [[26*z2 + 24, 21*z2 + 16, 0, 1], [25*z2 + 2, 16*z2 + 10, 28*z2 + 23, 0, 1], [11*z2 + 5, 22*z2 + 15, 4*z2 + 4, 0, 0, 1]])
    self.assertEqual(flip_61(D), A)

    D = C34CurveDivisor(C_1009, [[200, 650, 1008, 669, 672, 162, 1], [32, 473, 781, 41, 737, 264, 0, 1], [655, 215, 415, 353, 343, 902, 0, 0, 1]])
    A = C34CurveDivisor(C_1009, [[348, 780, 0, 1], [134, 20, 713, 0, 1], [643, 553, 823, 0, 0, 1]])
    self.assertEqual(flip_61(D), A)

    # Test case where <f, g> =/= <f, g, h> =/= <f, h>
    D = C34CurveDivisor(C_2, [[0, 1, 1, 0, 1, 0, 1], [0, 0, 0, 0, 1, 1, 0, 1], [1, 1, 0, 0, 1, 0, 0, 0, 1]])
    A = C34CurveDivisor(C_2, [[0, 1, 0, 1], [0, 0, 1, 0, 1], [0, 1, 0, 0, 0, 1]])
    self.assertEqual(flip_61(D), A)

    D = C34CurveDivisor(C_2_4, [[z4^3 + z4^2 + 1, z4^3 + z4, z4^2, z4, z4^3 + z4^2 + z4 + 1, z4^2 + 1, 1], [z4^3, z4^2 + z4, z4^3, z4^3 + 1, z4, z4^2 + z4, 0, 1], [z4^3 + z4^2 + z4, z4^3 + z4, z4^3 + z4^2, z4^3 + z4 + 1, z4^3 + z4 + 1, z4^3 + z4^2 + z4, 0, 0, 1]])
    A = C34CurveDivisor(C_2_4, [[z4^3 + z4, z4^2 + z4 + 1, 0, 1], [z4^3 + z4^2 + z4 + 1, z4 + 1, z4^2 + 1, 0, 1], [z4^3 + z4^2 + z4, z4, z4^2 + 1, 0, 0, 1]])
    self.assertEqual(flip_61(D), A)

    D = C34CurveDivisor(C_3, [[0, 2, 2, 2, 0, 2, 1], [0, 2, 0, 1, 0, 0, 0, 1], [0, 0, 2, 1, 0, 2, 0, 0, 1]])
    A = C34CurveDivisor(C_3, [[2, 0, 0, 1], [1, 2, 2, 0, 1], [1, 0, 1, 0, 0, 1]])
    self.assertEqual(flip_61(D), A)

    D = C34CurveDivisor(C_3_3, [[2*z3^2 + 2, 2*z3^2 + 1, 2*z3^2, z3^2 + 2*z3, z3^2 + 2, 2*z3 + 2, 1], [2*z3^2 + z3, 2*z3^2 + 2*z3, z3 + 2, z3^2 + z3 + 2, z3^2 + z3, 2*z3^2 + z3 + 2, 0, 1], [z3^2 + 1, z3^2 + z3 + 1, z3^2 + 1, 1, 2*z3^2 + z3 + 1, z3^2 + z3 + 2, 0, 0, 1]])
    A = C34CurveDivisor(C_3_3, [[z3 + 2, 1, 0, 1], [2, z3^2 + 2*z3 + 2, z3^2 + z3 + 1, 0, 1], [z3^2, 2*z3^2 + 2*z3 + 1, 2*z3^2, 0, 0, 1]])
    self.assertEqual(flip_61(D), A)

    D = C34CurveDivisor(C_31, [[9, 14, 27, 24, 26, 0, 1], [16, 29, 5, 9, 17, 26, 0, 1], [18, 3, 0, 21, 16, 6, 0, 0, 1]])
    A = C34CurveDivisor(C_31, [[21, 10, 0, 1], [26, 17, 7, 0, 1], [21, 5, 21, 0, 0, 1]])
    self.assertEqual(flip_61(D), A)

    D = C34CurveDivisor(C_31_2, [[23*z2 + 3, 22*z2 + 27, 21*z2 + 12, 18*z2 + 30, 3*z2 + 1, 14*z2 + 27, 1], [12, 6*z2 + 30, 9*z2 + 9, 21*z2 + 2, 13*z2 + 1, 17*z2 + 15, 0, 1], [18*z2 + 26, 3*z2 + 21, 12*z2 + 15, 18*z2, 16*z2 + 2, 2*z2 + 14, 0, 0, 1]])
    A = C34CurveDivisor(C_31_2, [[8*z2 + 28, 15*z2 + 7, 0, 1], [29*z2 + 30, 4*z2 + 26, 7*z2 + 2, 0, 1], [29*z2 + 28, 28*z2 + 28, 20*z2 + 12, 0, 0, 1]])
    self.assertEqual(flip_61(D), A)

    D = C34CurveDivisor(C_1009, [[400, 243, 940, 734, 368, 690, 1], [205, 314, 698, 13, 475, 495, 0, 1], [100, 639, 724, 960, 184, 964, 0, 0, 1]])
    A = C34CurveDivisor(C_1009, [[92, 496, 0, 1], [881, 298, 521, 0, 1], [546, 590, 92, 0, 0, 1]])
    self.assertEqual(flip_61(D), A)



  def test_flip_62(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    D = C34CurveDivisor(C_2, [[1, 1, 1, 1, 1, 1], [0, 0, 0, 0, 0, 0, 1], []])
    A = C34CurveDivisor(C_2, [[1, 1], [1, 0, 0, 0, 0, 1], []])
    self.assertEqual(flip_62(D), A)

    D = C34CurveDivisor(C_2_4, [[z4^3 + z4, 1, z4 + 1, 1, 1, 1], [z4^2 + z4 + 1, z4^3 + z4^2, z4^3 + z4^2 + z4, z4^2 + 1, z4^2 + z4 + 1, 0, 1], []])
    A = C34CurveDivisor(C_2_4, [[z4^2 + z4 + 1, 1], [z4^3 + z4 + 1, 0, z4^2, 0, 0, 1], []])
    self.assertEqual(flip_62(D), A)

    D = C34CurveDivisor(C_3, [[1, 1, 2, 2, 1, 1], [1, 0, 0, 0, 0, 0, 1], []])
    A = C34CurveDivisor(C_3, [[0, 1], [1, 0, 2, 0, 0, 1], []])
    self.assertEqual(flip_62(D), A)

    D = C34CurveDivisor(C_3_3, [[z3^2, z3^2 + z3 + 1, 2*z3^2 + 2, z3 + 2, 1, 1], [2*z3, z3^2 + z3 + 1, z3 + 1, z3, 2*z3^2 + 2*z3 + 1, 0, 1], []])
    A = C34CurveDivisor(C_3_3, [[2*z3^2 + 2*z3 + 1, 1], [2*z3^2 + z3 + 1, 0, z3 + 1, 0, 0, 1], []])
    self.assertEqual(flip_62(D), A)

    D = C34CurveDivisor(C_31, [[7, 12, 25, 13, 1, 1], [24, 18, 10, 13, 15, 0, 1], []])
    A = C34CurveDivisor(C_31, [[13, 1], [2, 0, 12, 0, 0, 1], []])
    self.assertEqual(flip_62(D), A)

    D = C34CurveDivisor(C_31_2, [[23*z2 + 30, 2*z2 + 17, 25*z2, 13*z2 + 26, 1, 1], [30*z2 + 18, 4*z2, 14*z2 + 25, z2 + 20, 23*z2 + 11, 0, 1], []])
    A = C34CurveDivisor(C_31_2, [[30*z2 + 10, 1], [8*z2 + 9, 0, 26*z2 + 21, 0, 0, 1], []])
    self.assertEqual(flip_62(D), A)

    D = C34CurveDivisor(C_1009, [[392, 351, 279, 867, 1, 1], [541, 164, 403, 467, 764, 0, 1], []])
    A = C34CurveDivisor(C_1009, [[714, 1], [692, 0, 574, 0, 0, 1], []])
    self.assertEqual(flip_62(D), A)



  def test_flip_63(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    D = C34CurveDivisor(C_2, [[0, 1, 1, 0, 1, 1], [1, 1, 1, 1, 1, 0, 0, 1], []])
    A = C34CurveDivisor(C_2, [[0, 1, 1], [0, 1, 0, 1], []])
    self.assertEqual(flip_63(D), A)

    D = C34CurveDivisor(C_2_4, [[z4^3 + z4^2 + z4 + 1, z4^3 + z4, z4^3, z4, z4^3 + z4^2 + 1, 1], [z4^3 + 1, 0, z4^2 + z4 + 1, 0, z4^2, 0, z4^2 + 1, 1], []])
    A = C34CurveDivisor(C_2_4, [[z4^3 + z4^2, z4^3, 1], [z4^3, z4^2, 0, 1], []])
    self.assertEqual(flip_63(D), A)

    D = C34CurveDivisor(C_3, [[0, 1, 0, 2, 2, 1], [1, 1, 0, 1, 0, 0, 2, 1], []])
    A = C34CurveDivisor(C_3, [[0, 0, 1], [0, 2, 0, 1], []])
    self.assertEqual(flip_63(D), A)

    D = C34CurveDivisor(C_3_3, [[2, 2, 2*z3^2 + z3, z3^2 + 1, 2*z3 + 1, 1], [z3^2 + 2*z3 + 2, 2*z3 + 1, z3^2 + 1, z3, z3^2 + z3 + 2, 0, z3, 1], []])
    A = C34CurveDivisor(C_3_3, [[z3^2 + z3 + 1, z3 + 1, 1], [z3^2 + 2*z3, z3^2 + 1, 0, 1], []])
    self.assertEqual(flip_63(D), A)

    D = C34CurveDivisor(C_31, [[14, 12, 27, 14, 10, 1], [18, 0, 21, 1, 22, 0, 14, 1], []])
    A = C34CurveDivisor(C_31, [[3, 27, 1], [16, 16, 0, 1], []])
    self.assertEqual(flip_63(D), A)

    D = C34CurveDivisor(C_31_2, [[z2 + 26, z2 + 13, 12*z2 + 3, 16*z2 + 26, z2 + 16, 1], [26*z2 + 15, 16*z2 + 2, z2 + 15, 17*z2 + 24, 30*z2 + 7, 0, 5*z2 + 3, 1], []])
    A = C34CurveDivisor(C_31_2, [[9*z2 + 1, 27*z2 + 13, 1], [11*z2 + 4, 13*z2 + 10, 0, 1], []])
    self.assertEqual(flip_63(D), A)

    D = C34CurveDivisor(C_1009, [[346, 980, 871, 950, 515, 1], [150, 685, 106, 273, 241, 0, 340, 1], []])
    A = C34CurveDivisor(C_1009, [[961, 175, 1], [214, 189, 0, 1], []])
    self.assertEqual(flip_63(D), A)



  def test_flip_64(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    D = C34CurveDivisor(C_2, [[0, 1, 1, 1, 1], [1, 1, 0, 0, 0, 0, 1, 0, 0, 1], []])
    A = C34CurveDivisor(C_2, [[0, 1], [0, 0, 1], []])
    self.assertEqual(flip_64(D), A)

    D = C34CurveDivisor(C_2_4, [[z4^3, z4, z4^3 + 1, z4^3 + z4, 1], [z4^3 + z4^2 + z4 + 1, 1, z4^3 + z4, z4^3 + z4, 0, z4^3 + z4^2, 0, 0, 0, 1], []])
    A = C34CurveDivisor(C_2_4, [[z4^3 + z4^2 + z4, 1], [z4^3 + z4^2 + 1, 0, 1], []])
    self.assertEqual(flip_64(D), A)

    D = C34CurveDivisor(C_3, [[0, 0, 1, 0, 1], [2, 0, 0, 0, 0, 0, 0, 0, 0, 1], []])
    A = C34CurveDivisor(C_3, [[0, 1], [0, 0, 1], []])
    self.assertEqual(flip_64(D), A)

    D = C34CurveDivisor(C_3_3, [[z3^2 + 2*z3, 2*z3 + 1, 1, 2*z3, 1], [z3^2, 2*z3^2 + z3, 1, 2*z3^2 + 2*z3, 0, 2*z3, z3 + 2, 0, 0, 1], []])
    A = C34CurveDivisor(C_3_3, [[z3^2 + z3 + 2, 1], [z3^2 + 2*z3, 0, 1], []])
    self.assertEqual(flip_64(D), A)

    D = C34CurveDivisor(C_31, [[13, 5, 0, 22, 1], [18, 0, 20, 5, 0, 17, 27, 0, 0, 1], []])
    A = C34CurveDivisor(C_31, [[12, 1], [6, 0, 1], []])
    self.assertEqual(flip_64(D), A)

    D = C34CurveDivisor(C_31_2, [[29*z2 + 6, 6*z2 + 27, 25*z2 + 21, 14*z2 + 2, 1], [z2 + 27, 8*z2 + 1, 15*z2 + 20, 21*z2 + 23, 0, 21*z2 + 19, 4*z2 + 20, 0, 0, 1], []])
    A = C34CurveDivisor(C_31_2, [[19*z2 + 8, 1], [6*z2 + 7, 0, 1], []])
    self.assertEqual(flip_64(D), A)

    D = C34CurveDivisor(C_1009, [[473, 951, 963, 517, 1], [670, 714, 1003, 501, 0, 672, 801, 0, 0, 1], []])
    A = C34CurveDivisor(C_1009, [[517, 1], [276, 0, 1], []])
    self.assertEqual(flip_64(D), A)



  def test_flip_65(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    D = C34CurveDivisor(C_2, [[0, 1, 0, 1], [], []])
    A = C34CurveDivisor(C_2, [[1], [], []])
    self.assertEqual(flip_65(D), A)

    D = C34CurveDivisor(C_2_4, [[z4^3 + z4 + 1, z4^3 + z4^2 + 1, z4^3 + 1, 1], [], []])
    A = C34CurveDivisor(C_2_4, [[1], [], []])
    self.assertEqual(flip_65(D), A)

    D = C34CurveDivisor(C_3, [[1, 0, 1, 1], [], []])
    A = C34CurveDivisor(C_3, [[1], [], []])
    self.assertEqual(flip_65(D), A)

    D = C34CurveDivisor(C_3_3, [[2*z3 + 1, 2*z3^2 + 2*z3 + 2, z3 + 1, 1], [], []])
    A = C34CurveDivisor(C_3_3, [[1], [], []])
    self.assertEqual(flip_65(D), A)

    D = C34CurveDivisor(C_31, [[24, 23, 20, 1], [], []])
    A = C34CurveDivisor(C_31, [[1], [], []])
    self.assertEqual(flip_65(D), A)

    D = C34CurveDivisor(C_31_2, [[30*z2 + 29, 21*z2 + 9, 10*z2 + 25, 1], [], []])
    A = C34CurveDivisor(C_31_2, [[1], [], []])
    self.assertEqual(flip_65(D), A)

    D = C34CurveDivisor(C_1009, [[470, 927, 973, 1], [], []])
    A = C34CurveDivisor(C_1009, [[1], [], []])
    self.assertEqual(flip_65(D), A)

