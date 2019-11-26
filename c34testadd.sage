load("c34crv.sage")

class TestAdd(unittest.TestCase) :
  def setUp(self) :
    self.z2 = GF(31^2).gen()
    self.z3 = GF(3^3).gen()
    self.z4 = GF(2^4).gen()
    z2, z3, z4 = self.z2, self.z3, self.z4
    
    self.C_2 = C34Curve(GF(2), [0, 0, 1, 0, 1, 1, 1, 1, 1])
    self.C_3 = C34Curve(GF(3), [0, 2, 1, 1, 0, 2, 2, 0, 2])
    self.C_11 = C34Curve(GF(11), [6, 5, 1, 8, 6, 3, 7, 10, 4]) # The vertical tangent line at P = (9 : 6 : 1) intersects P with multiplicity 2
    self.C_31 = C34Curve(GF(31), [2, 29, 13, 7, 28, 25, 17, 13, 17])
    #self.C_41 = C34Curve(GF(41), [26, 22, 16, 37, 30, 22, 7, 22, 29]) # Vertical tangents at (31 : 33 : 1) and (28 : 22 : 1) intersect with m = 2
    self.C_2_4 = C34Curve(GF(2^4), [z4^3 + 1, z4, z4^3 + z4^2 + 1, z4^2, z4^3 + z4, z4^3 + z4 + 1, z4^3 + z4, z4^3 + z4^2 + z4, z4^3 + 1])
    self.C_3_3 = C34Curve(GF(3^3), [z3 + 2, z3 + 1, z3^2 + z3 + 2, 2*z3^2 + z3 + 1, 2*z3^2 + z3, z3^2 + 2, 2*z3^2 + z3, 2*z3^2 + 1, 2])
    self.C_31_2 = C34Curve(GF(31^2), [9*z2 + 11, 11*z2 + 28, 17*z2 + 2, 24*z2 + 22, 6*z2 + 29, 21*z2 + 8, 18*z2 + 5, 18*z2 + 15, 13*z2 + 10])
    self.C_1009 = C34Curve(GF(1009), [715, 286, 748, 7, 92, 446, 122, 935, 314])
  
  def test_add_11_11(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4
    
    # Test most common case, type(lcm(D1, D2)) = 21, over several fields
    D1 = C34CurveDivisor(C_2, [[0, 1], [0, 0, 1], []])
    D2 = C34CurveDivisor(C_2, [[0, 1], [0, 0, 1], []])
    D3 = C34CurveDivisor(C_2, [[0, 0, 1], [0, 0, 0, 1], []])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_2_4, [[z4^3 + z4 + 1, 1], [z4^2 + z4, 0, 1], []])
    D2 = C34CurveDivisor(C_2_4, [[z4^3, 1], [1, 0, 1], []])
    D3 = C34CurveDivisor(C_2_4, [[z4^3 + z4 + 1, z4^3 + z4^2, 1],
                                 [z4^2 + z4 + 1, z4 + 1, 0, 1], []])
    self.assertEqual(D1 + D2, D3)
    
    D1 = C34CurveDivisor(C_3, [[1, 1], [2, 0, 1], []])
    D2 = C34CurveDivisor(C_3, [[0, 1], [0, 0, 1], []])
    D3 = C34CurveDivisor(C_3, [[0, 1, 1], [0, 1, 0, 1], []])
    self.assertEqual(D1 + D2, D3)
    
    D1 = C34CurveDivisor(C_3_3, [[2*z3^2 + 2*z3 + 1, 1], [z3 + 2, 0, 1], []])
    D2 = C34CurveDivisor(C_3_3, [[z3^2 + z3 + 1, 1], [z3 + 2, 0, 1], []])
    D3 = C34CurveDivisor(C_3_3, [[z3 + 2, 0, 1], [z3^2 + 2*z3, 2, 0, 1], []])
    self.assertEqual(D1 + D2, D3)
    
    D1 = C34CurveDivisor(C_31, [[22, 1], [13, 0, 1], []])
    D2 = C34CurveDivisor(C_31, [[28, 1], [17, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[19, 20, 1], [27, 19, 0, 1], []])
    self.assertEqual(D1 + D2, D3)
    
    D1 = C34CurveDivisor(C_31_2, [[17*z2 + 25, 1], [13*z2 + 9, 0, 1], []])
    D2 = C34CurveDivisor(C_31_2, [[16*z2 + 4, 1], [14*z2 + 18, 0, 1], []])
    D3 = C34CurveDivisor(C_31_2, [[13*z2 + 4, 5*z2 + 10, 1], [20*z2 + 28, 2*z2 + 29, 0, 1], []])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_1009, [[202, 1], [648, 0, 1], []])
    D2 = C34CurveDivisor(C_1009, [[861, 1], [413, 0, 1], []])
    D3 = C34CurveDivisor(C_1009, [[743, 475, 1], [374, 54, 0, 1], []])
    self.assertEqual(D1 + D2, D3)
    self.assertEqual(add_11_11(D1, D2), D3)

    # Test rarer case, type(lcm(D1, D2)) = 22
    D1 = C34CurveDivisor(C_2, [[1, 1], [0, 0, 1], []])
    D2 = C34CurveDivisor(C_2, [[1, 1], [1, 0, 1], []])
    D3 = C34CurveDivisor(C_2, [[1, 1], [0, 0, 1, 0, 0, 1], []])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_2_4, [[z4^2 + 1, 1], [z4^2 + 1, 0, 1], []])
    D2 = C34CurveDivisor(C_2_4, [[z4^2 + 1, 1], [1, 0, 1], []])
    D3 = C34CurveDivisor(C_2_4, [[z4^2 + 1, 1], [z4^2 + 1, 0, z4^2, 0, 0, 1], []])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_3, [[2, 1], [2, 0, 1], []])
    D2 = C34CurveDivisor(C_3, [[2, 1], [0, 0, 1], []])
    D3 = C34CurveDivisor(C_3, [[2, 1], [0, 0, 2, 0, 0, 1], []])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_3_3, [[z3^2 + 2*z3 + 2, 1], [2, 0, 1], []])
    D2 = C34CurveDivisor(C_3_3, [[z3^2 + 2*z3 + 2, 1], [2*z3 + 2, 0, 1], []])
    D3 = C34CurveDivisor(C_3_3, [[z3^2 + 2*z3 + 2, 1], [z3 + 1, 0, 2*z3 + 1, 0, 0, 1], []])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_31, [[11, 1], [26, 0, 1], []])
    D2 = C34CurveDivisor(C_31, [[11, 1], [28, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[11, 1], [15, 0, 23, 0, 0, 1], []])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_31_2, [[27*z2 + 9, 1], [26*z2 + 7, 0, 1], []])
    D2 = C34CurveDivisor(C_31_2, [[27*z2 + 9, 1], [30*z2 + 6, 0, 1], []])
    D3 = C34CurveDivisor(C_31_2, [[27*z2 + 9, 1], [4*z2 + 27, 0, 25*z2 + 13, 0, 0, 1], []])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_1009, [[70, 1], [605, 0, 1], []])
    D2 = C34CurveDivisor(C_1009, [[70, 1], [605, 0, 1], []])
    D3 = C34CurveDivisor(C_1009, [[70, 1], [767, 0, 201, 0, 0, 1], []])
    self.assertEqual(D1 + D2, D3)



  def test_add_21_11(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4
    
    # Test most common case, type(lcm(D1, D2)) = 31, over several fields
    D1 = C34CurveDivisor(C_2, [[1, 1, 1], [1, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_2, [[0, 1], [0, 0, 1], []])
    D3 = C34CurveDivisor(C_2, [[0, 1, 1, 1], [0, 0, 1, 0, 1], [0, 0, 0, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)
    
    D1 = C34CurveDivisor(C_2_4, [[z4^3 + z4 + 1, z4^3 + z4^2, 1],
                                 [z4^2 + z4 + 1, z4 + 1, 0, 1], []])
    D2 = C34CurveDivisor(C_2_4, [[z4^3 + z4 + 1, 1], [z4^2 + z4, 0, 1], []])
    D3 = C34CurveDivisor(C_2_4, [[z4 + 1, 1, z4^2 + z4 + 1, 1],
                                 [z4^3 + z4^2 + z4, z4^3 + z4, z4^3 + 1, 0, 1],
                                 [z4^3 + z4 + 1, z4^3 + z4^2 + z4 + 1, z4^3 + z4 + 1, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)
    
    D1 = C34CurveDivisor(C_3, [[2, 0, 1], [2, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_3, [[2, 1], [2, 0, 1], []])
    D3 = C34CurveDivisor(C_3, [[2, 0, 0, 1], [1, 2, 2, 0, 1], [1, 0, 1, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)
    
    D1 = C34CurveDivisor(C_3_3, [[2*z3 + 1, 2*z3^2, 1], [2*z3^2 + 1, 2*z3^2 + z3, 0, 1], []])
    D2 = C34CurveDivisor(C_3_3, [[z3^2 + 2*z3 + 1, 1], [2*z3^2 + 2, 0, 1], []])
    D3 = C34CurveDivisor(C_3_3, [[z3^2 + 2*z3 + 1, 2*z3^2 + 2, z3^2 + 2*z3 + 2, 1],
                                 [z3^2 + z3, 2*z3^2 + 2*z3, z3^2 + 2, 0, 1],
                                 [2*z3, z3^2 + 2*z3 + 1, 2*z3^2 + z3, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)
    
    D1 = C34CurveDivisor(C_31, [[1, 10, 1], [14, 24, 0, 1], []])
    D2 = C34CurveDivisor(C_31, [[30, 1], [25, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[19, 12, 5, 1], [26, 26, 11, 0, 1], [13, 21, 9, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)
    
    D1 = C34CurveDivisor(C_31_2, [[22*z2 + 11, 18*z2 + 9, 1], [29*z2 + 2, 14*z2 + 29, 0, 1], []])
    D2 = C34CurveDivisor(C_31_2, [[28*z2 + 30, 1], [5*z2 + 24, 0, 1], []])
    D3 = C34CurveDivisor(C_31_2, [[16*z2 + 28, 9*z2 + 8, 28*z2 + 21, 1], [15*z2 + 24, 2*z2 + 20, 2*z2 + 20, 0, 1], [6*z2 + 1, 21*z2 + 29, 11*z2 + 25, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)
    
    D1 = C34CurveDivisor(C_1009, [[118, 576, 1], [929, 39, 0, 1], []])
    D2 = C34CurveDivisor(C_1009, [[958, 1], [670, 0, 1], []])
    D3 = C34CurveDivisor(C_1009, [[930, 232, 419, 1],
                                  [135, 568, 765, 0, 1],
                                  [291, 230, 72, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    # TODO: Test case where type(lcm(D1, D2)) = 32
    # ...
    # TODO: Test case where type(lcm(D1, D2)) = 21
    # ...

  
  
  def test_add_21_21(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4
    
    # Test case where type(lcm(D1, D2)) = 41 under various fields
    D1 = C34CurveDivisor(C_2, [[1, 1, 1], [1, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_2, [[0, 1, 1], [0, 1, 0, 1], []])
    D3 = C34CurveDivisor(C_2, [[0, 0, 1, 1], [1, 1, 0, 0, 1], [0, 1, 1, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_2_4, [[0, z4^2 + z4, 1], [z4^3 + z4^2 + z4, z4^2, 0, 1], []])
    D2 = C34CurveDivisor(C_2_4, [[z4^2 + 1, z4^3 + z4^2 + z4 + 1, 1], [z4 + 1, z4^2, 0, 1], []])
    D3 = C34CurveDivisor(C_2_4, [[z4^3 + 1, z4^3 + z4 + 1, z4^3 + z4^2 + 1, 1],
                                 [0, z4^2 + 1, z4^3 + 1, 0, 1],
                                 [z4^3 + z4, z4^3 + z4, z4 + 1, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_3, [[1, 2, 1], [0, 1, 0, 1], []])
    D2 = C34CurveDivisor(C_3, [[0, 2, 1], [0, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_3, [[1, 0, 1, 1], [0, 0, 2, 0, 1], [0, 0, 2, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_3_3, [[z3^2 + 2*z3, z3^2 + 2*z3 + 1, 1], [z3 + 1, z3^2, 0, 1], []])
    D2 = C34CurveDivisor(C_3_3, [[z3^2 + z3 + 1, 2*z3^2 + z3 + 1, 1],
                                 [1, z3^2 + z3 + 1, 0, 1], []])
    D3 = C34CurveDivisor(C_3_3, [[z3, 2*z3^2 + 2*z3 + 1, 2*z3 + 2, 1],
                                 [2*z3, 2*z3^2 + 1, z3^2, 0, 1],
                                 [z3, z3^2 + 2*z3 + 1, 2*z3^2 + 2*z3 + 2, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_31, [[10, 17, 1], [15, 14, 0, 1], []])
    D2 = C34CurveDivisor(C_31, [[29, 25, 1], [8, 17, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[30, 13, 29, 1], [18, 11, 13, 0, 1], [21, 15, 27, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_31_2, [[9*z2 + 19, 17*z2 + 9, 1], [14*z2 + 22, 28*z2 + 18, 0, 1], []])
    D2 = C34CurveDivisor(C_31_2, [[27, 13*z2 + 27, 1], [6*z2 + 24, 25*z2 + 10, 0, 1], []])
    D3 = C34CurveDivisor(C_31_2, [[6*z2 + 13, 3*z2 + 13, 11*z2 + 6, 1],
                                  [20*z2 + 15, 28*z2 + 17, 3*z2 + 1, 0, 1],
                                  [16*z2 + 16, 5*z2 + 22, 19*z2 + 9, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_1009, [[688, 179, 1], [360, 620, 0, 1], []])
    D2 = C34CurveDivisor(C_1009, [[13, 218, 1], [131, 369, 0, 1], []])
    D3 = C34CurveDivisor(C_1009, [[754, 638, 275, 1],
                                  [904, 369, 1007, 0, 1],
                                  [104, 959, 281, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    
    # TODO: Test case where type(lcm(D1, D2)) = 42
    # ...
    # TODO: Test case where type(lcm(D1, D2)) = 43
    # ...
    # TODO: Test case where type(lcm(D1, D2)) = 44
    # ...
    # TODO: Test case where type(lcm(D1, D2)) = 31
    # ...
    # TODO: Test case where type(lcm(D1, D2)) = 32
    # ...
    # TODO: Test case where type(lcm(D1, D2)) = 33 (impossible?)
    # ...



  def test_add_21_22(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    # Test case where type(lcm(D1, D2)) = 41 under various fields
    D1 = C34CurveDivisor(C_2, [[1, 1, 1], [1, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_2, [[0, 1], [1, 0, 1, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_2, [[0, 1, 1, 1], [1, 0, 1, 0, 1], [0, 1, 0, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_2_4, [[z4^3 + z4^2 + 1, z4^2 + z4, 1],
                                 [z4^3 + z4 + 1, z4^2 + z4 + 1, 0, 1], []])
    D2 = C34CurveDivisor(C_2_4, [[z4^3, 1], [z4^3, 0, z4^3 + z4^2 + z4, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_2_4, [[z4^2 + z4, z4^3 + z4^2 + z4, z4 + 1, 1],
                                 [z4^3 + z4^2 + z4, z4^3 + z4^2 + z4, z4^2, 0, 1],
                                 [z4^3 + z4^2 + 1, 1, z4^2 + z4 + 1, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_3, [[2, 0, 1], [2, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_3, [[0, 1], [1, 0, 2, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_3, [[0, 2, 2, 1], [2, 1, 1, 0, 1], [2, 1, 2, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_3_3, [[2*z3^2 + 2*z3, 2*z3^2 + 2*z3 + 1, 1],
                                 [z3^2, z3^2 + 2, 0, 1], []])
    D2 = C34CurveDivisor(C_3_3, [[z3^2 + z3 + 1, 1], [z3^2 + 2*z3, 0, 2*z3 + 2, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_3_3, [[z3^2 + 2*z3 + 1, z3 + 1, 2*z3 + 1, 1],
                                 [z3, z3 + 1, 2*z3^2 + z3 + 1, 0, 1],
                                 [z3^2 + z3, 2*z3^2 + 2*z3, 0, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_31, [[7, 23, 1], [8, 5, 0, 1], []])
    D2 = C34CurveDivisor(C_31, [[15, 1], [6, 0, 2, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[0, 18, 3, 1], [15, 20, 21, 0, 1], [15, 11, 10, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_31_2, [[28*z2 + 9, 24*z2 + 5, 1], [27*z2 + 3, 29*z2 + 3, 0, 1], []])
    D2 = C34CurveDivisor(C_31_2, [[16*z2 + 14, 1], [7*z2 + 15, 0, 29*z2 + 27, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_31_2, [[z2 + 29, 29*z2 + 3, 8*z2 + 13, 1],
                                  [3*z2 + 10, 2*z2 + 20, 13*z2 + 1, 0, 1],
                                  [18*z2 + 13, 14*z2 + 23, 3*z2 + 29, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_1009, [[35, 881, 1], [1006, 1007, 0, 1], []])
    D2 = C34CurveDivisor(C_1009, [[858, 1], [369, 0, 444, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_1009, [[809, 604, 751, 1],
                                  [535, 20, 774, 0, 1],
                                  [101, 384, 817, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)
    
    # TODO: Test case where type(lcm(D1, D2)) = 42
    # ...
    # TODO: Test case where type(lcm(D1, D2)) = 43
    # ...
    # TODO: Test case where type(lcm(D1, D2)) = 44
    # ...
    # TODO: Test case where type(lcm(D1, D2)) = 31
    # ...
    # TODO: Test case where type(lcm(D1, D2)) = 32
    # ...
    # TODO: Test case where type(lcm(D1, D2)) = 33
    # ...



  def test_add_22_11(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4
    
    # Test case where type(lcm(D1, D2)) = 31 under various fields
    D1 = C34CurveDivisor(C_2, [[1, 1], [1, 0, 0, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_2, [[1, 1], [1, 0, 1], []])
    D3 = C34CurveDivisor(C_2, [[1, 0, 0, 1], [1, 1, 1, 0, 1], [0, 1, 0, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_2_4, [[z4^2 + z4, 1], [z4^3 + z4, 0, z4^2 + z4, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_2_4, [[z4^3 + z4^2 + z4 + 1, 1], [z4^3 + z4 + 1, 0, 1], []])
    D3 = C34CurveDivisor(C_2_4, [[z4^2, z4^3 + 1, 0, 1],
                                 [z4^3 + z4^2 + z4 + 1, z4^3 + z4 + 1, z4^2 + z4, 0, 1],
                                 [z4^2 + 1, z4^3 + z4 + 1, z4^2 + z4, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_3, [[0, 1], [0, 0, 1, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_3, [[1, 1], [2, 0, 1], []])
    D3 = C34CurveDivisor(C_3, [[0, 1, 0, 1], [0, 2, 0, 0, 1], [0, 2, 1, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_3_3, [[z3^2 + 2*z3 + 2, 1], [z3^2, 0, 2*z3^2 + 2, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_3_3, [[z3^2 + 2*z3 + 2, 1], [2, 0, 1], []])
    D3 = C34CurveDivisor(C_3_3, [[2*z3, 2*z3^2 + z3 + 1, 0, 1],
                                 [2*z3^2 + z3 + 1, 2, z3^2 + 2*z3 + 2, 0, 1],
                                 [2*z3^2 + z3, z3 + 2, 2*z3^2 + 2, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_31, [[30, 1], [1, 0, 17, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_31, [[19, 1], [18, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[12, 18, 0, 1], [13, 18, 30, 0, 1], [14, 18, 17, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_31_2, [[5*z2 + 18, 1], [22*z2 + 30, 0, 30*z2 + 3, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_31_2, [[24*z2 + 6, 1], [9*z2 + 25, 0, 1], []])
    D3 = C34CurveDivisor(C_31_2, [[20*z2 + 27, 29*z2 + 24, 0, 1],
                                  [5*z2 + 5, 9*z2 + 25, 5*z2 + 18, 0, 1],
                                  [17*z2 + 6, 28*z2 + 22, 30*z2 + 3, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_1009, [[341, 1], [435, 0, 354, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_1009, [[566, 1], [794, 0, 1], []])
    D3 = C34CurveDivisor(C_1009, [[287, 907, 0, 1],
                                  [342, 794, 341, 0, 1],
                                  [485, 299, 354, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    # TODO: Test case where type(lcm(D1, D2)) = 32
    # ...
    # TODO: Test case where type(lcm(D1, D2)) = 33
    # ...
    # TODO: Test case where type(lcm(D1, D2)) = 22




  def test_add_22_22(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4
    
    # Test most common case where type(lcm(D1, D2)) = 43 under various fields
    D1 = C34CurveDivisor(C_2, [[0, 1], [1, 0, 1, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_2, [[0, 1], [1, 0, 1, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_2, [[0, 0, 1], [0, 1, 0, 1], []])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_2_4, [[z4 + 1, 1], [z4^2, 0, z4^3 + z4 + 1, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_2_4, [[z4^2 + z4 + 1, 1], [z4^3 + z4 + 1, 0, z4^2, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_2_4, [[z4^3 + 1, z4^3 + z4^2 + z4, 1],
                                 [z4^3 + z4^2 + 1, 1, 0, 1], []])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_3, [[2, 1], [0, 0, 2, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_3, [[0, 1], [1, 0, 2, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_3, [[0, 2, 1], [0, 0, 0, 1], []])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_3_3, [[z3, 1], [z3^2 + 2, 0, z3^2 + 2*z3 + 1, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_3_3, [[1, 1], [0, 0, 2*z3 + 2, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_3_3, [[z3^2 + 2*z3 + 1, z3, 1],
                                 [z3^2 + 2*z3 + 1, z3^2 + 2*z3 + 2, 0, 1], []])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_31, [[19, 1], [23, 0, 21, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_31, [[25, 1], [28, 0, 14, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[18, 21, 1], [30, 13, 0, 1], []])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_31_2, [[17*z2 + 18, 1], [21*z2 + 14, 0, 10*z2 + 9, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_31_2, [[21*z2 + 1, 1], [3*z2 + 27, 0, 21*z2 + 26, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_31_2, [[27*z2 + 10, 27*z2 + 21, 1],
                                  [26*z2 + 9, 11*z2 + 22, 0, 1], []])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_1009, [[67, 1], [81, 0, 579, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_1009, [[697, 1], [870, 0, 233, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_1009, [[836, 464, 1], [424, 51, 0, 1], []])
    self.assertEqual(D1 + D2, D3)

    # TODO: Test other cases (42, ...)
    # ...



  def test_add_31_11(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    # Test most common case where type(lcm(D1, D2)) = 41 under various fields
    D1 = C34CurveDivisor(C_2, [[0, 0, 1, 1], [0, 0, 1, 0, 1], [0, 0, 1, 0, 0, 1]])
    D2 = C34CurveDivisor(C_2, [[1, 1], [0, 0, 1], []])
    D3 = C34CurveDivisor(C_2, [[0, 1, 0, 1], [0, 1, 0, 0, 1], [1, 1, 1, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_2_4, [[z4^3 + z4, z4 + 1, z4^2 + 1, 1],
                                 [z4^3 + z4^2 + z4, z4^2 + 1, z4^2 + z4, 0, 1],
                                 [z4^2, z4^3 + z4^2 + z4, 1, 0, 0, 1]])
    D2 = C34CurveDivisor(C_2_4, [[z4^3 + z4^2 + z4 + 1, 1], [1, 0, 1], []])
    D3 = C34CurveDivisor(C_2_4, [[1, z4^3 + z4^2 + z4 + 1, 1, 1],
                                 [z4^3 + z4^2 + 1, z4^3, z4^3 + z4 + 1, 0, 1],
                                 [z4^3 + 1, z4^3 + z4, z4 + 1, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_3, [[0, 2, 2, 1], [2, 1, 1, 0, 1], [2, 1, 2, 0, 0, 1]])
    D2 = C34CurveDivisor(C_3, [[1, 1], [2, 0, 1], []])
    D3 = C34CurveDivisor(C_3, [[0, 2, 0, 1], [0, 0, 0, 0, 1], [0, 0, 1, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_3_3, [[z3 + 1, 2*z3^2 + 2*z3, z3^2 + 2, 1],
                                 [2*z3^2 + 2*z3 + 2, 2*z3^2 + 2*z3, z3^2 + 2, 0, 1],
                                 [z3, z3^2 + 1, 2*z3 + 2, 0, 0, 1]])
    D2 = C34CurveDivisor(C_3_3, [[z3^2 + z3 + 2, 1], [2*z3 + 1, 0, 1], []])
    D3 = C34CurveDivisor(C_3_3, [[2*z3^2 + 2, 2*z3^2 + 1, 2*z3^2 + z3 + 2, 1],
                                 [z3^2 + z3, 2, z3^2 + z3, 0, 1],
                                 [z3^2 + 1, 2*z3^2 + z3 + 2, z3^2 + 2*z3 + 1, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_31, [[15, 4, 17, 1], [23, 4, 29, 0, 1], [10, 0, 22, 0, 0, 1]])
    D2 = C34CurveDivisor(C_31, [[15, 1], [16, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[18, 13, 7, 1], [29, 4, 14, 0, 1], [10, 26, 13, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_31_2, [[26*z2 + 15, 2*z2 + 10, 9*z2 + 26, 1],
                                  [11*z2 + 17, 27*z2 + 14, 14*z2 + 11, 0, 1],
                                  [z2 + 28, 28*z2 + 11, 16*z2 + 14, 0, 0, 1]])
    D2 = C34CurveDivisor(C_31_2, [[2*z2 + 30, 1], [14*z2 + 28, 0, 1], []])
    D3 = C34CurveDivisor(C_31_2, [[4*z2 + 24, 28*z2 + 13, 20*z2 + 8, 1],
                                  [10, 24*z2 + 30, 29*z2 + 7, 0, 1],
                                  [18*z2 + 11, 12*z2 + 9, 10*z2 + 1, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_1009, [[221, 507, 55, 1],
                                  [112, 772, 678, 0, 1],
                                  [38, 507, 389, 0, 0, 1]])
    D2 = C34CurveDivisor(C_1009, [[448, 1], [82, 0, 1], []])
    D3 = C34CurveDivisor(C_1009, [[296, 653, 266, 1],
                                  [802, 145, 938, 0, 1],
                                  [482, 789, 878, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    # TODO: Test other cases
    # ...


    
  def test_add_31_21(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    # Test most common case where type(lcm(D1, D2)) = 51 under various fields
    D1 = C34CurveDivisor(C_2, [[0, 1, 0, 1], [0, 0, 1, 0, 1], [0, 1, 0, 0, 0, 1]])
    D2 = C34CurveDivisor(C_2, [[0, 1, 1], [0, 1, 0, 1], []])
    D3 = C34CurveDivisor(C_2, [[1, 0, 1, 1], [1, 1, 0, 0, 1], [1, 1, 0, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_2_4, [[z4^3 + z4^2 + 1, z4^3 + z4 + 1, z4^3 + z4 + 1, 1], [z4^3 + z4^2, z4^3 + z4 + 1, z4^2, 0, 1], [z4 + 1, z4^3 + z4^2 + 1, z4^3 + z4^2, 0, 0, 1]])
    D2 = C34CurveDivisor(C_2_4, [[z4^3, z4^3 + z4^2, 1], [z4^3 + z4^2 + z4 + 1, z4^3 + z4^2 + z4 + 1, 0, 1], []])
    D3 = C34CurveDivisor(C_2_4, [[z4^2 + z4, z4^3 + z4^2 + z4, z4^2 + 1, 1], [z4^3 + z4^2 + 1, 0, z4^3 + z4 + 1, 0, 1], [z4^3 + z4^2 + 1, z4^2 + z4, z4^2, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_3, [[0, 2, 0, 1], [0, 2, 0, 0, 1], [1, 2, 2, 0, 0, 1]])
    D2 = C34CurveDivisor(C_3, [[1, 1, 1], [0, 2, 0, 1], []])
    D3 = C34CurveDivisor(C_3, [[1, 0, 1, 1], [0, 0, 0, 0, 1], [0, 0, 1, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_3_3, [[z3^2 + z3, z3^2 + 1, z3^2 + z3 + 1, 1], [2*z3^2 + 2*z3 + 2, 2*z3^2 + 2*z3 + 2, 2*z3^2 + 2*z3 + 1, 0, 1], [z3^2, z3^2 + z3, 1, 0, 0, 1]])
    D2 = C34CurveDivisor(C_3_3, [[z3^2 + z3 + 2, z3^2 + z3, 1], [z3^2 + 2*z3, 2*z3^2 + 2*z3 + 2, 0, 1], []])
    D3 = C34CurveDivisor(C_3_3, [[z3^2, 0, z3^2 + z3 + 2, 1], [2*z3^2 + 2*z3 + 1, z3 + 1, z3^2 + 2, 0, 1], [2*z3^2 + 2, z3^2 + z3 + 1, 2*z3^2 + 2*z3, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_31, [[6, 20, 7, 1], [6, 19, 16, 0, 1], [4, 16, 24, 0, 0, 1]])
    D2 = C34CurveDivisor(C_31, [[23, 26, 1], [8, 4, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[16, 0, 20, 1], [11, 18, 23, 0, 1], [10, 0, 22, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_31_2, [[3*z2 + 19, 9*z2 + 22, 6*z2 + 9, 1], [10*z2 + 18, 10*z2 + 6, 24*z2 + 3, 0, 1], [22*z2, 4*z2 + 21, 22*z2 + 16, 0, 0, 1]])
    D2 = C34CurveDivisor(C_31_2, [[7*z2 + 10, 13*z2 + 19, 1], [2*z2 + 9, 28*z2 + 10, 0, 1], []])
    D3 = C34CurveDivisor(C_31_2, [[25*z2 + 14, 11*z2 + 22, 10*z2 + 21, 1], [z2 + 15, 13*z2 + 27, 8*z2 + 28, 0, 1], [20*z2 + 27, 26*z2, 11*z2 + 20, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_1009, [[505, 854, 752, 1], [933, 33, 629, 0, 1], [167, 606, 248, 0, 0, 1]])
    D2 = C34CurveDivisor(C_1009, [[847, 396, 1], [816, 514, 0, 1], []])
    D3 = C34CurveDivisor(C_1009, [[32, 831, 986, 1], [91, 1000, 91, 0, 1], [703, 566, 373, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    # TODO: Test other cases
    # ...



  def test_add_31_22(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4
    
    # Test most common case where type(lcm(D1, D2)) = 51 under various fields
    D1 = C34CurveDivisor(C_2, [[1, 0, 1, 1], [1, 1, 0, 0, 1], [1, 1, 0, 0, 0, 1]])
    D2 = C34CurveDivisor(C_2, [[1, 1], [0, 0, 1, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_2, [[0, 1, 1, 1], [1, 0, 1, 0, 1], [0, 1, 0, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)
    D1 = C34CurveDivisor(C_2_4, [[z4^3, 1, z4^3 + z4^2 + z4, 1], [1, z4^3 + 1, 1, 0, 1], [z4^3 + z4^2, z4^3 + z4 + 1, z4, 0, 0, 1]])
    D2 = C34CurveDivisor(C_2_4, [[z4^3 + z4 + 1, 1], [z4^3 + z4^2 + z4, 0, 1, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_2_4, [[0, z4^3 + z4^2 + z4, z4^2 + z4, 1], [z4^3 + z4^2 + z4, z4^3 + z4^2 + 1, z4^2, 0, 1], [1, z4^3 + z4 + 1, 1, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)
    D1 = C34CurveDivisor(C_3, [[0, 2, 1, 1], [1, 0, 0, 0, 1], [1, 2, 0, 0, 0, 1]])
    D2 = C34CurveDivisor(C_3, [[2, 1], [0, 0, 2, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_3, [[0, 2, 0, 1], [0, 0, 2, 0, 1], [0, 1, 1, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)
    D1 = C34CurveDivisor(C_3_3, [[2*z3^2 + z3, 0, 0, 1], [z3^2 + 2*z3 + 1, z3 + 2, 2*z3^2, 0, 1], [2*z3^2 + z3, 2*z3^2 + 1, z3^2 + z3 + 1, 0, 0, 1]])
    D2 = C34CurveDivisor(C_3_3, [[z3, 1], [z3^2 + 2, 0, z3^2 + 2*z3 + 1, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_3_3, [[2*z3^2 + 2*z3 + 2, 2*z3^2 + 2*z3 + 2, 2*z3^2 + z3, 1], [z3^2 + 2, 0, z3^2 + z3 + 1, 0, 1], [2*z3, 2*z3^2 + 2, 2, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)
    D1 = C34CurveDivisor(C_31, [[13, 20, 30, 1], [29, 9, 0, 0, 1], [29, 29, 27, 0, 0, 1]])
    D2 = C34CurveDivisor(C_31, [[11, 1], [26, 0, 27, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[26, 11, 12, 1], [18, 11, 4, 0, 1], [3, 28, 16, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)
    D1 = C34CurveDivisor(C_31_2, [[11*z2 + 28, 20*z2 + 27, 18*z2 + 27, 1], [18*z2 + 21, 3*z2 + 28, 27*z2 + 23, 0, 1], [z2 + 16, 28*z2 + 23, 9*z2 + 4, 0, 0, 1]])
    D2 = C34CurveDivisor(C_31_2, [[24*z2 + 28, 1], [9*z2 + 12, 0, 7*z2 + 1, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_31_2, [[10*z2 + 22, 29*z2 + 13, 7*z2 + 6, 1], [11, 22*z2 + 20, 25*z2 + 9, 0, 1], [19*z2 + 7, 17*z2 + 9, 16*z2 + 26, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)
    D1 = C34CurveDivisor(C_1009, [[35, 840, 550, 1], [118, 47, 846, 0, 1], [340, 960, 370, 0, 0, 1]])
    D2 = C34CurveDivisor(C_1009, [[436, 1], [327, 0, 276, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_1009, [[529, 488, 1, 1], [538, 172, 205, 0, 1], [283, 416, 199, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    # TODO: Test other cases
    # ...



  def test_add_31_31(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4
    
    # Test most common case where type(lcm(D1, D2)) = 61 under various fields
    D1 = C34CurveDivisor(C_2, [[0, 1, 0, 1], [0, 1, 0, 0, 1], [1, 1, 1, 0, 0, 1]])
    D2 = C34CurveDivisor(C_2, [[0, 1, 1, 1], [0, 0, 1, 0, 1], [0, 0, 0, 0, 0, 1]])
    D3 = C34CurveDivisor(C_2, [[1, 0, 0, 1], [0, 0, 1, 0, 1], [1, 1, 1, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)
    D1 = C34CurveDivisor(C_2_4, [[z4^3 + z4^2 + z4 + 1, z4^3 + z4^2 + z4 + 1, 1, 1], [z4^3, z4^2 + z4 + 1, z4^3 + z4^2, 0, 1], [0, z4^3 + z4, z4^3 + z4^2 + z4 + 1, 0, 0, 1]])
    D2 = C34CurveDivisor(C_2_4, [[z4^3 + z4 + 1, z4^3 + z4, z4^2 + z4 + 1, 1], [z4^3 + z4^2 + z4 + 1, z4^3 + 1, z4^2 + z4 + 1, 0, 1], [z4^3 + z4^2 + z4 + 1, z4^3 + z4^2 + 1, z4^3 + z4 + 1, 0, 0, 1]])
    D3 = C34CurveDivisor(C_2_4, [[z4^2 + z4, z4^3 + z4^2 + z4 + 1, z4 + 1, 1], [z4^2 + 1, z4 + 1, z4, 0, 1], [z4, 1, z4^2 + z4, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)
    D1 = C34CurveDivisor(C_3, [[0, 1, 0, 1], [0, 2, 0, 0, 1], [1, 1, 2, 0, 0, 1]])
    D2 = C34CurveDivisor(C_3, [[1, 1, 2, 1], [2, 1, 0, 0, 1], [1, 2, 0, 0, 0, 1]])
    D3 = C34CurveDivisor(C_3, [[0, 2, 0, 1], [0, 0, 2, 0, 1], [0, 0, 2, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)
    D1 = C34CurveDivisor(C_3_3, [[2*z3^2, 2*z3^2 + 1, z3^2, 1], [2*z3, 2*z3^2 + 2*z3 + 2, 2*z3^2, 0, 1], [z3 + 2, z3 + 2, 2*z3^2 + 2*z3 + 2, 0, 0, 1]])
    D2 = C34CurveDivisor(C_3_3, [[z3^2, z3^2 + 2, 0, 1], [z3^2 + 2, z3^2 + 2*z3 + 2, 2*z3^2 + z3 + 2, 0, 1], [z3^2 + 2, z3 + 1, z3^2, 0, 0, 1]])
    D3 = C34CurveDivisor(C_3_3, [[2*z3^2 + z3, 2*z3^2 + 2*z3 + 2, z3^2 + 1, 1], [2*z3^2 + 2*z3, 2*z3^2 + 2, z3^2 + 2*z3 + 1, 0, 1], [2*z3^2, z3^2 + z3 + 2, z3^2, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)
    D1 = C34CurveDivisor(C_31, [[26, 4, 2, 1], [17, 24, 22, 0, 1], [0, 23, 18, 0, 0, 1]])
    D2 = C34CurveDivisor(C_31, [[9, 21, 4, 1], [0, 2, 4, 0, 1], [20, 2, 26, 0, 0, 1]])
    D3 = C34CurveDivisor(C_31, [[28, 9, 21, 1], [15, 20, 21, 0, 1], [19, 6, 23, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)
    D1 = C34CurveDivisor(C_31_2, [[5*z2 + 18, 10*z2 + 19, 9*z2 + 13, 1], [19*z2 + 30, 11*z2 + 3, 12*z2 + 16, 0, 1], [6*z2 + 24, 10*z2 + 6, 10*z2 + 5, 0, 0, 1]])
    D2 = C34CurveDivisor(C_31_2, [[24*z2 + 26, z2 + 29, 30*z2 + 6, 1], [21*z2 + 2, 14*z2 + 28, 8, 0, 1], [5*z2 + 6, 25*z2 + 9, 10*z2 + 27, 0, 0, 1]])
    D3 = C34CurveDivisor(C_31_2, [[7*z2 + 1, 5, 11*z2 + 3, 1], [16*z2 + 26, 14*z2 + 9, 12*z2 + 26, 0, 1], [16*z2 + 9, 13*z2 + 3, z2 + 18, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)
    D1 = C34CurveDivisor(C_1009, [[351, 509, 773, 1], [395, 206, 757, 0, 1], [450, 72, 961, 0, 0, 1]])
    D2 = C34CurveDivisor(C_1009, [[1007, 845, 679, 1], [452, 771, 469, 0, 1], [165, 413, 223, 0, 0, 1]])
    D3 = C34CurveDivisor(C_1009, [[877, 43, 18, 1], [191, 561, 937, 0, 1], [1000, 380, 341, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    # TODO: Test other cases
    # ...
