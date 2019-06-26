
"""
def run_tests() :
  def test_create() :
    print("Creating degree 1 divisors.")
    for i in range(N_TRIALS) :
      if not test_create1() :
        return False
    print("OK!")

    print("Creating degree 2 divisors.")
    for i in range(N_TRIALS) :
      if not test_create2() :
        return False
    print("OK!")

    print("Creating degree 3 divisors.")
    for i in range(N_TRIALS) :
      if not test_create3() :
        return False
    print("OK!")
    return True

  if not test_create() :
    return
  print("All tests OK!")
"""


def add_formal_sums(sum1, sum2) :
  """
    A formal sum is represented by a list of pairs.
    The first element of each pair is a C34CrvPt and the second element is an integer.
    For example, if P = (3 : 2 : 1) and Q = (5 : 0 : 1), the formal sum 2P + Q is represented by
    
      [((3 : 2 : 1), 2), ((5 : 0 : 1), 1)]
    
    This function adds two formal sums.
    E.g. add_formal_sums(2P + Q, P + Q + R) returns the list representation of 3P + 2Q + R
    
    Lists are assumed to be sorted according to an order on C34CrvPt's
  """
  ret = []
  i, j = 0, 0
  while (i < len(sum1)) or (j < len(sum2)) :
    if i == len(sum1) :
      ret = ret + [sum2[j]]
      j = j + 1
    elif j == len(sum2) :
      ret = ret + [sum1[i]]
      i = i + 1
    elif sum1[i][0] == sum2[j][0] :
      ret = ret + [(sum1[i][0], sum1[i][1] + sum2[j][1])]
      i = i + 1
      j = j + 1
    elif sum1[i][0] < sum2[j][0] :
      ret = ret + [sum1[i]]
      i = i + 1
    else :
      ret = ret + [sum2[j]]
      j = j + 1
  return ret

def random_curve(q) :
  base = FiniteField(q).prime_subfield()
  degree = FiniteField(q).degree()
  Kbar = base.algebraic_closure()
  K = Kbar.subfield(degree)[0]
  
  # Pick a random curve
  c = [K.random_element() for i in range(9)]
  C = C34Crv(K, c)
  return C



def test_create1() :
  """
    Tests creation of divisors of degree 1.
  """
  # Pick a random curve
  C = random_curve()

  # Pick a random point on the curve
  P = C.random_point()
  
  # Create a divisor
  D = C34CrvDiv(C, [P])
  
  # Verify
  V = D.variety()
  x, y = C.R.gens()
  gave = [[P[0],P[1]]]
  got = [[t[x], t[y]] for t in V]
  gave.sort()
  got.sort()
  
  if (gave == got) :
    return True

  print("p = {}".format(p))
  print("K = FiniteField(p)")
  print("c = {}".format(c))
  print("C = C34Crv(K, c)")
  print("P = C.point({}, {})".format(P[0], P[1]))
  print("D = C34CrvDiv(C, [P])")
  return False



def test_create2() :
  """
    Tests creation of divisors of degree 2.
  """
  # Pick a random curve
  C = random_curve()

  # Pick two random points on the curve
  P1 = C.random_point()
  P2 = C.random_point()
  
  # Create a divisor
  D = C34CrvDiv(C, [P1, P2])
  
  # Verify
  V = D.variety()
  x, y = C.R.gens()
  gave = [[P1[0],P1[1]], [P2[0],P2[1]]]
  got = [[t[x], t[y]] for t in V]
  gave.sort()
  got.sort()
  
  if (gave == got) :
    return True

  print("p = {}".format(p))
  print("K = FiniteField(p)")
  print("c = {}".format(c))
  print("C = C34Crv(K, c)")
  print("P1 = C.point({}, {})".format(P1[0], P1[1]))
  print("P2 = C.point({}, {})".format(P2[0], P2[1]))
  print("D = C34CrvDiv(C, [P1, P2])")
  return False



def test_create3() :
  """
    Tests creation of divisors of degree 3.
  """
  # Pick a random curve
  C = random_curve()

  # Pick two random points on the curve
  P1 = C.random_point()
  P2 = C.random_point()
  P3 = C.random_point()
  
  # Create a divisor
  D = C34CrvDiv(C, [P1, P2, P3])
  
  # Verify
  V = D.variety()
  x, y = C.R.gens()
  gave = [[P1[0],P1[1]], [P2[0],P2[1]], [P3[0],P3[1]]]
  got = [[t[x], t[y]] for t in V]
  gave.sort()
  got.sort()
  
  if (gave == got) :
    return True

  print("p = {}".format(p))
  print("K = FiniteField(p)")
  print("c = {}".format(c))
  print("C = C34Crv(K, c)")
  print("P1 = C.point({}, {})".format(P1[0], P1[1]))
  print("P2 = C.point({}, {})".format(P2[0], P2[1]))
  print("P2 = C.point({}, {})".format(P3[0], P3[1]))
  print("D = C34CrvDiv(C, [P1, P2, P3])")
  
  return False



def test_double() :
  return test_double1()

def test_double1() :
  # For each curve, if possible, find a point whose tangent line is vertical
  # This happens where the partial derivative w.r.t. y of the curve equation is 0
  print "Testing doubling of degree 1 divisors..."

  # Test over a curve of characteristic 2
  def test_case(q) :
    try :
      C = random_curve(q)
    except ValueError as e:
      # If curve was singular
      print e
      return test_case(q)
    R = C.R
    K = C.K
    x, y = R.gens()

    # Find points whose tangent line is vertical
    x, y = C.R.gens()
    F = C.poly()
    Fy = C.partial_y(x,y)
    V = R.ideal(F, Fy).variety(K)
    for d in V :
      x0, y0 = d[x], d[y]
      P = C.point(x0, y0)
      D = C.divisor([P])
      DD = double1(D)
      FS1 = D.formal_sum()
      expected = add_formal_sums(FS1, FS1)
      got = DD.formal_sum()
      if expected != got :
        return C, D

    # Test some random points
    n_points = 10
    for i in range(n_points) :
      P = C.random_point()
      D = C.divisor([P])
      DD = double1(D)
      FS1 = D.formal_sum()
      expected = add_formal_sums(FS1, FS1)
      got = DD.formal_sum()
      if expected != got :
        return C, D
    return

  for q in ORDERS :
    print("Testing over field of order {} (characteristic {}).".format(q, radical(q)))
    res = test_case(q)
    if res != None :
      print ("   Fail. Returning curve and divisor.")
      return res
  print "All tests passed!"
  print
  return

