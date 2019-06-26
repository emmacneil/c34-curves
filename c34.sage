import random, timeit

load("c34add.sage")
load("c34crv.sage")
load("c34crvdiv.sage")
load("c34crvpt.sage")
load("c34double.sage")
load("c34flip.sage")
load("c34slow.sage")
load("c34test.sage")
load("c34triple.sage")
load("c34util.sage")

#suite = unittest.TestLoader().loadTestsFromTestCase(TestAdd)
#unittest.TextTestRunner(verbosity=2).run(suite)
suite = unittest.TestLoader().loadTestsFromTestCase(TestDouble)
unittest.TextTestRunner(verbosity=2).run(suite)
#suite = unittest.TestLoader().loadTestsFromTestCase(TestFlip)
#unittest.TextTestRunner(verbosity=2).run(suite)

q = 1009
Kbar = GF(q).prime_subfield().algebraic_closure()
K = Kbar.subfield(GF(q).degree())[0]

R = PolynomialRing(K, 2, "x,y", order="m(3,4,0,1)")
x,y = R.gens()
#c = [K.random_element() for t in range(9)]
c = [K(t) for t in [105, 392, 370, 162, 517, 810, 26, 260, 148]]
C = C34Crv(K, c)
F = C.poly()
print C

AA = C.ambient_space()
X = C.scheme() # Curve as an affine scheme.



def find_case() : 
  max_time = 30 # seconds
  t0 = timeit.default_timer()
  ret = 0, 0
  while timeit.default_timer() - t0 < max_time : 
    D1 = find_reduced_divisor(C_41, 21)
    D2 = sage_compose(D1, D1)
    if D1.type == 42 :
      ret = D1, D2
      print("D1 = C34CrvDiv(C_41, {})".format([D1.f, D1.g, D1.h]))
      print("D2 = C34CrvDiv(C_41, {})".format([D2.f, D2.g, D2.h]))
      print("self.assertEqual(double_21(D1), D2)")
      break
  return ret

def find_atypicals(C) :
  Pts = C.points()
  for P in Pts :       
    for Q in Pts :       
      D = C.divisor([P, Q])   
      E = sage_compose(D, D)
      if E.type == 42 :
      #if not E.typical :    
        print("Type {} : {}".format(E.type, E))
        print E.formal_sum()
                  

def find_type(C, t, typical) :
  """
    Finds a divisor of type t on the curve C.
    Only works for divisors of degree > 1
  """
  pt_set = C.points()
  deg = t // 10
  deg2 = deg // 2
  deg1 = deg - deg2
  
  while True :
    D1 = C.divisor([ random.choice(pt_set) for i in range(deg1) ])
    D2 = C.divisor([ random.choice(pt_set) for i in range(deg2) ])
    
    if (D1.type not in [11, 21, 22, 31]) :
      continue
    if (D2.type not in [11, 21, 22, 31]) :
      continue
    
    D3 = C34CrvDiv(C, [[K.one()], [], []])
    case = (D1.type, D2.type)
    if case == (11, 11) :
      D3 = add11(D1, D2)
    elif case == (21, 11) :
      D3 = add2a1(D1, D2)
    elif case == (21, 21) :
      D3 = add2a2a(D1, D2)
    elif case == (21, 22) :
      D3 = add2b2a(D2, D1)
    elif case == (31, 11) :
      D3 = add31(D1, D2)
    elif case == (31, 21) :
      D3 = add32a(D1, D2)
    elif case == (31, 22) :
      #D3 = add32b(D1, D2)
      continue
    elif case == (31, 31) :
      D3 = add33(D1, D2)

    if (D3.type == t) and (D3.typical == typical) :
      return D3

def test_script() :
  N_CURVES = 100
  N_TRIALS_PER_CURVE = 100
  MAX_PRIME = 100
  prime_list = list(primes(MAX_PRIME))
  
  for j in range(N_CURVES) :
    # Get a random prime-order finite field, K
    q = random.choice(prime_list)
    Kbar = GF(q).prime_subfield().algebraic_closure()
    K = Kbar.subfield(GF(q).degree())[0]

    # Get a random curve over K
    R = PolynomialRing(K, 2, "x,y", order="m(3,4,0,1)")
    c = [K.random_element() for t in range(9)]
    try :
      C = C34Crv(K, c)
    except :
      continue
    print C
    pt_set = C.points()
    if len(pt_set) == 0:
      continue

    for i in range(N_TRIALS_PER_CURVE) :
      deg1 = randint(1, 2)
      deg2 = randint(1, 2)
      D1 = C.divisor([ random.choice(pt_set) for t in range(deg1) ])
      D2 = C.divisor([ random.choice(pt_set) for t in range(deg2) ])
      E = sage_add(D1, D2)
      try :
        G = D1 + D2
        if (E != G) :
          return D1, D2
      except :
        return D1, D2
  return 0, 0

def find_fail_case(verify = True) :
  def int_to_bin(n) :
    ret = [0]*9
    # Assume 0 <= n < 512
    for i in range(9) :
      ret[i] = n % 5
      n = n // 5
    return ret
  
  pt_set = C.points()

  # Find a semi-typical type 31 divisor 

  ctr = 0
  while ctr < 1000 :
    ctr = ctr + 1
    print("Test {}".format(ctr))
    #P1, P2, P3, P4, P5, P6 = [ random.choice(pt_set) for t in range(6) ]
    #shuffle(pt_set)
    #P1, P2, P3, P4, P5 = pt_set[0:5]
    #D1 = C.divisor([P1, P2])
    #D2 = C.divisor([P4, P5])
    #if (D1.type == 22) or (D2.type == 22):
    #  continue
    #try :
    #  D3 = D1 + D2
    #except :
    #  return D1, D2
    #if verify :
    #  expected = sage_add(D1, D2)
    #  if (expected != D3) :
    #    return D1, D2
    D = find_type(C, 31, False)
    expected = sage_flip(D)
    got = flip(D)
    if (expected != got) :
      return D
  return 0

def subtract(D1, D2) :
  """ 
    Returns the divisor D = D1 - D2.
    I.e. returns a divisor D such that D2 + D = D1.
    Assumes D1 and D2 are reduced divisors of degrees 3 and 1, respectively.
  """
  K = D1.K
  f, g, h = D1.f, D1.g, D1.h
  F, H = D2.f, D2.g
  new_f, new_g, new_h = [], [], []
  
  if f[2] != 0 :
    s1 = (F[0] - g[2])/f[2]
    s0 = g[1] + s1*(f[1] - F[0])
    t1 = f[1] - s1*f[2]
    t0 = f[0] - s0*f[2]
    new_f = [ s0, s1, K.one() ]
    new_g = [ t0, t1, K.zero(), K.one() ]
  elif g[1] != -G[0] :
    s1 = h[1]/(g[1] + G[0])
    s0 = s1*g[2] - G[0]
    t1 = f[1] - s1*f[2]
    t0 = f[0] - s0*f[2]
    new_f = [ s0, s1, K.one() ]
    new_g = [ t0, t1, K.zero(), K.one() ]
  else :
    s0 = g[2]
    t2 = h[2]
    t0 = h[0] - h[1]*g[2]
    new_f = [ s0, K.one() ]
    new_g = [ t0, K.zero(), t2, K.zero(), K.zero(), K.one() ]
    
  return C34CrvDiv(D1.C, [new_f, new_g, new_h])


