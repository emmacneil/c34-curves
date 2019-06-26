

def sato_tate(f) :
  def count(f, p) :
    """
      Count the number of F_p-rational points on the curve defined by f.
    """
    K = GF(p)
    R.<x,y> = K[]
    A.<x,y> = AffineSpace(K, 2)
    X = A.subscheme(f)
    return X.count_points(1)
  
  for p in primes(2, 10) :
    print ("p = {}, Np = {}".format(p, count(f,p)))
