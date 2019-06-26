def sage_add(D1, D2) :
  D3 = sage_compose(D1, D2)
  return D3 if D3.reduced else sage_flip(sage_flip(D3))

def sage_compose(D1, D2) :
  C = D1.C
  R = C.R
  F = C.poly()
  x, y = R.gens()
  
  I = R.ideal(D1.polys())*R.ideal(D2.polys()) + R.ideal(F)
  J = R.ideal(I.groebner_basis())
  polys = J.gens()[0:]
  f = polys[0]
  if f.lm() == y^3 :
    polys = polys[1:]
    
  return C34CrvDiv(C, polys)

def sage_flip(D) :
  C = D.C
  R = C.R
  F = C.poly()
  x, y = R.gens()
  
  I = D.ideal()
  f = D.polys()[0]
  J = R.ideal(f, F)
  Q = R.ideal(J.quotient(I).groebner_basis())
  polys = Q.gens()[0:]
  f = polys[0]
  if f.lm() == y^3 :
    polys = polys[1:]
  
  return C34CrvDiv(C, polys)

def sage_gcd(D1, D2) :
  C = D1.C
  R = C.R
  F = C.poly()
  x, y = R.gens()
  
  I = R.ideal(D1.polys()) + R.ideal(D2.polys()) + R.ideal(F)
  J = R.ideal(I.groebner_basis())
  polys = J.gens()[0:]
  f = polys[0]
  if f.lm() == y^3 :
    polys = polys[1:]
    
  return C34CrvDiv(C, polys)

def sage_scale(D, n) :
  if (n < 0) :
    return sage_scale(sage_flip(D), -n)

  C = D.C
  R = C.R
  F = C.poly()
  x, y = R.gens()

  ret = C.zero_divisor()
  base = D
  while (n > 0) :
    if (n & 1 == 1) :
      ret = sage_add(ret, base)
    n = n >> 1
    base = sage_add(base, base)
  return ret

