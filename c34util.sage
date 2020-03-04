def differentials(D) :
  """
    Given a type 31 divisor D = div(f, g, h), finds polynomials r, s, t, r',
    s', and t' (with leading monomials y, x, 1, x^2, 1, y, resp.) satisfying

       rf +  sg +  th = 0, and
      r'f + s'g + t'h = F,

    computes

      df := st' - ts',
      dg := tr' - rt', and
      dh := rs' - sr',

    and returns (df, dg, dh). The divisor A = div(df, dg, dh) has type 41 and
    is linearly equivalent to D.
  """
  if (D.type != 31) :
    raise ValueError("D must be of type 31.")

  c0, c1, c2, c3, c4, c5, c6, c7, c8 = D.parent_curve().coefficients()
  f0, f1, f2 = D.f[0:3]
  g0, g1, g2 = D.g[0:3]
  h0, h1, h2 = D.h[0:3]

  r0 = g1
  s0 = g2 - f1
  t0 = -f2

  T1 = c8
  R2 = c7 - f2
  R1 = c6 - f1
  T0 = c5 - h2 - f2*R2
  S0 = c4 - h2*T1 - h1 - f2*R1 - f1*R2
  R0 = c3 - h1*T1 - f1*R1 - f0

  x, y = D.parent_curve().polynomial_ring().gens()
  r = y + r0
  s = -x + s0
  t = t0
  rr = x*x + R2*y + R1*x + R0
  ss = S0
  tt = y + T1*x + T0

  f, g, h = D.groebner_basis()
  assert f*r + s*g + t*h == 0, "rf + sg + th != 0"
  assert f*rr + g*ss + h*tt == C.defining_polynomial(), "r'f + s'g + t'h != F"

  df = t*ss - s*tt
  dg = r*tt - t*rr
  dh = s*rr - r*ss
  return df, dg, dh
