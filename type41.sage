var("f0,f1,f2,f3,g0,g1,g2,g3,c0,c1,c2,c3,c4,c5,c6,c7,c8")

s0 = f0 - f3*c3
s1 =    - f3*c4
s2 =    - f3*c5
s3 = f1 - f3*c6
s4 = f2 - f3*c7
s5 =    - f3*c8

s6 =    - c5*s3 - c2
s7 = s0 - c6*s3
s8 = s1 - c7*s3 - c3
s9 = s2 - c8*s3 - c4
s10 = -s3 - c5
s11 =  s4 - c6
s12 =  s5 - c7
s13 = -f3 - c8

j0 = f2
j1 = -f3^2

k0 =    - f3*s2
k1 =    - f3*s3
k2 = f1 - f3*s4
k3 = f2 - f3*s5
k4 =      f3*f3

n0 = s6  - f0*s13
n1 = s7
n2 = s8
n3 = s9  - f1*s13
n4 = s10 - f2*s13
n5 = s11
n6 = s12 - f3*s13

r4 = g3 + n6
r5 =      n2 - r4*k2 - n5*s4
r6 = g1 + n3 - r4*k3 - n5*s5

a1 = g2
a3 = -g3*f3

a2 = g0 + n0 - r4*k0 - n5*s2 - r6*f2
a4 =      n1 - r4*k1 - n5*s3 + r6*f3^2 - r5*f3
a5 = g2 + n4 - r4*k4 + n5*f3

print("a4 = {}".format(a4.subs(g3 = -f3^2).expand()))
print("a5 = {}".format(a5.subs(g3 = -f3^2).expand()))
