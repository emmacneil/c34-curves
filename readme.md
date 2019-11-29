# C<sub>3,4</sub> Curves

Last updated 2019-Nov-28.

## Getting Started

To begin, run Sage. At the prompt, load the file c34.sage.
```
sage: load("c34.sage")
```
This file loads all of the other required files for C34 curve arithmetic.



## Constructing C<sub>3,4</sub> Curves

To construct a C34 curve, one must specify its base field and a `list` of coefficients of its defining polynomial. This `list` must be of length 9. The elements of the `list` represent, in order of first to last, the coefficients of 1, x, y, x², xy, y², x³, x²y, and xy².

```
sage: C = C34Curve(GF(97), [1,2,3,4,5,6,7,8,9]); C
C34 curve defined by y^3 + x^4 + 9*x*y^2 + 8*x^2*y + 7*x^3 + 6*y^2 + 5*x*y + 4*x^2 + 3*y + 2*x + 1 over Finite Field of size 97
```

If the given field and coefficients would result in a singular curve, the constructor throws a`ValueError`.

When the field is of a characteristic other than 2 or 3, an invertible change of coordinates puts the curve's defining polynomial in short form, where the coefficients of y², x³, and xy² are all zero. The `short_form` method returns an equivalent curve in short form.
```
sage: C = C34Curve(GF(97), [1,2,3,4,5,6,7,8,9])
sage: C = C.short_form(); C
C34 curve defined by y^3 + x^4 - 19*x^2*y - 19*x*y - 8*x^2 + 4*y + 26*x - 21 over Finite Field of size 97
```
In characteristic 2, the short form has zero for the coefficients of y² and xy². In characteristic 3, the short form of the curve has zero for the coefficient of x³.
```
sage: C = C34Curve(GF(2), [1,1,1,1,1,0,0,1,1])
sage: C.short_form()
C34 curve defined by y^3 + x^4 + x^3 + x*y + y + 1 over Finite Field of size 2
sage: C = C34Curve(GF(3), [1,1,1,2,2,2,1,1,1])
sage: C.short_form()
C34 curve defined by y^3 + x^4 + x*y^2 + x^2*y + y^2 - x^2 - x - 1 over Finite Field of size 3
```

It is possible to generate curves randomly, specifying only the base field, via the static method `random_curve`. The method returns a curve, guaranteed to be non-singular and in short form.
```
sage: C = C34Curve.random_curve(GF(95803)); C
C34 curve defined by y^3 + x^4 - 31498*x^2*y + 12556*x*y - 31969*x^2 + 15891*y + 15938*x + 40981 over Finite Field of size 95803
```
This curve is selected by choosing each coefficient of the curve's defining polynomial according to the `random_element` method of the base field.

Curves may be constructed over finite fields of orders equal to a prime power, e.g.
```
sage: C = C34Curve.random_curve(GF(25)); C
C34 curve defined by y^3 + x^4 + x^2*y + (-2*z2)*x*y - x^2 + (2*z2 - 1)*y + (2*z2)*x + (-z2 - 2) over Finite Field in z2 of size 5^2
```
Curves may also be constructed over infinite fields as well, although functionality has not been tested over infinite fields.
```
sage: C = C34Curve.random_curve(QQ); C
C34 curve defined by y^3 + x^4 + 28*x*y + 607/24*x^2 - 2351/132*y - 8101/264*x + 21078379/2356992 over Rational Field
```



## C34Curve Methods

The `C34Curve` class contains several helpful methods. After constructing a `C34Curve`, one may want to prepare the environment by storing the curve's defining polynomial, base field, and polynomial ring in variables.
```
sage: C = C34Curve(GF(97), [1,2,3,4,5,6,7,8,9]).short_form(); C
C34 curve defined by y^3 + x^4 - 19*x^2*y - 19*x*y - 8*x^2 + 4*y + 26*x - 21 over Finite Field of size 97
sage: F = C.defining_polynomial(); F
y^3 + x^4 - 19*x^2*y - 19*x*y - 8*x^2 + 4*y + 26*x - 21
sage: R = C.polynomial_ring(); R
Multivariate Polynomial Ring in x, y over Finite Field of size 97
sage: K = C.base_field(); K
Finite Field of size 97
```

A list of rational points on the curve (points with coordinates in the base field) may be obtained via the `rational_points` method. The elements of this list are instances of a `C34CurvePoint` class.
```
sage: C = C34Curve(GF(5), [1,2,3,4,0,1,2,3,4]).short_form(); C
C34 curve defined by y^3 + x^4 + x^2*y - 2*x*y + 2*x^2 - 2*y + 2 over Finite Field of size 5
sage: C.rational_points()

[(0 : 1 : 0),
 (1 : 0 : 1),
 (2 : 1 : 1),
 (2 : 2 : 1),
 (4 : 0 : 1),
 (4 : 2 : 1),
 (4 : 3 : 1)]
```
This list always begins with the point at infinity, (0:1:0). All other rational points are sorted by their x-coordinates, then by their y-coordinates. This method requires computing roots `q` polynomials, where `q` is the order of the (finite) base field, and so is very slow over large finite fields. This method will not work over infinite fields.

Rather than computing an entire list of points, one may request a single random rational point via the method `random_rational_point`.
```
sage: C = C34Curve.random_curve(GF(next_prime(10^6))); C
C34 curve defined by y^3 + x^4 + 433720*x^2*y + 236806*x*y + 358651*x^2 + 304614*y - 135338*x + 18781 over Finite Field of size 1000003
sage: C.random_rational_point()
(661666 : 736333 : 1)
sage: P, Q, R = [C.random_rational_point() for i in range(3)]
sage: P;Q;R
(599528 : 714161 : 1)
(631944 : 550594 : 1)
(305327 : 829729 : 1)
```

The point at infinity on the curve may be requested via the `point_at_infinity` method. Other points may be requested via the `point` method.
```
sage: C = C34Curve(GF(5), [1,2,3,4,0,1,2,3,4]).short_form()
sage: P_inf = C.point_at_infinity()
sage: P = C.point(2,2)
sage: P_inf; P
(0 : 1 : 0)
(2 : 2 : 1)
```
If the arguments given to `point` are not coordinates of a rational point on the curve, a `ValueError` is thrown.



## Constructing Divisors

There are several ways to construct divisors on C34 curves.
One way is via the `divisor` method of the `C34Curve` class. This method takes a `list` of `C34CurvePoints` and returns the reduced divisor equivalent to the sum of those points. However, this method assumes that all points in the `list` are rational.
```
sage: C = C34Curve(GF(97), [1,2,3,4,5,6,7,8,9]).short_form()
sage: P, Q, R = [C.random_rational_point() for i in range(3)]
sage: P;Q;R
(35 : 48 : 1)
(90 : 50 : 1)
(12 : 38 : 1)
sage: D = C.divisor([P,Q,R]); D
<x^2 - 4*y + 18*x - 14, x*y + 39*y - 28*x + 47, y^2 - 11*y + 18*x + 19>
sage: Pts = [C.random_rational_point() for i in range(100)]
sage: D = C.divisor(Pts); D
<x^2 - 37*y - 31*x + 4, x*y - 20*y + 42*x - 28, y^2 + 19*y - 20*x + 30>
```
This method is not appropriate for very large lists of points, otherwise it may exceed Python's maximum recursion depth.

Rather than a list of points, one may specify a list of polynomials describing a Gröbner basis for the divisor's ideal. These polynomials must be in the same polynomial ring as the divisor's parent curve.

```
sage: R = C.polynomial_ring()
sage: x, y = R.gens()
sage: f = x^2 - 9*y + 25*x - 16
sage: g = x*y + 28*y + 14*x + 45
sage: h = y^2 + 28*y - 17*x - 44
sage: D = C.divisor([f,g,h]); D
<x^2 - 9*y + 25*x - 16, x*y + 28*y + 14*x + 45, y^2 + 28*y - 17*x - 44>
sage: D.formal_sum()

[((46*z3^2 + 18*z3 + 32 : 67*z3^2 + 64*z3 : 1), 1),
 ((15*z3^2 + 39*z3 + 40 : 22*z3^2 + 89*z3 + 21 : 1), 1),
 ((36*z3^2 + 40*z3 + 69 : 8*z3^2 + 41*z3 + 34 : 1), 1)]
```
It is up to the user to ensure that the polynomials given do in fact form a Gröbner basis. If they do not, any arithmetic done with this divisor is unlikely to be correct. Garbage in, garbage out.

Over a *finite* field, divisors may be produced randomly using the `random_divisor` method of the `C34Curve` class. Unlike the first method above, it may produce a divisor whose support consists of irrational points. Unlike the second method, it will always be correctly given by Gröbner basis.
```
sage: D = C.random_divisor(); D
<x^2 - 47*y + 6*x - 45, x*y - 10*y - 30*x + 11, y^2 - 18*y + 31*x + 39>
sage: D = C.random_divisor(); D
<x^2 - 7*y + 13*x - 33, x*y + 39*y + 16*x - 7, y^2 - 41*y - 7*x + 46>
```
The randomly generated divisor will always be reduced, though may occasionally be atypical. If a divisor of a particular type is needed, one may use the `random_divisor_of_type` method.
```
sage: D = C.random_divisor_of_type(21); D
<y - 36*x - 45, x^2 + 19*x + 44>
sage: D = C.random_divisor_of_type(22); D
<x + 45, y^2 + 21*y - 24>```
```
For semi-typical divisors (divisors of type 31, 41, 51, or 61), one may specify whether the randomly generated divisor should be typical or atypical. If unspecified, it may be either, though will be typical with probability approximately `1/q`, where `q` is the order of the finite field.
```
sage: D = C.random_divisor_of_type(31); D
<x^2 + 38*y + 6*x + 11, x*y + 18*y - 39*x - 18, y^2 + 41*y - 18*x + 6>
sage: D = C.random_divisor_of_type(31, typical = False); D
<x^2 + 9*x + 14, x*y + 7*y + 7*x - 48, y^2 - 47*y + 10*x + 30>
sage: D = C.random_divisor_of_type(31, typical = True); D
<x^2 - 33*y + 18*x - 20, x*y + 9*y - x + 43, y^2 + 5*y - 19*x + 17>
```
If so desired, one may also generate unreduced divisors.
```
sage: D = C.random_divisor_of_type(44); D
<y + 22*x - 18>
sage: D = C.random_divisor_of_type(53); D
<x*y + 18*x^2 - 13*y - 47*x - 24, x^3 - 24*y^2 - 39*x^2 - 21*y - 46*x - 46>
sage: D = C.random_divisor_of_type(61, typical = False); D
<x^3 - 4*y^2 + 29*x^2 + 9*y + 31*x + 42, x^2*y + 16*y^2 + 46*x*y + 40*x^2 + 26*y + 37*x - 32, x*y^2 - 40*y^2 + 48*x*y + 40*x^2 - 25*y - 45*x + 42>
```



## Divisor Arithmetic