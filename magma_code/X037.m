/*
    X037.m

    The function `ComputePreimages` below carries out the computation of pulling back
    the rational points on E(K) = E(Q) to determine that X0(37)(K) = X0(37)(Q),
    and hence dealing with the (N=37, K=Qsqrt213) case of Proposition 6.1.

*/

ComputePreimages := function(C, phi)

G := AutomorphismGroup(C,[phi]);
E,m := CurveQuotient(G);
assert Type(E) eq CrvEll;
assert RankBound(E) eq 0;

GG, map := MordellWeilGroup(E);

S := [map(g) : g in GG];

preimages := {};

for P in S do
    preimageofP := P @@ m;
    preimages := preimages join Points(preimageofP);
 end for;


return preimages;
end function;

_<x> := PolynomialRing(Rationals());
K := QuadraticField(213);
_<x> := PolynomialRing(K);
C := HyperellipticCurve(-x^6-9*x^4-11*x^2+37);
F<x,y> := FunctionField(C);
phi := iso< C -> C | [-x,-y,1],[-x,-y,1]>;

ComputePreimages(C,phi);