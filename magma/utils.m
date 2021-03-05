/*
    utils.m

    Useful functions

*/


// Returns values of d such that X_0(d) has genus 1
function Genus1values(n)
S:=[11,14,15,17,19,20,21,24,27,32,36,49];
A:={};
_<x>:=PolynomialRing(Rationals());
for i:=1 to #S do
d:=S[i];
C:=SmallModularCurve(d);
r:=DescentInformation(QuadraticTwist(C,n));
if r[2] ge 1 then A:=A join {d}; end if;
end for;
return A;
end function;

// Code for 65
E := EllipticCurve([1, 0, 0, -1, 0]);
DescentInformation(QuadraticTwist(E,5));


// Code for 125:
R<x> := PolynomialRing(Rationals());
C := HyperellipticCurve(R![0, 0, 1, 2, 2, 2], R![1, 1, 0, 1]);
X,pi:= SimplifiedModel(C);