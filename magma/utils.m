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

// Code for 50:
C:=SmallModularCurve(50);
s:=AtkinLehnerInvolution(C,50,2);
G:=AutomorphismGroup(C,[s]);
E,f:=CurveQuotient(G);
f;
Conductor(E);
jInvariant(E);
Degree(f);
DescentInformation(E);
g,m:=TorsionSubgroup(E);
P:=m(g.1);
P;
K<w>:=QuadraticField(5);
A:=ChangeRing(P@@f, QuadraticField(5));
Points(A);
A:=ChangeRing((2*P)@@f, QuadraticField(5));
Points(A);
A:=ChangeRing((3*P)@@f, QuadraticField(5));
Points(A);
C1:=ChangeRing(C,K);
Points(C1,-1);
j1:=9845745509376*w + 22015749613248;
j2:=-9845745509376*w + 22015749613248;
j1:=jInvariant(A[1],50);
j2:=jInvariant(A[2],50);
f:=ClassicalModularPolynomial(25);
Evaluate(f,[j1,j1]);
f:=ClassicalModularPolynomial(2);
Evaluate(f,[j1,j2]);
f:=ClassicalModularPolynomial(5);
Evaluate(f,[1728,j2]);
Evaluate(f,[1728,j1]);

// Code for 91:
R<x> := PolynomialRing(Rationals());
