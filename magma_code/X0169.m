//This file is part of the code for the paper B. S. Banwait, F. Najman, O. Padurariu: Cyclic isogenies of elliptic curves over a fixed quadratic field. It determines all the quadratic points on X0(169)

load "X0p_NiceModel.m";
load "Chabauty_MWSieve_131.m";

Curve_and_Map := function(X,d);
	R := AmbientSpace(X);
	RR<[u]> := CoordinateRing(R);
	n := Dimension(AmbientSpace(X));
	P := ProjectiveSpace(Rationals(), d - 1);
	proj := map<R -> P|[u[i] : i in [1..d]]>;
	Xwd := proj(X);
	mp := map<X -> Xwd|[u[i] : i in [1..d]]>;
	return Xwd, mp;
end function;


//we find models for X and X/w169

C := CuspForms(169);
"Dimension of CuspForms(169) is: ", Dimension(C);

AL := AtkinLehnerOperator(C, 169);
N := Nullspace(Matrix(AL - 1));

"Dimesion of eigenspace lambda = 1 for w169 is: ", Dimension(N);

Nc := Nullspace(Matrix(AL + 1));

"Dimesion of eigenspace lambda = -1 for w169 is: ", Dimension(Nc);
"";

B := [&+[(Integers()!(2*Eltseq(Basis(N)[i])[j]))*C.j : j in [1..Dimension(C)]] : i in [1..Dimension(N)]];
Bc := [&+[(Integers()!(2*Eltseq(Basis(Nc)[i])[j]))*C.j : j in [1..Dimension(C)]] : i in [1..Dimension(Nc)]];

X := modformeqns(Bc cat B, 169, 500, 13);
"Nice model for X0(169) is:";
X;
"";
RR<[u]> := CoordinateRing(AmbientSpace(X));
n := Dimension(AmbientSpace(X));

H := Matrix(RationalField(), 8, 8, [1,0,0,0,0,0,0,0, 0,1,0,0,0,0,0,0, 0,0,1,0,0,0,0,0, 0,0,0,1,0,0,0,0, 0,0,0,0,1,0,0,0, 0,0,0,0,0,-1,0,0, 0,0,0,0,0,0,-1,0, 0,0,0,0,0,0,0,-1]);
rows := [[&+[RowSequence(H)[i][j]*u[j] : j in [1..n+1]] : i in [1..n+1]]];
w := iso<X -> X | rows, rows>;
"w on X is given by:";
w;

Xw, quotMap := Curve_and_Map(X, 5);
RR<[v]> := CoordinateRing(AmbientSpace(Xw));
"Genus of X0(169) is ", Genus(X);
"Genus of X0(169)/w169 is ", Genus(Xw);
"";
ptsXw:=PointSearch(Xw,100);
P1:=X![1,0,0,0,0,1,0,0];
P2:=w(P1);
deg2pb:=[Pullback(quotMap,Place(p)):p in ptsXw];
pls1 := [Place(P1), Place(P2)];
deg2:=[];
for i in [1..#pls1] do
	for j in [i..#pls1] do
		deg2 := Append(deg2, 1*pls1[i] + 1*pls1[j]);
		if w(RepresentativePoint(pls1[i])) eq RepresentativePoint(pls1[j]) then
			deg2pb := Append(deg2pb, 1*pls1[i] + 1*pls1[j]);
		end if;
	end for;
end for;
#deg2;
Dtor:=Divisor(P1)-Divisor(P2);

//We now show that the torsion is Z/7Z

S:={};
for i:=2 to 4 do
d:=NthPrime(i);
S:=S join {#TorsionSubgroup(ClassGroup(ChangeRing(X,GF(d))))};
end for;
S;
GCD(S);
IsPrincipal(Dtor);
IsPrincipal(7*Dtor);
primes:=[3,5];
A:=AbelianGroup([7]);
genusC:=3;
wMatrix:=Matrix(w);
gens:=[Dtor];
B0, iA0 := sub<A | A.1>;
W0 := {0*A.1};
bp:=deg2pb[4];
assert Degree(ResidueClassField(Decomposition(bp)[1,1])) eq 1;

B,iA,W:= MWSieve(X,wMatrix,genusC,primes, A, gens, bp, B0,iA0,W0,deg2);
B;
W;


//It remains to determine all the quadratic points
g:=Genus(X);
for i in [1..#deg2pb] do
	if Degree(ResidueClassField(Decomposition(deg2pb[i])[1,1])) gt 1 then
		K1<z>:=ResidueClassField(Decomposition(deg2pb[i])[1,1]);
		d:=SquarefreeFactorization(Discriminant(K1));
		K<w>:=QuadraticField(d);
		tr,f:=IsIsomorphic(K1,K);
		assert tr;
		P:=RepresentativePoint(Decomposition(deg2pb[i])[1,1]);
		Pm:=ChangeRing(X,K)![f(P[i]): i in [1..g]];
		w^2,Pm;
	end if;
end for;