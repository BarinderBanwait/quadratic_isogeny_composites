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
C := HyperellipticCurve(x^6 + 2*x^5 - x^4 - 8*x^3 - x^2 + 2*x +1);  // X0(91)+
J := Jacobian(C);

load "ozmansiksek.m";

ModCrvQuot:=function(N,wlist,remwlist)
C:=CuspForms(N);
ws:=[AtkinLehnerOperator(C,n) : n in wlist];
remws:=[AtkinLehnerOperator(C,n) : n in remwlist];
W:=AtkinLehnerOperator(C,N);
NN:=&meet([Nullspace(Matrix(w-1)) : w in ws] cat [Nullspace(1-W^2)]);
dim:=Dimension(NN);
seqs:=[[Coordinates(NN,Basis(NN)[i]*Matrix(w)) : i in [1..dim]] : w in remws];
BB:=[&+[(Integers()!(2*Eltseq(Basis(NN)[i])[j]))*C.j : j in [1..Dimension(C)]] : i in [1..dim]];
prec:=500;
L<q>:=LaurentSeriesRing(Rationals(),prec);
R<[x]>:=PolynomialRing(Rationals(),dim);
Bexp:=[L!qExpansion(BB[i],prec) : i in [1..dim]];
eqns:=[R | ];
	d:=1;
	tf:=false;
	while tf eq false do
		d:=d+1;
		mons:=MonomialsOfDegree(R,d);
		monsq:=[Evaluate(mon,Bexp) : mon in mons];
		V:=VectorSpace(Rationals(),#mons);
		W:=VectorSpace(Rationals(),prec-10);
		h:=hom<V->W | [W![Coefficient(monsq[i],j) : j in [1..(prec-10)]] : i in [1..#mons]]>;
		K:=Kernel(h);
		eqns:=eqns cat [ &+[Eltseq(V!k)[j]*mons[j] : j in [1..#mons] ] : k in Basis(K)  ];
        I:=Radical(ideal<R | eqns>);
		X:=Scheme(ProjectiveSpace(R),I);
		if Dimension(X) eq 1 then
			if IsSingular(X) eq false then
				X:=Curve(ProjectiveSpace(R),eqns);
				if Genus(X) eq dim then
					tf:=true;
				end if;
			end if;
		end if;
	end while;
	eqns:=GroebnerBasis(ideal<R | eqns>); // Simplifying the equations.
	tf:=true;
	repeat
		t:=#eqns;
		tf:=(eqns[t] in ideal<R | eqns[1..(t-1)]>);
		if tf then
			Exclude(~eqns,eqns[t]);
		end if;
	until tf eq false;
	t:=0;
	repeat
		t:=t+1;
		tf:=(eqns[t] in ideal<R | Exclude(eqns,eqns[t])>);
		if tf then
			Exclude(~eqns,eqns[t]);
			t:=0;
		end if;
	until tf eq false and t eq #eqns;
	X:=Curve(ProjectiveSpace(R),eqns); // Our model for X_0(N) discovered via the canonical embedding.
	assert Genus(X) eq dim;
	assert IsSingular(X) eq false;
    indexGam:=N*&*[Rationals() | 1+1/p : p in PrimeDivisors(N)];
	indexGam:=Integers()!indexGam; // Index of Gamma_0(N) in SL_2(Z)
	for eqn in eqns do
		eqnScaled:=LCM([Denominator(c) : c in Coefficients(eqn)])*eqn;
		wt:=2*Degree(eqn); // Weight of eqn as a cuspform.
		hecke:=Ceiling(indexGam*wt/12);  // Hecke=Sturm bound.
										// See Stein's book, Thm 9.18.
		Bexp1:=[qExpansion(BB[i],hecke+10) : i in [1..dim]]; // q-expansions
                        // of basis for S
                        // up to precision hecke+10.
		assert Valuation(Evaluate(eqnScaled,Bexp1)) gt hecke+1;
	end for; // We have now checked the correctness of the equations for X.
seqlist:=[[&+[seq[i][j]*x[j] : j in [1..dim]] : i in [1..dim]] : seq in seqs];
wmaplist:=[iso<X->X | seq,seq> : seq in seqlist];
return X,wmaplist,seqs;
end function;


X,ws:=ModCrvQuot(91,[],[91]); //Just a few seconds.
w:=ws[1];
assert Genus(X) eq 7;

//This function computes J_X(F_p).
JacobianFp:=function(X)
CC,phii,psii:=ClassGroup(X); //Algorithm of Hess
Z:=FreeAbelianGroup(1);
degr:=hom<CC->Z | [ Degree(phii(a))*Z.1 : a in OrderedGenerators(CC)]>;
JFp:=Kernel(degr); // This is isomorphic to J_X(\F_p).
return JFp,phii,psii;
end function;

cusps:=PointSearch(X,500);
numcusps := #cusps;
Dtors:=[Place(cusps[i])-Place(cusps[1]) : i in [2..4]];

