//B. Vukorepa: This code was written by Ozman and Siksek and was also used by Box.
//We just adapted it so that we can select our own basis for cuspforms as described in the article

//B is a basis for Cuspforms(N)
//N is N for which we want the model for X0(N)
//divN is a divisor of N, used only for jMap, divN = 1 is allowed

modformeqns:=function(B,N,prec,divN)
dim:=#B;
L<q>:=LaurentSeriesRing(Rationals(),prec);
R<[x]>:=PolynomialRing(Rationals(),dim);
Bexp:=[L!qExpansion(B[i],prec) : i in [1..dim]];
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
			if IsIrreducible(X) then
				X:=Curve(ProjectiveSpace(R),eqns);
				if Genus(X) eq dim then
					tf:=true;
				end if;
			end if;
		end if;
	end while;

	//We commented out this part because it is slow and only potentially simplifies the equations
	/*eqns:=GroebnerBasis(ideal<R | eqns>); // Simplifying the equations.
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
	until tf eq false and t eq #eqns;*/

	X:=Curve(ProjectiveSpace(R),eqns); // Our model for X_0(N) discovered via the canonical embedding.
	assert Genus(X) eq dim;
    	indexGam:=N*&*[Rationals() | 1+1/p : p in PrimeDivisors(N)];
	indexGam:=Integers()!indexGam; // Index of Gamma_0(N) in SL_2(Z)

	for eqn in eqns do
		eqnScaled:=LCM([Denominator(c) : c in Coefficients(eqn)])*eqn;
		wt:=2*Degree(eqn); // Weight of eqn as a cuspform.
		hecke:=Ceiling(indexGam*wt/12);  // Hecke=Sturm bound.
						 // See Stein's book, Thm 9.18.
		Bexp1:=[qExpansion(B[i],hecke+10) : i in [1..dim]]; // q-expansions
                        					    // of basis for S
                        					    // up to precision hecke+10.
		assert Valuation(Evaluate(eqnScaled,Bexp1)) gt hecke+1;
	end for; // We have now checked the correctness of the equations for X.

	return X;
end function;


// The function takes a curve X and computes a quotient X/w, where w acts like -1 on the first d coordinates and trivially on the rest

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