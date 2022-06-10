/*
    NonHyperellipticVerifs.m

    This carries out the conputations in Section 7.5, regarding N = 163.

*/


// Here we define X0^+(163)
X := X0NQuotient(163,[163]);
K := QuadraticField(213);


D := Decomposition(JZero(163));
J := JZero(163);
w163 := AtkinLehnerOperator(J,163);
W163 := Kernel(-Matrix(w163)+1);
w := Kernel(Matrix(w163)+1);


_,Jplus :=DefinesAbelianSubvariety(J,W163);
TorsionMultiple(Jplus);
Decomposition(Jplus);

// This shows that the torsion does not grow
TorsionMultiple(ChangeRing(Jplus,K)); //1 1


// Checking that the rank does not grow was already done as part of the
// sage `search_convenient_d` function; here we sanity check this in Magma

p := 163;
d := 213;
M := ModularSymbols(p);
S := CuspidalSubspace(M);
S_plus := AtkinLehnerDecomposition(S)[1];
assert Dimension(S_plus) eq 12;  // checks we got the plus part
chi := BaseExtend(KroneckerCharacter(d),RationalField());
tw := TwistedWindingElement(S_plus, 1, chi);
twmap := RationalMapping(S_plus)(tw);
twmap eq Parent(twmap)!0; // false, so rank of twist is zero
