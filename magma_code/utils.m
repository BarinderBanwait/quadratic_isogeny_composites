/*
    utils.m

    Useful functions

    ====================================================================

    This file is part of Quadratic Kenku Solver.

    Copyright (C) 2022 Barinder S. Banwait, Filip Najman, and Oana
    Padurariu

    Quadratic Kenku Solver is free software: you can redistribute it
    and/or modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation, either version 3 of
    the License, or any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <https://www.gnu.org/licenses/>.

    ====================================================================

*/

// Globals

Genus1Values:=[11,14,15,17,19,20,21,24,27,32,36,49];
R<x> := PolynomialRing(Rationals());


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

procedure EllipticJInvs(d)
	f := R![-d,0,1];
	K<a> := NumberField(f);
	for N in Genus1Values do
		X0N := SmallModularCurve(N);
		X0Ntw := QuadraticTwist(X0N,d);
		if Rank(X0Ntw) eq 0 then
			j_invs := {};
			X0NK := BaseExtend(X0N,K);
			Tors, m := TorsionSubgroup(X0NK);
			num_ers:=0;
			for P in Tors do
				if P ne Tors!0 then
					try
						a_j_inv := jInvariant(m(P),N);
						j_invs := j_invs join {a_j_inv};
					catch e
						// means the point is a cusp
						num_ers:=num_ers+1;
					end try;
				end if;
			end for;
			print N, j_invs, num_ers, #Cusps(Gamma0(N));
		end if;
	end for;
end procedure;

function EllipticCount(d)
	f := R![-d,0,1];
	K<a> := NumberField(f);
	finalDict := [];
	for N in Genus1Values do
		X0N := SmallModularCurve(N);
		X0Ntw := QuadraticTwist(X0N,d);
		if Rank(X0Ntw) eq 0 then
			j_invs := {};
			X0NK := BaseExtend(X0N,K);
			Tors, m := TorsionSubgroup(X0NK);
			num_ers:=0;
			for P in Tors do
				if P ne Tors!0 then
					try
						a_j_inv := jInvariant(m(P),N);
						j_invs := j_invs join {a_j_inv};
					catch e
						// means the point is a cusp
						num_ers:=num_ers+1;
					end try;
				end if;
			end for;
			// print N, j_invs, num_ers, #Cusps(Gamma0(N));
			if #j_invs ne 0 then
				Append(~finalDict, <N,#j_invs>);
			end if;
		end if;
	end for;
	return finalDict;
end function;

procedure CheckTorsionGrowth(d)
	f := R![-d,0,1];
	K<a> := NumberField(f);
	for N in Genus1Values do
		X0N := SmallModularCurve(N);
		X0NK := BaseExtend(X0N,K);
		Tors, m := TorsionSubgroup(X0N);
		TorsK, mK := TorsionSubgroup(X0NK);
		print N, #Tors, #TorsK;
	end for;
end procedure;


function Ranks(d)
	f := R![-d,0,1];
	K<a> := NumberField(f);
	MyVals := [];
	for N in Genus1Values do
		X0N := SmallModularCurve(N);
		X0Ntw := QuadraticTwist(X0N,d);
		if AnalyticRank(X0Ntw) eq 0 then
			MyVals := Append(MyVals, N);
		end if;
	end for;
	return MyVals;
end function;

// Code to generate rank data for subsequent sage computation
for d in [-10000..10000] do
	if d ne 1 then
		if d ne 0 then
			if IsSquarefree(d) then
				MyVals:= Ranks(d);
				Write("RankData.txt", Sprintf("%o: %o", d, MyVals));
			end if;
		end if;
	end if;
end for;

