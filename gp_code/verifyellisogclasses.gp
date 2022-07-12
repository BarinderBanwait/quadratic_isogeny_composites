/*
This GP code can be used to verify the computation of isogeny classes of
elliptic curves over number fields, in case the solver reports that
it failed in Sage, which it does sometimes.

*/
D = 9066;  \\ change to desired value
K = nfinit(a ^ 2 - D);
myJ = Mod(54000, a ^ 2 - D);  \\ change to desired value
v = ellfromj(myJ);
E = ellinit(v, K);
[L, M] = ellisomat(E, 1);  // inspect M for any unrecorded isogeny degrees

/*
Here is the same code in a convenient one-liner, in case your version of
gp also suffers from a readline issue where parentheses get stacked at the end.

D = 4569;K = nfinit(a ^ 2 - D);myJ = Mod(-140625/8, a ^ 2 - D);v = ellfromj(myJ);E = ellinit(v, K);[L, M] = ellisomat(E, 1);M

*/