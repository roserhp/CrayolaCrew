needsPackage "TensorComplexes"
needsPackage "FourTiTwo"
--ring with variables for concentration and covariance matrix
--n = number vertices
--c = number colors
createRing = method()
createRing(ZZ) := (n) -> (
    c :=binomial(n+1,2);
    M := multiSubsets(toList(1..n),2);
    R := QQ[for k in M list s_(toSequence(k)),t_1..t_c,for k in M list p_(toSequence(k)) ];
    R
    );

sourceRing = method()
sourceRing(ZZ) := (n) -> (
    c :=binomial(n+1,2);
    M := multiSubsets(toList(1..n),2);
    R := QQ[for k in M list s_(toSequence(k))];
    R
    )



--creates ideal of inverse linear space
--R should be createRing 
--n = number vertices
--colors is a list containing lists of edges with the same color (including loops)
--each edge should be given as a list so this is really a list of list of lists....sorry
--inverseIdeal = method()
--inverseIdeal(ZZ, Ring, List) := (n, R, colors) -> (
  --  Sig := map(R^n, R^n, {});
  --  M := local M;
  --  k := local k;
  --  i := local i;
  --  l := local l;
  --  shift := {1,1};
  --  for i from 0 to length(colors)-1 do (
--	M = map(R^n, R^n, join(for k in colors_i list toSequence(k-shift)=>t_(i+1),for k in colors_i list toSequence(reverse(k-shift))=>t_(i+1)));
--	Sig = Sig + M;
--	);
  --  N :=  multiSubsets(toList(1..n),2);
  --  N = for l in N list toSequence(l-shift);
  --  L1 := for l in N list l=>s_(toSequence(toList(l) + shift));
  --  L2 := for l in N list reverse(l)=>s_(toSequence(toList(l) + shift));
  --  L := join(L1,L2);
  --  K := map(R^n,R^n,L);
  --  I := minors(1,K*Sig - id_(R^n));
  --  I = eliminate(I, toList(t_1..t_(length(colors))));
  --  I
  -- )

inverseIdeal = method()
inverseIdeal(ZZ,Ring,List) := (n,R,colors) -> (
    M := local M;
    k := local k;
    i := local i;
    l := local l;
    S := local S;
    shift := {1,1};
    S = QQ[for k from 0 to length(colors)-1 list t_k];
    K := map(S^n,S^n,{});
    for i from 0 to length(colors)-1 do(
	M = map(S^n, S^n, join(for k in colors_i list toSequence(k-shift)=>t_(i),for k in colors_i list toSequence(reverse(k-shift))=>t_(i)));
	K = K + M;
	);
    f := map(S,R,for x in  multiSubsets(toList(1..n),2) list determinant(submatrix'(K,{x_0-1}, {x_1-1})));
    kernel f
    )



n = 4
R = sourceRing n
colors = {{{1,1},{3,3}},{{2,2},{4,4}},{{1,3}},{{1,2},{2,3},{3,4},{1,4}}}
inverseIdeal(n,R,colors)


n=4
R = sourceRing(n)
colors = {{{1,1}},{{2,2}},{{3,3}},{{4,4}},{{1,2}},{{2,3}},{{3,4}},{{1,4}}}
I = inverseIdeal(n,R,colors)
mingens I