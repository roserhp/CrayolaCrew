needsPackage "TensorComplexes"

createRing = method() -- n= number of vertices, a= number of vertex colors, b= number of edge colors
createRing(ZZ,List,List) := (n,V,E) -> (
    M := multiSubsets(toList(1..n),2);
    R := QQ[for k in M list s_(toSequence(k)),v_1..v_(#V),e_1..e_(#E)];
    R
    );


createMatrix = method() 
createMatrix(Ring,ZZ) :=(R,n) -> (
    M := mutableMatrix(R,n,n);
    for i from 0 to n-1 do(
    for j from i to n-1 do(
	M_(i,j)=s_(i+1,j+1);
	M_(j,i)=s_(i+1,j+1);
	);
    );
matrix M
);


createSufStats = method()
createSufStats(Ring,List,List) := (R,V,E) -> (
    L := {};
    for i from 0 to (#V-1) do (
	L=append(L,v_(i+1)-sum(for k in V_i list s_(k,k)));
	    );
    for j from 0 to (#E-1) do (
	L=append(L,e_(j+1)-2*sum(for k in E_j list s_(k_0,k_1)));
	    );
   ideal L
   );
   
    
createSufIdeal = method()
createSufIdeal(Ring,ZZ,List,List) := (R,m,V,E) -> (
    M := multiSubsets(toList(1..n),2);
    eliminate(minors(m,createMatrix(R,n))+createSufStats(R,V,E), for k in M list s_(toSequence(k)))
);

	      

-------

--Example 
n=3 -- number of vertices
m=2 -- m\times m minors
V={{1,3},{2}} -- vertex class colors 
E={{{1,2},{2,3}},{{1,3}}} -- edge class colors
R=createRing(n,V,E)
createSufIdeal(R,m,V,E) -- ideal of sufficient statistics for rank m-1 matrices. If this is zero, then MLT=m-1. 


createMatrix(R,n)
createSufStats(R,V,E)
