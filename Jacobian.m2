needsPackage "TensorComplexes"

createRing = method() -- n= number of vertices, a= number of vertex colors, b= number of edge colors
createRing(ZZ,List,List) := (n,V,E) -> (
    M := multiSubsets(toList(1..n),2);
    R := QQ[for k in M list s_(toSequence(k)),v_1..v_(#V),e_1..e_(#E)];
    R
    );
createRingJ = method() -- n= number of vertices, a= number of vertex colors, b= number of edge colors
createRingJ(ZZ,List,List) := (n,V,E) -> (
     L := {};
     for i from 0 to n-1 do(
	 for j from 0 to n-1 do( 
	     L=append(L,x_(i+1,j+1));
	     );
	 );
    R := QQ[L];
    R
    );

createJMatrix = method() 
createJMatrix(Ring,List,List,List) :=(R,L,V,E) -> (
l=L#0;
n=L#1;
    M := mutableMatrix(R,#V+#E,n);
    for i from 0 to #V-1 do(
    for j from 0 to n-1 do(
	for k from 0 to #(V#i)-1 do(
	    if V#i#k==j+1 then M_(i,j)=x_(j+1,l);
	);
	);
    );

 for i from #V to #V+#E-1 do(
    for j from 0 to n-1 do(
	for k from 0 to #(E#(i-#V))-1 do(
	    if E#(i-#V)#k#0==j+1 then M_(i,j)=x_(E#(i-#V)#k#1,l);
	    if E#(i-#V)#k#1==j+1 then M_(i,j)=x_(E#(i-#V)#k#0,l);
	);
	);
    );
matrix M
);

Jdet = method()
Jdet(Ring,ZZ,List,List) :=(R,n,V,E) -> (
    K = createJMatrix(R,{1,n},V,E);
    Kbig=K;
    for i from 0 to n-1 do(

    	Kbig=K*transpose(K);
	if det Kbig==0 then break i+1;
      	K= K|createJMatrix(R,{i+1,n},V,E);
);
n
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
n=4 -- number of vertices
V={{1,2,3,4}} -- vertex class colors 
E={{{1,2}},{{1,3}},{{1,4}},{{2,3}},{{2,4}},{{3,4}}} -- edge class colors
R=createRingJ(n,V,E)
Jdet(R,n,V,E) --calculates the minium number of data points for the jacobian to have non-zero det 





--Example 
n=6 -- number of vertices
m=6 -- m\times m minors
V={{1,2,3,4}} -- vertex class colors 
E={{{1,2}},{{1,3}},{{1,4}},{{2,3}},{{2,4}},{{3,4}}} -- edge class colors
R=createRing(n,V,E)
createSufIdeal(R,m,V,E) -- ideal of sufficient statistics for rank m-1 matrices. If this is zero, then MLT=m-1. 

--Example 
n=4 -- number of vertices
m=2 -- m\times m minors
V={{1},{2},{3},{4}} -- vertex class colors 
E={{{1,2}},{{2,3}},{{3,4}}} -- edge class colors
R=createRing(n,V,E)
createSufIdeal(R,m,V,E) -- ideal of sufficient statistics for rank m-1 matrices. If this is zero, then MLT=m-1. 

--Example 
n=4 -- number of vertices
m=2 -- m\times m minors
V={{1,2},{3},{4}} -- vertex class colors 
E={{{1,2}},{{2,3}},{{3,4}}} -- edge class colors
R=createRing(n,V,E)
createSufIdeal(R,m,V,E) -- ideal of sufficient statistics 

--Example 
n=4 -- number of vertices
m=2 -- m\times m minors
V={{1,2},{3,4}} -- vertex class colors 
E={{{1,2}},{{2,3}},{{3,4}}} -- edge class colors
R=createRing(n,V,E)
createSufIdeal(R,m,V,E) -- ideal of sufficient statistics

--Example 
n=4 -- number of vertices
m=2 -- m\times m minors
V={{1,3},{2},{4}} -- vertex class colors 
E={{{1,2}},{{2,3}},{{3,4}}} -- edge class colors
R=createRing(n,V,E)
createSufIdeal(R,m,V,E) -- ideal of sufficient statistics 
tex oo


--Example 
n=4 -- number of vertices
m=2 -- m\times m minors
V={{1,3},{2,4}} -- vertex class colors 
E={{{1,2}},{{2,3}},{{3,4}}} -- edge class colors
R=createRing(n,V,E)
createSufIdeal(R,m,V,E) -- ideal of sufficient statistics 
tex oo


--Example 
n=4 -- number of vertices
m=2 -- m\times m minors
V={{1,4},{2},{3}} -- vertex class colors 
E={{{1,2}},{{2,3}},{{3,4}}} -- edge class colors
R=createRing(n,V,E)
createSufIdeal(R,m,V,E) -- ideal of sufficient statistics 
tex oo

--Example 
n=4 -- number of vertices
m=2 -- m\times m minors
V={{1,4},{2,3}} -- vertex class colors 
E={{{1,2}},{{2,3}},{{3,4}}} -- edge class colors
R=createRing(n,V,E)
createSufIdeal(R,m,V,E) -- ideal of sufficient statistics 
tex oo

--Example 
n=4 -- number of vertices
m=2 -- m\times m minors
V={{1},{4},{2,3}} -- vertex class colors 
E={{{1,2}},{{2,3}},{{3,4}}} -- edge class colors
R=createRing(n,V,E)
createSufIdeal(R,m,V,E) -- ideal of sufficient statistics 
tex oo


--Example 
n=4 -- number of vertices
m=2 -- m\times m minors
V={{1},{2,3,4}} -- vertex class colors 
E={{{1,2}},{{2,3}},{{3,4}}} -- edge class colors
R=createRing(n,V,E)
createSufIdeal(R,m,V,E) -- ideal of sufficient statistics 
tex oo

--Example 
n=4 -- number of vertices
m=2 -- m\times m minors
V={{1,2,3,4}} -- vertex class colors 
E={{{1,2}},{{2,3}},{{3,4}}} -- edge class colors
R=createRing(n,V,E)
createSufIdeal(R,m,V,E) -- ideal of sufficient statistics 
tex oo


--Example 
n=4 -- number of vertices
m=4 -- m\times m minors
V={{1},{2},{3},{4}} -- vertex class colors 
E={{{1,2}},{{1,3}},{{1,4}},{{2,3}},{{2,4}},{{3,4}}} -- edge class colors
R=createRing(n,V,E)
createSufIdeal(R,m,V,E) -- ideal of sufficient statistics 
tex oo

--Example 
n=4 -- number of vertices
m=4 -- m\times m minors
V={{1,2},{3},{4}} -- vertex class colors 
E={{{1,2}},{{1,3}},{{1,4}},{{2,3}},{{2,4}},{{3,4}}} -- edge class colors
R=createRing(n,V,E)
createSufIdeal(R,m,V,E) -- ideal of sufficient statistics 
tex oo


--Example 
n=4 -- number of vertices
m=4 -- m\times m minors
V={{1,2,3},{4}} -- vertex class colors 
E={{{1,2}},{{1,3}},{{1,4}},{{2,3}},{{2,4}},{{3,4}}} -- edge class colors
R=createRing(n,V,E)
createSufIdeal(R,m,V,E) -- ideal of sufficient statistics 
tex oo

--Example 
n=4 -- number of vertices
m=3 -- m\times m minors
V={{1,2,3,4}} -- vertex class colors 
E={{{1,2}},{{1,3}},{{1,4}},{{2,3}},{{2,4}},{{3,4}}} -- edge class colors
R=createRing(n,V,E)
createSufIdeal(R,m,V,E) -- ideal of sufficient statistics 
tex oo


--Example 
n=4 -- number of vertices
m=4 -- m\times m minors
V={{1,2},{3,4}} -- vertex class colors 
E={{{1,2}},{{1,3}},{{1,4}},{{2,3}},{{2,4}},{{3,4}}} -- edge class colors
R=createRing(n,V,E)
createSufIdeal(R,m,V,E) -- ideal of sufficient statistics 
tex oo







