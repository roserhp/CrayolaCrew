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








needsPackage "SchurRings"
n=3
--number of vertices
m=2
--m\times m minors
Vertices={1,2,3}
Edges= {{1,2},{1,3},{2,3}};
--list of edges

permutationVertices= permutations Vertices;
permutationEdges= permutations Edges;
--- creates list of every order of the vertices for order in which vertices will be asigned colors (if all vertices are interchangable this doesnt matter)

IdealList= {};
--FinalList of suff Ideals
ColorList = {};
--List of the colorings producing those in the same order, its listed vertex then edge
PartV =partitions(n);

PartE = partitions(#Edges);
--List of all partitions of the size of vertices and edges


--for t from 0 to (#permutationVertices-1) do(
--Vertices= permutationVertices_t;

--- creates list of every order of the vertices for order in which vertices will be asigned colors (if all vertices are interchangable this doesnt matter)
-- comment out this for loop if you want just want one ordering

--for s from 0 to (#permutationEdges-1) do(
--Edges= permutationEdges_s;


--- creates list of every order of the edges for order in which vertices will be asigned colors (if all vertices are interchangable this doesnt matter)
-- comment out this for loop if you want just want one ordering


VertList= {};
EdgeList= {};

VertListfull = {};
EdgeListfull= {};
for p from 0 to (#PartV-1) do(
    V= toList PartV#p;
    curtotalvertex=0;
    
    for k from 0 to (#V-1) do(
	for i from    curtotalvertex to (curtotalvertex+V_k-1) do(   
	VertList= append (VertList, Vertices#i);
    );
   curtotalvertex=   curtotalvertex +(V_k);
   VertListfull= append (VertListfull, VertList);
   VertList= {};   
   );
  --For each partition adds vertices and edges into color classes in the order they are in

for q from 0 to (#PartE-1) do(
curtotaledge=0;


E= toList PartE#q;


    
for l from 0 to (#E-1) do( 
       

for j from  curtotaledge to (curtotaledge+E_l-1) do(  
    EdgeList= append (EdgeList, Edges#j); 
    );
  curtotaledge=   curtotaledge +(E_l);

EdgeListfull= append (EdgeListfull, EdgeList);
EdgeList= {};
);




R=createRing(n,VertListfull,EdgeListfull);
I=createSufIdeal(R,m,VertListfull,EdgeListfull);
IdealList= append (IdealList, I);
ColorList = append (ColorList,VertListfull);
ColorList = append (ColorList,EdgeListfull);
EdgeListfull= {};
);
VertListfull = {};
);



--);

--);







