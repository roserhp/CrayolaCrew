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




R=QQ[t_1,t_2,t_3,s_(1,1),s_(1,1),s_(1,2),s_(2,2),s_(1,3),s_(2,3),s_(3,3)]
M=mutableMatrix(R,3,3)
for i from 0 to 2 do(
    for j from i to 2 do(
	M_(i,j)=s_(i+1,j+1);
	M_(j,i)=s_(i+1,j+1);
	)
    )
--M_(1,2)=0;
--M_(2,1)=0;
S=matrix M
detI=minors(2,S)
I=ideal(t_1-s_(1,1), t_2-2*s_(1,2)-2*s_(1,3),t_3-s_(2,2)-s_(3,3))
eliminate(I+detI,{s_(1,1),s_(1,1),s_(1,2),s_(2,2),s_(1,3),s_(2,3),s_(3,3)})	


--cycle 2 colors
R=QQ[t_1,t_2,s_(1,1),s_(1,1),s_(1,2),s_(2,2),s_(1,3),s_(2,3),s_(3,3)]
M=mutableMatrix(R,3,3)
for i from 0 to 2 do(
    for j from i to 2 do(
	M_(i,j)=s_(i+1,j+1);
	M_(j,i)=s_(i+1,j+1);
	)
    )
--M_(1,2)=0;
--M_(2,1)=0;
S=matrix M
detI=minors(2,S)
I=ideal(t_1-s_(1,1)-s_(2,2)-s_(3,3), t_2-2*s_(1,2)-2*s_(1,3)-2*s_(2,3))
eliminate(I+detI,{s_(1,1),s_(1,1),s_(1,2),s_(2,2),s_(1,3),s_(2,3),s_(3,3)})	


--cycle 3 colors
R=QQ[t_1,t_2,t_3,t_4,s_(1,1),s_(1,1),s_(1,2),s_(2,2),s_(1,3),s_(2,3),s_(3,3)]
M=mutableMatrix(R,3,3)
for i from 0 to 2 do(
    for j from i to 2 do(
	M_(i,j)=s_(i+1,j+1);
	M_(j,i)=s_(i+1,j+1);
	)
    )
--M_(1,2)=0;
--M_(2,1)=0;
S=matrix M
detI=minors(3,S)
I=ideal(t_1-s_(1,1),t_2-s_(2,2)-s_(3,3), t_3-2*s_(1,2)-2*s_(1,3),t_4-2*s_(2,3))
eliminate(I+detI,{s_(1,1),s_(1,1),s_(1,2),s_(2,2),s_(1,3),s_(2,3),s_(3,3)})	



----- 


R=QQ[t_1,t_2,t_3,t_4,t_5,t_6,t_7,s_(1,1),s_(1,2),s_(2,2),s_(1,3),s_(2,3),s_(3,3)]
M=mutableMatrix(R,3,3)
for i from 0 to 2 do(
    for j from i to 2 do(
	M_(i,j)=s_(i+1,j+1);
	M_(j,i)=s_(i+1,j+1);
	)
    )
--M_(1,2)=0;
--M_(2,1)=0;
S=matrix M
detI=minors(3,S)
I=ideal(t_1-s_(1,1),t_2-s_(2,2)-s_(3,3), t_3-2*s_(1,2)-2*s_(1,3),t_4-2*s_(2,3))
eliminate(I+detI,{s_(1,1),s_(1,1),s_(1,2),s_(2,2),s_(1,3),s_(2,3),s_(3,3)})	


R=QQ[t_1,t_2,t_3,t_4,t_5,t_6,t_7,t_8,t_9,t_{10},t_{11},t_{12},t_{13},s_(1,1),s_(1,2),s_(2,2),s_(1,3),s_(2,3),s_(3,3),s_(1,4),s_(2,4),s_(3,4),s_(4,4),s_(1,5),s_(2,5),s_(3,5),s_(4,5),s_(5,5),s_(1,6),s_(2,6),s_(3,6),s_(4,6),s_(5,6),s_(6,6),s_(1,7),s_(2,7),s_(3,7),s_(4,7),s_(5,7),s_(6,7),s_(7,7)]
M=mutableMatrix(R,7,7)
for i from 0 to 6 do(
    for j from i to 6 do(
	M_(i,j)=s_(i+1,j+1);
	M_(j,i)=s_(i+1,j+1);
	)
    )
S=matrix M
detI=minors(2,S)



I=ideal(t_1-s_(1,1),t_2-s_(2,2)-s_(3,3), t_3-2*s_(1,7)-2*s_(2,7)-2*s_(6,7)-2*s_(5,6)-2*s_(4,6)-2*s_(3,6),t_4-s_(4,4)-s_(5,5),t_5-s_(6,6)-s_(7,7))
eliminate(I+detI,{s_(1,1),s_(1,2),s_(2,2),s_(1,3),s_(2,3),s_(3,3),s_(1,4),s_(2,4),s_(3,4),s_(4,4),s_(1,5),s_(2,5),s_(3,5),s_(4,5),s_(5,5),s_(1,6),s_(2,6),s_(3,6),s_(4,6),s_(5,6),s_(6,6),s_(1,7),s_(2,7),s_(3,7),s_(4,7),s_(5,7),s_(6,7),s_(7,7)})
elim = 0




I=ideal(t_1-s_(1,1)-s_(3,3)-s_(4,4),t_2-s_(2,2)-s_(5,5)-s_(7,7), t_4-s_(6,6), t_3-2*s_(1,7), t_5-2*s_(3,6)-2*s_(4,6), t_6-2*s_(2,7),t_7-2*s_(5,6)-2*s_(6,7))
eliminate(I+detI,{s_(1,1),s_(1,2),s_(2,2),s_(1,3),s_(2,3),s_(3,3),s_(1,4),s_(2,4),s_(3,4),s_(4,4),s_(1,5),s_(2,5),s_(3,5),s_(4,5),s_(5,5),s_(1,6),s_(2,6),s_(3,6),s_(4,6),s_(5,6),s_(6,6),s_(1,7),s_(2,7),s_(3,7),s_(4,7),s_(5,7),s_(6,7),s_(7,7)})
elim = 0

I=ideal(t_1-s_(1,1)-s_(3,3)-s_(4,4),t_2-s_(2,2)-s_(5,5)-s_(7,7), t_4-s_(6,6), t_5-2*s_(3,6)-2*s_(4,6)-2*s_(1,7),t_7-2*s_(5,6)-2*s_(6,7)-2*s_(2,7))
eliminate(I+detI,{s_(1,1),s_(1,2),s_(2,2),s_(1,3),s_(2,3),s_(3,3),s_(1,4),s_(2,4),s_(3,4),s_(4,4),s_(1,5),s_(2,5),s_(3,5),s_(4,5),s_(5,5),s_(1,6),s_(2,6),s_(3,6),s_(4,6),s_(5,6),s_(6,6),s_(1,7),s_(2,7),s_(3,7),s_(4,7),s_(5,7),s_(6,7),s_(7,7)})
elim = 0




%%%%%%%%
I=ideal(t_1-s_(1,1),t_2-s_(2,2), t_4-s_(3,3), t_3-2*s_(1,7)-2*s_(2,7)-2*s_(6,7)-2*s_(5,6)-2*s_(4,6)-2*s_(3,6),t_5-s_(4,4), t_6-s_(5,5),t_7-s_(6,6), t_8-s_(7,7))
eliminate(I+detI,{s_(1,1),s_(1,2),s_(2,2),s_(1,3),s_(2,3),s_(3,3),s_(1,4),s_(2,4),s_(3,4),s_(4,4),s_(1,5),s_(2,5),s_(3,5),s_(4,5),s_(5,5),s_(1,6),s_(2,6),s_(3,6),s_(4,6),s_(5,6),s_(6,6),s_(1,7),s_(2,7),s_(3,7),s_(4,7),s_(5,7),s_(6,7),s_(7,7)})
elim = too long

I=ideal(t_1-s_(1,1),t_2-s_(2,2)-s_(3,3), t_3-2*s_(1,7)-2*s_(2,7), t_6-2*s_(6,7)-2*s_(5,6)-2*s_(4,6),t_7-2*s_(3,6),t_4-s_(4,4)-s_(5,5),t_5-s_(6,6)-s_(7,7))
eliminate(I+detI,{s_(1,1),s_(1,2),s_(2,2),s_(1,3),s_(2,3),s_(3,3),s_(1,4),s_(2,4),s_(3,4),s_(4,4),s_(1,5),s_(2,5),s_(3,5),s_(4,5),s_(5,5),s_(1,6),s_(2,6),s_(3,6),s_(4,6),s_(5,6),s_(6,6),s_(1,7),s_(2,7),s_(3,7),s_(4,7),s_(5,7),s_(6,7),s_(7,7)})
elim = 

I=ideal(t_1-s_(1,1)-s_(2,2)-s_(5,5),t_2-s_(7,7)-s_(3,3),t_3-2*s_(1,7)-2*s_(6,7),t_6-2*s_(5,6)-2*s_(4,6),t_7-2*s_(3,6)-2*s_(2,7),t_4-s_(4,4),t_5-s_(6,6))
eliminate(I+detI,{s_(1,1),s_(1,2),s_(2,2),s_(1,3),s_(2,3),s_(3,3),s_(1,4),s_(2,4),s_(3,4),s_(4,4),s_(1,5),s_(2,5),s_(3,5),s_(4,5),s_(5,5),s_(1,6),s_(2,6),s_(3,6),s_(4,6),s_(5,6),s_(6,6),s_(1,7),s_(2,7),s_(3,7),s_(4,7),s_(5,7),s_(6,7),s_(7,7)})
elim =

I=ideal(t_1-s_(1,1),t_2-s_(2,2), t_6-s_(3,3), t_3-2*s_(1,7), t_7-2*s_(2,7),t_8-2*s_(6,7),t_9-2*s_(5,6),t_{10}-2*s_(4,6),t_{11}-2*s_(3,6),t_4-s_(4,4),t_{12}-s_(5,5),t_5-s_(6,6),t_{13}-s_(7,7))
eliminate(I+detI,{s_(1,1),s_(1,2),s_(2,2),s_(1,3),s_(2,3),s_(3,3),s_(1,4),s_(2,4),s_(3,4),s_(4,4),s_(1,5),s_(2,5),s_(3,5),s_(4,5),s_(5,5),s_(1,6),s_(2,6),s_(3,6),s_(4,6),s_(5,6),s_(6,6),s_(1,7),s_(2,7),s_(3,7),s_(4,7),s_(5,7),s_(6,7),s_(7,7)})
elim= insanely long

I=ideal(t_1-s_(1,1),t_2-s_(2,2), t_6-s_(3,3), t_3-2*s_(1,7), t_7-2*s_(2,7),t_8-2*s_(6,7),t_9-2*s_(5,6),t_{10}-2*s_(4,6),t_{11}-2*s_(3,6),t_4-s_(4,4),t_{12}-s_(5,5),t_5-s_(6,6)-s_(7,7))
eliminate(I+detI,{s_(1,1),s_(1,2),s_(2,2),s_(1,3),s_(2,3),s_(3,3),s_(1,4),s_(2,4),s_(3,4),s_(4,4),s_(1,5),s_(2,5),s_(3,5),s_(4,5),s_(5,5),s_(1,6),s_(2,6),s_(3,6),s_(4,6),s_(5,6),s_(6,6),s_(1,7),s_(2,7),s_(3,7),s_(4,7),s_(5,7),s_(6,7),s_(7,7)})
elim=  

I=ideal(t_1-s_(1,1),t_2-s_(2,2), t_6-s_(3,3), t_3-2*s_(1,7), t_7-2*s_(2,7),t_8-2*s_(6,7),t_9-2*s_(5,6),t_{10}-2*s_(4,6),t_{11}-2*s_(3,6),t_4-s_(4,4)-s_(5,5),t_5-s_(6,6),t_{13}-s_(7,7))
eliminate(I+detI,{s_(1,1),s_(1,2),s_(2,2),s_(1,3),s_(2,3),s_(3,3),s_(1,4),s_(2,4),s_(3,4),s_(4,4),s_(1,5),s_(2,5),s_(3,5),s_(4,5),s_(5,5),s_(1,6),s_(2,6),s_(3,6),s_(4,6),s_(5,6),s_(6,6),s_(1,7),s_(2,7),s_(3,7),s_(4,7),s_(5,7),s_(6,7),s_(7,7)})
elim = 

I=ideal(t_1-s_(1,1),t_2-s_(2,2), t_6-s_(3,3), t_3-2*s_(1,7)-2*s_(3,6), t_7-2*s_(2,7),t_8-2*s_(6,7),t_9-2*s_(5,6),t_{10}-2*s_(4,6),t_4-s_(4,4),t_{12}-s_(5,5),t_5-s_(6,6),t_{13}-s_(7,7))
eliminate(I+detI,{s_(1,1),s_(1,2),s_(2,2),s_(1,3),s_(2,3),s_(3,3),s_(1,4),s_(2,4),s_(3,4),s_(4,4),s_(1,5),s_(2,5),s_(3,5),s_(4,5),s_(5,5),s_(1,6),s_(2,6),s_(3,6),s_(4,6),s_(5,6),s_(6,6),s_(1,7),s_(2,7),s_(3,7),s_(4,7),s_(5,7),s_(6,7),s_(7,7)})
elim= 


time for paths%%%%%


R=QQ[t_1,t_2,t_3,t_4,t_5,s_(1,1),s_(1,2),s_(2,2),s_(1,3),s_(2,3),s_(3,3)]
M=mutableMatrix(R,3,3)
for i from 0 to 2 do(
    for j from i to 2 do(
	M_(i,j)=s_(i+1,j+1);
	M_(j,i)=s_(i+1,j+1);
	)
    )
S=matrix M
detI=minors(2,S)
I=ideal(t_1-s_(1,1),t_2-s_(2,2), t_3-s_(3,3), t_4-2*s_(1,2), t_5-2*s_(2,3))
eliminate(I+detI,{s_(1,1),s_(1,2),s_(2,2),s_(1,3),s_(2,3),s_(3,3)})

needsPackage "SchurRings"
P= partitions(set{s_(1,1),s_(2,2),s_(3,3)},{2,1})





%%%% Here is the path Calulator, n is the number of vertices %%%%
needsPackage "SchurRings"
n=6
d=2*n-1
NonZeroGuys= new MutableList; 
ZeroGuys= new MutableList; 
R=QQ[t_1..t_d,s_(1,1)..s_(n,n)]


PartE =partitions(n-1);
PartV = partitions(n);

SameS=new MutableList;
for r from 0 to n-1 do(
	SameS#r=s_(r+1,r+1);
    )
DiffS=new MutableList;
for t from 0 to n-2 do(
	DiffS#t=s_(t+1,t+2);
    )

Ilist= new MutableList; 

L= new MutableList;
for p from 0 to (#PartV-1) do(
    
    
for q from 0 to (#PartE-1) do(



V= partitions(set toList SameS,PartV#p);
E= partitions(set toList DiffS,PartE#q);
for k from 0 to (#V-1) do(
for l from 0 to (#E-1) do(    
for i from 0 to ((#toList V#k)-1) do(   
	L#i=t_(i+1)- sum(toList(toList V#k)#i);
    );
for j from  (#toList V#k) to ((#toList V#k +#toList E#l)-1) do(   
	L#j=t_(j+1)- 2*sum(toList(toList E#l)#(j-#toList V#k));
    );
Ilist= append (Ilist, ideal(toList L));
L= new MutableList;
);
)

);

)




FinalList = toList Ilist
#FinalList
M=mutableMatrix(R,n,n)
for i from 0 to n-1 do(
    for j from i to n-1 do(
	M_(i,j)=s_(i+1,j+1);
	M_(j,i)=s_(i+1,j+1);
	)
    )
S=matrix M
detI=minors(2,S)

For for i from 0 to (#FinalList-1) do(
    I=FinalList#i;
ElimI =eliminate(I+detI,toList(s_(1,1)..s_(n,n)));
if(ElimI == ideal(0*s_(1,1))) then (ZeroGuys= append (ZeroGuys, I)) else  NonZeroGuys= append (NonZeroGuys, I)
)
Z=toList (ZeroGuys);
NZ=toList (NonZeroGuys);
#Z
#NZ
