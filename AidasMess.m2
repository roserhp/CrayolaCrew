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
detI=minors(3,S)
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
detI=minors(2,S)
I=ideal(t_1-s_(1,1),t_2-s_(2,2)-s_(3,3), t_3-2*s_(1,2)-2*s_(1,3),t_4-2*s_(2,3))
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
detI=minors(2,S)
I=ideal(t_1-s_(1,1)-s_(2,2)-s_(3,3),t_2-2*s_(1,2),t_3-2*s_(1,3),t_4-2*s_(2,3))
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





R=QQ[t_1,t_2,t_3,t_4,t_5,s_(1,1),s_(1,2),s_(2,2),s_(1,3),s_(2,3),s_(3,3),s_(1,4),s_(2,4),s_(3,4),s_(4,4)]
M=mutableMatrix(R,4,4)
for i from 0 to 3 do(
    for j from i to 3 do(
	M_(i,j)=s_(i+1,j+1);
	M_(j,i)=s_(i+1,j+1);
	)
    )
--M_(1,2)=0;
--M_(2,1)=0;
S=matrix M
detI=minors(2,S)
I=ideal(t_1-s_(1,1),t_2-s_(2,2),t_3-s_(3,3)-s_(4,4), t_4-2*s_(1,2)-2*s_(3,4),t_5-2*s_(2,3))
eliminate(I+sub(detI,R),{s_(1,1),s_(1,1),s_(1,2),s_(2,2),s_(1,3),s_(2,3),s_(3,3),s_(1,4),s_(2,4),s_(3,4),s_(4,4)})	


R=QQ[t_1,t_2,t_3,t_4,t_5,s_(1,1),s_(1,2),s_(2,2),s_(1,3),s_(2,3),s_(3,3),s_(1,4),s_(2,4),s_(3,4),s_(4,4)]
M=mutableMatrix(R,4,4)
for i from 0 to 3 do(
    for j from i to 3 do(
	M_(i,j)=s_(i+1,j+1);
	M_(j,i)=s_(i+1,j+1);
	)
    )
--M_(1,2)=0;
--M_(2,1)=0;
S=matrix M
detI=minors(2,S)
I=ideal(t_1-s_(1,1),t_2-s_(2,2),t_3-s_(3,3)-s_(4,4), t_4-2*s_(1,2)-2*s_(3,4),t_5-2*s_(2,3))
eliminate(I+sub(detI,R),{s_(1,1),s_(1,1),s_(1,2),s_(2,2),s_(1,3),s_(2,3),s_(3,3),s_(1,4),s_(2,4),s_(3,4),s_(4,4)})	



        



	