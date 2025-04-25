-- K1 union K5 with all vertices equal
------------------------------------------------------------
restart
load "ColoredGraphicalModels_OK_RH.m2"
-- Ring with as many variables as the number of colors
R=QQ[l_1..l_11]
-- Concentration matrix
K=matrix{{l_1,0,0,0,0,0},{0,l_1,l_2,l_3,l_4,l_5},{0,l_2,l_1,l_6,l_7,l_8},{0,l_3,l_6,l_1,l_9,l_10},{0,l_4,l_7,l_9,l_1,l_11},{0,l_5,l_8,l_10,l_11,l_1}}
-- Setup: p=size of matrix, n=number of colors (dim of linear space)
-- Rtotal=aux ring with s and t, S=generic sample covariance matrix
(p,n,Rtotal,S)=coloredData(K)
-- ideal of sufficient statistics
stats=sub(suffStat(K),Rtotal)
netList stats_*
--Elimination ideals for rank n (i.e. minors of size n+1)
IG3=time sub(rankProjection(stats,3,S),ring(suffStat(K)));
IG3 --0
IG2=sub(rankProjection(stats,2,S),ring(suffStat(K)));
IG2 --not 0

K=matrix{{l_1,l_2,l_3,l_4,l_5},{l_2,l_1,l_6,l_7,l_8},{l_3,l_6,l_1,l_9,l_10},{l_4,l_7,l_9,l_1,l_11},{l_5,l_8,l_10,l_11,l_1}}
(p,n,Rtotal,S)=coloredData(K)
-- ideal of sufficient statistics
stats=sub(suffStat(K),Rtotal)
netList stats_*
--Elimination ideals for rank n (i.e. minors of size n+1)
IG3=time sub(rankProjection(stats,3,S),ring(suffStat(K)));
IG3 --0
IG2=sub(rankProjection(stats,2,S),ring(suffStat(K)));
IG2 --not 0


-- K1 union K3 with all vertices equal
------------------------------------------------------------
restart
load "ColoredGraphicalModels_OK_RH.m2"
-- Ring with as many variables as the number of colors
R=QQ[l_1..l_4]
-- Concentration matrix
K=matrix{{l_1,0,0,0},{0,l_1,l_2,l_3},{0,l_2,l_1,l_4},{0,l_3,l_4,l_1}}
-- Setup: p=size of matrix, n=number of colors (dim of linear space)
-- Rtotal=aux ring with s and t, S=generic sample covariance matrix
(p,n,Rtotal,S)=coloredData(K)
-- ideal of sufficient statistics
stats=sub(suffStat(K),Rtotal)
netList stats_*
--Elimination ideals for rank n (i.e. minors of size n+1)
IG1=time sub(rankProjection(stats,1,S),ring(suffStat(K)));
IG1 --0

K=matrix{{l_1,l_2,l_3},{l_2,l_1,l_4},{l_3,l_4,l_1}}
-- Setup: p=size of matrix, n=number of colors (dim of linear space)
-- Rtotal=aux ring with s and t, S=generic sample covariance matrix
(p,n,Rtotal,S)=coloredData(K)
-- ideal of sufficient statistics
stats=sub(suffStat(K),Rtotal)
netList stats_*
--Elimination ideals for rank n (i.e. minors of size n+1)
IG1=time sub(rankProjection(stats,1,S),ring(suffStat(K)));
IG1 --not 0

IG2=time sub(rankProjection(stats,2,S),ring(suffStat(K)));
IG2 --0

---!!! GCR(K3)=2 but GCR(K1 union K3)=1

-- K1 union K7 with all vertices equal
------------------------------------------------------------
restart
load "ColoredGraphicalModels_OK_RH.m2"
-- Ring with as many variables as the number of colors
R=QQ[l_1..l_22]
-- Concentration matrix
K=matrix{{l_1,0,0,0,0,0,0,0},
         {0,l_1,l_2,l_3,l_4,l_5,l_6,l_7},
	 {0,l_2,l_1,l_8,l_9,l_10,l_11,l_12},
	 {0,l_3,l_8,l_1,l_13,l_14,l_15,l_16},
	 {0,l_4,l_9,l_13,l_1,l_17,l_18,l_19},
	 {0,l_5,l_10,l_14,l_17,l_1,l_20,l_21},
	 {0,l_6,l_11,l_15,l_18,l_20,l_1,l_22},
	 {0,l_7,l_12,l_16,l_19,l_21,l_22,l_1}}
-- Setup: p=size of matrix, n=number of colors (dim of linear space)
-- Rtotal=aux ring with s and t, S=generic sample covariance matrix
(p,n,Rtotal,S)=coloredData(K)
-- ideal of sufficient statistics
stats=sub(suffStat(K),Rtotal)
netList stats_*
--Elimination ideals for rank n (i.e. minors of size n+1)
IG5=time sub(rankProjection(stats,5,S),ring(suffStat(K)));
IG5==ideal{0}
IG4=time sub(rankProjection(stats,3,S),ring(suffStat(K)));
IG4==ideal{0} 

K=matrix{{l_1,l_2,l_3,l_4,l_5,l_6,l_7},
	 {l_2,l_1,l_8,l_9,l_10,l_11,l_12},
	 {l_3,l_8,l_1,l_13,l_14,l_15,l_16},
	 {l_4,l_9,l_13,l_1,l_17,l_18,l_19},
	 {l_5,l_10,l_14,l_17,l_1,l_20,l_21},
	 {l_6,l_11,l_15,l_18,l_20,l_1,l_22},
	 {l_7,l_12,l_16,l_19,l_21,l_22,l_1}}
     
(p,n,Rtotal,S)=coloredData(K)
-- ideal of sufficient statistics
stats=sub(suffStat(K),Rtotal)
netList stats_*
--Elimination ideals for rank n (i.e. minors of size n+1)
IG5=time sub(rankProjection(stats,5,S),ring(suffStat(K)));
IG5==ideal{0}
IG4=time sub(rankProjection(stats,3,S),ring(suffStat(K)));
IG4==ideal{0} 


-- K1 union K6 with all vertices equal
------------------------------------------------------------
restart
load "ColoredGraphicalModels_OK_RH.m2"
-- Ring with as many variables as the number of colors
R=QQ[l_1..l_16]
-- Concentration matrix
K=matrix{{l_1,0,0,0,0,0,0},
         {0,l_1,l_2,l_3,l_4,l_5,l_6},
	 {0,l_2,l_1,l_7,l_8,l_9,l_10},
	 {0,l_3,l_7,l_1,l_11,l_12,l_13},
	 {0,l_4,l_8,l_11,l_1,l_14,l_15},
	 {0,l_5,l_9,l_12,l_14,l_1,l_16},
	 {0,l_6,l_10,l_13,l_15,l_16,l_1}}
-- Setup: p=size of matrix, n=number of colors (dim of linear space)
-- Rtotal=aux ring with s and t, S=generic sample covariance matrix
(p,n,Rtotal,S)=coloredData(K)
-- ideal of sufficient statistics
stats=sub(suffStat(K),Rtotal)
netList stats_*
--Elimination ideals for rank n (i.e. minors of size n+1)
IG4=time sub(rankProjection(stats,4,S),ring(suffStat(K)));
 -- used 14.4722 seconds
IG4 --0
IG3=time sub(rankProjection(stats,3,S),ring(suffStat(K)));
IG3==ideal{0} 

K=matrix{{l_1,l_2,l_3,l_4,l_5,l_6},
	 {l_2,l_1,l_7,l_8,l_9,l_10},
	 {l_3,l_7,l_1,l_11,l_12,l_13},
	 {l_4,l_8,l_11,l_1,l_14,l_15},
	 {l_5,l_9,l_12,l_14,l_1,l_16},
	 {l_6,l_10,l_13,l_15,l_16,l_1}}

(p,n,Rtotal,S)=coloredData(K)
-- ideal of sufficient statistics
stats=sub(suffStat(K),Rtotal)
netList stats_*
--Elimination ideals for rank n (i.e. minors of size n+1)
IG4=time sub(rankProjection(stats,4,S),ring(suffStat(K)));
IG4==ideal{0}
IG3=time sub(rankProjection(stats,3,S),ring(suffStat(K)));
IG3==ideal{0} 
