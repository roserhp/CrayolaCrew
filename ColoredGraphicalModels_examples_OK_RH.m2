-- 1-2-3-1 + 3-4
------------------------------------------------------------
restart
load "ColoredGraphicalModels_OK_RH.m2"
-- Ring with as many variables as the number of colors
R=QQ[l_1..l_6]
-- Concentration matrix
K=matrix{{l_1,l_4,l_5,0},{l_4,l_1,l_5,0},{l_5,l_5,l_2,l_6},{0,0,l_6,l_3}}
-- Setup: p=size of matrix, n=number of colors (dim of linear space)
-- Rtotal=aux ring with s and t, S=generic sample covariance matrix
(p,n,Rtotal,S)=coloredData(K)
-- ideal of sufficient statistics
stats=sub(suffStat(K),Rtotal)

--Elimination ideals for rank n (i.e. minors of size n+1)
IG3=sub(rankProjection(stats,3,S),ring(suffStat(K)))
IG2=sub(rankProjection(stats,2,S),ring(suffStat(K)))
IG1=sub(rankProjection(stats,1,S),ring(suffStat(K)))
codim IG1,dim IG1,degree IG1
netList (trim IG1)_* --2 quadrics

-- Empirical check of whether there are PD matrices that vanish for some
-- of the generators of the elimination ideal
use ring IG1
netList empiricalVanishingPolynomials(IG1,4,10000,p,n,stats) --second entry: rank of sample cov matrix
-- If no polynomials are returned here, it is likely that the MLE does not exist, so we do infeasibility test
--Infeasibility test (of the vanishing of a single generator of the elimination ideal)
restart
load "ColoredGraphicalModels_OK_RH.m2"
RL=QQ[t_1..t_6,u11,u12,u13,u14,u22,u23,u24,u33,u34,u44]
LT=matrix{{u11,0,0,0},{u12,u22,0,0},{u13,u23,u33,0},{u14,u24,u34,u44}}
PSD=LT*transpose(LT)
f=t_5^2-4*t_4*t_2-4*t_2*t_1 --does not work
f=t_6^2-4*t_3*t_2 --works for -fu
fu=sub(f,{t_1=>PSD_(0,0)+PSD_(1,1),
	t_2=>PSD_(2,2),
	t_3=>PSD_(3,3),
	t_4=>2*PSD_(0,1),
	t_5=>2*PSD_(0,2)+2*PSD_(1,2),
	t_6=>2*PSD_(2,3)})

loadPackage "SemidefiniteProgramming"
loadPackage "SumsOfSquares"
sol=solveSOS (-fu)
peek sol
s=sosPoly sol
peek s
coefs=s#coefficients
gene=s#generators
netList gene
prod=apply(coefs,gene,(i,j)->i*j^2)
sos=sum(prod)
-fu==sos --true



-- 3-cycle with all colored edges and vertices
------------------------------------------------------------
--Test ML degree for different ranks
restart
load "ColoredGraphicalModels_OK_RH.m2"
R=QQ[l_1,l_2]
K=matrix{{l_1,l_2,l_2},{l_2,l_1,l_2},{l_2,l_2,l_1}}
(p,n,Rtotal,S)=coloredData(K)
stats=sub(suffStat(K),Rtotal)

IG2=sub(rankProjection(stats,2,S),ring(suffStat(K)))
IG1=sub(rankProjection(stats,1,S),ring(suffStat(K)))


-- undirected 3-path
------------------------------------------------------------
--Test ML degree for different ranks
restart
load "ColoredGraphicalModels_OK_RH.m2"
R=QQ[l_1..l_4]
K=matrix{{l_1,l_3,0,0},{l_3,l_2,l_4,0},{0,l_4,l_2,l_3},{0,0,l_3,l_1}}
(p,n,Rtotal,S)=coloredData(K)
stats=sub(suffStat(K),Rtotal)

IG2=sub(rankProjection(stats,2,S),ring(suffStat(K)))
IG1=sub(rankProjection(stats,1,S),ring(suffStat(K)))


-----------cherry
restart
load "ColoredGraphicalModels_OK_RH.m2"
R=QQ[l_1..l_5]
K=matrix{{l_1,l_4,0},{l_4,l_2,l_5},{0,l_5,l_3}}
(p,n,Rtotal,S)=coloredData(K)
stats=sub(suffStat(K),Rtotal)

IG2=sub(rankProjection(stats,2,S),ring(suffStat(K)))
IG1=sub(rankProjection(stats,1,S),ring(suffStat(K)))


-----------colored cherry
restart
load "ColoredGraphicalModels_OK_RH.m2"
R=QQ[l_1..l_4]
K=matrix{{l_1,l_3,0},{l_3,l_2,l_4},{0,l_4,l_2}}
(p,n,Rtotal,S)=coloredData(K)
stats=sub(suffStat(K),Rtotal)

IG2=sub(rankProjection(stats,2,S),ring(suffStat(K)))
IG1=sub(rankProjection(stats,1,S),ring(suffStat(K)))

restart
load "ColoredGraphicalModels_OK_RH.m2"
RL=QQ[t_1..t_4,u11,u12,u13,u22,u23,u33]
LT=matrix{{u11,0,0},{u12,u22,0},{u13,u23,u33}}
PSD=LT*transpose(LT)
f=t_3^4-4*t_3^2*t_2*t_1+4*t_4^2*t_1^2
fu=sub(f,{t_1=>PSD_(0,0), 
	t_2=>PSD_(2,2)+PSD_(1,1),
	t_3=>2*PSD_(0,1),
	t_4=>2*PSD_(1,2)})

loadPackage "SemidefiniteProgramming"
loadPackage "SumsOfSquares"
sol=solveSOS (-fu)  ---no solution
peek sol
s=sosPoly sol
peek s
coefs=s#coefficients
gene=s#generators
netList gene
prod=apply(coefs,gene,(i,j)->i*j^2)
sos=sum(prod)
-fu==sos --true



-----------colored 5-star
restart
load "ColoredGraphicalModels_OK_RH.m2"
R=QQ[l_1..l_7]
K=matrix{{l_1,l_5,l_5,l_6,l_6,l_7},{l_5,l_2,0,0,0,0},{l_5,0,l_2,0,0,0},{l_6,0,0,l_3,0,0},{l_6,0,0,0,l_3,0},{l_7,0,0,0,0,l_4}}
(p,n,Rtotal,S)=coloredData(K)
stats=sub(suffStat(K),Rtotal)

IG2=sub(rankProjection(stats,2,S),ring(suffStat(K)))
IG1=sub(rankProjection(stats,1,S),ring(suffStat(K)))
--t_7^2-4t_4t_1 is a ppal determinant of a PD matrix: no MLE existence

-----------colored path 5 colors
restart
load "ColoredGraphicalModels_OK_RH.m2"
R=QQ[l_1..l_5]
K=matrix{{l_1,l_4,0,0},{l_4,l_2,l_5,0},{0,l_5,l_3,l_4},{0,0,l_4,l_3}}
(p,n,Rtotal,S)=coloredData(K)
stats=sub(suffStat(K),Rtotal)

IG2=sub(rankProjection(stats,2,S),ring(suffStat(K)))
IG1=sub(rankProjection(stats,1,S),ring(suffStat(K)))

use ring IG1
netList empiricalVanishingPolynomials(IG1,5,10000,p,n,stats)
-- It changes signs: evidence for existence 
--Check existence of MLE for several sample covariance matrix with a given rank: 
--empiricalMLEexistence(rank,number of cov matrix to consider,K)
empiricalMLEexistence(1,10,K) --some but not all  MLE exist
empiricalMLEexistence(2,10,K) --all exist
empiricalMLEexistence(3,10,K) --all exist


-- 4-cycle with all vertices with the same color

restart
load "ColoredGraphicalModels_OK_RH.m2"
R=QQ[l_1..l_5]
K=matrix{{l_1,l_2,0,l_5},{l_2,l_1,l_3,0},{0,l_3,l_1,l_4},{l_5,0,l_4,l_1}}
(p,n,Rtotal,S)=coloredData(K)
stats=sub(suffStat(K),Rtotal)
netList stats_*

IG2=sub(rankProjection(stats,2,S),ring(suffStat(K)))
IG1=sub(rankProjection(stats,1,S),ring(suffStat(K)))

use ring IG1
netList empiricalVanishingPolynomials(IG1,3,10000,p,n,stats)
empiricalMLEexistence(1,100,K) --seems like ALL MLE exist although the elimination ideal is no zero

      
-- 1-2-3 with some colors
------------------------------------------------------------
restart
load "ColoredGraphicalModels_OK_RH.m2"
R=QQ[l_1..l_3]
-- Concentration matrix
K=matrix{{l_1,l_3,0},{l_3,l_2,l_3},{0,l_3,l_1}}
-- Setup: p=size of matrix, n=number of colors (dim of linear space)
-- Rtotal=aux ring with s and t, S=generic sample covariance matrix
(p,n,Rtotal,S)=coloredData(K)
-- ideal of sufficient statistics
stats=sub(suffStat(K),Rtotal)
netList stats_*

--Elimination ideals for rank n (i.e. minors of size n+1)
IG2=sub(rankProjection(stats,2,S),ring(suffStat(K)))
IG1=sub(rankProjection(stats,1,S),ring(suffStat(K)))


restart
load "ColoredGraphicalModels_OK_RH.m2"
R=QQ[l_1..l_4]
-- Concentration matrix
K=matrix{{l_1,l_3,0},{l_3,l_1,l_4},{0,l_4,l_2}}
-- Setup: p=size of matrix, n=number of colors (dim of linear space)
-- Rtotal=aux ring with s and t, S=generic sample covariance matrix
(p,n,Rtotal,S)=coloredData(K)
-- ideal of sufficient statistics
stats=sub(suffStat(K),Rtotal)
netList stats_*

--Elimination ideals for rank n (i.e. minors of size n+1)
IG1=sub(rankProjection(stats,1,S),ring(suffStat(K)))


-- 1-2 with some colors
------------------------------------------------------------
restart
load "ColoredGraphicalModels_OK_RH.m2"
R=QQ[l_1..l_2]
-- Concentration matrix
K=matrix{{l_1,l_2},{l_2,l_1}}
-- Setup: p=size of matrix, n=number of colors (dim of linear space)
-- Rtotal=aux ring with s and t, S=generic sample covariance matrix
(p,n,Rtotal,S)=coloredData(K)
-- ideal of sufficient statistics
stats=sub(suffStat(K),Rtotal)
netList stats_*

--Elimination ideals for rank n (i.e. minors of size n+1)
IG1=sub(rankProjection(stats,1,S),ring(suffStat(K)))


-- 1-2 3 two vertices equal
------------------------------------------------------------
restart
load "ColoredGraphicalModels_OK_RH.m2"
-- Ring with as many variables as the number of colors
R=QQ[l_1..l_3]
-- Concentration matrix
K=matrix{{l_1,l_3,0},{l_3,l_1,0},{0,0,l_2}}
-- Setup: p=size of matrix, n=number of colors (dim of linear space)
-- Rtotal=aux ring with s and t, S=generic sample covariance matrix
(p,n,Rtotal,S)=coloredData(K)
-- ideal of sufficient statistics
stats=sub(suffStat(K),Rtotal)

--Elimination ideals for rank n (i.e. minors of size n+1)
IG2=sub(rankProjection(stats,2,S),ring(suffStat(K)))
IG1=sub(rankProjection(stats,1,S),ring(suffStat(K)))


-- 1-2 3 all vertices equal
------------------------------------------------------------
restart
load "ColoredGraphicalModels_OK_RH.m2"
-- Ring with as many variables as the number of colors
R=QQ[l_1,l_2]
-- Concentration matrix
K=matrix{{l_1,l_2,0},{l_2,l_1,0},{0,0,l_1}}
-- Setup: p=size of matrix, n=number of colors (dim of linear space)
-- Rtotal=aux ring with s and t, S=generic sample covariance matrix
(p,n,Rtotal,S)=coloredData(K)
-- ideal of sufficient statistics
stats=sub(suffStat(K),Rtotal)

--Elimination ideals for rank n (i.e. minors of size n+1)
IG2=sub(rankProjection(stats,2,S),ring(suffStat(K)))
IG1=sub(rankProjection(stats,1,S),ring(suffStat(K)))



-----------------------------------------------------------
-----------------------------------------------------------
-- Clique sums of colored 3-cliques
-----------------------------------------------------------
-----------------------------------------------------------

-------------------------------------------------------------------------------------------
--Example 1: clique sum over edge with equally colored vertices. No other colors elsewhere.
--- gcr(G_i)=2 but gcr(G)=3!!!
-------------------------------------------------------------------------------------------
restart
load "ColoredGraphicalModels_OK_RH.m2"
R=QQ[l_1..l_8]
-- Concentration matrix
K=matrix{{l_1,l_7,l_6,0},{l_7,l_2,l_4,l_8},{l_6,l_4,l_2,l_5},{0,l_8,l_5,l_3}}
-- Setup: p=size of matrix, n=number of colors (dim of linear space)
-- Rtotal=aux ring with s and t, S=generic sample covariance matrix
(p,n,Rtotal,S)=coloredData(K)
-- ideal of sufficient statistics
stats=sub(suffStat(K),Rtotal)
netList stats_*

--Elimination ideals for rank n (i.e. minors of size n+1)
IG3=sub(rankProjection(stats,3,S),ring(suffStat(K)))
IG2=sub(rankProjection(stats,2,S),ring(suffStat(K)))
IG1=sub(rankProjection(stats,1,S),ring(suffStat(K)))

---What about the MLE?
time empiricalMLEexistence(2,100,K) --62
 -- used 156.482 seconds
use ring IG2
netList empiricalVanishingPolynomials(IG2,4,10000,p,n,stats)
use ring IG1
netList empiricalVanishingPolynomials(IG1,4,10000,p,n,stats)--some but not all vanish

restart
load "ColoredGraphicalModels_OK_RH.m2"
RL=QQ[t_1..t_8,u11,u12,u13,u14,u22,u23,u24,u33,u34,u44]
LT=matrix{{u11,0,0,0},{u12,u22,0,0},{u13,u23,u33,0},{u14,u24,u34,u44}}
PSD=LT*transpose(LT)
f=t_7^2+t_6^2-4*t_2*t_1 
fu=sub(f,{t_1=>PSD_(0,0),
	t_2=>PSD_(1,1)+PSD_(2,2),
	t_3=>PSD_(3,3),
	t_4=>2*PSD_(1,2),
	t_5=>2*PSD_(2,3),
	t_6=>2*PSD_(0,2),
	t_7=>2*PSD_(0,1),
	t_8=>2*PSD_(1,3)})
--fu is already an SOS
--SUMMARY: wmlt=2, mlt=gcr=3

--double-check its components before clique-sum
restart
load "ColoredGraphicalModels_OK_RH.m2"
R=QQ[l_1,l_2,l_4,l_6,l_7]
-- Concentration matrix
K=matrix{{l_1,l_7,l_6},{l_7,l_2,l_4},{l_6,l_4,l_2}}
(p,n,Rtotal,S)=coloredData(K)
-- ideal of sufficient statistics
stats=sub(suffStat(K),Rtotal)
netList stats_*

--Elimination ideals for rank n (i.e. minors of size n+1)
IG2=sub(rankProjection(stats,2,S),ring(suffStat(K)))
IG1=sub(rankProjection(stats,1,S),ring(suffStat(K)))


-------------------------------------------------------------------------------------------
--Example 2: clique sum over edge with equally colored vertices. Remaining two edges of G_2 have same color
--- gcr(G_i)=gcr(G)=2 and I_G,1=I_G1,1+I_G2,1
-------------------------------------------------------------------------------------------
restart
load "ColoredGraphicalModels_OK_RH.m2"
R=QQ[l_1..l_7]
-- Concentration matrix
K=matrix{{l_1,l_7,l_6,0},{l_7,l_2,l_4,l_5},{l_6,l_4,l_2,l_5},{0,l_5,l_5,l_3}}
-- Setup: p=size of matrix, n=number of colors (dim of linear space)
-- Rtotal=aux ring with s and t, S=generic sample covariance matrix
(p,n,Rtotal,S)=coloredData(K)
-- ideal of sufficient statistics
stats=sub(suffStat(K),Rtotal)
netList stats_*

--Elimination ideals for rank n (i.e. minors of size n+1)
IG3=sub(rankProjection(stats,3,S),ring(suffStat(K)))
IG2=sub(rankProjection(stats,2,S),ring(suffStat(K)))
IG1=sub(rankProjection(stats,1,S),ring(suffStat(K)))
netList IG1_*
--ideal(t_5^2-4*t_4*t_3-4*t_3*t_2,t_7*t_6-2*t_4*t_1,t_7^2+t_6^2-4*t_2*t_1)
--I_G,1=I_G1,1+I_G2,1?

IG11=ideal(t_7*t_6-2*t_4*t_1,t_7^2+t_6^2-4*t_2*t_1)
IG21=ideal(t_5^2-4*t_4*t_3-4*t_3*t_2)

IG1==IG11+IG21 --true

--computations of elimination ideal of sum components
restart
load "ColoredGraphicalModels_OK_RH.m2"
R=QQ[l_1,l_2,l_4,l_6,l_7]
-- Concentration matrix
K=matrix{{l_1,l_7,l_6},{l_7,l_2,l_4},{l_6,l_4,l_2}}
(p,n,Rtotal,S)=coloredData(K)
-- ideal of sufficient statistics
stats=sub(suffStat(K),Rtotal)
netList stats_*

--Elimination ideals for rank n (i.e. minors of size n+1)
IG2=sub(rankProjection(stats,2,S),ring(suffStat(K)))
IG1=sub(rankProjection(stats,1,S),ring(suffStat(K)))
toString IG1
--ideal(t_5*t_4-2*t_3*t_1,t_5^2+t_4^2-4*t_2*t_1)
--t_3->t_4, t_4->t_6, t_5->t_7

restart
load "ColoredGraphicalModels_OK_RH.m2"
R=QQ[l_2..l_5]
-- Concentration matrix
K=matrix{{l_2,l_4,l_5},{l_4,l_2,l_5},{l_5,l_5,l_3}}
-- Setup: p=size of matrix, n=number of colors (dim of linear space)
-- Rtotal=aux ring with s and t, S=generic sample covariance matrix
(p,n,Rtotal,S)=coloredData(K)
-- ideal of sufficient statistics
stats=sub(suffStat(K),Rtotal)
netList stats_*
--Elimination ideals for rank n (i.e. minors of size n+1)
IG2=sub(rankProjection(stats,2,S),ring(suffStat(K)))
IG1=sub(rankProjection(stats,1,S),ring(suffStat(K)))
toString IG1
--ideal(t_4^2-4*t_3*t_2-4*t_2*t_1)
--t_1->t_2,t_2->t_3,t_3->t_4,t_4->t_5


-------------------------------------------------------------------------------------------
--Example 2: clique sum over edge with equally colored vertices. Additional edge of G_2
-- has same color as edge in the intersection
--- gcr(G_i)=gcr(G)=2 but I_G,1 not I_G1,1+I_G2,1
-------------------------------------------------------------------------------------------
restart
load "ColoredGraphicalModels_OK_RH.m2"
R=QQ[l_1..l_7]
-- Concentration matrix
K=matrix{{l_1,l_7,l_6,0},{l_7,l_2,l_4,l_4},{l_6,l_4,l_2,l_5},{0,l_4,l_5,l_3}}
-- Setup: p=size of matrix, n=number of colors (dim of linear space)
-- Rtotal=aux ring with s and t, S=generic sample covariance matrix
(p,n,Rtotal,S)=coloredData(K)
-- ideal of sufficient statistics
stats=sub(suffStat(K),Rtotal)
netList stats_*

--Elimination ideals for rank n (i.e. minors of size n+1)
IG3=sub(rankProjection(stats,3,S),ring(suffStat(K)))
IG2=sub(rankProjection(stats,2,S),ring(suffStat(K)))
IG1=sub(rankProjection(stats,1,S),ring(suffStat(K)))
netList IG1_*

IG11=ideal(t_7*t_6-2*t_4*t_1,t_7^2+t_6^2-4*t_2*t_1)
IG11_1 % gb IG1 -- 0
IG11_0 % gb IG1 -- not in ideal! comes from a 2-minor that involves the edge in the intersection
-- that has the same color as an edge in G2\intersection
IG12=ideal(t_5^4+4*t_5^3*t_3+4*t_5^2*t_3^2+4*t_4^2*t_3^2-4*t_5^2*t_3*t_2-16*t_5*t_3^2*t_2-16*t_3^3*t_2)

IG12_0 % gb IG1 --0


restart
load "ColoredGraphicalModels_OK_RH.m2"
R=QQ[l_2..l_5]
-- Concentration matrix
K=matrix{{l_2,l_4,l_4},{l_4,l_2,l_5},{l_4,l_5,l_3}}
-- Setup: p=size of matrix, n=number of colors (dim of linear space)
-- Rtotal=aux ring with s and t, S=generic sample covariance matrix
(p,n,Rtotal,S)=coloredData(K)
-- ideal of sufficient statistics
stats=sub(suffStat(K),Rtotal)
netList stats_*
--Elimination ideals for rank n (i.e. minors of size n+1)
IG2=sub(rankProjection(stats,2,S),ring(suffStat(K)))
IG1=sub(rankProjection(stats,1,S),ring(suffStat(K)))
toString IG1
--ideal(t_4^4+4*t_4^3*t_2+4*t_4^2*t_2^2+4*t_3^2*t_2^2-4*t_4^2*t_2*t_1-16*t_4*t_2^2*t_1-16*t_2^3*t_1)
--t_1->t_2,t_2->t_3,t_3->t_4,t_4->t_5




      
-- disjoint union of K_1 and K_3
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
IG3=sub(rankProjection(stats,3,S),ring(suffStat(K)))
IG2=sub(rankProjection(stats,2,S),ring(suffStat(K)))
IG1=sub(rankProjection(stats,1,S),ring(suffStat(K)))


-- disjoint union of K_1 and K_4
restart
load "ColoredGraphicalModels_OK_RH.m2"
-- Ring with as many variables as the number of colors
R=QQ[l_1..l_7]
-- Concentration matrix
K=matrix{{l_1,0,0,0,0},{0,l_1,l_2,l_5,l_7},{0,l_2,l_1,l_3,l_6},{0,l_5,l_3,l_1,l_4},{0,l_7,l_6,l_4,l_1}}
-- Setup: p=size of matrix, n=number of colors (dim of linear space)
-- Rtotal=aux ring with s and t, S=generic sample covariance matrix
(p,n,Rtotal,S)=coloredData(K)
-- ideal of sufficient statistics
stats=sub(suffStat(K),Rtotal)
netList stats_*

--Elimination ideals for rank n (i.e. minors of size n+1)
IG3=sub(rankProjection(stats,3,S),ring(suffStat(K)))
IG2=sub(rankProjection(stats,2,S),ring(suffStat(K)))
IG1=sub(rankProjection(stats,1,S),ring(suffStat(K)))

use ring IG1
time empiricalMLEexistence(1,1,K) --1
time empiricalMLEexistence(1,50,K) --50


-----------------------------------------------------------
-----------------------------------------------------------
-- Disjoint unions of K_1 and K_n
-----------------------------------------------------------
-----------------------------------------------------------


-- disjoint union of K_1 and K_5
restart
load "ColoredGraphicalModels_OK_RH.m2"
-- Ring with as many variables as the number of colors
R=QQ[l_1..l_11]
-- Concentration matrix
K=matrix{{l_1,0,0,0,0,0},{0,l_1,l_2,l_5,l_7,l_8},{0,l_2,l_1,l_3,l_6,l_9},{0,l_5,l_3,l_1,l_4,l_10},{0,l_7,l_6,l_4,l_1,l_11},{0,l_8,l_9,l_10,l_11,l_1}}
-- Setup: p=size of matrix, n=number of colors (dim of linear space)
-- Rtotal=aux ring with s and t, S=generic sample covariance matrix
(p,n,Rtotal,S)=coloredData(K)
-- ideal of sufficient statistics
stats=sub(suffStat(K),Rtotal)
netList stats_*

--Elimination ideals for rank n (i.e. minors of size n+1)
IG3=time sub(rankProjection(stats,3,S),ring(suffStat(K)));--0, very fast
IG2=time sub(rankProjection(stats,2,S),ring(suffStat(K)));--very fast
netList IG2_*
IG1=time sub(rankProjection(stats,1,S),ring(suffStat(K)));--very fast
netList IG1_*

-- K_5
restart
load "ColoredGraphicalModels_OK_RH.m2"
-- Ring with as many variables as the number of colors
R=QQ[l_1..l_11]
-- Concentration matrix
K=matrix{{l_1,l_2,l_5,l_7,l_8},{l_2,l_1,l_3,l_6,l_9},{l_5,l_3,l_1,l_4,l_10},{l_7,l_6,l_4,l_1,l_11},{l_8,l_9,l_10,l_11,l_1}}
-- Setup: p=size of matrix, n=number of colors (dim of linear space)
-- Rtotal=aux ring with s and t, S=generic sample covariance matrix
(p,n,Rtotal,S)=coloredData(K)
-- ideal of sufficient statistics
stats=sub(suffStat(K),Rtotal)
netList stats_*

--Elimination ideals for rank n (i.e. minors of size n+1)
IG3=time sub(rankProjection(stats,3,S),ring(suffStat(K)));--0, very fast
IG2=time sub(rankProjection(stats,2,S),ring(suffStat(K)));--very fast
netList IG2_* ---SOOOO LONG
betti IG2 --17 gens
IG1=time sub(rankProjection(stats,1,S),ring(suffStat(K)));--very fast
betti IG1 --25 gens
netList IG1_*

use ring IG1
time empiricalMLEexistence(1,1,K) --0 --> MLE>=2
 -- used 2813.17 seconds
length empiricalVanishingPolynomials(IG1,5,1000,p,n,stats) --25
use ring IG2
length empiricalVanishingPolynomials(IG2,5,1000,p,n,stats) --17
time empiricalMLEexistence(2,1,K) --INTERRUPTED PROCESS



-----------------------------------------------------------
-----------------------------------------------------------
-- Trees
-----------------------------------------------------------
-----------------------------------------------------------


-- 7-tree with colors
------------------------------------------------------------
restart
load "ColoredGraphicalModels_OK_RH.m2"
-- Ring with as many variables as the number of colors
R=QQ[l_1..l_4]
-- Concentration matrix
K=matrix{{l_1,0,0,0,0,0,l_4},{0,l_1,0,0,0,0,l_3},{0,0,l_1,0,0,0,l_3},{0,0,0,l_1,0,0,l_3},{0,0,0,0,l_1,0,l_3},{0,0,0,0,0,l_1,l_3},{l_4,l_3,l_3,l_3,l_3,l_3,l_2}}
-- Setup: p=size of matrix, n=number of colors (dim of linear space)
-- Rtotal=aux ring with s and t, S=generic sample covariance matrix
(p,n,Rtotal,S)=coloredData(K)
-- ideal of sufficient statistics
stats=sub(suffStat(K),Rtotal)
netList stats_*
--Elimination ideals for rank n (i.e. minors of size n+1)
IG3=time sub(rankProjection(stats,3,S),ring(suffStat(K)));
IG2=time sub(rankProjection(stats,2,S),ring(suffStat(K)));
IG1=time sub(rankProjection(stats,1,S),ring(suffStat(K)));



time empiricalMLEexistence(1,100,K) --100


restart
load "ColoredGraphicalModels_OK_RH.m2"
-- Ring with as many variables as the number of colors
R=QQ[l_1..l_4]
-- Concentration matrix
K=matrix{{l_1,0,0,0,l_4},{0,l_1,0,0,l_3},{0,0,l_1,0,l_3},{0,0,0,l_1,l_3},{l_4,l_3,l_3,l_3,l_2}}
-- Setup: p=size of matrix, n=number of colors (dim of linear space)
-- Rtotal=aux ring with s and t, S=generic sample covariance matrix
(p,n,Rtotal,S)=coloredData(K)
-- ideal of sufficient statistics
stats=sub(suffStat(K),Rtotal)
netList stats_*
--Elimination ideals for rank n (i.e. minors of size n+1)
IG3=time sub(rankProjection(stats,3,S),ring(suffStat(K)));
IG2=time sub(rankProjection(stats,2,S),ring(suffStat(K)));
use Rtotal
IG1=time sub(rankProjection(stats,1,S),ring(suffStat(K)));
IG1 --0


time empiricalMLEexistence(1,100,K) --100

