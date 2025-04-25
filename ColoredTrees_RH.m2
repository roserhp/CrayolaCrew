---------------------------------------------------
-- TREES 
--------------------------------------------------

restart
load "ColoredGraphicalModels_OK_RH.m2" 
R=QQ[l_1..l_4]
K=matrix{{l_1,l_4,0},{l_4,l_2,l_4},{0,l_4,l_3}}
(p,n,Rtotal,S)=coloredData(K)
stats=sub(suffStat(K),Rtotal)

-------------------------------------------------------------------------------
-- ELIMINATION CRITERION
-------------------------------------------------------------------------------
IG2=sub(rankProjection(stats,2,S),ring(suffStat(K))) --0
IG1=sub(rankProjection(stats,1,S),ring(suffStat(K))) --not 0
--Principal ideal: t_4^4-8*t_4^2*t_3*t_2+16*t_3^2*t_2^2-8*t_4^2*t_2*t_1-32*t_3*t_2^2*t_1+16*t_2^2*t_1^2
f=(sub(IG1,t_4=>0))_0
factor f -- f=t_2^2

-- CONCLUSION: gcr(G)=2

-------------------------------------------------------------------------------
-- EMPIRICAL CHECK 1
-------------------------------------------------------------------------------
--Empirical check: are there PD matrices (i.e. rank 3 sample covariance matrices)
-- that vanish at the generators of IG1? 
-- What it actually checks is whether you can find positive and negative evaluations of the generator at a PD matrix
-- empiricalVanishingPolynomials(elimination ideal, rank of cov matrix, max number of points we try,
-- size of cov matrix, dimension of linear space of concentration matrix, sufficient stats)
netList empiricalVanishingPolynomials(IG1,3,10000,p,n,stats)

-- CONCLUSION: the fact that the single generator changes signs
-- provides evidence for existence of the MLE with some probability (i.e. wmlt(G)=1)
-- At least this ensures that we cannot use a positivity certificate (such as SOS decomposition)
-- to prove that the MLE NEVER exists (which would imply wmlt(G)=mlt(G)=2)


-------------------------------------------------------------------------------
-- EMPIRICAL CHECK 2
-------------------------------------------------------------------------------
--Empirical check: does the MLE for n observations exist (i.e. sample cov matrix of rank n)?
--Here we use score equations 
--empiricalMLEexistence(rank,number of cov matrix to consider,K)
empiricalMLEexistence(1,10000,K) --none exist
empiricalMLEexistence(2,10000,K) --all exist

-- CONCLUSION: computations suggest that wmlt(G)=mlt(G)=2  
-- WARNING: It seems that the problem is the random function in M2 that we are using
-- to generate samples: it always gives positive values. See next section

-------------------------------------------------------------------------------
-- FIND A RANK 1 COVARIANCE MATRIX FOR WHICH THE MLE EXISTS
-------------------------------------------------------------------------------

--Start with easy examples
--Assume t_4=0
use ring(IG1)
f=(sub(IG1,t_4=>0))_0
factor f -- f=16t_2^2(t_3-t_1)^2
--If t_2=0, then Sigma_(1,1)=0 and it cannot be a PD matrix, 
-- so it should be t_3=-t_1 which would mean that x_3=-x_1
--Note that it's an if and only if: If x_3=-x_1, then SCov_(0,1)=-SCov(1,2) and hence t_4=0,
-- namely we are in a submodel with l_4=0, i.e. no edges

--Example:
X=sub(transpose matrix{{1,-2,-1}},QQ)
Semp=X*transpose(X)
I=ideal{jacobian(matrix{{det K}})-det(K)*jacobian(matrix{{trace(Semp*K)}})};
J=saturate(I,det K);
dim J,degree J --(0,1) MLdegree(G)=1
criticalPoints=zeroDimSolve(J) --for the concentration matrix K
criticalMatrices=genListMatrix(criticalPoints,K) --for the concentration matrix K
Sigma=inverse criticalMatrices_0 --for the covariance matrix Sigma
--note that Sigma and S have the same sufficient statistics
--but Sigma lies in a submodel where l_4=0!!!

--Look for more general samples: find a PD matrix which after perturbing entries (1,2) and (2,3) 
-- in a rank 1 matrix such that preserves the suff stats and with any (1,3)

restart
Rx=QQ[x_1..x_3,e,a]
M=mutableMatrix (transpose matrix{{x_1,x_2,x_3}}*matrix{{x_1,x_2,x_3}})
S=matrix M
det S==0
M_(0,1)=M_(0,1)+e
M_(1,2)=M_(1,2)-e
M_(1,0)=M_(1,0)+e
M_(2,1)=M_(2,1)-e
M_(0,2)=a
M_(2,0)=a
Se=matrix M
--Se and S have same sufficient statistics

--Joe's argument in the document
det Se_{0,1}^{0,1}
det Se_{1,2}^{1,2}

-- CONCLUSION: For a single observation (x_1,x_2,x_3) we have x_1x_2<-e/2 and and x_2x_3>e/2
-- A necessary condition for the MLE to exist for a single obs
-- is that either x_1>0,x_2<0,x_3<0 or x_1<0,x_2>0,x_3>0
-- In particular, this proves that MLT(G)=2

--To prove that in the previous open set the MLE always exists we need to prove
--that we can always choose a such that det Se is strictly positive
det Se
factor (det Se) 

netList terms det Se
--NOT CLEAR TO ME WHETHER THIS IS TRUE. Looking at the signs of the terms is not enough

--Example: x_1>0,x_2<0,x_3<0 
restart
load "ColoredGraphicalModels_OK_RH.m2" 
R=QQ[l_1..l_4]
K=matrix{{l_1,l_4,0},{l_4,l_2,l_4},{0,l_4,l_3}}
(p,n,Rtotal,S)=coloredData(K)

X=sub(transpose matrix{{4,-2,-3}},QQ)
Semp=X*transpose(X)
I=ideal{jacobian(matrix{{det K}})-det(K)*jacobian(matrix{{trace(Semp*K)}})};
J=saturate(I,det K);
dim J,degree J --(0,2) MLdegree(G)=2
criticalPoints=zeroDimSolve(J) --for the concentration matrix K
criticalMatrices=genListMatrix(criticalPoints,K) --for the concentration matrix K
MLEK=checkPD(criticalMatrices)
Sigma=inverse MLEK_0 --for the covariance matrix Sigma


--Example:x_1<0,x_2>0,x_3>0
X=sub(transpose matrix{{-4,10,5}},QQ)
Semp=X*transpose(X)
I=ideal{jacobian(matrix{{det K}})-det(K)*jacobian(matrix{{trace(Semp*K)}})};
J=saturate(I,det K);
dim J,degree J --(0,2) MLdegree(G)=2
criticalPoints=zeroDimSolve(J) --for the concentration matrix K
criticalMatrices=genListMatrix(criticalPoints,K) --for the concentration matrix K
MLEK=checkPD(criticalMatrices)
Sigma=inverse MLEK_0 --for the covariance matrix Sigma

--I tried a couple of examples and it seems to work!


-- CONCLUSION: wmlt(G)=1 and mlt(G)=gcr(G)=1


















----------------------------
-- DON'T LOOK AT THIS
----------------------------

Rx=QQ[x_1..x_3,support S]
Sigma=sub(S,Rx)
SCov=transpose matrix{{x_1,x_2,x_3}}*matrix{{x_1,x_2,x_3}}
eq=trim ideal{Sigma_(0,0)-SCov_(0,0),Sigma_(1,1)-SCov_(1,1),Sigma_(2,2)-SCov_(2,2),Sigma_(0,1)+Sigma_(1,2)-SCov_(0,1)-SCov_(1,2)}
netList eq_*
dim eq,degree eq --(5,8)
isPrime eq --true
sateq=time saturate(eq,ideal{prodPrincipalMinors(Sigma)});
sateq==eq --true
--Nothing useful from studying these equations like this

Rxu=QQ[x_1..x_3,u_(1,1),u_(1,2),u_(1,3),u_(2,2),u_(2,3),u_(3,3)]
L=matrix {{u_(1,1), 0,0}, {u_(1,2), u_(2,2),0},{u_(1,3), u_(2,3), u_(3,3)}}
Sigma=L*transpose L --Cholesky decomposition of Sigma (with u_(i,i)>0)
SCov=transpose matrix{{x_1,x_2,x_3}}*matrix{{x_1,x_2,x_3}} --Sample covariance matrix of rank 1
eq=trim ideal{Sigma_(0,0)-SCov_(0,0),Sigma_(1,1)-SCov_(1,1),Sigma_(2,2)-SCov_(2,2),Sigma_(0,1)+Sigma_(1,2)-SCov_(0,1)-SCov_(1,2)}
netList eq_*
dim eq,degree eq --(5,16)
sateq=time saturate(eq,ideal{u_(1,1)*u_(2,2)*u_(3,3)});
sateq==eq --true
PD=time minimalPrimes eq;
apply(PD,i->netList i_*)

eliminate(support L,eq)
