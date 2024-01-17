loadPackage "Graphs"
--EXAMPLE binary tree on 4 leaves using Digraph
n = 4
bfl = digraph ({(0,5), (5,6), (5,7), (6,1), (6,2), (7,3), (7,4)})
Rbfl = PRing(n, ZZ/101, p)
bflThetaLabels = {t_0, a, b, t_1, t_2, t_3, t_4}
Sbfl = getThetaRing(ZZ/101, bflThetaLabels);
bflPhi = getPtoThetaMap(bfl, n, Rbfl, Sbfl);
K = getConcMatrixP(n, p);
D = det K;
bflThetaD = bflPhi(D);
bflFactors = factor bflThetaD;
bflFactors

bflJac = ideal jacobian ideal bflThetaD
bflI = ideal bflThetaD
bflIdeal = bflJac + bflI
degree bflIdeal
dim bflIdeal
mingens bflIdeal
bflIdealFactors = for P in bflIdeal_* list factor P;
netList bflIdealFactors
netList primaryDecomposition bflIdeal
netList minimalPrimes bflIdeal
saturate(bflJac, thetaD)
isSubset(bflI, bflJac)
netList bflJac_*
thetaD
t_0*bflJac_0 + t_1*bflJac_1 + t_2*bflJac_2 + t_3*bflJac_3 == 6*thetaD

X4 = random(QQ^4,QQ^4)
S4 = X4*transpose(X4)
S4 = sub(S4, Sbfl)
T4 = trace(S4*bflPhi(K))

bflLtilde = ideal (jacobian ideal bflThetaD - bflThetaD*(jacobian ideal T4))
degree bflLtilde
dim bflLtilde
netList bflLtilde_*
bflLtildeSat = saturate(bflLtilde, bflThetaD)
degree bflLtildeSat
dim bflLtildeSat
netList bflLtildeSat_*
decompI = primaryDecomposition bflLtildeSat
netList minimalPrimes bflLtildeSat
isPrime(bflLtildeSat)