-----------------------------------------------------------
--3-cycle uncolored
-----------------------------------------------------------
restart
load "ColoredGraphicalModels_OK_RH.m2"
R=QQ[l_1..l_6]
K=matrix{{l_1,l_4,l_5},{l_4,l_2,l_6},{l_5,l_6,l_3}}
(p,n,Rtotal,S)=coloredData(K)
stats=sub(suffStat(K),Rtotal)
J2=ideal{det S}+stats;
netList J2_*

describe Rtotal
Relim=QQ[gens Rtotal,MonomialOrder => Eliminate 6]
describe Relim
JG2=trim sub(J2,Relim);
netList JG2_*
netList flatten entries gens gb JG2

use Rtotal
J1=(trim minors(2,S)+stats);
netList J1_*
JG1=trim sub(J1,Relim);
netList JG1_*
netList flatten entries gens gb JG1


-----------------------------------------------------------
--3-cycle with 2 vertices equal
-----------------------------------------------------------
restart
load "ColoredGraphicalModels_OK_RH.m2"
R=QQ[l_1..l_5]
K=matrix{{l_1,l_3,l_4},{l_3,l_1,l_5},{l_4,l_5,l_2}}
(p,n,Rtotal,S)=coloredData(K)
stats=sub(suffStat(K),Rtotal)
J2=ideal{det S}+stats;
netList J2_*

Relim=QQ[gens Rtotal,MonomialOrder => Eliminate 6]
JG2=trim sub(J2,Relim);
netList JG2_*
netList flatten entries gens gb JG2

use Rtotal
J1=(trim minors(2,S)+stats);
netList J1_*
JG1=trim sub(J1,Relim);
netList JG1_*
netList flatten entries gens gb JG1

-----------------------------------------------------------
--3-cycle with 2 edges equal
-----------------------------------------------------------
restart
load "ColoredGraphicalModels_OK_RH.m2"
R=QQ[l_1..l_5]
K=matrix{{l_1,l_4,l_4},{l_4,l_2,l_5},{l_4,l_5,l_3}}
(p,n,Rtotal,S)=coloredData(K)
stats=sub(suffStat(K),Rtotal)
J2=ideal{det S}+stats;
netList J2_*

Relim=QQ[gens Rtotal,MonomialOrder => Eliminate 6]
JG2=trim sub(J2,Relim);
netList JG2_*
netList flatten entries gens gb JG2

use Rtotal
J1=(trim minors(2,S)+stats);
netList J1_*
JG1=trim sub(J1,Relim);
netList JG1_*
netList flatten entries gens gb JG1


-----------------------------------------------------------
--complete 4-graph uncolored
-----------------------------------------------------------
restart
load "ColoredGraphicalModels_OK_RH.m2"
R=QQ[l_1..l_10]
K=matrix{{l_1,l_5,l_6,l_7},{l_5,l_2,l_8,l_9},{l_6,l_8,l_3,l_10},{l_7,l_9,l_10,l_4}}
(p,n,Rtotal,S)=coloredData(K)
stats=sub(suffStat(K),Rtotal)
J3=ideal{det S}+stats;
netList J3_*
netList flatten entries gens gb J3

Relim=QQ[gens Rtotal,MonomialOrder => Eliminate 10]
JG3=trim sub(J3,Relim);
netList JG3_*
netList flatten entries gens gb JG3

use Rtotal
J2=(trim minors(3,S)+stats);
netList J2_*
JG2=trim sub(J2,Relim);
netList JG2_*
netList flatten entries gens gb JG2


-----------------------------------------------------------
--complete 4-graph with 2 vertices equal
-----------------------------------------------------------
restart
load "ColoredGraphicalModels_OK_RH.m2"
R=QQ[l_1..l_9]
K=matrix{{l_1,l_4,l_5,l_6},{l_4,l_1,l_7,l_8},{l_5,l_7,l_2,l_9},{l_6,l_8,l_9,l_3}}
(p,n,Rtotal,S)=coloredData(K)
stats=sub(suffStat(K),Rtotal)
J3=ideal{det S}+stats;
netList J3_*

Relim=QQ[gens Rtotal,MonomialOrder => Eliminate 10]
JG3=trim sub(J3,Relim);
netList JG3_*
netList flatten entries gens gb JG3

use Rtotal
J2=(trim minors(3,S)+stats);
netList J2_*
JG2=trim sub(J2,Relim);
netList JG2_*
netList flatten entries gens gb JG2


-----------------------------------------------------------
--complete 4-graph with 2 edges equal
-----------------------------------------------------------
restart
load "ColoredGraphicalModels_OK_RH.m2"
R=QQ[l_1..l_9]
K=matrix{{l_1,l_5,l_5,l_6},{l_5,l_2,l_7,l_8},{l_5,l_7,l_3,l_9},{l_6,l_8,l_9,l_4}}
(p,n,Rtotal,S)=coloredData(K)
stats=sub(suffStat(K),Rtotal)
J3=ideal{det S}+stats;
netList J3_*

Relim=QQ[gens Rtotal,MonomialOrder => Eliminate 10]
JG3=trim sub(J3,Relim);
netList JG3_*
netList flatten entries gens gb JG3

use Rtotal
J2=(trim minors(3,S)+stats);
netList J2_*
JG2=trim sub(J2,Relim);
netList JG2_*
netList flatten entries gens gb JG2
