load"AidaMLTcolored.m2"

--example 1 
n=3 -- number of vertices
m=2 -- m\times m minors
V={{1,3},{2}} -- vertex class colors 
E={{{1,2},{2,3}},{{1,3}}} -- edge class colors
R=createRing(n,V,E)
createSufIdeal(R,m,V,E) -- ideal of sufficient statistics for rank m-1 matrices. If this is zero, then MLT=m-1. 


--example 2
n=3
m=2
V={{1,2,3}}
E={{{1,2},{2,3}},{{1,3}}} 
R=createRing(n,V,E)
createSufIdeal(R,m,V,E) 

--example 3
n=3 
m=2 
V={{1},{2},{3}} 
E={{{1,2},{2,3}},{{1,3}}} 
R=createRing(n,V,E)
createSufIdeal(R,m,V,E)


--example 4
n=3 
m=2
V={{1},{2},{3}} 
E={{{1,2},{2,3},{1,3}}} 
R=createRing(n,V,E)
I=createSufIdeal(R,m,V,E)


