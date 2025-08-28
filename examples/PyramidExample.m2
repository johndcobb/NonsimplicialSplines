load "../ChowMap.m2";

V = {{0,0,-1}, {-1,-1,1}, {1,-1,1},  {-1,1,1}, {1,1,1}};
F = {{0,1,2}, {0,1,3}, {0,2,4}, {0,3,4}, {1,2,3,4}};
Sigma = fan(V,F);
Splines = splineModule(Sigma,0, Homogenize => false)
R = ring Splines;

splines = splineList(Splines, V, F, R)
-*
c0 = chowMap(splines, 0) -- 0
operationalChowGroup(Sigma, 0) --- image (1,1,1,1,1)
c1 = chowMap(splines, 1) -- (1,1,1,1,0,0,0,0)
chowGroup(Sigma, 1)
Hom(coker transpose chowGroup(Sigma, 1),ZZ)
operationalChowGroup(Sigma, 1)

prune coker transpose chowGroup(Sigma, 1)
Hom(prune cokernel chowGroup(Sigma, 1), ZZ)

operationalChowGroup(Sigma, 1) --- image (1,1,1,1,1,1,1)
c2 = chowMap(splines, 2) 
chowGroup(Sigma, 2)
operationalChowGroup(Sigma, 2)
c3 = chowMap(splines, 3) -- ZZ --> ZZ given by multiplication by 2
operationalChowGroup(Sigma,3) -- ZZ
*-




-- Fulton Example
V = {{1,1,-1},{1,-1,-1},{-1,1,-1},{-1,-1,-1},{1,1,1},{1,-1,1},{-1,1,1},{-1,-1,1}};
F = {{4,5,7,6},{0,1,3,2},{0,2,6,4},{0,1,5,4},{1,3,7,5},{3,2,6,7}};
Sigma = fan(V,F);
Splines = splineModule(Sigma,0, Homogenize => false)
R = ring Splines;

splines = splineList(Splines, V, F, R)
rank operationalChowGroup(Sigma, 0) -- 1
c1 = chowMap(splines, 1)
operationalChowGroup(Sigma, 1) -- rank 0

--- my chow map maps from A^k(X)/M*A^(k-1)(X)
c2 = chowMap(splines, 2)
operationalChowGroup(Sigma, 2)

ker c2
prune coker c2

operationalChowGroup(Sigma, 2) -- rank 5
c3 = chowMap(splines, 3)
operationalChowGroup(Sigma, 3) -- rank 1


needsPackage "gfanInterface"
V = {{-1,-1,-1}, {1,0,0}, {0,1,0}, {0,0,1}};

maxCones gfanSecondaryFan(V)

