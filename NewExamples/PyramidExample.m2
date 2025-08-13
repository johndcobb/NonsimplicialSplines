load "../ChowMap.m2";

V = {{0,0,-1}, {-1,-1,1}, {1,-1,1},  {-1,1,1}, {1,1,1}, {0,0,0}};
F = {{0,1,2,5}, {0,1,3,5}, {0,2,4,5}, {0,3,4,5}, {1,2,3,4,5}};
Splines = splineModule(V,F,0, Homogenize => false)
R = ring Splines; 
minimalPresentation Splines -- Spline module is free: R + R(-1) + R^(-2)^2 + R(-3)

Sigma = fan(V,F)
sigmaCones = facesAsCones(0,Sigma)
sigmaCones / rays
topCone = sigmaCones_4
liftedCone = simplicialLift(topCone)
rays liftedCone
SigmaLift = simplicialLift(Sigma)
facesAsCones(0,SigmaLift) / rays

maxCones(SigmaLift)
V = entries transpose rays SigmaLift
(for maxCone in maxCones(SigmaLift) list (V_maxCone)) / transpose / matrix / coneFromVData

VBad = entries transpose rays Sigma
Fbad = maxCones(Sigma)

isSimplicial topCone
disc = numRays topCone - dim topCone
heightFunction = table(numRays topCone, disc, (i,j) -> 0)
liftedCone = coneFromVData(rays topCone || transpose matrix heightFunction)
V := entries transpose rays Sigma;
F := maxCones(Sigma);
maxDiscrepancy = max apply(F, face -> length face - dim Sigma)
heightFunction = table(length V, maxDiscrepancy, (i,j) -> random(maxDiscrepancy+2)) -- the max height is chosen mostly arbitrarily
liftedV := entries ( matrix V |  matrix heightFunction);

splines = splineList(Splines, V, F, R)
c0 = chowMap(splines, 0)
c1 = chowMap(splines, 1)
c2 = chowMap(splines, 2)
c3 = chowMap(splines, 3)

ker c0
ker c1

origin = coneFromVData transpose matrix{{0,0,0}}
f1 = splines_1
chowMap(f1, origin)