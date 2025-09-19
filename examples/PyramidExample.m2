load "../code/ChowMap.m2";

V = {{0,0,-1}, {-1,-1,1}, {1,-1,1},  {-1,1,1}, {1,1,1}};
F = {{0,1,2}, {0,1,3}, {0,2,4}, {0,3,4}, {1,2,3,4}};
Sigma = fan(V,F);
Splines = splineModule(Sigma,0, Homogenize => false)
R = ring Splines;


linearSpline = spline(entries Splines_1,Sigma, R);
describe linearSpline
c = chowMap(linearSpline)
tau = coneFromVData transpose matrix{{1,-1,1}, {0,0,-1}}
weight(c, tau)

sigma1 = coneFromVData transpose matrix{{0,0,-1}, {1,-1,1}, {1,1,1}}
equivariantMultiplicity(sigma1, tau, R)
sigma2 = coneFromVData transpose matrix{{0,0,-1}, {1,-1,1}, {-1,-1,1}}
equivariantMultiplicity(sigma1, tau, R)
isUnimodular sigma1
unimodularTriangulation = findUnimodularTriangulation(sigma1)
unimodularTriangulation / rays
trianglesContainingTau := select(unimodularTriangulation, triangle -> (rays intersect(triangle, tau)) != map(ZZ^(ambDim tau),0, 0));



splines = splineList(Splines, Sigma, R)
chowMap(splines, 1)
operationalChowGroup(Sigma, 1)
chowMap(splines,2)
chowMap(splines, 3)

sigma1 = coneFromVData transpose matrix{{0,0,-1}, {1,0,1}, {1,1,1}}
sigma2 = coneFromVData transpose matrix{{0,0,-1}, {1,-1,1}, {1,0,1}}
rays dualCone sigma1
rays dualCone sigma2


V = {{-1,-3}, {-1,0},{0,1}, {1,0}}
F = {{0,1}, {1,2}, {2,3}, {3,0}}
Sigma = fan(V,F);
Splines = splineModule(Sigma,0, Homogenize => false)
R = ring Splines;

splineObjs = splineList(Splines, Sigma, R);
chowMap(splineObjs, 1)
facesAsCones(1,Sigma) / rays

entries Splines_2
operationalChowGroup(Sigma, 1)

linearSpline = spline({t_2, 0, 0, t_2}, Sigma, R)
c = chowMap(linearSpline)
mat c



sigma3 = coneFromVData transpose matrix{{1,0,1}, {0,0,1}, {1,1,1}}
sigma4 = coneFromVData transpose matrix{{1,0,1}, {0,0,1}, {1,-1,1}}
tau = coneFromVData transpose matrix{{1,-1,1}, {1,1,1}}

equivariantMultiplicity(sigma1, tau, R)
equivariantMultiplicity(sigma2, tau, R)
equivariantMultiplicity(sigma3, tau, R)
equivariantMultiplicity(sigma4, tau, R)



-- Fulton Example
V = {{1,1,-1},{1,-1,-1},{-1,1,-1},{-1,-1,-1},{1,2,3},{1,-1,1},{-1,1,1},{-1,-1,1}};
F = {{4,5,7,6},{0,1,3,2},{0,2,6,4},{0,1,5,4},{1,3,7,5},{3,2,6,7}};
Sigma = fan(V,F);
Splines = splineModule(Sigma,0, Homogenize => false)
R = ring Splines;

splines = splineList(Splines, Sigma, R)
rank operationalChowGroup(Sigma, 0) -- 1
c1 = chowMap(splines, 1)
prune coker transpose chowGroup(Sigma, 1) -- matches Example 1.3 in Fulton-Sturmfels.
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

