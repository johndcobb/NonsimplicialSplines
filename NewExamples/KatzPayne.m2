load "../ChowMap.m2";

-- This is Example 4.1 from [Katz,Payne]
V = {{1,1},{1,-1}, {-1,-1}, {-1,1}, {0,0}};
F = {{0,1,4},{1,2,4},{2,3,4}, {3,0,4}};
Splines = splineModule(V,F,0, Homogenize => false)
minimalPresentation Splines -- R + R(-1)^2 + R(-2)

f = Splines_3 

sigma1 = coneFromVData transpose matrix{{1,1},{1,-1}}
unimodularTriangulation1 = {coneFromVData transpose matrix{{1,1},{1,0}}, coneFromVData transpose matrix{{1,0},{1,-1}}}

sigma2 = coneFromVData transpose matrix{{1,-1},{-1,-1}}
unimodularTriangulation2 = {coneFromVData transpose matrix{{1,-1},{0,-1}}, coneFromVData transpose matrix{{0,-1},{-1,-1}}}

sigma3 = coneFromVData transpose matrix{{-1,-1},{-1,1}}
unimodularTriangulation3 = {coneFromVData transpose matrix{{-1,-1},{-1,0}}, coneFromVData transpose matrix{{-1,0},{-1,1}}}

sigma4 = coneFromVData transpose matrix{{-1,1},{1,1}}
unimodularTriangulation4 = {coneFromVData transpose matrix{{-1,1},{0,1}}, coneFromVData transpose matrix{{0,1},{1,1}}}

R = ZZ/32003[x,y]
origin = coneFromVData transpose matrix{{0,0}}
equivariantMultiplicity(unimodularTriangulation1, origin, R)
equivariantMultiplicity(unimodularTriangulation2, origin, R)
equivariantMultiplicity(unimodularTriangulation3, origin, R)
equivariantMultiplicity(unimodularTriangulation4, origin, R)


tau = coneFromVData transpose matrix{{1,1}}
equivariantMultiplicity(unimodularTriangulation1, tau, R)
equivariantMultiplicity(unimodularTriangulation2, tau, R)
equivariantMultiplicity(unimodularTriangulation3, tau, R)
equivariantMultiplicity(unimodularTriangulation4, tau, R)