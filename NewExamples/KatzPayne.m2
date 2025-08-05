load "../ChowMap.m2";

V = {{1,1},{1,-1}, {-1,-1}, {-1,1}, {0,0}};
F = {{0,1,4},{1,2,4},{2,3,4}, {3,0,4}};
Splines = splineModule(V,F,0, Homogenize => false)
minimalPresentation Splines -- R + R(-1)^2 + R(-2)
R = ring Splines;


entries Splines_3
f = spline(entries Splines_3, V, F);
describe f
origin = coneFromVData transpose matrix{{0,0}}
-- The following does the same thing as in Example 4.1 from [Katz,Payne]
chowMap(f, origin)

-- we can also compute the equivariant multiplicity directly:
sigma1 = coneFromVData transpose matrix{{1,1},{1,-1}};
equivariantMultiplicity(sigma1, origin, R)
assert(equivariantMultiplicity(sigma1, origin, R)*restriction(f,sigma1) == chowMap(f, origin))

-- Now, on a degree 1 spline, the chowMap should give integers on the codim 1 cones. That is, the rays of tha fan.
f = spline(entries Splines_1, V, F)
describe f
v1 = coneFromVData transpose matrix{{1,1}}
chowMap(f, v1)
equivariantMultiplicity(sigma1, v1, R)

v2 = coneFromVData transpose matrix{{1,-1}}
chowMap(f, v2)


polarFace(polyhedron v1, polyhedron sigma1)

sigma1 = coneFromVData transpose matrix{{1,1},{1,-1}}
unimodularTriangulation1 = {coneFromVData transpose matrix{{1,1},{1,0}}, coneFromVData transpose matrix{{1,0},{1,-1}}}

equivariantMultiplicity(sigma1, v2, R)
equivariantMultiplicity(unimodularTriangulation1, v2, R)

sigma2 = coneFromVData transpose matrix{{1,-1},{-1,-1}}
unimodularTriangulation2 = {coneFromVData transpose matrix{{1,-1},{0,-1}}, coneFromVData transpose matrix{{0,-1},{-1,-1}}}

sigma3 = coneFromVData transpose matrix{{-1,-1},{-1,1}}
unimodularTriangulation3 = {coneFromVData transpose matrix{{-1,-1},{-1,0}}, coneFromVData transpose matrix{{-1,0},{-1,1}}}

sigma4 = coneFromVData transpose matrix{{-1,1},{1,1}}
unimodularTriangulation4 = {coneFromVData transpose matrix{{-1,1},{0,1}}, coneFromVData transpose matrix{{0,1},{1,1}}}

equivariantMultiplicity(unimodularTriangulation1, origin, R)
equivariantMultiplicity(unimodularTriangulation2, origin, R)
assert(equivariantMultiplicity(unimodularTriangulation2, origin, R) == equivariantMultiplicity(sigma2, origin, R))
equivariantMultiplicity(unimodularTriangulation4, origin, R)


tau = coneFromVData transpose matrix{{1,1}}
equivariantMultiplicity(unimodularTriangulation1, tau, R)
equivariantMultiplicity(unimodularTriangulation2, tau, R)
equivariantMultiplicity(unimodularTriangulation3, tau, R)
equivariantMultiplicity(unimodularTriangulation4, tau, R)
assert(equivariantMultiplicity(unimodularTriangulation2, tau, R) == equivariantMultiplicity(sigma2, tau, R))