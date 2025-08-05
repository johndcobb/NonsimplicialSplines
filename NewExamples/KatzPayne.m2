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
sigma1 = coneFromVData transpose matrix{{1,1},{1,-1}}
equivariantMultiplicity(sigma1, origin, R)

f = spline(entries Splines_1, V, F)
describe f

-- Now, since f is degree 1, the chowMap should give integers on the codim 1 cones. That is, the rays of tha fan.
v1 = coneFromVData transpose matrix{{1,1}}
chowMap(f, v1)

-- but it doesnt.



sigma1 = coneFromVData transpose matrix{{1,1},{1,-1}}
unimodularTriangulation1 = {coneFromVData transpose matrix{{1,1},{1,0}}, coneFromVData transpose matrix{{1,0},{1,-1}}}

sigma2 = coneFromVData transpose matrix{{1,-1},{-1,-1}}
unimodularTriangulation2 = {coneFromVData transpose matrix{{1,-1},{0,-1}}, coneFromVData transpose matrix{{0,-1},{-1,-1}}}

sigma3 = coneFromVData transpose matrix{{-1,-1},{-1,1}}
unimodularTriangulation3 = {coneFromVData transpose matrix{{-1,-1},{-1,0}}, coneFromVData transpose matrix{{-1,0},{-1,1}}}

sigma4 = coneFromVData transpose matrix{{-1,1},{1,1}}
unimodularTriangulation4 = {coneFromVData transpose matrix{{-1,1},{0,1}}, coneFromVData transpose matrix{{0,1},{1,1}}}

origin = coneFromVData transpose matrix{{0,0}}
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



R = QQ[x,y,a,b,c,d]
f = x^3 + a * x^2 + b*x*y +  c
g = x - y + d
J = det (jacobian(ideal(f,g)))^{0,1}
I = ideal(f,g,J)
decompose I

f = x^2 + a*x + b
J = det (jacobian(ideal(f)))^{0}
I = ideal(f,J)
decompose I

radical(I) == I
decompose I
associatedPrimes I
toString decompose ideal(f,g,J)

h = eliminate({x,y},I)
decompose h
toString h
degree h
factor ((gens h)_0)_0
radical h == h