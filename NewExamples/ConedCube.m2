needsPackage "AlgebraicSplines"; needsPackage "NormalToricVarieties"; needsPackage "SimplicialComplexes"; needsPackage "Polyhedra";

V = {{1,1,-1},{1,-1,-1},{-1,1,-1},{-1,-1,-1},{1,1,1},{1,-1,1},{-1,1,1},{-1,-1,1},{0,0,0}};
F = {{0,1,3,2,8},{0,2,6,4,8},{0,1,5,4,8},{1,3,7,5,8},{3,2,6,7,8},{4,5,7,6,8}};

Splines = splineModule(V,F,0, Homogenize => false)
minimalPresentation Splines -- so its R + R(-1) + R(-2)^2 + R(-3) + R(-4)
reduceHilbert hilbertSeries Splines
SR = stanleyReisner(V,F)

simplicialComplex(S)


X = normalToricVariety(V,F)
alexDual = dual monomialIdeal(ideal X) 
S = simplicialComplex(alexDual) 
SRdual = minimalPresentation(ring(alexDual)/alexDual)
reduceHilbert hilbertSeries SRdual

V = {{1,1,-1},{1,-1,-1},{-1,1,-1},{-1,-1,-1},{1,1,1},{1,-1,1},{-1,1,1},{-1,-1,1}};
F = {{0,1,3,2},{0,2,6,4},{0,1,5,4},{1,3,7,5},{3,2,6,7},{4,5,7,6}};
Splines = splineModule(V,F,0, Homogenize => false)


V = {{1,1},{1,-1}, {-1,1}, {-1,-1}, {0,0}};
F = {{0,1,4}, {0,2,4}, {1,3,4}, {2,3,4}}
Splines = splineModule(V,F,0, Homogenize => false)
minimalPresentation Splines


V = {{1,1},{1,-1}, {-1,1}, {-1,-1}};
F = {{0,1}, {0,2}, {1,3}, {2,3}}
X = normalToricVariety(V,F)
alexDual = dual monomialIdeal(ideal X) 
reduceHilbert hilbertSeries minimalPresentation(ring(alexDual)/alexDual)


V = {{0,1},{1,0}, {-1,-1}, {0,0}}; -- we have our problem even for P^2!!! Something is wrong....
F = {{0,1,3}, {1,2,3}, {2,0,3}}
Splines = splineModule(V,F,0, Homogenize => false)
minimalPresentation Splines   