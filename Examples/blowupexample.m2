restart
load "../SplinesCode.m2"

V = {{0,0,-1}, {-1,-1,1}, {1,-1,1},  {-1,1,1}, {1,1,1}};
F = {{0,1,2}, {0,1,3}, {0,2,4}, {0,3,4}, {1,2,3,4}};

Sigma = polyhedralComplex(V,F)
Pfaces = applyPairs(faces Sigma, (i,lst) -> ((dim Sigma)-i-1,apply(lst,first)))
R = QQ[x_0,x_1,x_2]

B = billeraComplex(Sigma, R, 1)
Splines = minimalPresentation HH_2 B
reduceHilbert hilbertSeries Splines

X = normalToricVariety(V,F)
BV = dual monomialIdeal(ideal X) -- (x_0*x_2*x_3, x_0*x_1*x_4)
simplicialComplex(BV)
stanReisner = coker gens BV


---------------------------------------

Fsimplicial = flatten(listToPolyhedra(maxPolyhedra(Sigma),Sigma) / barycentricTriangulation) / convexHull
SigmaSimplicial = polyhedralComplex(Fsimplicial)
SigmaSimplicialFaces = applyPairs(faces SigmaSimplicial, (i,lst) -> ((dim SigmaSimplicial)-i-1,apply(lst,first)))

BSimplicial = billeraComplex(SigmaSimplicial, R, 1)
SplinesSimplicial = minimalPresentation HH_2 BSimplicial
reduceHilbert hilbertSeries SplinesSimplicial

XSimplicial = normalToricVariety(fan(SigmaSimplicial))
BVSimplicial = dual monomialIdeal(ideal XSimplicial)
SSimplicial = simplicialComplex(BVSimplicial)
stanReisnerSimplicial = (module ring BVSimplicial)/BVSimplicial



i = map(BSimplicial_0, B_0, promote(matrix {{1,0,0,0,0}, {0,0,0,0,0}, {0,1,0,0,0}, {0,0,1,0,0}, {0,0,0,1,0}, {0,0,0,0,1}}, R) )
phi = extend(BSimplicial,B, i)
--- this fails

R = QQ[t,z,y,x, MonomialOrder => Lex]
I = ideal ( x-t, y-t^2, z-t^3)
G = gens gb I
