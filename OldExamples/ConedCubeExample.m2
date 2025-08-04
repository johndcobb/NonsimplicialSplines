restart
load "../SplinesCode.m2"
needsPackage "AlgebraicSplines"
debug needsPackage "Complexes"

V = {{-1,-1,1}, {-1,1,1}, {1,-1,1}, {1,1,1}}
F = {{0,1}, {0,2}, {1,3}, {2,3}, {0,1,2,3}}

Sigma = polyhedralComplex(V,F)
Sigmafaces = applyPairs(faces Sigma, (i,lst) -> ((dim Sigma)-i-1,apply(lst,first))) 
-*
                0 => {{0}, {1}, {2}, {3}}
               1 => {{0, 1}, {0, 2}, {1, 3}, {2, 3}}
               2 => {{0, 1, 2, 3}}
*-
isSimplicial(Sigma) -- true

B = billeraComplex(V,F,1)
Splines = minimalPresentation HH_2 B
reduceHilbert hilbertSeries Splines 

X = normalToricVariety(V,F)
isSimplicial(X) -- false

alexDual = dual monomialIdeal(ideal X) -- == (1)
S = simplicialComplex(alexDual)


-------------------------------
V = {{0,0,0},{-1,-1,1}, {-1,1,1}, {1,-1,1}, {1,1,1}}
F = {{0,1,2}, {0,1,3}, {0,2,4}, {0,3,4}, {1,2,3,4}} -- if I remove the zero from the last face, then simplicial complex has the two diagonals as the cones...



Sigma = polyhedralComplex(V,F)
Sigmafaces = applyPairs(faces Sigma, (i,lst) -> ((dim Sigma)-i-1,apply(lst,first)))
-* 
                0 => {{0}, {1}, {2}, {3}, {4}}
                1 => {{0, 1}, {0, 2}, {0, 3}, {0, 4}, {1, 2}, {1, 3}, {2, 4}, {3, 4}}
                2 => {{0, 1, 2}, {0, 1, 3}, {0, 2, 4}, {0, 3, 4}, {1, 2, 3, 4}}
*-
isSimplicial(Sigma) -- true

X = normalToricVariety(entries rays fan(Sigma),F)
isSimplicial(X) -- false
B = billeraComplex(V,F,1)
Splines = minimalPresentation HH_2 B
reduceHilbert hilbertSeries Splines 

alexDual = dual monomialIdeal(ideal X) 
S = simplicialComplex(alexDual) -- (x_1*x_4, x_2*x_3)



-------------------------
V = {{0,0,-1}, {-1,-1,1}, {-1,1,1}, {1,-1,1}, {1,1,1}}
F = {{0,1,2},{0,2,4},{0,3,4},{0,1,3},{1,2,3,4}}

Sigma = polyhedralComplex(V,F)
B = billeraComplex(Sigma, 1)
Splines = minimalPresentation HH_2 B
hilbertPolynomial Splines -- 5*P_0 - 8*P_1 + 5*P_2
reduceHilbert hilbertSeries Splines -- 1 + T + 2T^2 + T^3



X = normalToricVariety(V,F)
alexDual = dual monomialIdeal(ideal X) 
S = simplicialComplex(alexDual) 
fVector(S)
SRdual = minimalPresentation(ring(alexDual)/alexDual)

----------
needsPackage "AlgebraicSplines"
V = {{0,0,-1}, {-1,-1,1}, {1,-1,1}, {-1,1,1}, {1,1,1}, {0,0,0}}
F = {{1,2,3,4,5},{0,1,2,5},{0,1,3,5},{0,2,4,5},{0,3,4,5}}
phi = stanleyReisner(V,F, Homogenize => false)
ker phi
minimalPresentation splineModule(V,F,0)
stanleyReisner(V,F)
