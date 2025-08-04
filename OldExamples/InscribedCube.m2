restart
load "../SplinesCode.m2"

V = {{0,0,-1}, {-1,-1,1}, {-1,1,1}, {1,-1,1}, {1,1,1}, {-1/2, -1/2, 1}, {-1/2, 1/2, 1}, {1/2, -1/2, 1}, {1/2, 1/2, 1}}
F = {{0,1,2},{0,2,4},{0,3,4},{0,1,3},{5,6,7,8}, {1,2,5,6}, {2,4,6,7}, {3,4,7,8}, {1,3,5,8}}

Sigma = polyhedralComplex(V,F)
B = billeraComplex(Sigma, 1)
Splines = minimalPresentation HH_2 B
hilbertPolynomial Splines -- 14*P_0 - 20*P_1 + 9*P_2
reduceHilbert hilbertSeries Splines -- T + 5T^2 + 3*T^3

X = normalToricVariety(V,F)
alexDual = dual monomialIdeal(ideal X) 
S = simplicialComplex(alexDual) 
fVector(S)

-----------------
V = {{0,0,-1}, {-1,-1,1}, {-1,1,1}, {1,-1,1}, {1,1,1}, {-1/2, -1/2, 1}, {-1/2, 1/2, 1}, {1/2, 50/101, 1}, {1/2, -1/2, 1}}
F = {{0,1,2},{0,2,4},{0,3,4},{0,1,3},{5,6,7,8}, {1,2,5,6}, {2,4,6,7}, {3,4,7,8}, {1,3,5,8}}

Sigma = polyhedralComplex(V,F)
B = billeraComplex(Sigma, 1)
Splines = minimalPresentation HH_2 B
hilbertPolynomial Splines -- 9*P_0 - 16*P_1 + 9 * P_2
reduceHilbert hilbertSeries Splines -- 1+T+6T^2+T^3

X = normalToricVariety(V,F)
alexDual = dual monomialIdeal(ideal X) 
S = simplicialComplex(alexDual) 


C = normalToricVariety({{1,0}, {1,2}}, {{0,1}})
isSmooth C
degrees ring C
