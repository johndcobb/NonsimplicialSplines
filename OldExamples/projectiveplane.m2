restart
load "../SplinesCode.m2"
needsPackage "AlgebraicSplines"

V = {{1,0}, {0,1}, {-1,-1}}
F = {{0,1},{1,2}, {2,0}}

R = QQ[x,y]
Sigma = polyhedralComplex(V,F)
B = billeraComplex(Sigma,R, 1)
Splines = minimalPresentation HH_1 B

needsPackage "AlgebraicSplines"
--- use the algebraic splines code and kill off the extra variables.
V = {{0,0}, {1,0}, {0,1}, {-1,1}}
F = {{0,1,2}, {0,2,3}, {0,1,3}}

phi = stanleyReisner(V,F)
splinesModule(V,F,0)


-- In order to get map from R[x_0,x_1,x_2] to the splines module R[x,y] + R[x,y](-1) + R[x,y](-2)
-- first, add an all 1's row to the vertices.
V = {{0,1,1},{1,0,1}, {-1,-1,1}}
varlist = {x,y,1_R}
M = transpose matrix(R,V)

S = QQ[x,y,p_0,p_1,p_2]

phi = matrix({{y, 0, x}, {0, -y, x-y}, {y-x,-x,0}})*matrix({{p_0},{p_1},{p_2}})


linForms = ideal flatten entries (vars R * vertices Sigma)

minimalPresentation(Splines/linForms)

-- P^2 is smooth so, this should be giving the cohomology ring k[x]/(x^3), but it is not. 