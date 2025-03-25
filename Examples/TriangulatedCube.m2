restart 
load "SplinesCode.m2"

R = QQ[x_0,x_1,x_2];

V = {{1,1,-1},{1,-1,-1},{-1,1,-1},{-1,-1,-1},{1,1,1},{1,-1,1},{-1,1,1},{-1,-1,1}};
F = {{0,1,3,2},{0,2,6,4},{0,1,5,4},{1,3,7,5},{3,2,6,7},{4,5,7,6}};

Sigma = polyhedralComplex(V,F)
B = billeraComplex(Sigma, R, 1)
Splines = minimalPresentation HH_2 B
reduceHilbert hilbertSeries Splines  -- 1 + T + 2T^2 + T^3 + T^4

X = normalToricVariety(V,F);
BV = dual monomialIdeal(ideal X) 
S = simplicialComplex(BV) -- the simplicial complex corresponding to the stanley reisner ideal 

linForms = ideal flatten entries ((vars ring S)*(matrix V))
minimalPresentation ((ring S)/(BV + linForms))


-*-------------------------------
- Now, to take a barycentric subdivision of each face
*--------------------------------
VT = {{1,1,-1},{1,-1,-1},{-1,1,-1},{-1,-1,-1},{1,1,1},{1,-1,1},{-1,1,1},{-1,-1,1},
     {0,0,-1}, {1,0,0}, {0,1,0}, {0,-1,0}, {-1,0,0}, {0,0, 1}}
FT = {{4,6,13},{6,7,13},{7,5,13}, {5,4,13},
     {3,7,12},{7,6,12},{6,2,12},{2,3,12},
     {7,3,11},{3,1,11},{1,5,11},{5,7,11},
     {0,2,10},{4,0,10},{2,6,10},{6,4,10},
     {0,4,9},{4,5,9},{5,1,9},{1,0,9},
     {0,1,8},{1,3,8},{3,2,8},{2,0,8}} 

SigmaT = polyhedralComplex(VT,FT)
BT = billeraComplex(SigmaT, R,  1)
SplinesT = minimalPresentation HH_2 BT --This is equivariant cohomology ring due to payne

reduceHilbert hilbertSeries SplinesT -- 1 + 11T + 11T^2 + T^3

