restart
load "SplinesCode.m2"

V = {{0,0,-1}, {-1,-1,1}, {1,-1,1},  {1,1,1}, {-1,1,1}};
F = {{0,1,2}, {0,2,3}, {0,3,4}, {0,1,4}, {1,2,3,4}};
r=2
Sigma = polyhedralComplex(V,F)
listToPolyhedra(polyhedra(1, Sigma), Sigma)

P = convexHull transpose matrix {V_1, V_4}


X = normalToricVariety(myRays, myCones);
isSimplicial X --false
isComplete X --true
B = ideal X
BV = dual monomialIdeal(ideal X) --When X is simplicial, this is the Stanley-Riesner ring of S. Its still the Stanley Reisner ring of SOMETHING.
R = ring BV
simplicialComplex(BV) -- the simplicial complex corresponding to the stanley reisner ideal 
vars R
linForms = ideal flatten entries ((vars R)*(matrix myRays))
I = BV + linForms
S = R/I
minimalPresentation(S)
reduceHilbert hilbertSeries S -- should be artinian and gorenstein, which is visible in the hilbert series
-- even if its the cohomology ring of a blowup, it should be simpler
-- how to tell if variety is smooth from determinant of rays matrix

-- smooth this thing by subdividing the fan. 


--make a funtion to construct the Billera complex

------------------------
myRays = {{0,0,-1}, {-1,-1,1}, {1,-1,1},  {1,1,1}, {-1,1,1}};
myCones = {{0,1,2}, {0,2,3}, {0,3,4}, {0,1,4}, {1,2,3},{1,3,4}};
X = normalToricVariety(myRays, myCones);
isSimplicial X --false
isComplete X --true
B = ideal X
BV = dual monomialIdeal(ideal X) --When X is simplicial, this is the Stanley-Riesner ring of S. Its still the Stanley Reisner ring of SOMETHING.
R = ring BV
simplicialComplex(BV) -- the simplicial complex corresponding to the stanley reisner ideal 
vars R
linForms = ideal flatten entries ((vars R)*(matrix myRays))
I = BV + linForms
S = R/I
minimalPresentation(S)
reduceHilbert hilbertSeries S -- sh


debug needsPackage "AlgebraicSplines"
splineMatrix(myRays, myCones, 0)

V = {{-1,0},{0,1},{1,0},{0,-1},{0,0}}
F = {{0,1,4},{1,2,4},{2,3,4},{0,3,4}}
splineMatrix(V,F,0)
splineModule(V,F,0)
stanleyReisner(V,F)
stanleyReisnerPresentation(V,F,0)

---------------------------

V = {{0,0,-1}, {-1,-1,1}, {1,-1,1},  {1,1,1}, {-1,1,1}};
F = {{0,1,2}, {0,2,3}, {0,3,4}, {0,1,4}, {1,2,3,4}}; -- i don't know why this is failing, and it may require me to reconstruct the entire thing.

V={{0,1},{-1,-1},{1,-1},{0,10},{-2,-2},{2,-2}}
F={{0,1,2},{0,1,3,4},{0,2,3,5},{1,2,4,5}};
E = getCodim1Intersections(F)
r = 0
d := #(first V);
V = apply(V, v-> append(v,1))
S := createSplineRing(d,opts)
varlist = (vars S)|(matrix {{sub(1,S)}})
varCol := transpose varlist
M := (transpose(matrix(S,V)))
mM := numrows M
minorList := apply(E, e-> gens gb minors(mM,matrix(M_e)|varCol))
if any(minorList, I-> ideal I === ideal 1) then (
    error "Some vertices on entered face are not in codimension 1 face."
    );
flatten apply(minorList, m -> (m_(0,0))^(r+1))

---------------------


myRays = {{0,0,-1},{0,-4,1}, {3,3,1},{3,-3,1}, {0,-1,1}, {1,1,1},{1,-1,1}};
myCones = {{1,2,4,5},{2,3,5,6},{1,3,4,6},{4,5,6}, {0,1,2}, {0,1,3}, {0,2,3}}
X = normalToricVariety(myRays, myCones);
isSimplicial X --false
isComplete X
B = ideal X
BV = dual monomialIdeal(ideal X)
simplicialComplex(BV)

