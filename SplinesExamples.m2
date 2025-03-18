restart
load "SplinesCode.m2"

V = {{0,0,-1}, {-1,-1,1}, {1,-1,1},  {1,1,1}, {-1,1,1}};
F = {{0,1,2}, {0,2,3}, {0,3,4}, {0,1,4}, {1,2,3,4}};
r = 1


X = normalToricVariety(V,F);
isSimplicial X --false
isComplete X --true
B = ideal X
BV = dual monomialIdeal(ideal X) --When X is simplicial, this is the Stanley-Riesner ring of S. Its still the Stanley Reisner ring of SOMETHING.
R = ring BV
simplicialComplex(BV) -- the simplicial complex corresponding to the stanley reisner ideal 
vars R
linForms = ideal flatten entries ((vars R)*(matrix V))
I = BV + linForms
S = R/I
minimalPresentation(S)
reduceHilbert hilbertSeries S -- should be artinian and gorenstein, which is visible in the hilbert series
-- even if its the cohomology ring of a blowup, it should be simpler.

C = billeraComplex(V,F,r)

minimalPresentation HH_2 C

-- smooth this thing by subdividing the fan. 


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