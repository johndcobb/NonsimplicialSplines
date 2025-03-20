-*---------------------------------------

-- Do stellar subdivision of the cube, get simplicial thing. Compare to the non-simplicial cube.
-- Hope: Get map between them. This should be P^1 x P^1 x P^1, blown up at the 8 T-fixed points.


*----------------------------------------

restart
load "SplinesCode.m2"

V = {{0,0,-1}, {-1,-1,1}, {1,-1,1},  {-1,1,1}, {1,1,1}}; -- this defines the vertices 
-- Note that regardless of how you arrange V, the ordering upon creating the polyhedral complex will follow reverse lexicographic order.
F = {{0,2,4}, {0,3,4}, {0,1,3}, {0,1,2}, {1,2,3,4}}; -- this defines the maximal cones
Sigma = polyhedralComplex(V,F)
B = billeraComplex(V,F,1)
vertices Sigma
maxPolyhedra Sigma
Vcells = {{},{},{},{},{}} / newCell
--Fcells = apply(F, f -> newCell(apply(f, v -> Vcells_v)))

R = QQ[x_0,x_1,x_2]

cellD = cellComplex(R, Sigma) 
boundaries = (maxCells cellD)#2 / boundary
complexD = chainComplex(cellD)
boundaryMap(2, cellD)

cellC = cellComplex(R, Sigma)
complexC = chainComplex(cellC)
 -- cellC and cellD have different 

boundaryMap(2, cellD)
boundaryMap(2, cellC)


X = normalToricVariety(V,F);
isSimplicial X --false
isComplete X --true
B = ideal X
BV = dual monomialIdeal(ideal X) --When X is simplicial, this is the Stanley-Riesner ring of S. Its still the Stanley Reisner ring of SOMETHING.
R = ring BV
simplicialComplex(BV) -- the simplicial complex corresponding to the stanley reisner ideal 
linForms = ideal flatten entries ((vars R)*(matrix V)) -- this is what we should mod out by to get cohomology ring from equivariant cohomology ring
I = BV + linForms
S = R/I
minimalPresentation(S)
reduceHilbert hilbertSeries S -- should be artinian and gorenstein, which is visible in the hilbert series
-- even if its the cohomology ring of a blowup, it should be simpler.

R = QQ[x_0,x_1,x_2]
SigmaD = polyhedralComplex(V,F)
SigmaC = polyhedralComplex(V,F)

cellD = cellComplex(R, SigmaD) 
cellC = cellComplex(R, SigmaC)

D = billeraComplex(V,F,R, 1)
C = billeraComplex(V,F,R,1)

D == C

reduceHilbert (hilbertSeries (HH_2 D)) == reduceHilbert hilbertSeries HH_2 C

-- smooth this thing by subdividing the fan. 
-----------------------
-- Look at the cube centered at origin, 6 maximal cones for each face.

restart 
load "SplinesCode.m2"
V = {{1,1,-1},{1,-1,-1},{-1,1,-1},{-1,-1,-1},{1,1,1},{1,-1,1},{-1,1,1},{-1,-1,1}};
F = {{0,1,3,2},{0,2,6,4},{0,1,5,4},{1,3,7,5},{3,2,6,7},{4,5,7,6}};
r = 1

vertices polyhedralComplex(V,F)

d = #(V_0)
R = QQ[x_0..x_(d-1)]
C = billeraComplex(V,F,R,r)
Splines = minimalPresentation HH_2 C
betti Splines
L = ideal random(R^3, R^{-1}) -- three random linear forms
reduceHilbert hilbertSeries(Splines ** R/L) -- artinian reduction. 

X = normalToricVariety(V,F);
B = ideal X
EChow = minimalPresentation HH_2 C -- this is the equivariant cohomology ring by Payne 
linForms = ideal flatten entries ((vars ring B )*(matrix V)) 


hilbertPolynomial HH_2 C
hilbertPolynomial HH_1 C
hilbertPolynomial HH_0 C

B = ideal X
BV = dual monomialIdeal(ideal X) --When X is simplicial, this is the Stanley-Riesner ring of S. Its still the Stanley Reisner ring of SOMETHING.
R = ring BV
simplicialComplex(BV) -- the simplicial complex corresponding to the stanley reisner ideal 

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


Vmat := vertices Sigma;
VList := for numcol in 0..(numColumns(Vmat)-1) list Vmat_numcol;
for M in L list (
      for numcol in 0..numColumns(M)-1 list (
      position(VList, v -> lift(v, ZZ) == M_numcol )
      )
)