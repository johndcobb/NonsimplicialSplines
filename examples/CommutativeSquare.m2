load "../code/ChowMap.m2"
------------------------------------------------------
-* Explanation
Let X' --> X be a simplicialization of a complete toric variety X. We have the following maps:
(1) From Katz-Payne, we have a map (chowMap) i^k : A_T^k(X) --> A^k(X) and A_T^k(X') --> A^k(X')
(2) from Proposition 2.7 in Fulton-Sturmfels: A^k(X) --> A^k(X')
(3) We have a trivial inclusion of splines A_T^k(X) --> A_T^k(X').

Maps (1-4) fit into a commutative square
A_T^k(X)  -----------------> A_T^k(X')
   |                             |
   |                             |
   v                             v
A^k(X) --------------------> A^k(X')

I want to compute each of these maps....
*- 
------------------------------------------------------



------------------------------------------------------
-*
The following is a pyramid with a hexagonal base
*-
------------------------------------------------------
V={{0,0,-1},{-1,-1,1},{0,-2,1},{1,-1,1},{1,2,3/2}, {0,2,1},{-1,1,1}}
F={{1,2,3,4,5,6},{0,1,2},{0,2,3},{0,3,4},{0,4,5}, {0,5,6},{0,6,1}}
fileName = "hexagonalBase.m2"
fileName' = "hexagonalBaseSimplicial.m2" 


V = {{0,0,-1}, {-2,-1,1}, {-1,-2,1}, {1,-2,1}, {2,-1,1}, {2,1,1}, {1,2,1}, {-1,2,1}, {-2,1,1}}
F = {{1,2,3,4,5,6,7,8}, {0,1,2}, {0,2,3}, {0,3,4}, {0,4,5}, {0,5,6}, {0,6,7}, {0,7,8}, {0,8,1}}
fileName = "octagonalBase.m2"

V = {{-1/3,-1/3,-1/3}, {1,0,0}, {0,1,0}, {0,0,1}, {2/3,1/6,1/6}, {1/6, 2/3, 1/6}, {1/6, 1/6, 2/3}}
F = {{1,2,4,5}, {2,3,5,6}, {3,1,4,6}, {4,5,6}, {0,1,2}, {0,2,3}, {0,3,1}}
--- hal's homology triangle


--- 3D homology tetrahedron
V = {{0,0,-3}, {-3,-2,2}, {3,-2,2}, {0,3,2}, {-1,-1,2}, {1,-1,2}, {0,1,2}, {0,-3,-1}, {1,-3,1}, {-1,-3,1}, {-3,1,-1}, {-4,0,1}, {-3,3,1},{3,1,-1}, {3,3,1}, {4,0,1}}
F = {{1,2,4,5}, {2,3,5,6}, {3,6,1,4}, {4,5,6}, {0,1,7,9}, {0,2,7,8}, {1,2,8,9}, {7,8,9}, {0,1,10,11}, {1,3,11,12}, {0,3,10,12}, {10,11,12}, {0,3,13,14}, {3,2,14,15}, {0,2,13,15}, {13,14,15}}


V = {{0,0,-1}, {-2,-2,1}, {2,-2,1}, {2,2,1}, {-2,2,1}, {-1,-1,1}, {1,-1,1}, {1,1,1}, {-1,1,1}}
F = {{1,2,5,6}, {2,3,6,7}, {3,4,7,8}, {1,4,5,8}, {5,6,7,8}, {0,1,2}, {0,2,3}, {0,3,4}, {0,4,1}}
--- hal's homology square

-- the following is our standard pyramid example

V = {{0,0,-1}, {-1,-1,2}, {1,-1,1},  {-1,1,1}, {1,1,1}};
F = {{0,1,2}, {0,1,3}, {0,2,4}, {0,3,4}, {1,2,3,4}};
fileName = "pyramid.m2"
fileName' = "pyramidSimplicial.m2"

-- The following is the cube example
V = {{1,1,-1},{1,-1,-1},{-1,1,-1},{-1,-1,-1},{1,1,1},{1,-1,1},{-1,1,1},{-1,-1,1}};
F = {{4,5,7,6},{0,1,3,2},{0,2,6,4},{0,1,5,4},{1,3,7,5},{3,2,6,7}};
fileName = "cube.m2"
fileName' = "cubeSimplicial.m2"

-- The following is Fultons examples
V = {{1,1,-1},{1,-1,-1},{-1,1,-1},{-1,-1,-1},{1,2,3},{1,-1,1},{-1,1,1},{-1,-1,1}};
F = {{4,5,7,6},{0,1,3,2},{0,2,6,4},{0,1,5,4},{1,3,7,5},{3,2,6,7}};
fileName = "fulton.m2"
fileName' = "fultonSimplicial.m2"

--- This is a "more generic Fulton". Its example 4.4 in Katz-Payne.
V = {{1,1,-1},{1,-1,-1},{-1,1,-1},{-1,-1,-1},{1,1,2},{1,-1,2},{-1,1,1},{-1,-1,1}};
F = {{4,5,7,6},{0,1,3,2},{0,2,6,4},{0,1,5,4},{1,3,7,5},{3,2,6,7}};
fileName = "fulton44.m2"
fileName' = "fultonSimplicial44.m2"

-- The following is the dodecahedron (platonic solid)
V = {{1,1,1}, {1,1,-1}, {1,-1,1}, {1,-1,-1}, {-1,1,1}, {-1,1,-1}, {-1,-1,1}, {-1,-1,-1}, {0,144/89,89/144}, {0,144/89,-89/144}, {0,-144/89,89/144}, {0,-144/89,-89/144}, {89/144,0,144/89}, {89/144,0,-144/89}, {-89/144,0,144/89}, {-89/144,0,-144/89}, {144/89,89/144,0}, {144/89,-89/144,0}, {-144/89,89/144,0}, {-144/89,-89/144,0}};
F = {{0,8,4,18,16}, {0,16,1,9,8}, {1,13,5,19,17}, {1,17,3,11,9}, {2,12,6,14,10}, {2,10,7,15,14}, {3,11,7,15,13}, {3,17,1,9,11}, {4,18,5,14,12}, {5,19,6,14,18}, {6,15,7,11,19}, {7,13,5,19,11}};
fileName = "dodecahedron.m2"
fileName' = "dodecahedronSimplicial.m2"
-- haven't computed.

--EgyptPentagon
V={{0,0,0},{0,2,1},{1,1,1},{1,-1,1},{-1,-1,1},{-1,1,1},{0,0,-1}}
F={{0,1,2,3,4,5},{0,6,1,2},{0,6,3,4},{0,6,4,5},{0,6,5,1},{0,6,2,3}}
fileName = "pentagonalBase.m2"
fileName' = "pentagonalBaseSimplicial.m2"


-- EgyptPentagonalPrism
V={{0,0,0},{0,2,1},{1,1,1},{1,-1,1},{-1,-1,1},{-1,1,1},{0,2,-1},{1,1,-1},{1,-1,-1},{-1,-1,-1},{-1,1,-1}}
F={{0,1,2,6,7},{0,2,3,7,8},{0,3,4,8,9},{0,4,5,9,10},{0,5,1,10,6},{0,1,2,3,4,5},{0,6,7,8,9,10}}
fileName = "pentagonalPrism.m2"
fileName' = "pentagonalPrismSimplicial.m2"

-- triangular prism
V = {{0,0,0},{1/2,3/2,1}, {1/2,-1/2,1},{-3/2,-1/2,1},{1/2,3/2,-1}, {1/2,-1/2,-1},{-3/2,-1/2,-1}}
F={{0,1,2,3},{0,4,5,6},{0,1,2,4,5}, {0,3,2,6,5}, {0,1,4,6,3}}
fileName = "triangularPrism.m2"
fileName' = "triangularPrismSimplicial.m2"

V = {{0,0,0},{1/2,3/2,1}, {1/2,-1/2,1},{-3/2,-1/2,1},{1/2,3/2,-1}, {1/2,-1/2,-1},{-3/2,-1/2,-1}}
F={{0,1,2,3},{0,4,5,6},{0,1,2,4,5}, {0,3,2,6,5}, {0,1,4,6,3}}
--- hexagonal cone



-*
To save a triangulation for later, you can use the saveTriangulation and loadTriangulation functions in helpers.m2. The pwd may need to be changed for where you want to store these.
Here is how it works:
1) Create Delta = fan(V,F)
2) triangulation = findUnimodularTriangulation(Delta)
3) saveTriangulation(triangulation, filename, pwd) 
*-

Delta = fan(V,F);
isComplete Delta
C = splineComplex(Delta, 0, Homogenize => false)
S = ring(C)
for i from 0 to 3 list prune HH_i C 
prune HH_2 C
hasExtraKernel(Delta, 3)
triangulation = findUnimodularTriangulation(Delta)
saveTriangulation(triangulation, fileName)

loadTriangulation(fileName, Delta) -- this will save the precomputed unimodularTriangulation into Delta 
C = splineComplex(Delta, 0, Homogenize => false)
Splines = splineModule(Delta,0, Homogenize => false)
R = ring Splines;

------ The right side of the square deals with the simplicialization of Delta
Delta' = simplicialization(Delta);
--triangulation = findUnimodularTriangulation(Delta')
--saveTriangulation(triangulation, fileName)

loadTriangulation(fileName', Delta') -- this will save the precomputed unimodularTriangulation into Delta'
Splines' = splineModule(Delta',0, Homogenize => false) 
splineObjs' = splineList(Splines', Delta', R);
--- The following is the map of lattices mapping Delta' --> Delta
psi = map(Delta,Delta', map(ZZ^(ambDim Delta), ZZ^(ambDim Delta'), 1))


--- This checks that the diagram commutes for all the splines:
--splineObjs = splineList(Splines, Delta, R);
--for testSpline in splineObjs list (
--  mat pullback(psi, chowMap(testSpline)) == mat chowMap(pullback(psi, testSpline))
--


-*
Our simplicialization map psi: Delta' --> Delta induces a dominant map X' --> X by Proposition 2.7 in Fulton-Sturmfels. Thus we should be able to pullback Minkowski weights in A^k(X') to A^k(X). 

Let c be a Minkowski weight of codimension k on Delta'. Let tau in Delta be a cone of codimension k. Let tau' be the smallest cone of Delta' that contains psi(tau)
*-

k=3


---- Lets compute the left side of the square.
kSplines = image super basis(k, Splines) -- this is the elements of A_T^*(X) that generate A^k_T(X).
prune kSplines -- these have a bunch of trivial relations
generators kSplines

-- Its expressed as a 5 x 11 matrix: There are 11 generators and each one has 5 parts, one for each maximal cone of Delta.
kSplinesObjs = splineList(kSplines, Delta, R)
ik = chowMap(kSplinesObjs, k) -- this is a matrix with the images of the 11 generators of A^k_T(X) in A^k(X) 

Ak = operationalChowGroup(Delta, k)
chowGroup(Delta, k) 
Ak
prune Ak -- this a presentation of ZZ^2

--- I want to express ik as a map from kSplines to A^k(X). 
-- ik maps into Ak, so each column of ik can be rewritten as a Z-linear combination of the generators of Ak.
-- The following matrix writes ik in the basis given by the generators of Ak.
ikMat = transpose matrix for col from 0 to numcols ik - 1 list entries solve(generators Ak, lift(ik_col,ZZ))

-- so (generators Ak) * M == ik.

--- this models the map A^k_T(X) --> A^k(X) as a map from a free module ZZ^(numcols ik) to Ak....
ikMap = map(Ak, ZZ^(numgens kSplines), ikMat)
tex ikMap
prune ker ikMap
prune coker ikMap
rank ikMap


-- Here, we are getting the kernel if ik as a module of splines and seeing if its the same as M*A^(k-1)_T(X)
kerModule = prune image( (gens kSplines)*(generators ker ikMap)) --- the kernel as a module of splines
guessModule = prune image((gens image super basis(k-1, Splines))**(vars R))
kerModule**QQ
guessModule**QQ
kerModule**QQ == guessModule**QQ -- true!




------ The right side of the square deals with the simplicialization of Delta
kSplines' = image super basis(k, Splines') -- this is the elements of A_T^*(X) that generate A^k_T(X).
prune kSplines' -- these have a bunch of trivial relations
generators kSplines'


kSplinesObjs' = splineList(kSplines', Delta', R)
ik' = chowMap(kSplinesObjs', k) ---- this will take a long time because there are 18 splines...

Ak' = operationalChowGroup(Delta', k)
prune Ak'

--- I want to express ik as a map from kSplines to A^k(X). 
-- ik maps into Ak, so each column of ik can be rewritten as a Z-linear combination of the generators of Ak.
-- The following matrix writes ik in the basis given by the generators of Ak.
ikMat' = transpose matrix for col from 0 to numcols ik' - 1 list entries solve(generators Ak', lift(ik'_col,ZZ))

-- so (generators Ak) * M == ik.

--- this models the map A^k_T(X) --> A^k(X) as a map from a free module ZZ^(numcols ik) to Ak....
ikMap' = map(Ak', ZZ^(numgens kSplines'), ikMat')
prune ker ikMap'
prune coker ikMap'



-*
Our simplicialization map psi: Delta' --> Delta induces a dominant map X --> X' by Proposition 2.7 in Fulton-Sturmfels. Thus we should be able to pullback Minkowski weights in A^k(X') to A^k(X). 

Let c be a Minkowski weight of codimension k on Delta'. Let tau in Delta be a cone of codimension k. Let tau' be the smallest cone of Delta' that contains psi(tau)
*-


----- The vertical maps in the square are implemented as pullback and pushforward
--- What about the inclusion maps?

kerSplines = (generators kSplines)*(transpose transpose promote(generators ker ikMap, R)) -- this is the basis composed with ikMat.
kerSplinesPullback = for s in splineList(image kerSplines, Delta, R) list pullback(psi, s)
kerSplinesPullback / mat

sourceMat = generators kSplines
iota = for s in kSplinesObjs list pullback(psi, s)
targetMat = transpose matrix(iota / mat / entries / flatten)

generators(kSplines' ** R)
targetMat ** R


kSplinesObjs = splineList(kSplines, Delta, R)
splineMats = for s in kSplinesObjs list pullback(psi, s)
iota = fold((a,b) -> a | b, splineMats / mat / entries / matrix)


prune ker ikMap'

kSplinePullback = fold( (a,b) -> a | b, apply(kSplinesObjs, s -> mat pullback(psi, s)))
generators (pullback(psi,kSplines**R))