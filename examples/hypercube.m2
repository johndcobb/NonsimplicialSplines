load "../code/ChowMap.m2"

V = {
 {1,1,-1,-1},{1,-1,-1,-1},{-1,1,-1,-1},{-1,-1,-1,-1},
 {1,1,1,-1},{1,-1,1,-1},{-1,1,1,-1},{-1,-1,1,-1},
 {1,1,-1,1},{1,-1,-1,1},{-1,1,-1,1},{-1,-1,-1,1},
 {1,1,1,1},{1,-1,1,1},{-1,1,1,1},{-1,-1,1,1}
};
F = {
 {0,1,2,3,4,5,6,7},      -- w = -1
 {8,9,10,11,12,13,14,15},-- w = +1
 {0,1,4,5,8,9,12,13},    -- x = +1
 {2,3,6,7,10,11,14,15},  -- x = -1
 {0,2,4,6,8,10,12,14},   -- y = +1
 {1,3,5,7,9,11,13,15},   -- y = -1
 {0,1,2,3,8,9,10,11},    -- z = -1
 {4,5,6,7,12,13,14,15}   -- z = +1
};
fileName = "hypercube.m2"
fileName' = "hypercubeSimplicial.m2" 

Delta = fan(V,F);
triangulation = findUnimodularTriangulation(Delta)
saveTriangulation(triangulation, fileName)

loadTriangulation(fileName, Delta) -- this will save the precomputed unimodularTriangulation into Delta 
Splines = splineModule(Delta,0, Homogenize => false)
R = ring Splines;
C = splineComplex(Delta, 0)
Splines = prune HH C -- gives a gradedModule
prune Splines
betti Splines

Delta' = simplicialization(Delta);
triangulation' = findUnimodularTriangulation(Delta')
saveTriangulation(triangulation', fileName')
Splines' = splineModule(Delta',0, Homogenize => false)**R
splineObjs' = splineList(Splines', Delta', R);
--- The following is the map of lattices mapping Delta' --> Delta
psi = map(Delta,Delta', map(ZZ^(ambDim Delta), ZZ^(ambDim Delta'), 1))