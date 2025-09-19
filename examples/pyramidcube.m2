load "../code/ChowMap.m2"

---- Here, we will define a sequence of blowups. The fan Delta will be a cube,
--- The fan Delta' will be the cube with 4 pyramids along the sides
--- the fan Delta'' will be the cube with 6 pyramids, one on each face.

-- cube
V = {{2,2,-2},{2,-2,-2},{-2,2,-2},{-2,-2,-2},{2,2,2},{2,-2,2},{-2,2,2},{-2,-2,2}};
F = {{4,5,7,6},{0,1,3,2},{0,2,6,4},{0,1,5,4},{1,3,7,5},{3,2,6,7}};
Delta = fan(V,F);
Splines = splineModule(Delta,0, Homogenize => false)
R = ring Splines;
minimalPresentation Splines --- 0^1, 1^1, 2^2, 3^1, 4^1  (total: 6)

--leave top and bottom faces of cube as squares
V'={{2,-2,-2},{-2,-2,-2},{2,2,-2},{-2,2,-2},{2,-2,2},{-2,-2,2},{2,2,2},{-2,2,2},{3,0,0},{0,3,0},{-3,0,0},{0,-3,0},{0,0,0}}
F'={
{4,5,6,7,12}, {12,0,1,11},{12,1,11,5},{12,11,4,5},{12,11,0,4},
{12,0,2,8},{12,0,4,8},{12,8,4,6},{12,8,2,6},
{12,1,3,10},{12,3,10,7},{12,10,5,7},{12,10,1,5},
{12,2,3,9},{12,2,6,9},{12,9,6,7},{12,9,3,7}, {0,1,2,3,12}}
Delta' = fan(V',F');
Splines'= splineModule(Delta', 0, Homogenize => false)
minimalPresentation Splines'
--Free, generator degrees and #’s: 0^1,1^7,2^9, 3^1 (total: 18)

--cube with 6 faces turned into Egyptian pyramids
V''={{2,-2,-2},{-2,-2,-2},{2,2,-2},{-2,2,-2},{2,-2,2},{-2,-2,2},{2,2,2},{-2,2,2},{3,0,0},{0,3,0},{-3,0,0},{0,-3,0},{0,0,3},{0,0,-3},{0,0,0}}
F''={
{14,0,1,13},{14,0,2,13},{14,13,2,3},{14,13,1,3},
{14,0,1,11},{14,1,11,5},{14,11,4,5},{14,11,0,4},
{14,0,2,8},{14,0,4,8},{14,8,4,6},{14,8,2,6},
{14,1,3,10},{14,3,10,7},{14,10,5,7},{14,10,1,5},
{14,2,3,9},{14,2,6,9},{14,9,6,7},{14,9,3,7},
{14,4,5,12},{14,5,12,7},{14,12,6,7},{14,12,4,6}}
Delta'' = fan(V'',F'');
Splines''= splineModule(Delta'',0, Homogenize => false)
minimalPresentation Splines'' -- 0^1,1^11,2^11, 3^1 (total: 24)

--- So we have a sequence of blowups Delta'' --> Delta' --> Delta
psi1 = map(Delta,Delta', map(ZZ^(ambDim Delta), ZZ^(ambDim Delta'), 1))
psi2 = map(Delta',Delta'', map(ZZ^(ambDim Delta'), ZZ^(ambDim Delta''), 1))

--- by the way, Delta'' is the same thing as simplicialization(Delta).


--- Lets investigate the chowMap for Delta
k=2
kSplines = image super basis(k, Splines) -- these are the elements of A_T^*(X) that generate A^k_T(X).
generators kSplines
prune kSplines -- these have a bunch of trivial (koszul) relations


kSplinesObjs = splineList(kSplines, Delta, R)
ik = chowMap(kSplinesObjs, k) -- this is a matrix with the images of the 11 generators of A^k_T(X) in A^k(X) 
-- This outputs a 8 x 11 matrix, but Ak is actually rank 5:

Ak = operationalChowGroup(Delta, k)
prune Ak -- this a presentation of ZZ^2

--- We ik as a map from kSplines to A^k(X) by rewriting each columns as a Z-linear combination of the generators of Ak.
ikMat = transpose matrix for col from 0 to numcols ik - 1 list entries solve(generators Ak, lift(ik_col,ZZ))
-- now its 5 x 11 and (generators Ak) * M == ik.
ikMap = map(Ak, ZZ^(numgens kSplines), ikMat)
prune ker ikMap
prune coker ikMap

Mdot = (vars R)**(super basis(k-1, Splines))
prune image Mdot

---- --- Lets investigate the chowMap for Delta'
kSplines' = image super basis(k, Splines') -- this is the elements of
generators kSplines'
prune kSplines' -- these have a bunch of trivial relations
kSplinesObjs' = splineList(kSplines', Delta', R)
ik' = chowMap(kSplinesObjs', k) ---- this will take a long time
-- This outputs a 12 x 36 matrix, but Ak' is actually rank 9:
Ak' = operationalChowGroup(Delta', k)
prune Ak' -- this a presentation of ZZ^9
--- We ik as a map from kSplines to A^k(X) by rewriting each columns as a Z-linear combination of the generators of Ak.
ikMat' = transpose matrix for col from 0 to numcols ik' - 1
    list entries solve(generators Ak', lift(ik'_col,ZZ))
-- now its 9 x 36 and (generators Ak) * M == ik.
ikMap' = map(Ak', ZZ^(numgens kSplines'), ikMat')
prune ker ikMap' -- 27!
prune coker ikMap' -- Z/2 x Z/2!

Mdot' = (vars ring Splines')**(super basis(k-1, Splines'))
prune image Mdot' -- also a quotient of 27.


--- Lets investigate the chowMap for Delta''



