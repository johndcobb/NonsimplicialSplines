load "../code/ChowMap.m2"
------------------------------------------------------
-* Explanation
Let X' --> X be a simplicialization of a complete toric variety X. We have the following maps:
(1) From Katz-Payne, we have a map (chowMap) i^k : A_T^k(X) --> A^k(X) 
(2) from Proposition 2.7 in Fulton-Sturmfels: A^k(X) --> A^k(X')
(3) From Brion, we have a map A_T^k(X') --> A^k(X'). Should just be the same thing as the chowMap for X', but is much simpler since X' is simplicial.
UPDATE: WAIT. Brion is for the chow GROUP. We really should just do the same chow map for X'.
(4) The trivial inclusion of splines A_T^k(X) --> A_T^k(X')

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
The following is the cube example
*-
------------------------------------------------------
V = {{1,1,-1},{1,-1,-1},{-1,1,-1},{-1,-1,-1},{1,1,1},{1,-1,1},{-1,1,1},{-1,-1,1}};
F = {{4,5,7,6},{0,1,3,2},{0,2,6,4},{0,1,5,4},{1,3,7,5},{3,2,6,7}};
Delta = fan(V,F);
Splines = splineModule(Delta,0, Homogenize => false)
R = ring Splines;

-- Lets compute the square for k = 2
k = 2

---- Lets compute the left side of the square.
kSplines = image super basis(k, Splines) -- this is the elements of A_T^*(X) that generate A^k_T(X).
-- Its expressed as a 6 x 11 matrix: There are 11 generators and each one has 6 parts, one for each maximal cone of Delta.
kSplinesObjs = splineList(kSplines, Delta, R)
ik = chowMap(kSplinesObjs, k) -- this is a matrix with the images of the 11 generators of A^k_T(X) in A^k(X) 


Ak = operationalChowGroup(Delta, k)

--- I want to express ik as a map from kSplines to A^k(X). 
-- ik maps into Ak, so each column of ik can be rewritten as a Z-linear combination of the generators of Ak.
-- The following matrix writes ik in the basis given by the generators of Ak.
ikMat = transpose matrix for col from 0 to numcols ik - 1 list (
   entries solve(generators Ak, lift(ik_col,ZZ))
)
-- so (generators Ak) * M == ik.

--- this models the map A^k_T(X) --> A^k(X) as a map from a free module ZZ^(numcols ik) to Ak....
ikMap = map(Ak, ZZ^(numcols ik), ikMat)


--- If I just apply Chow map to the generators of splines, then I get the map A^k_T(X)/M*A^(k-1)(X) --> A^k(X).
splineObjs = splineList(Splines, Delta, R);
ik = chowMap(splineObjs, k)


------ The right side of the square deals with the simplicialization of Delta
Delta' = simplicialization(Delta);
Splines' = splineModule(Delta',0, Homogenize => false) -- here, I tensor with R so that they both are in the same ring.
splineObjs' = splineList(Splines', Delta', R);

f = map(Delta,Delta', map(ZZ^(ambDim Delta), ZZ^(ambDim Delta'), 1))
linearSpline = spline(entries Splines_1,Delta, R)
describe linearSpline
c = chowMap(linearSpline)

for tau in cones c list (
   imagePhi := affineImage(mat f, tau);
   -- Here i should find tau, the smallest cone of Delta' containing imagePhi
   apply()
)


pullback(f, MinkowskiWeight) := MinkowskiWeight => opts -> c -> (
   for tau cones c
   affineImage(f, )
)
all(apply(maxFacesAsCones(Sigma1), sigma -> (
        imagePhi := affineImage(phi, sigma);
        any(apply(maxFacesAsCones(Sigma2), tau -> contains(tau, imagePhi)), bool -> bool)
            )), bool -> bool)
    )












------------------------------------------------------
-*
The following is the pyramid example
*-
------------------------------------------------------
V = {{0,0,-1}, {-1,-1,1}, {1,-1,1}, {-1,1,1}, {1,1,1}}
F = {{0,1,2}, {0,1,3}, {0,2,4}, {0,3,4}, {1,2,3,4}};
Delta = fan(V,F);
Splines = splineModule(Delta,0, Homogenize => false)
R = ring Splines;

-- Lets compute the square for k = 2
k = 2

splineObjs = splineList(Splines, Delta, R);
ik = chowMap(splineObjs, k)


Delta' = simplicialization(Delta);
Splines' = splineModule(Delta',0, Homogenize => false) -- here, I tensor with R so that they both are in the same ring.
splineObjs' = splineList(Splines', Delta', R);

ik' = chowMap(splineObjs', k)