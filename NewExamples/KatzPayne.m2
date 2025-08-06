load "../ChowMap.m2";

-- Here is Example 4.1 from Katz and Payne, the "mirror dual" of P^1 x P^1 
V = {{1,1},{1,-1}, {-1,-1}, {-1,1}, {0,0}};
F = {{0,1,4},{1,2,4},{2,3,4}, {3,0,4}};
Splines = splineModule(V,F, 0, Homogenize => false)
R = ring Splines; 
minimalPresentation Splines -- Spline module is free: R + R(-1)^2 + R(-2)

assert(entries Splines_3 == {t_1^2 - t_2^2, 0, 0, 0})
f = spline(entries Splines_3, V, F , R);
describe f
origin = coneFromVData transpose matrix{{0,0}}
-- The following does the same thing as in Example 4.1 from [Katz,Payne]
assert(chowMap(f, origin)==2)
-- we can also compute the equivariant multiplicity directly:
assert(chowMap(f, origin) == equivariantMultiplicity(sigma1, origin, R)*restriction(f, sigma1))

-- On a degree 1 spline, the chowMap should give integers on the codim 1 cones. That is, the rays of tha fan.
f = spline(entries Splines_1, V, F, R)
describe f
(v1, v2, v3, v4) = (coneFromVData transpose matrix{{1,1}}, coneFromVData transpose matrix{{1,-1}}, coneFromVData transpose matrix{{-1,-1}}, coneFromVData transpose matrix{{-1,1}});

assert(chowMap(f, v1) == 1)
assert(chowMap(f, v2) == 0)
assert(chowMap(f, v3) == 1)
assert(chowMap(f, v4) == 0)

-- So the minkowski weight corresponding to f is the rays v1 and v3. Is this reasonable?

-- Again, we could do this by hand:
assert(chowMap(f,v1) == equivariantMultiplicity(sigma1, v1, R) * restriction(f, sigma1))


---------- Cube Example
V = {{1,1,-1},{1,-1,-1},{-1,1,-1},{-1,-1,-1},{1,1,1},{1,-1,1},{-1,1,1},{-1,-1,1},{0,0,0}};
F = {{4,5,7,6,8},{0,1,3,2,8},{0,2,6,4,8},{0,1,5,4,8},{1,3,7,5,8},{3,2,6,7,8}};
SP=splineModule(V, F,0, Homogenize => false)
minimalPresentation SP -- so its R + R(-1) + R(-2)^2 + R(-3) + R(-4)
R = ring SP
entries SP_5

f4 = spline(entries SP_5, V, F, ring SP)
origin = coneFromVData transpose matrix{{0,0,0}}
chowMap(f4, origin) -- Since this is a degree 4 spline, it should give a homogeneous polynomial of degree 4 + 0 - 3 = 1 per Proposition 1.2

f3 = spline(entries SP_4,V,F,ring SP) -- this is a degree 3 spline, so  when we compute the chow map at the origin it should get a integer, and it is.
chowMap(f3, origin)

f2 = spline(entries SP_2, V, F, ring SP) -- this is a degree 2 spline, so when we compute the chow map at the origin we need to get zero since deg f + dim origin - 3 < 0
chowMap(f2, origin) -- this should be zero, and it is.
-- but at a ray, we should get an integer
tau = coneFromVData transpose matrix{{1,1,1}}
chowMap(f2, tau) -- and we do
apply(facesAsCones(2,fan(f2)), sigma -> chowMap(f2, sigma))

f1 = spline(entries SP_1, V, F, ring SP) -- this is a degree 1 spline
apply(facesAsCones(1, fan(f1)), sigma -> chowMap(f1, sigma)) -- this should give integers at the rays, and it does.

--- Okay, now to try and get the Table.
-- i = 0
f0 = spline(entries SP_0, V, F, ring SP) -- deg 0 means i need a codim 0 cone, so a dim 3 cone.
apply(facesAsCones(0,fan(f0)), sigma -> chowMap(f0, sigma))

------
V = {{0,1}, {1,0}, {-1,-1}, {0,0}}
F = {{0,1,3},{1,2,3}, {2,0,3}}
Splines = splineModule(V,F,0, Homogenize => false)

f = spline(entries Splines_1, V, F, ring Splines)
g = spline(entries Splines_2, V, F, ring Splines)

apply(facesAsCones(degree f, fan(f)), sigma -> chowMap(f, sigma))
apply(facesAsCones(degree g, fan(g)), sigma -> chowMap(g, sigma))

---------
------ Fulton Example 4.3

FultonVerts={{0,0,0},{-1,-1,1},{1,-1,1},{-1,1,1},{1,2,3},{-1,-1,-1},{1,-1,-1},{-1,1,-1},{1,1,-1}}
Cube={{0,1,2,3,4},{0,1,2,5,6},{0,1,3,5,7},{0,2,4,6,8},{0,3,4,7,8},{0,5,6,7,8}}
SP=splineModule(FultonVerts, Cube,0, Homogenize => false)
minimalPresentation SP


entries SP_2
f = spline(entries SP_2, FultonVerts, Cube);


