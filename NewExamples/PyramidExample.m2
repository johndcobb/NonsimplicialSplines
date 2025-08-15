load "../ChowMap.m2";

V = {{0,0,-1}, {-1,-1,1}, {1,-1,1},  {-1,1,1}, {1,1,1}};
F = {{0,1,2}, {0,1,3}, {0,2,4}, {0,3,4}, {1,2,3,4}};
Sigma = fan(V,F);
Splines = splineModule(Sigma,0, Homogenize => false)
R = ring Splines;

splines = splineList(Splines, V, F, R)
c0 = chowMap(splines, 0)
operationalChowGroup(Sigma, 0)
c1 = chowMap(splines, 1)
operationalChowGroup(Sigma, 1)
c2 = chowMap(splines, 2)
operationalChowGroup(Sigma, 2)
c3 = chowMap(splines, 3)
operationalChowGroup(Sigma,)