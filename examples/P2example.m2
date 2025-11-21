load "../code/ChowMap.m2"

V = {{-1,-3},{0,1},{1,0}};
F = {{0,1},{1,2},{2,0}};
Delta = fan(V,F);
Splines = splineModule(Delta,0, Homogenize => false)
R = ring Splines;

splineObjs = splineList(Splines, Delta, R)
chowMap(splineObjs,1) 
facesAsCones(0, Delta) / rays