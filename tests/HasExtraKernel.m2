needs "code/ChowMap.m2"

triPath = currentDirectory() | "examples/triangulations/";

cubeV = {{1,1,-1},{1,-1,-1},{-1,1,-1},{-1,-1,-1},{1,1,1},{1,-1,1},{-1,1,1},{-1,-1,1}};
cubeF = {{4,5,7,6},{0,1,3,2},{0,2,6,4},{0,1,5,4},{1,3,7,5},{3,2,6,7}};
cube = fan(cubeV, cubeF);
loadTriangulation("cube.m2", triPath, cube);

assert(not hasExtraKernel(cube, 2));
