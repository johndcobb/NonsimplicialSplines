needs "code/ChowMap.m2"

zeroCone3 = coneFromVData transpose matrix{{0,0,0}};
assert(dim zeroCone3 == 0);
assert(isUnimodular zeroCone3);

smoothCone3In4 = coneFromVData transpose matrix{{1,0,0,0},{0,1,0,0},{0,0,1,0}};
assert(dim smoothCone3In4 == 3);
assert(ambDim smoothCone3In4 == 4);
assert(isSimplicial smoothCone3In4);
assert(isUnimodular smoothCone3In4);

indexTwoCone2In3 = coneFromVData transpose matrix{{1,0,0},{1,2,0}};
indexTwoTriangulation = findUnimodularTriangulation indexTwoCone2In3;
assert(length indexTwoTriangulation == 2);
assert(all(indexTwoTriangulation, sigma -> dim sigma == 2 and ambDim sigma == 3 and isSimplicial sigma and isUnimodular sigma));

squareCone3In4 = coneFromVData transpose matrix{{1,0,1,0},{0,1,1,0},{-1,0,1,0},{0,-1,1,0}};
squareTriangulation = findUnimodularTriangulation squareCone3In4;
assert(length squareTriangulation == 4);
assert(all(squareTriangulation, sigma -> dim sigma == 3 and ambDim sigma == 4 and isSimplicial sigma and isUnimodular sigma));

fullCone4 = coneFromVData transpose matrix{{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
lowerSquareCone3In4 = coneFromVData transpose matrix{{-1,1,0,0},{-1,0,1,0},{-1,-1,0,0},{-1,0,-1,0}};
mixedDimensionalFan = fan{fullCone4, lowerSquareCone3In4};
mixedSimplicialization = simplicialization mixedDimensionalFan;
assert(isSimplicial mixedSimplicialization);
assert(length maxFacesAsCones mixedSimplicialization == 5);
mixedTriangulation = findUnimodularTriangulation mixedDimensionalFan;
assert(length mixedTriangulation == 2);
assert(all(flatten mixedTriangulation, sigma -> isSimplicial sigma and isUnimodular sigma));
assert(any(mixedTriangulation, triangulation -> length triangulation == 4 and all(triangulation, sigma -> dim sigma == 3 and ambDim sigma == 4)));
