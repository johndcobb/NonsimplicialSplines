Something is broken about the code in `SplinesCode.m2`. Luckily, we can make the `AlgebraicSplines` package work. 

Before, we were writing our vertices and faces like this:

`V = {{0,0,-1}, {-1,-1,1}, {1,-1,1},  {-1,1,1}, {1,1,1}};`
`F = {{0,1,2}, {0,1,3}, {0,2,4}, {0,3,4}, {1,2,3,4}};`

In order to do the same thing using the AlgebraicSplines package, we should add in the origin as a point. This will lift everything  by a dimension: 
`V = {{0,0,-1}, {-1,-1,1}, {1,-1,1},  {-1,1,1}, {1,1,1}, {0,0,0}};`
`F = {{0,1,2,5}, {0,1,3,5}, {0,2,4,5}, {0,3,4,5}, {1,2,3,4,5}};`
`Splines = splineModule(V,F,0)`




