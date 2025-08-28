-*
Is there a way so compute spline dimensions for some graph that’s not necessarily arising as the dual graph of a simplicial complex? (For what we were talking about today my dual graph would be the octahedron given as the span of plus and minus the standard basis vectors ei, but I want something more general).

I'm hoping that if we take the induced subgraph of the Cayley graph of Sn (with edge label ti-tj on the edge w—(i,j)w) given by the “noncrossing permutations” obtained by taking a noncrossing partition of {1,…,n} and putting a backwards cycle on each part, then the r-spline ring is the image of r-quasisymmetric polynomials of Hivert in the quotient ring Q[x1,…,xn]/(positive degree quasisymmetric polynomials).

If so this would be the first and only interesting thing ever proved about r quasisymmetric polynomials in 20 years.

It’s true for r=1!
*- 
restart
needsPackage "AlgebraicSplines";

-- Lets try S_3. 
-- I've computed the subgraph described by hunter by hand, can code it up later.

-- The function takes in
-- E = a list of edges, with each edge represented as a list with two elements.
-- I = a list of ideals in a ring

S = QQ[t_1..t_3];

r= 4
E = {{0,1},{0,2},{0,3},{1,4},{2,4},{3,4}}
I = {(t_3-t_2)^r, (t_3-t_1)^r, (t_2-t_1)^r, (t_3-t_1)^r, (t_2-t_1)^r, (t_3-t_2)^r}
rSplines = generalizedSplines(E,I)
minimalPresentation rSplines

R = QQ[x_1..x_3]
4:3

-- The conjecture:
-- The r-spline ring is the image of r-quasisymmetric polynomials of Hivert in the quotient ring Q[x1,…,xn]/(positive degree quasisymmetric polynomials).