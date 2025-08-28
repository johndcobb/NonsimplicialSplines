-* So it looks like we don’t need to translate the P monomials to monomials in the lattice M, can just use them straight up. Here is the way to embed X_\Sigma where \Sigma comes from coning over faces of symmetric cube, and using the divisor {0,0,0,2,2,2,2}.
*-
CubeRays={{1,-1,-1},{-1,-1,-1},{1,1,-1},{-1,1,-1},{1,-1,1},{-1,-1,1},{1,1,1},{-1,1,1}}
CubeCones={{0,1,2,3},{0,1,4,5},{0,2,4,6},{1,3,5,7},{2,3,6,7},{4,5,6,7}}
Cube=normalToricVariety(CubeRays, CubeCones)
classGroup(Cube)
picardGroup(Cube)
fromPictoCl(Cube)
--gives the right divisor, with 7 entries not 8 because Class group is (Z/2)^2+Z^5.
D=toricDivisor({0,0,0,0,2,2,2,2},Cube)
isAmple D --true
--DD = OO D --line bundle
--get the Cox ring
Cox=ring(Cube)
---use Cox_D to map to projective space
R2=QQ[t_0..t_6] --there are seven sections.
Csections=super basis({0,0,0,2,2,2,2}, Cox)
I=kernel map(Cox,R2, Csections)
hilbertPolynomial(coker gens I)
IS= singularLocus I
ISRad=primaryDecomposition(radical ideal presentation IS)