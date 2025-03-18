needsPackage "NormalToricVarieties"; needsPackage "SimplicialComplexes"; needsPackage "Polyhedra";
needsPackage "CellularResolutions"

billeraComplex = method()
billeraComplex(List, List, ZZ) := ChainComplex => (V, F, r) -> (
    x := symbol x;
    Sigma := polyhedralComplex(V,F); -- save data as a PolyhedralComplex
    d := #(V_0);
    R := QQ[x_0..x_(d-1)];
    -- This collects the modules that appears in each spot of the Billera complex
    B := append(for i in 1..d-1 list directSum(for P in listToPolyhedra(polyhedra(i, Sigma),Sigma) list (module R)/J(P, Sigma, R, r)) , (module R)^(#maxPolyhedra(Sigma))); -- i=1 is actually checking the vertices.
    -- then i have to create the maps between these modules....
    C := cellComplex(R, Sigma);
    maps := for i in 0..d-2 list map(B_i, B_(i+1), boundaryMap(i+1,C));
    chainComplex(maps)
)

polyhedralComplex(List, List) := PolyhedralComplex => (V,F) -> (
    polyhedralComplex(apply(F, C -> transpose matrix apply(C, idx -> V_idx)) / convexHull)
)

J = method()
J(Polyhedron, PolyhedralComplex, Ring, ZZ) := Ideal => (P, Sigma, R, r) -> (
    --V := vertices Sigma; 
    D := listToPolyhedra(polyhedra(dim Sigma, Sigma), Sigma); -- this lists codim 1 polyhedra
    adjacenttoD := select(D, d -> contains(d, P)); -- this lists the polyhedra containing P
    trim sum(apply(adjacenttoD, I -> (getHyperplaneEquation(I, R))^r)) 
)

getHyperplaneEquation = method()
getHyperplaneEquation(Polyhedron, Ring) := Ideal => (P,R) -> (
    d := ambDim P;
    A := transpose( vertices P | matrix(toList(d:{0}))); -- this is a matrix of points (including the origin) that lie upon a hyperplane.
    -- The kernel of this matrix is the normal vector of the plane containing these points. 
    ideal( (vars R)*(gens kernel A) ) 
)
getHyperplaneEquation(Polyhedron) := Ideal => (P) -> (
    d := ambDim P;
    R = QQ[x_0..x_(d-1)];
    getHyperplaneEquation(P, R)
)

listToPolyhedra = method()
listToPolyhedra(List, PolyhedralComplex) := List => (L,Sigma) -> (
    V := vertices Sigma;
    apply(L, l -> V_(l_0)) / convexHull
)