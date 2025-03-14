needsPackage "NormalToricVarieties"; needsPackage "SimplicialComplexes"; 
debug needsPackage "AlgebraicSplines"; needsPackage "Polyhedra";

billeraComplex = method()
billeraComplex(List, List, ZZ) := ChainComplex => (V, F, r) -> (
    Sigma := polyhedralComplex(V,F); -- save data as a PolyhedralComplex
    d := #(V_0);
    R := QQ[x_0..x_(d-1)];
    -- This collects the modules that appears in each spot of the Billera complex
    B := for i in 0..d-1 list directSum(for P in listToPolyhedra(polyhedra(i, Sigma),Sigma) list (module R)/J(P, Sigma, R, r)) -- there is some weird behavior for i = 0, and I'm missing i = d.
)

polyhedralComplex(List, List) := PolyhedralComplex => (V,F) -> (
    polyhedralComplex(apply(F, C -> transpose matrix apply(C, idx -> V_idx)) / convexHull)
)

J = method()
J(Polyhedron, PolyhedralComplex, Ring, ZZ) := Ideal => (P, Sigma, R, r) -> (
    V := vertices Sigma; 
    D := listToPolyhedra(polyhedra(dim Sigma, Sigma), Sigma);
    adjacenttoD := select(D, d -> contains(d, P));
    sum(apply(adjacenttoD, I -> (getHyperplaneEquation(I, R))^r))
)

getHyperplaneEquation = method()
getHyperplaneEquation(Polyhedron, Ring) := Ideal => (P,R) -> (
    d := ambDim P;
    A := transpose( vertices P | matrix(toList(d:{0})));
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