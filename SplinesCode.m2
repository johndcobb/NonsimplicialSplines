needsPackage "NormalToricVarieties"; needsPackage "SimplicialComplexes"; needsPackage "Polyhedra";

billeraComplex = method()
billeraComplex(List, List, Ring, ZZ) := ChainComplex => (V, F, R, r) -> (
    Sigma := polyhedralComplex(V,F); -- save data as a PolyhedralComplex
    d := dim Sigma;
    -- This collects the modules that appears in each spot of the Billera complex
    B := append(for i in 1..d list directSum(for P in listToPolyhedra(polyhedra(i, Sigma),Sigma) list (module R)/J(P, Sigma, R, r)) , (module R)^(#maxPolyhedra(Sigma))); -- i=1 is actually checking the vertices.
    -- then i have to create the maps between these modules....
    maps := for i in 1..d list map(B_(i-1), B_i, promote(boundaryMaps(i,Sigma), R));
    chainComplex(maps)
)
billeraComplex(List, List, ZZ) := ChainComplex => (V, F, r) -> (
    x := symbol x;
    d := #(V_0);
    R := QQ[x_0..x_(d-1)];
    billeraComplex(V,F,R,r)
)

-- This is the i^th boundary map in the Billera complex
boundaryMaps = method()
boundaryMaps(ZZ, PolyhedralComplex) := Matrix => (r, Sigma) -> (
    --Vmat := vertices Sigma;
    --if r == dim Sigma then (
    --    sourceCells := faces(0,fan(Sigma));
    --    sourcePolyhedra := apply(sourceCells, C -> Vmat_C) / convexHull;
    --) else (
    --    sourceCells := apply(polyhedra(r+1,Sigma), P -> P_0);
    --    sourcePolyhedra := listToPolyhedra(polyhedra(r+1,Sigma), Sigma);
    --);
    --targetCells := apply(polyhedra(r,Sigma), P -> P_0);

    Pfaces := applyPairs(faces Sigma, (i,lst) -> ((dim Sigma)-i-1,apply(lst,first)));
    sourceCells := Pfaces#r;
    targetCells := Pfaces#(r-1);
    boundaries := new MutableHashTable;
    for i from 0 to dim Sigma do (
        for face in Pfaces#i do (
            boundaries#face = for f in Pfaces#(i-1) list (if isSubset(f,face) then f else continue);
        );
    );

    --boundariesAsCones := apply(sourcePolyhedra, P -> facesAsCones(1, cone(P)));
   -- boundariesAsMats := apply(boundariesAsCones, C -> apply(C, v -> submatrix'(rays v, {0}, )));
   -- boundariesAsLists := apply(boundariesAsMats, L -> matsToIndex(L, Sigma));

    M := mutableMatrix map(ZZ^#targetCells, ZZ^#sourceCells, 0);
    for i in 0..#sourceCells-1 do (
        sourceBdd := boundaries#(sourceCells_i);
        targetIdxs := for i in sourceBdd list (position(targetCells, targets -> targets == i));
        orientation := 0;
        for targetIdx in reverse(targetIdxs) do (
            M_(targetIdx, i) = (-1)^(orientation);
            orientation += 1;
        );
    );
    matrix M
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
    x := symbol x; 
    d := ambDim P;
    R := QQ[x_0..x_(d-1)];
    getHyperplaneEquation(P, R)
)

listToPolyhedra = method()
listToPolyhedra(List, PolyhedralComplex) := List => (L,Sigma) -> (
    Vmat := vertices Sigma;
    apply(L, l -> Vmat_(l_0)) / convexHull
)

matsToIndex = method()
matsToIndex(List, PolyhedralComplex) := List => (L, Sigma) -> (
    Vmat := vertices Sigma;
    VList := for numcol in 0..(numColumns(Vmat)-1) list Vmat_numcol;
    for M in L list (
      for numcol in 0..numColumns(M)-1 list (
      position(VList, v -> lift(v, ZZ) == M_numcol )
        )
    )
)
matsToIndex(Matrix, PolyhedralComplex) := List => (M, Sigma) -> (
    Vmat := vertices Sigma;
    VList := for numcol in 0..(numColumns(Vmat)-1) list Vmat_numcol;
    for numcol in 0..numColumns(M)-1 list (
    position(VList, v -> lift(v, ZZ) == M_numcol )
    )
)