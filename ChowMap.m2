needsPackage "AlgebraicSplines"; needsPackage "NormalToricVarieties"; needsPackage "Polyhedra";
topLevelMode = Standard


equivariantMultiplicity = method()
equivariantMultiplicity(List, Cone, Ring) := ZZ => (unimodularTriangulation, tau, R) -> (
    -- Only need to compute multiplicities if the triangle contains tau
    trianglesContainingTau := select(unimodularTriangulation, triangle -> contains(triangle, tau));

    sum for triangle in trianglesContainingTau list (
        tauidxs := positions(entries transpose rays triangle, k -> k == flatten entries transpose rays tau);
        ej := flatten entries(vars R * rays dualCone triangle);
        -- drop the entries corresponding to the rays of tau
        -- at least, I think the column ordering respects the order of the rays
        -- I hope
        if not instance(tauidxs, Nothing) then (
            tauidxsComplement := toList(set(0.. numColumns rays triangle - 1)  - set(tauidxs));
            ej = ej_tauidxsComplement
        );
        product apply(ej, i -> 1/i)
    ) 
)
equivariantMultiplicity(Cone, Cone, Ring) := ZZ => (sigma, tau, R) -> (
    if isUnimodular(sigma) then equivariantMultiplicity({sigma}, tau, R) else (
        error "Need to generate a unimodular triangulation of the cone sigma, and then pass it as the third argument to equivariantMultiplicity."
    )
)

isUnimodular = method()
isUnimodular(Cone) := Boolean => (sigma) -> (abs(det rays sigma) == 1)