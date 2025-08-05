needsPackage "AlgebraicSplines"; needsPackage "NormalToricVarieties"; needsPackage "Polyhedra";
topLevelMode = Standard


equivariantMultiplicity = method()
equivariantMultiplicity(Cone, Cone, Ring, List) := ZZ => (sigma, tau, R, unimodularTriangulation) -> (

    -- Only need to compute multiplicities if the triangle contains tau
    trianglesContainingTau := for triangle in unimodularTriangulation when contains(triangle, tau) list triangle;

    sum for triangle in trianglesContainingTau list (
        tauidxs = positions(entries transpose rays triangle, k -> k == flatten entries transpose rays tau);
        ej = flatten entries(vars R * rays dualCone triangle);
        -- drop the entries corresponding to the rays of tau
        -- at least, I think the column ordering respects the order of the rays
        -- I hope
        if not instance(tauidxs, Nothing) then (
            for tauidx in tauidxs do (
                ej = drop(ej, {tauidx,tauidx});
            );
        );
        product apply(ej, i -> 1/i)
    ) 
)
equivariantMultiplicity(Cone, Cone, Ring) := ZZ => (sigma, tau, R) -> (
    if isUnimodular(sigma) then equivariantMultiplicity(sigma, tau, R, {sigma}) else (
        error "Need to generate a unimodular triangulation of the cone sigma, and then pass it as the third argument to equivariantMultiplicity."
    )
)

isUnimodular = method()
isUnimodular(Cone) := Boolean => (sigma) -> (abs(det rays sigma) == 1)