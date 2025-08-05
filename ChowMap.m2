needsPackage "AlgebraicSplines"; needsPackage "NormalToricVarieties"; needsPackage "Polyhedra";
topLevelMode = Standard


equivariantMultiplicity = method()
-- This computes the equivariant multiplicity of a unimodular triangulation at a cone tau
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
        equivariantMultiplicity(findUnimodularTriangulation(sigma), tau, R)
    )
)

isUnimodular = method()
isUnimodular(Cone) := Boolean => (sigma) -> (abs(det rays sigma) == 1)

findUnimodularTriangulation = method()
findUnimodularTriangulation(Cone) := List => (sigma) -> (
    if isUnimodular(sigma) then {sigma} else (
        numRays := numColumns rays sigma;
        barycenter := transpose matrix{(sum entries transpose rays sigma)/numRays};
        subdivision := facesAsCones(0,stellarSubdivision(fan sigma, barycenter));
        unimodularPartition := partition(triangle -> isUnimodular(triangle), subdivision);
        
        fixedTriangles := {};
        if isMember(false, keys unimodularPartition) then (
            fixedTriangles = join(unimodularPartition#false / findUnimodularTriangulation);
        );
        unimodularPartition#true | fixedTriangles
    )
)