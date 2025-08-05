needsPackage "AlgebraicSplines"; needsPackage "NormalToricVarieties"; needsPackage "Polyhedra";
topLevelMode = Standard


--------- Methods for Splines

Spline = new Type of MutableHashTable

expression Spline := X -> (
    if hasAttribute (X, ReverseDictionary) 
    then expression getAttribute (X, ReverseDictionary) else 
    (describe X)
)
describe Spline := f -> (
    for facet in facesAsCones(0,fan(f)) do (
        << "Facet: " | net rays facet | " Value: " | net f.splineFunction(facet)
        << endl << endl
    );
)

fan Spline := Fan => (f) -> (
    fan(f.cache#Vertices, f.cache#Facets)
)
ring Spline := Ring => (f) -> (
    f.cache#Ring
)
vertices Spline := List => (f) -> (
    f.cache#Vertices
)
facets Spline := List => (f) -> (
    f.cache#Facets
)

spline = method(
    TypicalValue => Spline
)
spline(List, List, List) := Spline => (f, V, F) -> (
    k := length(V_0);
    VFixed := delete(toList(k:0),V);
    zeroIdx := position(V, v -> v == toList(k:0));
    FFixed := apply(F, face -> delete(zeroIdx, face));

    coneHash := hashTable(for face in FFixed list face => coneFromVData transpose matrix apply(face, idx -> VFixed_idx));

    fCone := inputCone -> (
        coneVertices := (keys selectValues(coneHash, k -> k == inputCone))_0;
        coneNum := position(FFixed, face -> face == coneVertices);
        f_coneNum
    );
    return new Spline from {splineFunction => fCone, cache => hashTable{Ring => R, Vertices => VFixed, Facets => FFixed}};
)

restriction = method()
restriction(Spline, Cone) := RingElement => (f, sigma) -> (
    f.splineFunction(sigma)
)


-- This takes in V and F required to compute splines, and then creates the fan after removing the zero vector.
fan(List, List) := Fan => (V,F) -> (
    k := length(V_0);
    VFixed := delete(toList(k:0),V);
    zeroIdx := position(V, v -> v == toList(k:0));
    FFixed := apply(F, face -> delete(zeroIdx, face));
    fan(apply(FFixed, C -> transpose matrix apply(C, idx -> VFixed_idx)) / coneFromVData)
)


--------------------- Now methods for the chow map

chowMap = method() 
chowMap(Spline, Cone) := RingElement => (f, tau) -> (
    R := ring(f);
    Sigma := fan(f);
    SigmaFaces := facesAsCones(0, Sigma);
    tauFaces := select(SigmaFaces, face -> contains(face, tau));
    sum for tauFace in tauFaces list restriction(f,tauFace)*equivariantMultiplicity(tauFace, tau, R)
)
chowMap(List, List, List, Cone) := RingElement => (f, V, F, tau) -> chowMap(spline(f,V,F), tau)


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