needsPackage "AlgebraicSplines"; needsPackage "NormalToricVarieties"; needsPackage "Polyhedra"; needsPackage "Normaliz";
topLevelMode = Standard


--------- Methods for Splines

Spline = new Type of MutableHashTable

-*isWellDefined(Spline) := Boolean => (f) -> (
    -- Get all maximal cones of the fan
    maximalCones := cones(f);

    pairwiseIntersections := apply(subsets(maximalCones, 2), L -> (intersect(L_0, L_1), restriction(f, L_0) - restriction(f, L_1)));

    all(apply(pairwiseIntersections, (face, splineVal) -> (
        if dim face == 0 then true else splineVal % minors(1,vars R * rays face) == 0
    )))
)*-

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

degree Spline := ZZ => (f) -> f.cache#Degree
cones Spline := List => (f) -> f.cache#Cones
fan Spline := Fan => (f) -> f.cache#Fan
ring Spline := Ring => (f) -> f.cache#Ring
vertices Spline := List => (f) -> f.cache#Vertices
facets Spline := List => (f) -> f.cache#Facets

spline = method(
    TypicalValue => Spline
)
spline(List, List, List, Ring) := Spline => (f, V, F, R) -> (
    (VFixed, FFixed):= removeOrigin(V,F);

    coneHash := hashTable(for face in FFixed list face => coneFromVData transpose matrix apply(face, idx -> VFixed_idx));

    fCone := inputCone -> (
        coneVertices := (keys selectValues(coneHash, k -> k == inputCone))_0;
        coneNum := position(FFixed, face -> face == coneVertices);
        f_coneNum
    );

    Sigma := fan(V,F);

    return new Spline from {splineFunction => fCone, cache => hashTable{Ring => R, Vertices => VFixed, Facets => FFixed, Fan => Sigma, Cones => facesAsCones(0, Sigma), Degree => first max(f / degree)}};
)

restriction = method()
restriction(Spline, Cone) := RingElement => (f, sigma) -> (
    f.splineFunction(sigma)
)

removeOrigin = method()
removeOrigin(List, List) := Sequence => (V,F) -> (
    k := length(V_0);
    zeroIdx := position(V, v -> v == toList(k:0));
    VFixed := delete(toList(k:0),V);
    FFixed := apply(F, face -> apply(delete(zeroIdx, face), idx -> if idx < zeroIdx then idx else idx-1));
    (VFixed, FFixed)
)


-- This takes in V and F required to compute splines, and then creates the fan after removing the zero vector.
fan(List, List) := Fan => (V,F) -> (
    (VFixed, FFixed) := removeOrigin(V,F);
    fan(apply(FFixed, C -> transpose matrix apply(C, idx -> VFixed_idx)) / coneFromVData)
)

polyhedralComplex(List, List) := PolyhedralComplex => (V,F) -> (
    (VFixed, FFixed) := removeOrigin(V,F);
    polyhedralComplex(apply(FFixed, C -> transpose matrix apply(C, idx -> VFixed_idx)) / convexHull)
)

--------------------- Now methods for the chow map

chowMap = method() 
chowMap(Spline, Cone) := RingElement => (f, tau) -> (
    R := ring(f);
    SigmaFaces := cones(f);
    tauFaces := select(SigmaFaces, face -> contains(face, tau));
    sum for tauFace in tauFaces list restriction(f,tauFace)*equivariantMultiplicity(tauFace, tau, R)
)
chowMap(List, List, List, Cone) := RingElement => (f, V, F, tau) -> chowMap(spline(f,V,F), tau)
-- TODO: Make a chowMap(Spline) method that returns a new datatype MinkowskiWeight

equivariantMultiplicity = method()
-- This computes the equivariant multiplicity of a unimodular triangulation at a cone tau

equivariantMultiplicity(List, Cone, Ring) := ZZ => (unimodularTriangulation, tau, R) -> (
    -- Only need to compute multiplicities if the triangle contains tau
    trianglesContainingTau := select(unimodularTriangulation, triangle -> contains(triangle, tau));

    sum for triangle in trianglesContainingTau list (
        dualMatrix := rays dualCone triangle;
        tauidxs := select(toList(0..numColumns dualMatrix - 1), i -> transpose rays tau * dualMatrix_i == 0);
        ej := flatten entries(vars R * dualMatrix);
        -- drop the entries corresponding to the rays of tau
        -- at least, I think the column ordering respects the order of the rays
        -- I hope
        if length(tauidxs) > 0 then (
            -- comment out Complement if you want integers, but this doesn't align with code.
            --tauidxsComplement := toList(set(0.. numColumns rays triangle - 1)  - set(tauidxs));
            ej = ej_tauidxs
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
isUnimodular(Cone) := Boolean => (sigma) -> (
    sigmaRays = rays sigma;
    if numColumns sigmaRays != numRows sigmaRays then
        false
    else (
    abs(det rays sigma) == 1)
)


--- wait doing barycenter may not get unimodular triangulation
-- or wait, try rescaling barycenter to be on integer lattice at each step!!
findUnimodularTriangulation = method()

findUnimodularTriangulation(Cone) := List => (sigma) -> (
    if isMember(UnimodularTriangulation,keys sigma.cache) then return (sigma.cache)#UnimodularTriangulation else (
    numRays := numColumns rays sigma;
    if not isSimplicial(sigma) then (
        -- need to compute a barycentric subdivision and then find unimodular triangles in those!
        barycenter := (transpose matrix{(sum entries transpose rays sigma)})/numRays;
        subdivision := facesAsCones(0,stellarSubdivision(fan sigma, barycenter));
        result1 := flatten join(subdivision / findUnimodularTriangulation);
        (sigma.cache)#UnimodularTriangulation = result1;
        result1
        )
    else (
        if isUnimodular(sigma) then {sigma} else (
                hilbBasis := entries ((normaliz(transpose rays sigma, "integral_closure"))#"gen") ;
                hilbRays := apply(hilbBasis, b -> coneFromVData transpose matrix{b});
                interiorHilbRays := select(hilbRays, r -> not contains(facesAsCones(ambDim sigma- 1, sigma), r)) / rays;
                if length(interiorHilbRays) < 0 then error "No interior rays found and not unimodular, should implement higher translate.";
                -- That is, I need to get a ray from a higher translate of the cone.
                subdivisionHilb := facesAsCones(0,stellarSubdivision(fan sigma, interiorHilbRays_0));
                unimodularPartition := partition(triangle -> isUnimodular(triangle), subdivisionHilb);
                
                fixedTriangles := {};
                if isMember(false, keys unimodularPartition) then (
                    fixedTriangles = join(unimodularPartition#false / findUnimodularTriangulation);
                );
                result2 := unimodularPartition#true | fixedTriangles;
                (sigma.cache)#UnimodularTriangulation = result2;
                result2
            )
        )
    )
)
