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
    for facet in maxFacesAsCones(fan(f)) do (
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

    if length gens R != length first VFixed then error "The number of generators in the ring must match the dimension of fan.";

    coneHash := hashTable(for face in FFixed list face => coneFromVData transpose matrix apply(face, idx -> VFixed_idx));

    fCone := inputCone -> (
        coneVertices := (keys selectValues(coneHash, k -> k == inputCone))_0;
        coneNum := position(FFixed, face -> face == coneVertices);
        f_coneNum
    );

    Sigma := fan(V,F);

    return new Spline from {splineFunction => fCone, cache => new MutableHashTable from {Ring => R, Vertices => VFixed, Facets => FFixed, Fan => Sigma, Cones => facesAsCones(0, Sigma), Degree => first max(f / degree)}};
)

restriction = method()
restriction(Spline, Cone) := RingElement => (f, sigma) -> (
    f.splineFunction(sigma)
)

removeOrigin = method()
removeOrigin(List, List) := Sequence => (V,F) -> (
    k := length(V_0);
    zeroIdx := position(V, v -> v == toList(k:0));
    if zeroIdx === null then (V,F) else (
    VFixed := delete(toList(k:0),V);
    FFixed := apply(F, face -> apply(delete(zeroIdx, face), idx -> if idx < zeroIdx then idx else idx-1));
    (VFixed, FFixed))
)

addOrigin = method()
addOrigin(List, List) := Sequence => (V,F) -> (
    zeroVec := toList(ambDim Sigma:0);
    if isMember(zeroVec, V) then (V,F) else (
    VWithZero := append(V, zeroVec);
    FWithZero := apply(F, face -> append(face, numColumns rays Sigma));
    (VWithZero, FWithZero)
    )
)

-* Getting V and F from a fan Sigma:
    V := entries transpose rays Sigma;
    F := maxCones(Sigma);
*-


-- This takes in V and F required to compute splines, and then creates the fan after removing the zero vector.
fan(List, List) := Fan => (V,F) -> (
    (VFixed, FFixed) := removeOrigin(V,F);
    fan(apply(FFixed, C -> transpose matrix apply(C, idx -> VFixed_idx)) / coneFromVData)
)

polyhedralComplex(List, List) := PolyhedralComplex => (V,F) -> (
    (VFixed, FFixed) := removeOrigin(V,F);
    polyhedralComplex(apply(FFixed, C -> transpose matrix apply(C, idx -> VFixed_idx)) / convexHull)
)

splineModule(Fan, ZZ) := Matrix => opts -> (Sigma, r) -> (
    V := entries transpose rays Sigma;
    F := maxCones(Sigma);
    (VWithZero, FWithZero) := addOrigin(V,F);
    splineModule(VWithZero, FWithZero, r, opts)
)

-------------
MinkowskiWeight = new Type of MutableHashTable

minkowskiWeight = method(
    TypicalValue => MinkowskiWeight
)

cones MinkowskiWeight := List => (mw) -> (mw.cache)#Cones

minkowskiMatrix = method()
minkowskiMatrix(MinkowskiWeight) := List => (mw) -> (
    if isMember(Mat, keys mw.cache) then return (mw.cache)#Mat else (
        result := transpose matrix{apply(cones mw, c -> weight(mw, c))};
        (mw.cache)#Mat = result;
        result
    )
)

fan MinkowskiWeight := Fan => (mw) -> mw.cache#Fan
dim(MinkowskiWeight) := ZZ => (mw) -> dim fan mw - (mw.cache)#Codimension

weight = method()
weight(MinkowskiWeight, Cone) := ZZ => (mw, tau) -> (
    mw.weightFunction(tau)
)

minkowskiWeight(Fan, FunctionClosure, ZZ) := MinkowskiWeight => (Sigma, chowmap, coneCodim) -> (
    -- weightFunction should take in a cone and return an integer
    -- We should check the balancing condition here.
    

    new MinkowskiWeight from {weightFunction => chowmap, cache => new MutableHashTable from {Fan => Sigma, Cones => facesAsCones(coneCodim,Sigma), Codimension => coneCodim}}
)

--------------------- Now methods for the chow map
splineList = (Splines, V, F, R) -> for splineCol in entries transpose generators Splines list spline(splineCol, V, F, R)

chowMap = method() 
chowMap(List, ZZ) := RingElement => (splineList,d) -> (
    splineListd := select(splineList, f -> degree f == d);
    transpose matrix (splineListd / chowMap / minkowskiMatrix / entries / flatten)
)

chowMap(Spline) := MinkowskiWeight => f -> (
    Sigma := fan(f);
    if degree f < dim Sigma - ambDim Sigma then error "Degree of spline too small to produce a nontrivial Minkowski weight.";
    d := degree(f);
    ambDimSigma := ambDim Sigma;
    weightFunction := tau -> (
        if d + dim tau - ambDimSigma < 0 then 0 else chowMap(f, tau)
    );
    minkowskiWeight(Sigma, weightFunction, d)
)
chowMap(Spline, Cone) := RingElement => (f, tau) -> (
    R := ring(f);
    SigmaFaces := cones(f);
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
        dualMatrix := rays dualCone triangle;
        tauidxs := select(toList(0..numColumns dualMatrix - 1), i -> transpose rays tau * dualMatrix_i == 0);
        ej := flatten entries(vars R * dualMatrix);
        -- drop the entries corresponding to the rays of tau
        -- at least, I think the column ordering respects the order of the rays
        -- I hope
        if length(tauidxs) > 0 then (
            -- comment out Complement if you want integers, but this doesn't align with code.
            --tauidxsComplement := toList(set(0.. numColumns rays triangle - 1)  - set(tauidxs));
            ej = ej_tauidxs;
            product apply(ej, i -> 1/i)
        ) else (
            1
        )
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

maxFacesAsCones = method()
maxFacesAsCones(Fan) := List => (Sigma) -> (
    V := entries transpose rays Sigma;
    (for maxCone in maxCones(Sigma) list (V_maxCone)) / transpose / matrix / coneFromVData
)
maxFacesAsCones(Cone) := List => (sigma) -> (
    V := entries transpose rays sigma;
    (for maxCone in maxCones(sigma) list (V_maxCone)) / transpose / matrix / coneFromVData
)


findUnimodularTriangulation = method()

findUnimodularTriangulation(Cone) := List => (sigma) -> (
    if isMember(UnimodularTriangulation,keys sigma.cache) then return (sigma.cache)#UnimodularTriangulation else (
    if not isSimplicial(sigma) then (
        -- need to compute a barycentric subdivision and then find unimodular triangles in those!
        barycenter := (transpose matrix{(sum entries transpose rays sigma)})/(numColumns rays sigma);
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

numRays := (sigma) -> length entries transpose rays sigma

simplicialLift = method()
simplicialLift(Cone) := Cone => (sigma) -> (
    if isSimplicial(sigma) then sigma else (
        disc := numRays sigma - dim sigma;
        heightFunction := table(numRays sigma, disc, (i,j) -> random(1,10));
        liftedCone := coneFromVData(rays sigma || transpose matrix heightFunction);

        -- if the lifted cone does not have the correct dimension, try again.
        if (numRays liftedCone - dim liftedCone) != 0 then (
            simplicialLift(sigma)
        )
    )
)
simplicialLift(Fan) := Fan => (Sigma) -> (
    if isSimplicial Sigma then Sigma else (
        V := entries transpose rays Sigma;
        F := maxCones(Sigma);
        maxDiscrepancy := max apply(F, face -> length face - dim Sigma);
        heightFunction := table(length V, maxDiscrepancy, (i,j) -> random(maxDiscrepancy+2)); -- the max height is chosen mostly arbitrarily
        liftedV := entries ( matrix V |  matrix heightFunction);
        fan(liftedV, F)
    )
)
