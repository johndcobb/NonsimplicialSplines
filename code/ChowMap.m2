needsPackage "AlgebraicSplines"; needsPackage "NormalToricVarieties"; needsPackage "Polyhedra"; needsPackage "Normaliz"; needsPackage "LLLBases";
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
spline(List, Fan, Ring) := Spline => (f, Sigma, R) -> (
    (V,F) := (entries transpose rays Sigma, maxCones Sigma);

    if length gens R != length first V then error "The number of generators in the ring must match the dimension of fan.";

    coneHash := hashTable(for face in F list face => coneFromVData transpose matrix apply(face, idx -> V_idx));

    fCone := inputCone -> (
        coneVertices := (keys selectValues(coneHash, k -> k == inputCone))_0;
        coneNum := position(F, face -> face == coneVertices);
        promote(f_coneNum, R)
    );

    return new Spline from {splineFunction => fCone, cache => new MutableHashTable from {Ring => R, Vertices => VFixed, Facets => FFixed, Fan => Sigma, Cones => maxFacesAsCones(Sigma), Degree => first max(f / degree)}};
)
spline(List, List, List, Ring) := Spline => (f, V, F, R) -> (
    (VFixed, FFixed):= removeOrigin(V,F);
    spline(f, fan(VFixed,FFixed), R)
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
    k := length(V_0);
    zeroVec := toList(k:0);
    if isMember(zeroVec, V) then (V,F) else (
    VWithZero := append(V, zeroVec);
    FWithZero := apply(F, face -> append(face, length V));
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

mat = method()
mat(MinkowskiWeight) := List => (mw) -> (
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
splineList = method()
splineList(Module, Fan, Ring) := (Splines, Sigma, R) -> for splineCol in entries transpose generators (Splines**R) list spline(splineCol, Sigma, R)
splineList(Module, List, List, Ring) := (Splines, V, F, R) -> (
    Sigma = fan(V,F);
    splineList(Splines, Sigma, R)
)

chowMap = method() 
chowMap(List, ZZ) := RingElement => (splineList,d) -> (
    splineListd := select(splineList, f -> degree f == d);
    transpose matrix (splineListd / chowMap / mat / entries / flatten)
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
-*
maxFacesAsCones(Cone) := List => (sigma) -> (
    V := entries transpose rays sigma;
    (for maxCone in maxCones(sigma) list (V_maxCone)) / transpose / matrix / coneFromVData
)*-


faceRing = method()
faceRing(Fan) := Ring => Sigma -> (
    X := normalToricVariety(Sigma);
    B := ideal X;
    quotient dual monomialIdeal B
)
faceRing(List, List) := Ring => (V, F) -> faceRing(fan(V,F))


findUnimodularTriangulation = method()

findUnimodularTriangulation(Cone) := List => (sigma) -> (
    if isMember(UnimodularTriangulation,keys sigma.cache) then return (sigma.cache)#UnimodularTriangulation else (
    if not isSimplicial(sigma) then (
        -- need to compute a barycentric subdivision and then find unimodular triangles in those!
        result1 := flatten join(simplicialization(sigma) / findUnimodularTriangulation);
        (sigma.cache)#UnimodularTriangulation = result1;
        result1
        )
    else (
        if isUnimodular(sigma) then {sigma} else (
                hilbRays := getHilbRays(sigma);
                interiorHilbRays := select(hilbRays, r -> not contains(facesAsCones(ambDim sigma- 1, sigma), r)) / rays;
                if length(interiorHilbRays) < 0 then error "No interior rays found and not unimodular, should implement higher translate.";
                -- That is, I need to get a ray from a higher translate of the cone.
                subdivisionHilb := facesAsCones(0,stellarSubdivision(fan sigma, interiorHilbRays_0));
                unimodularPartition := partition(triangle -> isUnimodular(triangle), subdivisionHilb);
                
                correctTriangles := if isMember(true, keys unimodularPartition) then unimodularPartition#true else {};
                fixedTriangles := {};
                if isMember(false, keys unimodularPartition) then (
                    fixedTriangles = flatten join(unimodularPartition#false / findUnimodularTriangulation);
                );
                result2 := join(correctTriangles,fixedTriangles);
                (sigma.cache)#UnimodularTriangulation = result2;
                result2
            )
        )
    )
)
findUnimodularTriangulation(Fan) := Fan => (Sigma) -> (
    unimodularCones := flatten apply(maxFacesAsCones(Sigma), sigma -> findUnimodularTriangulation(sigma));
    fan unimodularCones
)

simplicialization = method()
simplicialization(Fan) := Fan => (Sigma) ->  (
    if isSimplicial Sigma then Sigma else (
        fan(flatten join(maxFacesAsCones(Sigma) / simplicialization))
    )
)
simplicialization(Cone) := List => (sigma) -> (
    if isSimplicial sigma then {sigma} else (
        barycenter := (transpose matrix{(sum entries transpose rays sigma)})/(numColumns rays sigma);
        subdivision := facesAsCones(0,stellarSubdivision(fan sigma, barycenter));
        subdivision
    )
)

numRays = (sigma) -> length entries transpose rays sigma

simplicialLift = method()
simplicialLift(Fan) := Fan => (Sigma) -> (
    if isSimplicial Sigma then Sigma else (
        V := entries transpose rays Sigma;
        F := maxCones(Sigma);
        maxDiscrepancy := max apply(F, face -> length face - dim Sigma);
        heightFunction := table(length V, maxDiscrepancy, (i,j) -> random(-maxDiscrepancy-2,maxDiscrepancy+2)); -- the max height is chosen mostly arbitrarily
        liftedV := entries ( matrix V |  matrix heightFunction);
        liftedFan := fan(liftedV, F);
        liftedFan
        -*
        if (isSimplicial liftedFan) and isComplete liftedFan then (
            liftedFan
        ) else (
            simplicialLift(Sigma) -- try again if not simplicial or not complete
        )
        *-
    )
)

lawrenceLift = method()
lawrenceLift(List, List, ZZ) := Sequence => (V,F,r) -> (
    n := length V;
    d := length V_0;
    lawrenceV := (flatten table(n, r, (i,j) -> ((basis ZZ^r)_j**(transpose matrix V)_i) || (basis ZZ^n)_i))/entries;
    lawrenceF := apply(F, face -> flatten apply(face, i -> toList(r*i..r*i+(r-1)) ));
    (lawrenceV, lawrenceF)
)

-- Finds the chow group A_k(X) of any toric variety X with fan Sigma.
chowGroup = method()

-- This is a matrix of integers which gives the relations between the cones of codimension k in Sigma. The ordering is given by facesAsCones(k, Sigma). So to compute the chowGroup from ths presentation, you just take the cokernel of the transpose of this matrix.
chowGroup(Fan, ZZ) := Matrix => (Sigma, k) -> (
    -- Collect cones of codimension k+1 (taus) and codimension k (sigmas).
    -- We will produce one relation row per tau and one column per sigma.
    tauClasses := facesAsCones(k+1, Sigma);
    sigmaClasses := facesAsCones(k, Sigma);
    numClasses := length sigmaClasses;

    -- If there are no (k+1)-cones then there are no relations: return zero matrix.
    if tauClasses == {} then (
        return map(ZZ^numClasses, ZZ^numClasses, 0)
    ) else (
    -- For each tau build a row of the relation matrix.  Each entry is obtained
    -- by evaluating basis elements of M_tau on the primitive vector that extends
    -- tau to an adjacent sigma, and then dividing by the lattice index to
    -- account for non-primitive extensions.
    matrix flatten for tau in tauClasses list (
            -- Mtau is a Z-basis for the lattice M restricted to the sublattice
            -- orthogonal to tau; these are the integer linear functionals we use
            -- to evaluate extension vectors v.
            Mtau := entries transpose kernelLLL ((transpose rays tau));
            Ntau := if dim tau != 0 then rays tau else map(ZZ^(dim Sigma), ZZ^1, 0);
            table(Mtau, sigmaClasses, (u, sigma) -> (
                if contains(sigma, tau) then (
                    -- Nsigma and Ntau are the ray matrices (columns are generators)
                    Nsigma := rays sigma;

                    -- Find a generator v for the quotient lattice N_sigma / N_tau.
                    -- kernelLLL(Nsigma | Ntau) returns integer relations; we pick
                    -- the appropriate column (vidx) and use it to form v in the
                    -- ambient coordinates of Nsigma.
                    
                    sigmaRays := set entries transpose Nsigma;
                    tauRays := set entries transpose Ntau;
                    vRays := toList(sigmaRays - tauRays);

                    if #vRays == 0 and dim sigma > dim tau then error "Could not find an extending ray from tau to sigma.";
                    if #vRays == 0 then return 0; -- Case where sigma == tau
                    v := transpose matrix{vRays_0};
                    

                    -- Compute the lattice index [N_tau + Z v : N_tau] by applying
                    -- Smith normal form to (Ntau | v).  The diagonal entry at
                    -- (rank Ntau, rank Ntau) equals this index (1 if v is primitive).
                    D := smithNormalForm(Ntau | v, ChangeMatrix => {false, false});
                    latticeIdx := D_(rank Ntau, rank Ntau);
        
                    -- Evaluate u on v and divide by the lattice index.  This
                    -- gives the integer relation coefficient for this (tau,sigma)
                    -- pair.  Note: orientation/sign conventions may affect the
                    -- global sign; keep that in mind if results disagree by -1.

                    (flatten entries (matrix{u}*v))_0 // abs latticeIdx
                ) else (
                    0
                )
            ))
        )
    )
)
-*chowGroup(Fan, ZZ) := Sequence => (Sigma, k) -> (
    symbol x;
    tauClasses := facesAsCones(k+1, Sigma);
    sigmaClasses := facesAsCones(k, Sigma);
    R := QQ[x_1..x_(length sigmaClasses)];
    if tauClasses == {} then (
        return R
    ) else (
    var := sigma -> (
        idx := position(sigmaClasses, sigmaclass -> sigmaclass == sigma);
        if idx === null then error "Cone not found in codimkClasses" else R_(idx)
    );
    relationList := flatten for tau in tauClasses list (
        tauDualGens := if dim tau == 0 then (
            apply(entries map(ZZ^(ambDim Sigma), ZZ^(ambDim Sigma), 1), column -> coneFromVData transpose matrix{column}) 
            ) else (
            getHilbRays(dualCone tau)
        );
        conesContainingTau := select(sigmaClasses, sigma -> contains(sigma, tau));
        flatten for u in tauDualGens list (
            sum for sigma in conesContainingTau list (
                dualMatrix := rays dualCone sigma;
                nidx := select(toList(0..numColumns dualMatrix - 1), i -> transpose rays tau * dualMatrix_i == 0);
                flatten entries ((transpose rays u)*dualMatrix_(nidx)*var(sigma))
            )
        )
    );
    (sigmaClasses, relationList))
)*-

operationalChowGroup = method()
operationalChowGroup(Fan, ZZ) := Module => (Sigma, k) -> (
    if isComplete Sigma then (
        Hom(coker transpose chowGroup(Sigma, k), ZZ) -- This is true by Theorem 2.1 in Fulton-Sturmfels
    ) else (
        error("Fan is not complete, need to implement Kimura's inductive method.")
    )
)


FanMap = new Type of MutableHashTable

map(Fan, Fan, Matrix) := FanMap => opts -> (Sigma2, Sigma1, phi) -> (
    (N, L) := (target phi, source phi);
    if ambDim Sigma1 != rank L then error "Lattice dimension of source does not match ambient dimension of source fan.";
    if ambDim Sigma2 != rank N then error "Lattice dimension of target does not match ambient dimension of target fan.";
    if not mapsConestoCones(Sigma2, Sigma1, phi) then error "Map does not send cones to cones.";
    new FanMap from {source => Sigma1, target => Sigma2, map => phi}
)

mat FanMap := Matrix => (f) -> f#map


getHilbRays = method()
getHilbRays(Cone) := List => sigma -> (
    hilbBasis := entries ((normaliz(transpose rays sigma, "integral_closure"))#"gen") ;
    hilbRays := apply(hilbBasis, b -> coneFromVData transpose matrix{b});
    hilbRays
)


mapsConestoCones = method()
mapsConestoCones(Fan, Fan, Matrix) := Boolean => (Sigma2, Sigma1, phi) -> (
    all(apply(maxFacesAsCones(Sigma1), sigma -> (
        imagePhi := affineImage(phi, sigma);
        any(apply(maxFacesAsCones(Sigma2), tau -> contains(tau, imagePhi)), bool -> bool)
            )), bool -> bool)
    )
