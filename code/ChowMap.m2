needsPackage "AlgebraicSplines"; needsPackage "NormalToricVarieties"; needsPackage "Polyhedra"; needsPackage "Normaliz"; needsPackage "LLLBases";
topLevelMode = Standard
needs "./helpers.m2"
needs "./Splines.m2"
needs "./MinkowskiWeights.m2"


--------------------- Now methods for the chow map

chowMap = method() 
chowMap(List, ZZ) := RingElement => (splineList,d) -> (
    splineListd := select(splineList, f -> degree f == d);
    transpose matrix (splineListd / chowMap / mat / entries / flatten)
)

chowMap(Spline) := MinkowskiWeight => f -> (
    if isMember(MinkowskiWeight, keys f.cache) then return (f.cache)#MinkowskiWeight else (
    Sigma := fan(f);
    if degree f < dim Sigma - ambDim Sigma then error "Degree of spline too small to produce a nontrivial Minkowski weight.";
    d := degree(f);
    ambDimSigma := ambDim Sigma;
    weightFunction := tau -> (
        if d + dim tau - ambDimSigma < 0 then 0 else chowMap(f, tau)
    );
    result := minkowskiWeight(Sigma, weightFunction, d);
    f.cache#MinkowskiWeight = result;
    result
    )
)
chowMap(Spline, Cone) := RingElement => (f, tau) -> (
    R := ring(f);
    SigmaFaces := cones(f);
    --tauFaces := select(SigmaFaces, face -> contains(face, tau));
    tauFaces := select(SigmaFaces, triangle -> isFace(tau, triangle));
    -*for tauFace in tauFaces do (
        << "Cone containing tau: " << rays tauFace << endl;
    );*-
    -- TODO: This calculations equivariant multiplicity even when spline is zero. Should fix this for speed.
    sum for tauFace in tauFaces list restriction(f,tauFace)*equivariantMultiplicity(tauFace, tau, R)
)
chowMap(List, List, List, Cone) := RingElement => (f, V, F, tau) -> chowMap(spline(f,V,F), tau)


equivariantMultiplicity = method()
-- This computes the equivariant multiplicity of a unimodular triangulation at a cone tau

equivariantMultiplicity(List, Cone, Ring) := ZZ => (unimodularTriangulation, tau, R) -> (
    -- Only need to compute multiplicities if the triangle contains tau
    trianglesIntersectingTau := select(unimodularTriangulation, triangle -> (dim intersect(triangle, tau)) == dim tau);
    trianglesContainingTau := select(unimodularTriangulation, triangle -> contains(triangle, tau));
    --- The below is trying to account for situations outside of Lemma 2.7 -- tau intersects sigma but does not contain it. In this case, this method overcounts....
    -- Honestly, a little bit confused about this but it seems to work. If there is a problem, check here.
    diffTriangles := length trianglesIntersectingTau - length trianglesContainingTau;
    latticeVolume := if diffTriangles != 0 then diffTriangles else 1;

    sum for triangle in trianglesIntersectingTau list (
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
            (product apply(ej, i -> 1/i))/ latticeVolume
        ) else (
            1 / latticeVolume
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
                if length(interiorHilbRays) == 0 then error "No interior rays found and not unimodular, should implement higher translate.";
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
    -*result := hashTable for k from 0 to dim Sigma list (
        k => for sigma in facesAsCones(1, Sigma) list if isMember(UnimodularTriangulation, keys sigma.cache) then sigma.cache#UnimodularTriangulation else findUnimodularTriangulation(sigma)
    );*-
    for sigma in maxFacesAsCones(Sigma) list (
        if isMember(UnimodularTriangulation, keys sigma.cache) then sigma.cache#UnimodularTriangulation else findUnimodularTriangulation(sigma)
    ) -- so this will save the unimodular triangulation of each maximal cone in the cache.
)

simplicialization = method()
simplicialization(Fan) := Fan => (Sigma) ->  (
    if isSimplicial Sigma then Sigma else (
        
        firstNonsimplicialCodim := first select(1, reverse toList(0..dim Sigma), k -> any(facesAsCones(k, Sigma) / isSimplicial, bool -> bool == false));
        simplicialization fan(flatten join(facesAsCones(firstNonsimplicialCodim, Sigma) / simplicialization))
    )
)
simplicialization(Cone) := List => (sigma) -> (
    if isSimplicial sigma then {sigma} else (
        barycenter := (transpose matrix{(sum entries transpose rays sigma)})/(numColumns rays sigma);
        subdivision := facesAsCones(0,stellarSubdivision(fan sigma, barycenter));
        subdivision
    )
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

