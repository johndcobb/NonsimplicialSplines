mat = method()


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

maxFacesAsCones = method()
maxFacesAsCones(Fan) := List => (Sigma) -> (
    if isMember(maxCones, keys Sigma.cache) then (Sigma.cache)#maxCones else (
        V := entries transpose rays Sigma;
        result := (for maxCone in maxCones(Sigma) list (V_maxCone)) / transpose / matrix / coneFromVData;
        (Sigma.cache)#maxCones = result;
        result
    )
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

splineModule(Fan, ZZ) := Matrix => opts -> (Sigma, r) -> (
    V := entries transpose rays Sigma;
    F := maxCones(Sigma);
    (VWithZero, FWithZero) := addOrigin(V,F);
    splineModule(VWithZero, FWithZero, r, opts)
)


faceRing = method()
faceRing(Fan) := Ring => Sigma -> (
    X := normalToricVariety(Sigma);
    B := ideal X;
    quotient dual monomialIdeal B
)
faceRing(List, List) := Ring => (V, F) -> faceRing(fan(V,F))

numRays = (sigma) -> length entries transpose rays sigma




FanMap = new Type of MutableHashTable


map(Fan, Fan, Matrix) := FanMap => opts -> (Sigma2, Sigma1, phi) -> (
    (N, L) := (target phi, source phi);
    if ambDim Sigma1 != rank L then error "Lattice dimension of source does not match ambient dimension of source fan.";
    if ambDim Sigma2 != rank N then error "Lattice dimension of target does not match ambient dimension of target fan.";
    if not mapsConestoCones(Sigma2, Sigma1, phi) then error "Map does not send cones to cones.";
    new FanMap from {source => Sigma1, target => Sigma2, map => phi}
)

--isStrict = method()
isStrict(FanMap) := Boolean => f -> all(maxFacesAsCones(f#source), sigma -> length imageCones(f,sigma) == 1)

mat = method()
mat FanMap := Matrix => (f) -> f#map

target FanMap := Fan => f -> f#target
source FanMap := Fan => f -> f#source

facesAsCones(Fan) := List => Sigma -> flatten for k from 0 to dim Sigma  list facesAsCones(k,Sigma)

--- Let phi: Delta' --> Delta. This collects the minimal cones tau such that phi(sigma) is contained in tau. If the morphism is strict, there will be exactly one such cone for each sigma. 
imageCones = method()
imageCones(FanMap, Cone) := Cone => (phi, sigma) -> (
    imagePhi := affineImage(phi#map, sigma);
   -- select(maxFacesAsCones(phi#target), tau -> contains(tau, imagePhi))
   -- gather every cone of the target by taking faces of each maximal cone

    -- all cones of the target fan
    allTargetCones := facesAsCones(phi#target);

    -- cones that contain the image
    candidates := select(allTargetCones, tau -> contains(tau, imagePhi));

    -- keep only those candidates that are minimal w.r.t. inclusion
    minimalCandidates := select(candidates, c -> not any(candidates, d -> (not (d == c)) and contains(c, d)));

    minimalCandidates
)

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





-*
Graveyard
*-

-*
This is a mostly broke method that was attempting to compute a simplicialization in a higher dimension 

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
        
        --if (isSimplicial liftedFan) and isComplete liftedFan then (
        --    liftedFan
        --) else (
        --    simplicialLift(Sigma) -- try again if not simplicial or not complete
        --)
        
    )
)
*-

---- This is unused but is slightly interesting
lawrenceLift = method()
lawrenceLift(List, List, ZZ) := Sequence => (V,F,r) -> (
    n := length V;
    d := length V_0;
    lawrenceV := (flatten table(n, r, (i,j) -> ((basis ZZ^r)_j**(transpose matrix V)_i) || (basis ZZ^n)_i))/entries;
    lawrenceF := apply(F, face -> flatten apply(face, i -> toList(r*i..r*i+(r-1)) ));
    (lawrenceV, lawrenceF)
)

saveTriangulation = method()
saveTriangulation(List, String, String) := () => (triangulation, filename, pwd) -> (
    (pwd | filename) << toExternalString(apply(triangulation, t -> t / rays)) << close;
)
saveTriangulation(List, String) := () => (triangulation, filename) -> (
    pwd := "/Users/John/Documents/GitHub/NonsimplicialSplines/triangulations/";
    saveTriangulation(triangulation, filename, pwd)
)

loadTriangulation = method()
loadTriangulation(String, String) := List => (filename, pwd) -> (
    apply(value get (pwd | filename), t -> t / coneFromVData)
)
loadTriangulation(String) := List => filename -> (
    pwd := "/Users/John/Documents/GitHub/NonsimplicialSplines/triangulations/";
    loadTriangulation(filename, pwd)
)

loadTriangulation(String, String, Fan) := List => (filename, pwd, Sigma) -> (
    triangulation := loadTriangulation(filename, pwd);
    for p in pairs maxFacesAsCones(Sigma) do (
        (p_1).cache#UnimodularTriangulation = triangulation_(p_0) 
    );
    triangulation
)
loadTriangulation(String, Fan) := List => (filename, Sigma) -> (
    pwd := "/Users/John/Documents/GitHub/NonsimplicialSplines/triangulations/";
    loadTriangulation(filename, pwd, Sigma)
)
