
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
cones Spline := List => (f) -> maxFacesAsCones(fan f)
fan Spline := Fan => (f) -> f.cache#Fan
ring Spline := Ring => (f) -> f.cache#Ring
vertices Spline := List => (f) -> f.cache#Vertices
facets Spline := List => (f) -> f.cache#Facets

mat(Spline) := Matrix => f -> (
    if isMember(Mat, keys f.cache) then return (f.cache)#Mat else (
        result := transpose matrix{apply(cones f, c -> restriction(f, c))};
        (f.cache)#Mat = result;
        result
    )
)

spline = method(
    TypicalValue => Spline
)
spline(List, Fan, Ring) := Spline => (f, Sigma, R) -> (
    (V,F) := (entries transpose rays Sigma, maxCones Sigma);

    if length gens R != length first V then error "The number of generators in the ring must match the dimension of fan.";

    coneList := maxFacesAsCones Sigma;
    if length f != length coneList then error "Spline data must have one entry for each maximal cone.";
    fValues := apply(f, g -> promote(g, R));
    restrictionTable := new MutableHashTable from (
        for i from 0 to length coneList - 1 list coneList_i => fValues_i
    );

    fCone := inputCone -> (
        if isMember(inputCone, keys restrictionTable) then restrictionTable#inputCone else (
            coneNum := position(coneList, sigma -> sigma == inputCone);
            if coneNum === null then error "Cone not found in spline fan.";
            restrictionTable#inputCone = fValues_coneNum;
            fValues_coneNum
        )
    );

    return new Spline from {splineFunction => fCone, cache => new MutableHashTable from {Ring => R, Vertices => V, Facets => F, Fan => Sigma, Degree => first max(f / degree), "RestrictionTable" => restrictionTable, "ConeValues" => fValues}};
)
spline(List, List, List, Ring) := Spline => (f, V, F, R) -> (
    (VFixed, FFixed):= removeOrigin(V,F);
    spline(f, fan(VFixed,FFixed), R)
)


restriction = method()
restriction(Spline, Cone) := RingElement => (f, sigma) -> (
    f.splineFunction(sigma)
)

pullback(FanMap, Spline) := Spline => {} >> opts -> (psi, inputSpline) -> (
    assert(isStrict(psi)); -- ensure every cone is sent to only one cone.
    assert(fan inputSpline == target psi); -- ensure the spline is on the source fan.
    pullbackSpline := apply(maxFacesAsCones psi#source, sigma -> restriction(inputSpline,(imageCones(psi, sigma))_0) );
    spline(pullbackSpline, psi#source, R)
)

splineList = method()
splineList(Module, Fan, Ring) := (Splines, Sigma, R) -> for splineCol in entries transpose generators (Splines**R) list spline(splineCol, Sigma, R)
splineList(Module, List, List, Ring) := (Splines, V, F, R) -> (
    Sigma = fan(V,F);
    splineList(Splines, Sigma, R)
)
