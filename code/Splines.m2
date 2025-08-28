
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

splineList = method()
splineList(Module, Fan, Ring) := (Splines, Sigma, R) -> for splineCol in entries transpose generators (Splines**R) list spline(splineCol, Sigma, R)
splineList(Module, List, List, Ring) := (Splines, V, F, R) -> (
    Sigma = fan(V,F);
    splineList(Splines, Sigma, R)
)