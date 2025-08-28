areIsomorphic = (X1, X2) -> (
	dim X1 == dim X2 and length rays X1 == length rays X2 and length max X1 == length max X2 and any(
		permutations length rays X1, 
		perm -> (all(max X1, kone -> isMember(set apply(kone, i -> perm#i), apply(max X2, set)))) and
		-- don't need to do other direction assuming toric varieties are well defined with no repeated cones
		-- because a permutation is bijective, so the map on cones is injective (no two cones can be permuted to the same cone),
		-- and because the number of cones for X1 is same as X2,
		-- it is also surjective, so it is a bijection on the sets of cones
		( -- check that there is a matrix sending v_i to w_sigma(i) for all i
			V = matrix rays X1;
			W = matrix rays X2;
			Wsig = (permMatrix perm) * W;
			T = solve(V, Wsig, MaximalRank=>true); -- 'solve' uses row reduction over ZZ. Already implemented in Macaulay2 (InvariantRing package)
			diffMatrix = V*T - Wsig;
			diffMatrix == 0 and (try inverse T then true else false) -- to make sure inverse is also defined over ZZ
		)
	)
)

genRow = (i, n) -> (
    l = {};
    for ind from 0 to n-1 do (
        if ind == i then l = append(l, 1) else l = append(l, 0);
    );
    l
)

permMatrix = perm -> (
    n = length perm;
    mat = {};
    for i in perm do (
        mat = append(mat, genRow(i, n));
    );
    matrix mat
)

ToricLinearSeries = new Type of HashTable
ToricLinearSeries.synonym = "toric linear series"

toricLinearSeries = method(TypicalValue => ToricLinearSeries)
toricLinearSeries ToricDivisor := d -> toricLinearSeries monomials d
toricLinearSeries List := ToricLinearSeries => m -> (
    if #m == 0 then error "toricLinearSeries expects a nonempty list";
    d := degree m#0;
    if not all(m, x -> d == degree x) then error "toricLinearSeries expects a list of monomials of the same degree";
    s := new ToricLinearSeries from {
        "monomials" => m,
        "degree" => d,
        "variety" => variety ring m#0,
        "projectiveSpace" => (if (#m - 1)>0 then toricProjectiveSpace (#m - 1) else null)
    };
    s
)

monomials ToricLinearSeries := List => o -> s -> (
    s#"monomials"
)

degree ToricLinearSeries := s-> (
    s#"degree"
)
-- the monomial map that defines a toric map
-- monomials ToricMap := List => o -> f -> first entries matrix inducedMap f

--isComplete = method(TypicalValue => Boolean)

isComplete ToricLinearSeries := linSeries -> (
    m := monomials linSeries;
    setM := set m;
    d := degree linSeries;
    if #m != #setM then return false;
    degDMonomials := flatten entries basis(d, ring m#0);
    setM == set degDMonomials
)

baseLocusIdeal = method(TypicalValue => Ideal)

baseLocusIdeal ToricLinearSeries := linSeries -> (
    ideal monomials linSeries
)

isBasepointFree = method(TypicalValue => Boolean)

isBasepointFree ToricLinearSeries := linSeries -> (
    m := monomials(linSeries);
    S := ring m#0;
    B := ideal(S.variety);
    I := radical (ideal m, Strategy => Monomial);
    isSubset(B,I)
) 

variety ToricLinearSeries :=
normalToricVariety ToricLinearSeries := linSeries -> (
    linSeries#"variety"
)


-- -- helper for listing monomials of given degree in the ring
-- -- TODO: move to Core
-- monomials(ZZ,   Ring) :=
-- monomials(List, Ring) := List => o -> (d, S) -> first entries basis(d, S)


-- getting map of tori from a divisor or linear series
map(NormalToricVariety, NormalToricVariety, ToricLinearSeries) :=
map(NormalToricVariety, NormalToricVariety, ToricDivisor) := ToricMap => opts -> (Y, X, D) -> map(Y, X, monomials D)

map(NormalToricVariety, NormalToricVariety, List) := ToricMap => opts -> (Y, X, L) -> (
    divmap := matrix transpose(first \ exponents \ L); -- map CDiv Y -> CDiv X    
    D := (divmap * matrix rays Y)//(matrix rays X);
    map(Y, X, transpose D) -- map    M_Y -> M_X
)

toricMap = method(TypicalValue => ToricMap)
toricMap ToricLinearSeries := linSeries -> (
    m := monomials linSeries;
    X := variety linSeries;
    if #m == 1 then error "toricMap expects a toric linear series of length at least 2";
    Pn := toricProjectiveSpace (#m - 1);
    map(Pn, X, m)
)

idealOfImage = method(Options => true)

idealOfImage ToricLinearSeries := {TargetRing => null} >> o -> linSeries ->(
    m := monomials linSeries;
    T := if o.TargetRing === null then ring linSeries#"projectiveSpace" else o.TargetRing;
    matrixOfMonomials := matrix transpose(first \ exponents \ m);
    K := kernelLLL matrixOfMonomials;
    returnPositive := n -> if n >= 0 then n else 0; -- returns n if n is positive
    I := ideal ((entries transpose K) / (c -> T_(c / returnPositive)  - T_((-c )/ returnPositive)));
    saturate(I, T_(m/(t->1))) -- saturate with respect to coordinate axes of Proj T
)

-- completeLinearSeriesFromDivisor Divisor := D ->(
-- )

-- if the source is not given, get the variety from the linear series or divisor
map(NormalToricVariety, Nothing, ToricLinearSeries) :=
map(NormalToricVariety, Nothing, ToricDivisor) := ToricMap => opts -> (Y, X, D) -> map(Y, variety D, monomials D)
map(NormalToricVariety, Nothing, List)         := ToricMap => opts -> (Y, X, L) -> map(Y, variety ring L#0, L)

-- if target is not given, use the appropriate projective space
map(Nothing, NormalToricVariety, ToricLinearSeries) :=
map(Nothing, NormalToricVariety, ToricDivisor) := ToricMap => opts -> (Y, X, D) -> map(, X, monomials D)
map(Nothing, NormalToricVariety, List)         := ToricMap => opts -> (Y, X, L) -> map(toricProjectiveSpace(#L - 1), X, L)

-- if neither source or target is given, deduce both of them!
map ToricLinearSeries :=
map ToricDivisor      := ToricMap => opts -> D -> map(, variety D, monomials D)
map List              := ToricMap => opts -> L -> map(, variety ring L#0, L)

needsPackage "NormalToricVarieties"

importFrom_Core "concatCols"
importFrom_NormalToricVarieties "cartierCoefficients"

-----------------------------------------------------------------------------
-- helpers

vecs = m -> entries transpose m

-- given a list return a permutation w such that L_w is increasingly ordered
-- e.g. perm {2,1,4} gives {1,0,2}
--perm = L -> last \ sort(reverse \ toList pairs L)

-- given a permutation and a list, permute elements of L according to w
-- e.g. move_{1,0,2} {2,0,1} gives {2,1,0}
move = (w, L) -> apply(L, i -> position(w, j -> i == j))

-- given the Cox ring of a projective bundle, give the indices of
-- rays corresponding to the base and the fiber of the bundle
-- TODO: this is not very robust and depends on a particular basis of Pic X
-- it should instead be determined by selecting the rays in the fan
baseVars  = S -> select(numgens S, i -> last degree S_i == 0)
fiberVars = S -> select(numgens S, i -> last degree S_i == 1)

-- ad-hoc check for two toric varieties being the same up to permutation of the rays
-- compare X and Y, possibly with a permutation of the coordinates given by perm
sameToricVariety  = (X, Y) -> X === Y or fan X == fan Y or any(permutations(d := dim X), perm -> sameToricVariety'(X, Y, id_(ZZ^d)_perm))
sameToricVariety' = (X, Y, m) -> try (f := inducedMap map(X, Y, m)) then f ideal X == ideal Y else false

-----------------------------------------------------------------------------
-- Projectivization of a sum of line bundles on a toric variety

protect Fiber
protect ProjectiveBundle

-- see CLS 7.3
PP = method()
PP CoherentSheaf := E -> E.cache.ProjectiveBundle ??= (
    X := variety E;
    r := rank E - 1; -- a rank r+1 sheaf gives a PP^r bundle on X
    if not isFreeModule module E then error "projectivization of arbitrary sheaves if not yet implemented";
    if r == 0 then return X;
    -- see CLS pp. 337
    p := inverse fromWDivToCl X * fromPicToCl X;
    -- TODO: should the degrees be sorted?
    a := apply(-degrees E, deg -> first vecs(p * transpose matrix {deg}));
    --
    d := dim X;
    N0 := ZZ^d; -- lattice of the base X
    N1 := ZZ^r; -- lattice of the fiber PP^r
    i0 := (id_N0 ++ 0 * id_N1)_{0 ..< d};
    i1 := (0 * id_N0 ++ id_N1)_{d ..< d+r};
    B := transpose matrix rays X;
    -- see CLS (7.3.2)
    Facets := affineImage_i1 \ facesAsCones(0, fan toricProjectiveSpace r);
    -- see CLS (7.3.3)
    Sigmas := apply(#rays X,
	rho -> i0 * B_{rho} + i1 * sum toList apply(1 .. r, k -> (a_k_rho - a_0_rho) * N1_{k-1}));
    -- This one is kind of out of order
    Y0 := normalToricVariety fan flatten table(max X, r + 1,
	(ell, i) -> coneFromVData(concatCols Sigmas_ell) + Facets_i);
    -- TODO: this part doesn't _always_ work
    eff0 := inverse nefGenerators Y0 * effGenerators Y0;
    -- FIXME: the columns are not sorted the way one might expect
    perm := sortColumns eff0;
    Y := normalToricVariety((rays Y0)_perm, move_perm \ max Y0,
	CoefficientRing => X.cache.CoefficientRing,
	Variable        => X.cache.Variable,
	WeilToClass     => eff0_perm);
    Y.cache.Base = X;
    Y.cache.Fiber = E;
    Y)


-- TODO: this should project the fan in some way
-- TODO: if X is not a fiber bundle, should this be an error instead?
base NormalToricVariety  := X -> if X.cache.?Base  then X.cache.Base  else X
-- TODO: technically not always correct
fiber = method()
fiber NormalToricVariety := X -> if X.cache.?Fiber then X.cache.Fiber else X.cache.CoefficientRing^(dim X + 1)

