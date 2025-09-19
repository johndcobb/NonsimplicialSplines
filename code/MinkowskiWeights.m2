

-------------
MinkowskiWeight = new Type of MutableHashTable

minkowskiWeight = method(
    TypicalValue => MinkowskiWeight
)

cones MinkowskiWeight := List => (mw) -> (mw.cache)#Cones

mat(MinkowskiWeight) := Matrix => (mw) -> (
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
    -- ensure weight table exists in cache
    if not isMember("Weights", keys mw.cache) then mw.cache#Weights = new MutableHashTable from {};
    wtTable := mw.cache#Weights;

    -- return cached value if available
    if isMember(tau, keys wtTable) then (
        wtTable#tau
    ) else (
        -- compute on demand via stored weightFunction closure and cache it
        val := mw.weightFunction(tau);
        wtTable#tau = val;
        val
    )
)

minkowskiWeight(Fan, Function, ZZ) := MinkowskiWeight => (Sigma, f, coneCodim) -> (
    -- weightFunction should take in a cone and return an integer
    -- We should check the balancing condition here.
    new MinkowskiWeight from {
        weightFunction => f, 
        cache => new MutableHashTable from 
            {Fan => Sigma, 
            Cones => facesAsCones(coneCodim, Sigma), Codimension => coneCodim}}
)
minkowskiWeight(Fan, Matrix, ZZ) := MinkowskiWeight => (Sigma, M, coneCodim) -> (
    if numrows M != #facesAsCones(coneCodim,Sigma) then error "Matrix must have a column for each cone of the specified codimension.";
    localHash := (hashTable for p in pairs facesAsCones(coneCodim,Sigma) list p_1 => M_(p_0,0));
    chowmap := (tau -> (values selectKeys(localHash, sigma -> sigma == tau))_0);
    minkowskiWeight(Sigma, chowmap, coneCodim)
)

pullback(FanMap, MinkowskiWeight) := MinkowskiWeight => {} >> opts -> (psi, mw) -> (
    k := dim(psi#source) - dim mw;
    assert(isStrict(psi));
    assert(fan mw == target psi); 
    weightFunction := tau -> (
        tau' := (imageCones(psi, tau))_0;
        if ambDim tau' - dim tau' != k then 0 else (
            --- here i need to compute the lattice index.
            phiMat := psi#map;
            -- columns generating psi(N') are the columns of phiMat
            NtauPrime := rays tau'; -- matrix whose columns are ray generators of tau' (in N)
            S := phiMat | NtauPrime; -- matrix with columns generating the sublattice
            -- compute smith normal form; request only D when supported
            snf := smithNormalForm(S, ChangeMatrix => {false, false});
            D := if class snf === List then snf#1 else snf;
            -- number of target lattice rows (ambient dimension)
            n := ambDim(psi#target);
            -- product of the first n diagonal entries (if D smaller, treat missing as 1)
            diagProduct := 1;
            for i from 0 to (n-1) do (
                entry := if i < numrows D and i < numcols D then D_(i,i) else 1;
                diagProduct = diagProduct * entry;
            );
            (abs diagProduct) * weight(mw, tau')
        )
    );
    minkowskiWeight(psi#source, weightFunction, k)
)

MinkowskiWeight + MinkowskiWeight := MinkowskiWeight => (mw1, mw2) -> (
    if fan mw1 != fan mw2 or dim mw1 != dim mw2 then error "Minkowski weights must be on the same fan and have the same codimension to be added.";
    new MinkowskiWeight from {weightFunction => (tau -> weight(mw1,tau) + weight(mw2,tau)), cache => new MutableHashTable from {Fan => fan mw1, Cones => cones mw1, Codimension => (mw1.cache)#Codimension}}
)