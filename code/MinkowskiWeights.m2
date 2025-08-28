

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