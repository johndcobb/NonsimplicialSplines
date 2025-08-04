restart
load "../SplinesCode.m2"

V = {{0,0,-1}, {-1,-1,1}, {1,-1,1},  {-1,1,1}, {1,1,1}};
F = {{0,1,2}, {0,1,3}, {0,2,4}, {0,3,4}, {1,2,3,4}};

Sigma = polyhedralComplex(V,F)
Pfaces = applyPairs(faces Sigma, (i,lst) -> ((dim Sigma)-i-1,apply(lst,first)))
R = QQ[x_0,x_1,x_2]

B = billeraComplex(Sigma, R, 1)
Splines = HH_2 B
reduceHilbert hilbertSeries Splines

X = normalToricVariety(V,F)
BV = dual monomialIdeal(ideal X) -- (x_0*x_2*x_3, x_0*x_1*x_4)
simplicialComplex(BV)
stanReisner = coker gens BV

-- how to do the same computation but with the algebraic splines package:
needsPackage "AlgebraicSplines"
V = {{0,0,-1}, {-1,-1,1}, {1,-1,1},  {-1,1,1}, {1,1,1}, {0,0,0}};
F = {{0,1,2,5}, {0,1,3,5}, {0,2,4,5}, {0,3,4,5}, {1,2,3,4,5}};
Splines = splineModule(V,F,0)
---------------------------------------

Fsimplicial = flatten(listToPolyhedra(maxPolyhedra(Sigma),Sigma) / barycentricTriangulation) / convexHull
SigmaSimplicial = polyhedralComplex(Fsimplicial)
SigmaSimplicialFaces = applyPairs(faces SigmaSimplicial, (i,lst) -> ((dim SigmaSimplicial)-i-1,apply(lst,first)))

VSimplicial = vertices SigmaSimplicial
maxConesOfSigma = for L in SigmaSimplicialFaces#2 list (
    for d in L list entries VSimplicial_d
) 
maxConesOfSigma = maxConesOfSigma / matrix
maxConesOfSigma = maxConesOfSigma / coneFromVData
dualConesOfSigma = maxConesOfSigma / dualCone
dualConesOfSigma / rays

BSimplicial = billeraComplex(SigmaSimplicial, R, 1)
SplinesSimplicial = minimalPresentation HH_2 BSimplicial
reduceHilbert hilbertSeries SplinesSimplicial

XSimplicial = normalToricVariety(fan(SigmaSimplicial))
BVSimplicial = dual monomialIdeal(ideal XSimplicial)
SSimplicial = simplicialComplex(BVSimplicial)
stanReisnerSimplicial = (module ring BVSimplicial)/BVSimplicial


------- Now, to lift one dimension higher
V = {{0,0,-1,0}, {-1,-1,1,1}, {1,-1,1,2},  {-1,1,1,3}, {1,1,1,4}};
F = {{0,1,2}, {0,1,3}, {0,2,4}, {0,3,4}, {1,2,3,4}};

SigmaLift = polyhedralComplex(V,F)
PfacesLift = applyPairs(faces SigmaLift, (i,lst) -> ((dim SigmaLift)-i-1,apply(lst,first)))
RLift = QQ[y_0,y_1,y_2,y_3]

BLift = billeraComplex(SigmaLift, RLift, 1)
SplinesLift = HH_2 BLift
minimalPresentation SplinesLift
reduceHilbert hilbertSeries SplinesLift


XLift = normalToricVariety(V,F)
BVLift = dual monomialIdeal(ideal XLift) -- (x_0*x_2*x_3, x_0*x_1*x_4)
simplicialComplex(BVLift)
stanReisnerLift = coker gens BVLift

pr = map(R,RLift,matrix{{x_0,x_1,x_2,0}})
pr(SplinesLift)
Splines
map(Splines, pr(SplinesLift), 1)

restart
debug needsPackage "AlgebraicSplines"
V = {{0,0,-1,0}, {-1,-1,1,1}, {1,-1,1,2},  {-1,1,1,3}, {1,1,1,4}, {0,0,0,0}};
F = {{0,1,2,5}, {0,1,3,5}, {0,2,4,5}, {0,3,4,5}, {1,2,3,4,5}};
SplinesLifted = splineModule(V,F,0) -- error
splineMatrix(V,F,0) -- error
E = getCodim1Intersections(F)
facetEdgeH := apply(#E, e-> positions(F, f-> all(E_e,v-> member(v,f))));
--Compute indices of interior edges, and replace edge list and 
--facet adjacencies to only include these interior edges:
indx := positions(facetEdgeH, i-> #i === 2);
E = E_indx;
facetEdgeH = facetEdgeH_indx;
--Compute top boundary map for complex:
BM := matrix apply(
    facetEdgeH, i-> apply(
    #F, j-> if (
        j === first i) then 1 else if (
        j===last i) then -1 else 0));
--List of forms defining interior codim one faces (raised to (r+1) power)
flist := formsList(V,E,0); -- the issue is here. 

d := #(first V);
V = apply(V, v-> append(v,1));
S = QQ[x_0..x_4]
if opts.Homogenize then (
varlist := vars S;
) else (
varlist = (vars S)|(matrix {{sub(1,S)}});
);
varCol := transpose varlist;
M := (transpose(matrix(S,V)));
mM := numrows M;
minorList := apply(E, e-> gens gb minors(mM,matrix(M_e)|varCol));
if any(minorList, I-> ideal I === ideal 1) then (
    error "Some vertices on entered face are not in codimension 1 face."
    );
flatten apply(minorList, m -> (m_(0,0))^(r+1))



T := diagonalMatrix(flist);
splineM := BM|T;
) else if opts.InputType === "ByLinearForms" then (
    print "Wrong inputs, put in lists of adjacent facets and linear forms and continuity r."
        );
splineM
--------- making E_1 page automatic
SigmaSimplicialFaces = applyPairs(faces SigmaSimplicial, (i,lst) -> ((dim SigmaSimplicial)-i-1,apply(lst,first)))
d = dim SigmaSimplicial
F = vertices SigmaSimplicial
C = matrix{toList((length SigmaSimplicialFaces#d):1)}
B = matrix table(length SigmaSimplicialFaces#(d-2), length SigmaSimplicialFaces#(d-1), (i,j) -> (
    sourceCone = SigmaSimplicialFaces#(d-1)_j;
    targetCone = SigmaSimplicialFaces#(d-2)_i;
    if isSubset(targetCone, sourceCone) then (
        (-1)^(position(sourceCone, i -> i == targetCone_0))
    ) else (
        0
    ))
)
A = matrix table(length SigmaSimplicialFaces#(d-1), length SigmaSimplicialFaces#d, (i,j) -> (
    sourceCone = SigmaSimplicialFaces#(d)_j;
    targetCone = SigmaSimplicialFaces#(d-1)_i;
    if isSubset(targetCone, sourceCone) then (
        1
    ) else (
        0
    ))
) -- need to add orientation

E = matrix table(d+1, 2*(length SigmaSimplicialFaces#(d-2)), (i,j) -> (
    sourceCone = SigmaSimplicialFaces#(d-2)_j;
    0
))

D = matrix table(2*(length SigmaSimplicialFaces#(d-2)), length SigmaSimplicialFaces#(d-1), (i,j) -> (
    sourceCone = SigmaSimplicialFaces#(d-1)_j;
    0
))


F = matrix{
    {1,-1,1,-1,1,-1,1,-1},
    {-1,-1,1,1,-1,-1,1,1},
    {-1,-1,-1,-1,1,1,1,1}
}

prune coker F
prune ker F

E = matrix{  
    {1,1,1,1,1,0,0,1,0,1,1,0,1,1,1,1},
    {0,1,0,-1,0,1,1,1,1,1,0,1,0,-1,0,1},
    {1,0,-1,0,1,1,1,0,1,0,1,1,-1,0,1,0}
}

for k from 0 to 7 list (transpose(matrix{E_(2*k)})*matrix{F_k},transpose(matrix{E_(2*k+1)})*matrix{F_k}) -- should all be zero

prune coker E



crossProduct = (u, v) -> (
    matrix{
        {u_(1)*v_(2) - u_(2)*v_(1)},
        {u_(2)*v_(0) - u_(0)*v_(2)},
        {u_(0)*v_(1) - u_(1)*v_(0)}
    }
)

n = 8

edgeList = {{0,1},{0,2},{0,4}, {1,3}, {1,5}, {2,6}, {2,3}, {3,7}, {4,5}, {4,6}, {5,7}, {6,7}}
DList = for edge in edgeList list (
    i := edge#0;
    j := edge#1;
    eijList := for k from 0 to n-1 list (
        normalij := crossProduct(F_(i), F_(j));
        if k == i then (
            solve(E_{2*i,2*i+1}, normalij)
        )
        else if k == j then (
            -solve(E_{2*j,2*j+1}, normalij) -- negative because of orientation
        ) else (
            matrix{{0},{0}}
        ));
    fold((a,b) -> a || b, eijList)
)
D = fold((a,b) -> a | b, DList)

E*D -- this is the zero map! (so its a complex)

prune ker D
prune (ker E/image(D)) 

