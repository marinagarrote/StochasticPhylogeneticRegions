
-- Round a complex number up to n decimals
roundCC = (z, n) -> (
    return round(n, realPart z) + round(n, imaginaryPart z)*ii;
);

-- Round a all coordinates of a point up to n decimals
roundPointCC = (p, n) -> (
    return point{for i to #p list(roundCC((coordinates p)_i, n))};
);

-- Check if an element e is contained in a list L
inList = (L, e) -> (
    for i to #L - 1 do(if areEqual(L_i, e) then return true;);
    return false;
);

-- Returns a list with partial derivatives of f defined over CC
diffList = (f) -> (
    variables := support f;
    S := CC[variables];
    DF := diff(matrix {variables}, f);
    DFList := for i from 0 to (numColumns DF - 1) list sub(DF_i_0, S);
    return DFList
);

-- Returns real solutions with coordinates in [-1/3, 1]
validSolutions = (sols) -> (
    realSols := realPoints sols;
    validSols := for i to #realSols-1 list (if not all (coordinates realSols_i, n -> n >= -1/3 and n <= 1) then continue; realSols_i);
    return validSols 
);

-- Returns a point in CC^5. Completes with the coordinates that are in the boundary
CompletePoint = (vars1, value1, vars2, value2, sol, gensR) -> (
    vars3 := toList (set gensR - set vars1 - set vars2);
    newPoint := new MutableList from gensR;
    for i to #vars1 - 1 do newPoint#(index vars1_i) = value1;
    for i to #vars2 - 1 do newPoint#(index vars2_i) = value2;
    for i to #vars3 - 1 do newPoint#(index vars3_i) = roundCC((coordinates sol)_i, 12);
    return toList newPoint
);

-- Returns a list with the border indices
BoundaryIndices = (variables) -> (
   BIaux = subsets variables;
   BI = for i from 0 to #BIaux - 1 list( {BIaux_i, subsets toList (set variables - set BIaux_i)});
   return BI;
);

-- Given a function f and a list of points L = {x_1..x_n} it returns the point
-- x_i such that f(x_i) is minimum and the value f(x_i)
FindingMinimum = (f, L) -> (
    p = null;
    minVal = 1e10;
    for i to #L -1 do (
	val = (evaluate(polySystem{f}, point{L_i}))_(0,0);
	if sub(val,CC) < minVal then (
	    minVal = val;
	    p = point{L_i};
	);
    );
    return (p, minVal);   
);






