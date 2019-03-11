restart
loadPackage "PHCpack";

localPath = currentFileDirectory;
load(concatenate{localPath, "functions.m2"});
tol = 1e-8;

-- Write here values for k and m:
k0 = 1/10;
m0 = 12/10;
-- 

---- PROBLEM DEFINITION ----

R = QQ[k,m][x_0..x_4];
f = 12*(x_0*x_1*x_2*x_3*x_4 - k^2*m)^2 + 9*(x_0*x_1*x_2*x_3 - k^2)^2 + 6*(x_0*x_1*x_2*x_4 - k^2*m)^2 + 6*(x_0*x_1*x_3*x_4 - k*m)^2 + 6*(x_0*x_2*x_3*x_4 - k^2*m)^2 + 6*(x_1*x_2*x_3*x_4 - k*m)^2 + 3*(k^2*m - x_0*x_2*x_4)^2 + 3*(x_1*x_2*x_4 - k*m)^2 + 3*(x_0*x_3*x_4 - k*m)^2 + 3*(x_1*x_3*x_4 - m)^2 + 3*(x_0*x_1 - k)^2 + 3*(x_2*x_3 - k)^2;
f = sub(f, {k => k0, m => m0});
I = diffList(f);

-- Uncomment the following 2 lines in case you want to 
-- compute the degree.

--J = saturate(ideal I, ideal(x_0*x_1*x_2*x_3*x_4));
--degJ = degree J;

degJ = 249;

---- NONSINGULAR CRITICAL POINTS AT THE INTERIOR ----

nonSingSolsList = new List;
L = new List;

while #nonSingSolsList < degJ do(
    sols = solveSystem(I);
    singularSols = unique flatten for i to 4 list(zeroFilter(sols, i, tol));
    
    -- Filter non-singular solutions
    nonSingularSols = select(sols, i-> not inList(singularSols, i));
    nonSingularSols = refineSolutions(I, nonSingularSols, 64);
    nonSingularSols = select(nonSingularSols, 
		             i -> all(flatten entries evaluate(polySystem I, i), n -> abs(n) < tol)); 


    -- Add new non-singular solutions
    newSols = select(nonSingularSols, i -> not inList(nonSingSolsList, i));
    nonSingSolsList = join(nonSingSolsList, newSols);
);

-- Filter solution at the interior of the feasible region
L = validSolutions(nonSingSolsList);


---- CRITICAL POINTS AT THE BOUNDARY ----

BI = BoundaryIndices(gens R);
for i from 0 to #BI - 1 do(
    S1 = BI_i_0; -- Define set S1
    for j from 0 to #BI_i_1 - 1 do(
	S2 = BI_i_1_j; -- Define set S2
	if #S1 == 0 and #S2 == 0 then continue; -- If both S1 and S2 are empty, then F = f
	F = sub(sub(f, for k to #S1-1 list (S1_k => 1)), 
	               for k to #S2-1 list (S2_k => -1/3));
	if #S1 + #S2 < length gens R then ( 
	    difF = diffList(F);
	    sols = solveSystem(difF);
	    sols = validSolutions(sols);
	    if #sols == 0 then continue;
	    feasibleBoundPoints = unique apply(sols, s -> CompletePoint(S1, 1, S2, -1/3, s, gens R));    
	)else feasibleBoundPoints = {CompletePoint(S1, 1, S2, -1/3, {}, gens R)};
    	L = join(L, feasibleBoundPoints);
    );
);
	
	
---- GLOBAL MINIMUM ----

(globMin, dist) = FindingMinimum(f, L);

print(concatenate("The global minimum is the point ", toString coordinates globMin,"\n", " and the distance to the Stochastic Phylogenetic Variety is dist = ", toString dist));

