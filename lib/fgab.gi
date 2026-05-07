#################################################################
# groups

# FGAb( m )
# It returns the FGAb.
InstallMethod( FGAb, "for a MatrixOrMatrixObj", [ IsMatrixOrMatrixObj ], 
	function(m)
		local fam, A, gens;

		m := Matrix(Integers, m);
		fam := NewFamily("FGAbElementsFam", IsFGAbElement);
		fam!.type := NewType(fam, IsFGAbElement and IsFGAbElementRep);
		fam!.matrix := m;
		fam!.snf := SmithNormalFormIntegerMatTransforms(m);
		A := Objectify( 
			NewType(CollectionsFamily( fam ), IsFGAb and IsFGAbRep),
			rec()
		);
		gens := List([1..DimensionsMat(m)[1]], 
			function(i)
				local v;
				
				v := ZeroVector(Integers, DimensionsMat(m)[1]);
				v[i] := 1;
				v := ObjByExtRep(fam, v);

				return v;
			end 
		);
		SetUnderlyingMatrix(A, m);
		SetZero( A, ObjByExtRep(fam, ZeroVector(Integers, DimensionsMat(m)[1])) );
		#SetGeneratorsOfAdditiveGroup(A, gens);

		return A;
	end 
);

# FGAb( v )
# It returns the FGAb.
# v is the diagonal entries.
InstallMethod( FGAb, "for a Vector", [ IsVector ], 
	function(v)
		local m;

		if Length(v) = 0 then
			m := ZeroMatrix(Integers, 0, 0);
		else
			m := Matrix(Integers, DiagonalMat(v));
		fi;

		return FGAb(m);
	end 
);

# FreeFGAb( n )
# It returns the FGAb.
# n is the dimension.
InstallMethod( FreeFGAb, "for a Int", [ IsInt ], 
	function(n)
		local m;

		if n < 0 then
			Error("FreeFGAb: n < 0.");
		fi;
		m := ZeroMatrix(Integers, n, 0);

		return FGAb(m);
	end 
);

#! @ChapterInfo `FGAb`, Additive Group
#! @Arguments A
#! @Returns a list of generators
InstallMethod( GeneratorsOfAdditiveGroup, "for a FGAb", [ IsFGAb and IsFGAbRep ], 
	function(A)
		local sA, dim_sA, from, fam_sA, gens_sA, gens_A;

		sA := SimplifiedFGAb(A);
		from := FromSimplifiedFGAb(sA);
		fam_sA := FamilyObj(Representative(sA));
		dim_sA := DimensionsMat(UnderlyingMatrix(sA))[1];

		gens_sA := List( [1..dim_sA], 
            function(i)
                local gen;

                gen := ZeroVector(Integers, dim_sA);
                gen[i] := 1;

                return ObjByExtRep(fam_sA, gen);
            end 
        );
        gens_A := List(gens_sA, x -> from(x));

		return gens_A;
	end 
);

InstallMethod( IndependentGeneratorsOfAdditiveGroup, "for a AdditiveGroup of FGAbElementCollection", [ IsAdditiveGroup and IsFGAbElementCollection ], 
    function(A)
        local gens, f, Aparent, gen, Asub, from, rst;

        gen := GeneratorsOfAdditiveGroup(A);
        if IsFGAb(A) then
            # If A is FGAb, then the generators are already independent.
            return gen;
        elif  HasFromSubFGAb(A) then
            # If A is SubFGAb
            f:= FromSubFGAb(A);

            return f(GeneratorsOfAdditiveGroup(Source(f)));
        elif HasParent(A) then
            Aparent := Parent(A);
        else
            Error("IndependentGeneratorsOfAdditiveGroup: A is not FGAb or does't have parent.");
        fi;

        Asub := SubFGAb(Aparent, gen);
        from := FromSubFGAb(Asub);
        rst := from(GeneratorsOfAdditiveGroup(Source(from)));

        return rst;
    end 
);


# It returns a FGAb B:B0 -> B1 that B is rectangular diagonal and injective, and the diagonal entries are all greater than 1. 
InstallMethod( SimplifiedFGAb, "for a FGAb", [ IsFGAb and IsFGAbRep ],
	function(A)
        local matA, n, m, snf, S, P, Q, invP, invQ, k1, k2, i, p, q,
              matD, D, mat_to0, to0, mat_to1, to1, mat_from0, from0, mat_from1, from1,
              sA;

        # 1. Get dimensions of the input MatrixObj A
        # In left-acting convention, A: Z^m -> Z^n means A is an n x m matrix.
		matA := UnderlyingMatrix(A);
        n := DimensionsMat(matA)[1]; # Number of rows (codomain Z^n)
        m := DimensionsMat(matA)[2]; # Number of columns (domain Z^m)
        
        k1 := 0; # Number of 1s on the diagonal
        k2 := 0; # Number of elements > 1 on the diagonal
        
        # 2. Compute Smith Normal Form if dimensions are non-zero
        if n > 0 and m > 0 then
            snf := FamilyObj(Representative(A))!.snf;
            S := snf.normal;
            P := snf.rowtrans;
            Q := snf.coltrans;
            invP := snf.invrowtrans;
            invQ := snf.invcoltrans;
            
            # Count k1 and k2
            for i in [1 .. Minimum(n, m)] do
                if S[i][i] = 1 then
                    k1 := k1 + 1;
                elif S[i][i] > 1 then
                    k2 := k2 + 1;
                else
                    break;
                fi;
            od;
        else
            # Handle cases where the matrix has 0 rows or 0 columns
            if n > 0 then
                P := IdentityMat(n);
                invP := P;
            fi;
            if m > 0 then
                Q := IdentityMat(m);
                invQ := Q;
            fi;
        fi;
        
        # 3. Dimensions for the simplified matrix D
        # In left-acting convention, D = ( \Lambda \\ 0 ) is a q x p matrix.
        p := k2;
        q := n - k1;
        
        # 4. Construct D (q x p)
        if q > 0 and p > 0 then
            matD := List([1..q], i -> List([1..p], j -> 0));
            for i in [1..p] do
                matD[i][i] := S[k1 + i][k1 + i];
            od;
            D := Matrix(Integers, matD);
        else
            D := ZeroMatrix(Integers, q, p);
        fi;
        
        # 5. Construct to1 (Left/Domain map, p x m): Submatrix of Q^-1 (rows for \Lambda, all columns)
        # Commutativity requires: D * to1 = to0 * A
        if p > 0 and m > 0 then
            mat_to1 := List([1..p], i -> List([1..m], j -> invQ[k1 + i][j]));
            to1 := Matrix(Integers, mat_to1);
        else
            to1 := ZeroMatrix(Integers, p, m);
        fi;
        
        # 6. Construct to0 (Right/Codomain map, q x n): Submatrix of P (rows for \Lambda and 0, all columns)
        if q > 0 and n > 0 then
            mat_to0 := List([1..q], i -> List([1..n], j -> P[k1 + i][j]));
            to0 := Matrix(Integers, mat_to0);
        else
            to0 := ZeroMatrix(Integers, q, n);
        fi;
        
        # 7. Construct from1 (Left/Domain map, m x p): Submatrix of Q (all rows, columns for \Lambda)
        # Commutativity requires: A * from1 = from0 * D
        if m > 0 and p > 0 then
            mat_from1 := List([1..m], i -> List([1..p], j -> Q[i][k1 + j]));
            from1 := Matrix(Integers, mat_from1);
        else
            from1 := ZeroMatrix(Integers, m, p);
        fi;
        
        # 8. Construct from0 (Right/Codomain map, n x q): Submatrix of P^-1 (all rows, columns for \Lambda and 0)
        if n > 0 and q > 0 then
            mat_from0 := List([1..n], i -> List([1..q], j -> invP[i][k1 + j]));
            from0 := Matrix(Integers, mat_from0);
        else
            from0 := ZeroMatrix(Integers, n, q);
        fi;

        # 9. Assemble the SimplifiedFGAb

		# The SimplifiedFGAb
		sA := FGAb(D);
		# Set the FGAbHomomorphism from A to sA
		SetToSimplifiedFGAb(sA, FGAbHomomorphismNC(A, sA, to0, to1));
		# Set the FGAbHomomorphism from sA to A
		SetFromSimplifiedFGAb(sA, FGAbHomomorphismNC(sA, A, from0, from1));

		return sA;
	end 
);

# FGAbElement( A, v )
InstallMethod( FGAbElement, "for a FGAb and a Vector", [ IsFGAb, IsVector ],
	function(A, v)
		if DimensionsMat(UnderlyingMatrix(A))[1] <> Length(v) then
			Error("FGAbElement: The length of v is not correct.");
		fi;

		v := Vector(Integers, v);

		return ObjByExtRep(FamilyObj(Representative(A)), v);
	end 
);

InstallMethod( PrintObj, "for a FGAb", [ IsFGAb and IsFGAbRep ], 
	function(A)
		Print("< FGAb ", UnderlyingMatrix(A), ">");
	end 
);

InstallMethod( ViewObj, "for a FGAb", [ IsFGAb and IsFGAbRep ], 
	function(A)
		Print("< FGAb ", UnderlyingMatrix(A), ">");
	end 
);

# FGAbDirectSum( A, B )
InstallMethod( FGAbDirectSum, "for two FGAb", [ IsFGAb and IsFGAbRep, IsFGAb and IsFGAbRep ], 
    function(A, B)
        local Amat, Bmat, rA, cA, rB, cB, rD, cD, mat, m, rst, i, j;

        Amat := UnderlyingMatrix(A);
        Bmat := UnderlyingMatrix(B);
        
        rA := DimensionsMat(Amat)[1]; cA := DimensionsMat(Amat)[2];
        rB := DimensionsMat(Bmat)[1]; cB := DimensionsMat(Bmat)[2];
        rD := rA + rB; cD := cA + cB;

        # Safely construct the block diagonal matrix even if some dimensions are 0
        if rD = 0 or cD = 0 then
            m := ZeroMatrix(Integers, rD, cD);
        else
            mat := NullMat(rD, cD);
            if rA > 0 and cA > 0 then
                for i in [1..rA] do
                    for j in [1..cA] do
                        mat[i][j] := Amat[i, j];
                    od;
                od;
            fi;
            if rB > 0 and cB > 0 then
                for i in [1..rB] do
                    for j in [1..cB] do
                        mat[rA + i][cA + j] := Bmat[i, j];
                    od;
                od;
            fi;
            m := Matrix(Integers, mat);
        fi;

        rst := FGAb(m);
        SetFGAbDirectSumInfo(rst, rec(FGAb := [A, B],
                                      embeddings := [],
                                      projections := []) );

        return rst;
    end 
);

# FGAbDirectSum( list )
InstallMethod( FGAbDirectSum, "for a list of FGAb", [ IsList ], 
    function(L)
        local r, c, rD, cD, mat, m, rst, i, j, offsetR, offsetC, A, Amat;
        
        if IsEmpty(L) then
            Error("List of summands cannot be empty.\n");
        fi;

        # Calculate total dimensions
        rD := Sum(L, A -> DimensionsMat(UnderlyingMatrix(A))[1]);
        cD := Sum(L, A -> DimensionsMat(UnderlyingMatrix(A))[2]);

        if rD = 0 or cD = 0 then
            m := ZeroMatrix(Integers, rD, cD);
        else
            mat := NullMat(rD, cD);
            offsetR := 0;
            offsetC := 0;
            for A in L do
                Amat := UnderlyingMatrix(A);
                r := DimensionsMat(Amat)[1];
                c := DimensionsMat(Amat)[2];
                if r > 0 and c > 0 then
                    for i in [1..r] do
                        for j in [1..c] do
                            mat[offsetR + i][offsetC + j] := Amat[i, j];
                        od;
                    od;
                fi;
                offsetR := offsetR + r;
                offsetC := offsetC + c;
            od;
            m := Matrix(Integers, mat);
        fi;

        rst := FGAb(m);
        SetFGAbDirectSumInfo(rst, rec(FGAb := L,
                                      embeddings := [],
                                      projections := []) );
        return rst;
    end 
);

# SubFGAb( A, gens )
InstallMethod( SubFGAb, "for a FGAb and a List", [ IsFGAb and IsFGAbRep, IsList ], 
	function(A, gens)
		local fam, free, rm, cm, m, f, cokerf, g, kerg, h, rst; 

		if ( not IsFGAbElementCollection(gens) ) and Length(gens) > 0 then
			if IsVectorCollection(gens) then
				fam := FamilyObj(Representative(A));
				gens := List(gens, v -> ObjByExtRep(fam, Vector(Integers, v)));
			else
				Error("SubFGAb: gens must be IsFGAbElementCollection, IsVectorCollection or [].");
			fi;
		fi;

		rm := Length(gens);
		cm := DimensionsMat(UnderlyingMatrix(A))[1];

		free := FreeFGAb(rm);
		if rm = 0 then
			m := ZeroMatrix(Integers, cm, rm);
			gens := [Zero(A)];
		else
			m := Matrix(Integers, List(gens, x -> ExtRepOfObj(x)));
			m := TransposedMatMutable(m);
		fi;
		# f:free -> A
		f := FGAbHomomorphism(free, A, m);
		# g:A -> cokerf
		cokerf := FGAbCokernel(f);
		g := ToFGAbCokernel(cokerf);
		# h:kerg -> A
		kerg := FGAbKernel(g);
		h := FromFGAbKernel(kerg);

		rst := SubadditiveGroup(A, gens);
		SetFromSubFGAb(rst, h);
		
		return rst;
	end 
);

#! @ChapterInfo `FGAb`, Additive Group
#! @Group `FGAbDirectSum`
#! @Arguments A, n
#! @Description `Embedding` can be used.
InstallMethod( Embedding, "FGAbDirectSum and integer",
    [ IsFGAb and IsFGAbRep and HasFGAbDirectSumInfo, IsPosInt ],
    function( D, i )
    local info, summands, rows, rowsD, rows_i, offset, k, j, mat, M, hom;

    info := FGAbDirectSumInfo( D );
    
    # Return cached embedding if it exists
    if IsBound( info.embeddings[i] ) then
        return info.embeddings[i];
    fi;

    summands := info.FGAb;
    if i < 1 or i > Length(summands) then
        Error("Index i is out of bounds for the direct sum.\n");
    fi;

    # Calculate the number of ROWS (generators) for each summand
    rows := List(summands, A -> DimensionsMat(UnderlyingMatrix(A))[1]);
    rowsD := Sum(rows);
    rows_i := rows[i];

    # Calculate the row offset for the i-th summand
    offset := 0;
    for k in [1 .. i-1] do
        offset := offset + rows[k];
    od;

    # Left action: Homomorphism A_i -> D is a (rowsD x rows_i) matrix
    if rowsD = 0 or rows_i = 0 then
        M := ZeroMatrix(Integers, rowsD, rows_i);
    else
        mat := NullMat(rowsD, rows_i);
        for j in [1..rows_i] do
            # Map the j-th generator of A_i to the (offset + j)-th generator of D
            mat[offset + j][j] := 1;
        od;
        M := Matrix(Integers, mat);
    fi;

    hom := FGAbHomomorphismNC(summands[i], D, M);
    
    # Store information
    info.embeddings[i] := hom;
    return hom;
end );

#! @ChapterInfo `FGAb`, Additive Group
#! @Group `FGAbDirectSum`
#! @Arguments A, n
#! @Description `Projection` can be used.
InstallMethod( Projection, "FGAbDirectSum and integer",
    [ IsFGAb and IsFGAbRep and HasFGAbDirectSumInfo, IsPosInt ],
    function( D, i )
    local info, summands, rows, rowsD, rows_i, offset, k, j, mat, M, hom;

    info := FGAbDirectSumInfo( D );
    
    # Return cached projection if it exists
    if IsBound( info.projections[i] ) then
        return info.projections[i];
    fi;

    summands := info.FGAb;
    if i < 1 or i > Length(summands) then
        Error("Index i is out of bounds for the direct sum.\n");
    fi;

    # Calculate the number of ROWS (generators) for each summand
    rows := List(summands, A -> DimensionsMat(UnderlyingMatrix(A))[1]);
    rowsD := Sum(rows);
    rows_i := rows[i];

    # Calculate the row offset for the i-th summand
    offset := 0;
    for k in [1 .. i-1] do
        offset := offset + rows[k];
    od;

    # Left action: Homomorphism D -> A_i is a (rows_i x rowsD) matrix
    if rows_i = 0 or rowsD = 0 then
        M := ZeroMatrix(Integers, rows_i, rowsD);
    else
        mat := NullMat(rows_i, rowsD);
        for j in [1..rows_i] do
            # Map the (offset + j)-th generator of D to the j-th generator of A_i
            mat[j][offset + j] := 1;
        od;
        M := Matrix(Integers, mat);
    fi;

    hom := FGAbHomomorphismNC(D, summands[i], M);
    
    # Store information
    info.projections[i] := hom;
    return hom;
end );

# Intersection2( A1, A2 )
# A1 and A2 are SubFGAb or FGAb.
#! @ChapterInfo `FGAb`, Additive Group
#! @Arguments A1, A2
#! @Returns a `SubFGAb`
#! @Description You can use just `Intersection`. `A1` and `A2` should both `HasParent`, or `Parent(A1)` is `A2`, or `Parent(A2)` is `A1`.
InstallMethod( Intersection2, "for two AdditiveGroup of FGAbElement", [ IsAdditiveGroup and IsFGAbElementCollection, IsAdditiveGroup and IsFGAbElementCollection ],
	function(A1, A2)
		local A1p, A2p, A, f1, f2, pb, f, rst;

		if HasParent(A1) then
			A1p := Parent(A1);
			# f1: A1 -> A
			f1 := FromSubFGAb(A1);
		else
			A1p := A1;
			f1 := IdentityMapping(A1);
		fi;
		if HasParent(A2) then
			A2p := Parent(A2);
			# f2: A2 -> A
			f2 := FromSubFGAb(A2);
		else
			A2p := A2;
			f2 := IdentityMapping(A2);
		fi;
		if not IsIdenticalObj(A1p, A2p) then
			Error("Intersection2: A1 and A2 must be in the same FGAb.");
		fi;
		A := A1p;

		pb := FGAbPullback(f1, f2);
		# f: pb -> A1 -> A
		f := CompositionMapping(f1, FromFGAbPullback(pb, 1));
		rst := SubFGAb(A, f(GeneratorsOfAdditiveGroup(pb)));

		return rst;
	end 
);

# \+( A1, A2 )
# A1 and A2 are SubFGAb or FGAb.
#! @ChapterInfo `FGAb`, Additive Group
#! @Arguments A1, A2
#! @Returns a `SubFGAb` `A1 + A2`
#! @Description `A1` and `A2` should both `HasParent`, or `Parent(A1)` is `A2`, or `Parent(A2)` is `A1`.
InstallMethod( \+, "for two AdditiveGroup of FGAbElement", [ IsAdditiveGroup and IsFGAbElementCollection, IsAdditiveGroup and IsFGAbElementCollection ],
	function(A1, A2)
		local A1p, A2p, A, f1, f2, gens, rst;

		if HasParent(A1) then
			A1p := Parent(A1);
			# f1: A1 -> A
			f1 := FromSubFGAb(A1);
		else
			A1p := A1;
			f1 := IdentityMapping(A1);
		fi;
		if HasParent(A2) then
			A2p := Parent(A2);
			# f2: A2 -> A
			f2 := FromSubFGAb(A2);
		else
			A2p := A2;
			f2 := IdentityMapping(A2);
		fi;
		if not IsIdenticalObj(A1p, A2p) then
			Error("Intersection2: A1 and A2 must be in the same FGAb.");
		fi;
		A := A1p;

        gens := Concatenation(GeneratorsOfAdditiveGroup(A1), GeneratorsOfAdditiveGroup(A2));
		rst := SubFGAb(A, gens);

		return rst;
	end 
);

#! @ChapterInfo `FGAb`, Additive Group
#! @Arguments A
#! @Description Determine whether A is a trivial additive group.
InstallMethod( IsTrivial, "for a FGAb", [ IsFGAb and IsFGAbRep ], 
	function(A)
		local m;

		m := UnderlyingMatrix(SimplifiedFGAb(A));

		return DimensionsMat(m)[1] = 0;
	end 
);

# A1 and A2 are SubFGAb or FGAb.
#! @ChapterInfo `FGAb`, Additive Group
#! @Arguments A1, A2
#! @Description Determine whether `A1 = A2`. `A1` and `A2` should both `HasParent`, or `Parent(A1)` is `A2`, or `Parent(A2)` is `A1`.
InstallMethod( \=, "for two AdditiveGroup of FGAbElement", [ IsAdditiveGroup and IsFGAbElementCollection, IsAdditiveGroup and IsFGAbElementCollection ],
	function(A1, A2)
		local A1p, A2p, A, f1, f2, pb, f, cokerf, g, cokerg, rst;

		if HasParent(A1) then
			A1p := Parent(A1);
			# f1: A1 -> A
			f1 := FromSubFGAb(A1);
		else
			A1p := A1;
			f1 := IdentityMapping(A1);
		fi;
		if HasParent(A2) then
			A2p := Parent(A2);
			# f2: A2 -> A
			f2 := FromSubFGAb(A2);
		else
			A2p := A2;
			f2 := IdentityMapping(A2);
		fi;
		if not IsIdenticalObj(A1p, A2p) then
			Error("\=: A1 and A2 must be in the same FGAb.");
		fi;
		A := A1p;

		pb := FGAbPullback(f1, f2);
		# f: pb = A1\cup A2 -> A1
		f := FromFGAbPullback(pb, 1);
		cokerf := FGAbCokernel(f);
		# g: pb = A1\cup A2 -> A2
		g := FromFGAbPullback(pb, 2);
		cokerg := FGAbCokernel(g);
		# If cokerf and cokerg are both trivial, then A1 = A2 in A.
		
		return IsTrivial(cokerf) and IsTrivial(cokerg);
	end 
);