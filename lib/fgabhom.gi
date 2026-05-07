#################################################################
# homomorphisms

# FGAbHomomorphismNC( A, B, m )
InstallMethod( FGAbHomomorphismNC, "for two FGAb and a MatrixOrMatrixObj", [ IsFGAb, IsFGAb, IsMatrixOrMatrixObj ],
	function(A, B, m)
		local dims_m, dim_A, dim_B, rst;

		dims_m := DimensionsMat(m);
		dim_A := DimensionsMat(UnderlyingMatrix(A))[1];
		dim_B := DimensionsMat(UnderlyingMatrix(B))[1];
		if dims_m[2] <> dim_A or dims_m[1] <> dim_B then
			Error("FGAbHomomorphismNC: The dimension of m doesn't match.");
		fi;

		m := Matrix(Integers, m);
		# SPGeneralMapping is structure repserving mapping.
		rst := Objectify( TypeOfDefaultGeneralMapping(A, B, IsSPGeneralMapping and IsFGAbHomomorphism and IsFGAbHomomorphismRep), rec());
		SetUnderlyingMatrix(rst, m);
		SetFilterObj(rst, IsMapping);

		return rst;
	end 
);

# FGAbHomomorphism( A, B, m )
# It checks if LiftingMatrix exists.
InstallMethod( FGAbHomomorphism, "for two FGAb and a MatrixOrMatrixObj", [ IsFGAb, IsFGAb, IsMatrixOrMatrixObj ],
	function(A, B, m)
		local rst;

		rst := FGAbHomomorphismNC(A, B, m);
		if LiftingMatrix(rst) = fail then
			Error("FGAbHomomorphism: m is not a FGAbHomomorphism.");
		fi;

		return rst;
	end 
);

# FGAbHomomorphismNC( A, B, m, n )
InstallMethod( FGAbHomomorphismNC, "for two FGAb and two MatrixOrMatrixObj", [ IsFGAb, IsFGAb, IsMatrixOrMatrixObj, IsMatrixOrMatrixObj ],
	function(A, B, m, n)
		local dims_m, dims_n, dims_A, dims_B, rst;

		dims_m := DimensionsMat(m);
		dims_n := DimensionsMat(n);
		dims_A := DimensionsMat(UnderlyingMatrix(A));
		dims_B := DimensionsMat(UnderlyingMatrix(B));
		if dims_m[2] <> dims_A[1] or dims_m[1] <> dims_B[1] or dims_n[2] <> dims_A[2] or dims_n[1] <> dims_B[2] then
			Error("FGAbHomomorphismNC: The dimension of m or n doesn't match.");
		fi;

		m := Matrix(Integers, m);
		n := Matrix(Integers, n);
		# SPGeneralMapping is structure repserving mapping.
		rst := Objectify( TypeOfDefaultGeneralMapping(A, B, IsSPGeneralMapping and IsFGAbHomomorphism and IsFGAbHomomorphismRep), rec());
		SetUnderlyingMatrix(rst, m);
		SetLiftingMatrix(rst, n);
		SetFilterObj(rst, IsMapping);

		return rst;
	end 
);

# FGAbHomomorphism( A, B, m, n )
# It checks if LiftingMatrix exists.
# n is the matrix A1 -> B1.
InstallMethod( FGAbHomomorphism, "for two FGAb and two MatrixOrMatrixObj", [ IsFGAb, IsFGAb, IsMatrixOrMatrixObj, IsMatrixOrMatrixObj ],
	function(A, B, m, n)
		local rst;

		rst := FGAbHomomorphismNC(A, B, m, n);
		A := UnderlyingMatrix(A);
		B := UnderlyingMatrix(B);
		if m * A <> B * n then
			Error("FGAbHomomorphism: m is not a FGAbHomomorphism.");
		fi;

		return rst;
	end 
);

InstallMethod( IdentityMapping, "for a FGAb", [ IsFGAb and IsFGAbRep ],
	function(A)
		local d, m, n;

		d := DimensionsMat(UnderlyingMatrix(A));
		m := IdentityMatrix(Integers, d[1]);
		n := IdentityMatrix(Integers, d[2]);

		return FGAbHomomorphismNC(A, A, m, n);
	end 
);

# It returns the lifting f1:A1 -> B1 of f0:A0 -> B0.
InstallMethod( LiftingMatrix, "for a FGAbHomomorphism", [ IsFGAbHomomorphism and IsFGAbHomomorphismRep ],
	function(map)
		local mat_f, A, B, dimA, nA, mA, dimB, nB, mB, 
			C_raw, BT_raw, BT, L_raw, j, i, c_list, c_vec, l_vec, l_list;

		# Extract matrices before overwriting the map variable
		A := UnderlyingMatrix(Source(map));
		B := UnderlyingMatrix(Range(map));
		mat_f := UnderlyingMatrix(map);

		dimA := DimensionsMat(A); nA := dimA[1]; mA := dimA[2];
		dimB := DimensionsMat(B); nB := dimB[1]; mB := dimB[2];

		# We want to find L (mB x mA) such that B * L = mat_f * A.
		# Since matrices act on the left (column vectors), we have:
		# B * L_{*, j} = (mat_f * A)_{*, j} for each column j.
		# GAP's SolutionIntMat solves x * M = v (right action).
		# Transposing our equation gives: (L_{*, j})^T * B^T = ((mat_f * A)_{*, j})^T.

		# 1. Compute C = mat_f * A (Raw lists)
		if nB > 0 and nA > 0 and mA > 0 then
			C_raw := Unpack(mat_f) * Unpack(A);
		else
			# Create a zero matrix of size nB x mA manually if any dimension is 0
			C_raw := List([1..nB], r -> List([1..mA], c -> 0));
		fi;

		# 2. Prepare B^T as a MatrixObj
		if nB > 0 and mB > 0 then
			BT_raw := TransposedMat(Unpack(B));
			BT := Matrix(Integers, BT_raw);
		else
			BT := ZeroMatrix(Integers, mB, nB);
		fi;

		# 3. Solve for each column of L
		L_raw := List([1..mB], r -> List([1..mA], c -> 0));

		for j in [1..mA] do
			# Extract j-th column of C
			if nB > 0 then
				c_list := List([1..nB], i -> C_raw[i][j]);
				c_vec := Vector(Integers, c_list);
			else
				c_vec := ZeroVector(Integers, 0);
			fi;

			# Solve (L_{*, j})^T * B^T = (C_{*, j})^T
			# Your extended SolutionIntMat takes MatrixObj and VectorObj, returns VectorObj
			l_vec := SolutionIntMat(BT, c_vec);

			if l_vec = fail then
				return fail; # No lifting exists
			fi;

			# Unpack the solution vector
			if mB > 0 then
				l_list := Unpack(l_vec);
			else
				l_list := [];
			fi;

			# Assign to the j-th column of L_raw
			for i in [1..mB] do
				L_raw[i][j] := l_list[i];
			od;
		od;

		# 4. Return L as a MatrixObj
		if mB > 0 and mA > 0 then
			return Matrix(Integers, L_raw);
		else
			return ZeroMatrix(Integers, mB, mA);
		fi;
	end 
);

InstallMethod( ImagesRepresentative, "for a FGAbHomomorphism and a FGAbElement", [ IsFGAbHomomorphism and IsFGAbHomomorphismRep, IsFGAbElement and IsFGAbElementRep ],
	function(f, v)
		local fam, mat, vec;

		fam := FamilyObj(Representative(Range(f)));
		mat := UnderlyingMatrix(f);
		vec := ExtRepOfObj(v);

		return ObjByExtRep(fam, mat * vec);
	end 
);

InstallMethod( ImagesElm, "for a FGAbHomomorphism and a FGAbElement", [ IsFGAbHomomorphism and IsFGAbHomomorphismRep, IsFGAbElement and IsFGAbElementRep ],
	function(f, v)
		local fam, mat, vec;

		fam := FamilyObj(Representative(Range(f)));
		mat := UnderlyingMatrix(f);
		vec := ExtRepOfObj(v);

		return [ObjByExtRep(fam, mat * vec)];
	end 
);

InstallMethod( ImagesSet, "for a FGAbHomomorphism and a AdditiveGroup of FGAbElement", [ IsFGAbHomomorphism and IsFGAbHomomorphismRep, IsAdditiveGroup and IsFGAbElementCollection ],
	function(f, elms)
		local A, B, gens;

        A := Source(f);
        B := Range(f);
        gens := GeneratorsOfAdditiveGroup(elms);
        
        return SubFGAb(B, f(gens));
	end 
);

InstallMethod( ImagesSource, "for a FGAbHomomorphism", [ IsFGAbHomomorphism and IsFGAbHomomorphismRep ],
	function(f)
		local A, B, gens;

        A := Source(f);
        B := Range(f);
        gens := GeneratorsOfAdditiveGroup(A);
        
        return SubFGAb(B, f(gens));
	end 
);

# Preimage of an element is not an additive group but a coset. So we can't use the methods in FGAb straightforwardly.
InstallMethod( PreImagesRepresentative, "for a FGAbHomomorphism and a FGAbElement", [ IsFGAbHomomorphism and IsFGAbHomomorphismRep, IsFGAbElement and IsFGAbElementRep ],
	function(f, elm)
		local A, B, gensC, C, g, gmat, pb, h1, h2, h2mat, vec, sol, gens, rst;

        # f: A -> B
		A := Source(f);
        B := Range(f);
        # g: C = Z -> B
        C := FreeFGAb(1);
        gmat := TransposedMatMutable(Matrix(Integers, [ExtRepOfObj(elm)]));
        g := FGAbHomomorphism(C, B, gmat);
        # h2: pb -> C h1: pb -> A
        pb := FGAbPullback(f, g);
        h1 := FromFGAbPullback(pb, 1);
        h2 := FromFGAbPullback(pb, 2);
        h2mat := UnderlyingMatrix(h2);
        vec := Vector(Integers, [1]);
        sol := SolutionIntMat(TransposedMat(h2mat), vec);
        if sol = fail then
            return fail;
        fi;
        rst := h1(FGAbElement(pb, sol));

		return rst;
	end 
);

# For the preimage of a subgroup, it is always nonempty, since it contains 0.
InstallMethod( PreImagesSet, "for a FGAbHomomorphism and a AdditiveGroup of FGAbElement", [ IsFGAbHomomorphism and IsFGAbHomomorphismRep, IsAdditiveGroup and IsFGAbElementCollection ],
	function(f, elms)
		local A, B, gensC, C, g, pb, h, gens, rst;

        # f: A -> B
		A := Source(f);
        B := Range(f);
        # g: C -> B
        gensC := GeneratorsOfAdditiveGroup(elms);
        C := SubFGAb(B, gensC);
        g := FromSubFGAb(C);
        # h: pb -> A
        pb := FGAbPullback(f, g);
        h := FromFGAbPullback(pb, 1);

        gens := h(GeneratorsOfAdditiveGroup(pb));
        rst := SubFGAb(A, gens);

		return rst;
	end 
);

InstallMethod( PreImagesRange, "for a FGAbHomomorphism", [ IsFGAbHomomorphism and IsFGAbHomomorphismRep ],
	function(f)
		return Source(f);
	end 
);

InstallMethod( CompositionMapping2,"for two FGAbHomomorphism", [ IsFGAbHomomorphism and IsFGAbHomomorphismRep, IsFGAbHomomorphism and IsFGAbHomomorphismRep ],
	function(g, f)
		local gmat, fmat, rst;

		if not IsIdenticalObj(Range(f), Source(g)) then
			Error("CompositionMapping2: Range(f) and Source(g) must be the same.");
		fi;

		gmat := UnderlyingMatrix(g);
		fmat := UnderlyingMatrix(f);

		if HasLiftingMatrix(f) or HasLiftingMatrix(g) then
			rst := FGAbHomomorphismNC(Source(f), Range(g), gmat * fmat, LiftingMatrix(g) * LiftingMatrix(f));
		else
			rst := FGAbHomomorphismNC(Source(f), Range(g), gmat * fmat);
		fi;

		return rst;
	end 
);

InstallMethod( PrintObj, "for a FGAbHomomorphism", [ IsFGAbHomomorphism and IsFGAbHomomorphismRep ], 
	function(f)
		Print("< FGAbHomomorphism ", UnderlyingMatrix(f), " from ", Source(f)," to ", Range(f), ">");
	end 
);

InstallMethod( ViewObj, "for a FGAbHomomorphism", [ IsFGAbHomomorphism and IsFGAbHomomorphismRep ], 
	function(f)
		Print("< FGAbHomomorphism ", UnderlyingMatrix(f), " from ", Source(f)," to ", Range(f), ">");
	end 
);

#! @ChapterInfo `FGAbHomomorphism`, `FGAbHomomorphism`
#! @Group `FGAbKernel`
#! @Arguments f
#! @Description You can use just `Kernel(f)` to get the same result as `KernelOfAdditiveGeneralMapping(f)`. `Kernel(f)` has an attribute `UnderlyingFGAbKernel(A)` which is a output of `FGAbKernel`.
#!
#! `Kernel(f)` is what you should use. Because `Kernel(f)` is actually a subset of `Range(f)`, in particular, a SubFGAb. However, `FGAbKernel(f)` is just a `FGAb` equipped with a canonical map. `FGAbKernel` should be treated as an inner-built function, which is used to get `kernel(f)`.
InstallMethod( KernelOfAdditiveGeneralMapping, "for a FGAbHomomorphism", [ IsFGAbHomomorphism and IsFGAbHomomorphismRep ],
	function(f)
		local ker, fromker, gensker, rst;
		
		ker := FGAbKernel(f);
		fromker := FromFGAbKernel(ker);
		gensker := GeneratorsOfAdditiveGroup(ker);
		
		rst := SubFGAb(Source(f), fromker(gensker));

		return rst;
	end 
);

# FGAbKernel( f )
InstallMethod( FGAbKernel, "for a FGAbHomomorphism", [ IsFGAbHomomorphism and IsFGAbHomomorphismRep ], 
    function(f)
        local A, B, sA, sB, fromsA, tosA, fromsB, tosB, sf, conemat, kermat, ker, fromker0, fromker1, fromker, kertoA, sker, fromsker, skertoA,
              matA, matB, sf0, sf1, dimA0, dimA1, dimB0, dimB1, d1c_list, d1c, snf, r, k, K_list, K, md2c_list, Qinv_md2c, kermat_list, fromker0_list, fromker1_list, i, j, row;

        A := Source(f);
        B := Range(f);
        sA := SimplifiedFGAb(A);
        sB := SimplifiedFGAb(B);
        # sA -> A
        fromsA := FromSimplifiedFGAb(sA);
        # A -> sA
        tosA := ToSimplifiedFGAb(sA);
        # sB -> B
        fromsB := FromSimplifiedFGAb(sB);
        # B -> sB
        tosB := ToSimplifiedFGAb(sB);
        # sA -> sB
        sf := CompositionMapping(tosB, f, fromsA);
        
        # Extract underlying matrices for sA, sB and sf
        matA := UnderlyingMatrix(sA);
        matB := UnderlyingMatrix(sB);
        sf0 := UnderlyingMatrix(sf);
        sf1 := LiftingMatrix(sf);
        
        dimA0 := DimensionsMat(matA)[1];
        dimA1 := DimensionsMat(matA)[2];
        dimB0 := DimensionsMat(matB)[1];
        dimB1 := DimensionsMat(matB)[2];
        
        # Construct d1c = ( d'_1, -f_0 )
        # This maps from B_1 + A_0 to B_0
        d1c_list := [];
        for i in [1..dimB0] do
            row := [];
            for j in [1..dimB1] do
                Add(row, matB[i,j]);
            od;
            for j in [1..dimA0] do
                Add(row, -sf0[i,j]);
            od;
            Add(d1c_list, row);
        od;
        if dimB0 = 0 then
            d1c := ZeroMatrix(Integers, 0, dimB1 + dimA0);
        else
            d1c := Matrix(Integers, d1c_list);
        fi;
        
        # Construct conemat = -d2c = ( f_1 \\ d_1 )
        # cone2 -> cone1 (which is A_1 -> B_1 + A_0)
        md2c_list := [];
        for i in [1..dimB1] do
            row := [];
            for j in [1..dimA1] do
                Add(row, sf1[i,j]);
            od;
            Add(md2c_list, row);
        od;
        for i in [1..dimA0] do
            row := [];
            for j in [1..dimA1] do
                Add(row, matA[i,j]);
            od;
            Add(md2c_list, row);
        od;
        if dimB1 + dimA0 = 0 or dimA1 = 0 then
            conemat := ZeroMatrix(Integers, dimB1 + dimA0, dimA1);
        else
            conemat := Matrix(Integers, md2c_list);
        fi;
        
        # Compute the kernel of d1c using Smith Normal Form
        # P * d1c * Q = S  =>  d1c * Q = P^-1 * S
        snf := SmithNormalFormIntegerMatTransforms(d1c);
        r := 0;
        for i in [1..Minimum(dimB0, dimB1 + dimA0)] do
            if snf.normal[i,i] <> 0 then
                r := r + 1;
            fi;
        od;
        k := dimB1 + dimA0 - r; # Dimension of the kernel
        
        # Extract K (the embedding matrix of ker(d1c) into B_1 + A_0)
        # K consists of the last k columns of Q (snf.coltrans)
        K_list := [];
        for i in [1..dimB1 + dimA0] do
            row := [];
            for j in [1..k] do
                Add(row, snf.coltrans[i, r + j]);
            od;
            Add(K_list, row);
        od;
        if dimB1 + dimA0 = 0 or k = 0 then
            K := ZeroMatrix(Integers, dimB1 + dimA0, k);
        else
            K := Matrix(Integers, K_list);
        fi;
        
        # cone2 = ker1 -> ker0
        # Compute kermat: K * kermat = conemat
        # Since Q^-1 * K = (0 \\ I), kermat is the bottom k rows of Q^-1 * conemat
        Qinv_md2c := snf.invcoltrans * conemat;
        kermat_list := [];
        for i in [1..k] do
            row := [];
            for j in [1..dimA1] do
                Add(row, Qinv_md2c[r + i, j]);
            od;
            Add(kermat_list, row);
        od;
        if k = 0 or dimA1 = 0 then
            kermat := ZeroMatrix(Integers, k, dimA1);
        else
            kermat := Matrix(Integers, kermat_list);
        fi;
        
        # ker
        ker := FGAb(kermat);
        
        # ker0 -> sA0
        # Compute fromker0: projection of K to A_0 (the bottom dimA0 rows of K)
        fromker0_list := [];
        for i in [1..dimA0] do
            row := [];
            for j in [1..k] do
                Add(row, K[dimB1 + i, j]);
            od;
            Add(fromker0_list, row);
        od;
        if dimA0 = 0 or k = 0 then
            fromker0 := ZeroMatrix(Integers, dimA0, k);
        else
            fromker0 := Matrix(Integers, fromker0_list);
        fi;
        
        # ker1 -> sA1
        # Compute fromker1: identity matrix on A_1
        if dimA1 = 0 then
            fromker1 := ZeroMatrix(Integers, 0, 0);
        else
            fromker1_list := List([1..dimA1], i -> List([1..dimA1], j -> 0));
            for i in [1..dimA1] do
                fromker1_list[i][i] := 1;
            od;
            fromker1 := Matrix(Integers, fromker1_list);
        fi;
        
        # ker -> sA
        fromker := FGAbHomomorphismNC(ker, sA, fromker0, fromker1);
        # ker -> A
        kertoA := CompositionMapping(fromsA, fromker);
        # sker
        sker := SimplifiedFGAb(ker);
        # sker -> A
        fromsker := FromSimplifiedFGAb(sker);
        skertoA := CompositionMapping(kertoA, fromsker);
        SetFromFGAbKernel(sker, skertoA);

        return sker;
    end 
);

# FGAbCokernel( f )
InstallMethod( FGAbCokernel, "for a FGAbHomomorphism", [ IsFGAbHomomorphism and IsFGAbHomomorphismRep ], 
    function(f)
        local A, B, sA, sB, fromsA, tosA, fromsB, tosB, sf,
              matA, matB, sf0, sf1, dimA0, dimA1, dimB0, dimB1, 
              d1c_list, d1c, coker, tocoker0_list, tocoker0, tocoker1_list, tocoker1, 
              tocoker, scoker, toscoker, BtosB, Bto_scoker, i, j, row;

        A := Source(f);
        B := Range(f);
        sA := SimplifiedFGAb(A);
        sB := SimplifiedFGAb(B);
        
        # sA -> A
        fromsA := FromSimplifiedFGAb(sA);
        # sB -> B
        fromsB := FromSimplifiedFGAb(sB);
        # B -> sB
        tosB := ToSimplifiedFGAb(sB);
        
        # sA -> sB
        sf := CompositionMapping(tosB, f, fromsA);
        
        # Extract underlying matrices for sA, sB and sf
        matA := UnderlyingMatrix(sA);
        matB := UnderlyingMatrix(sB);
        sf0 := UnderlyingMatrix(sf);
        sf1 := LiftingMatrix(sf);
        
        dimA0 := DimensionsMat(matA)[1];
        dimA1 := DimensionsMat(matA)[2];
        dimB0 := DimensionsMat(matB)[1];
        dimB1 := DimensionsMat(matB)[2];
        
        # Construct d1c = ( d'_1, -f_0 )
        # This maps from B_1 + A_0 to B_0. 
        # The cokernel is exactly presented by this matrix!
        d1c_list := [];
        for i in [1..dimB0] do
            row := [];
            for j in [1..dimB1] do
                Add(row, matB[i,j]);
            od;
            for j in [1..dimA0] do
                Add(row, -sf0[i,j]);
            od;
            Add(d1c_list, row);
        od;
        if dimB0 = 0 then
            d1c := ZeroMatrix(Integers, 0, dimB1 + dimA0);
        else
            d1c := Matrix(Integers, d1c_list);
        fi;
        
        # Create the unsimplified cokernel
        coker := FGAb(d1c);
        
        # Construct the canonical projection map from sB to coker
        # Degree 0 part: Identity matrix on B_0
        if dimB0 = 0 then
            tocoker0 := ZeroMatrix(Integers, 0, 0);
        else
            tocoker0_list := List([1..dimB0], i -> List([1..dimB0], j -> 0));
            for i in [1..dimB0] do
                tocoker0_list[i][i] := 1;
            od;
            tocoker0 := Matrix(Integers, tocoker0_list);
        fi;
        
        # Degree 1 part: Maps B_1 to the B_1 component of B_1 + A_0
        # It is a block matrix [ I_{dimB1} \\ 0_{dimA0 x dimB1} ]
        if dimB1 + dimA0 = 0 or dimB1 = 0 then
            tocoker1 := ZeroMatrix(Integers, dimB1 + dimA0, dimB1);
        else
            tocoker1_list := [];
            for i in [1..dimB1+dimA0] do
                row := [];
                for j in [1..dimB1] do
                    if i = j then
                        Add(row, 1);
                    else
                        Add(row, 0);
                    fi;
                od;
                Add(tocoker1_list, row);
            od;
            tocoker1 := Matrix(Integers, tocoker1_list);
        fi;
        
        # sB -> coker
        tocoker := FGAbHomomorphismNC(sB, coker, tocoker0, tocoker1);
        
        # Simplify the cokernel
        scoker := SimplifiedFGAb(coker);
        # coker -> scoker
        toscoker := ToSimplifiedFGAb(scoker);
        
        # B -> sB -> coker -> scoker
        Bto_scoker := CompositionMapping(toscoker, tocoker, tosB);
        
        # Set the attribute (assuming SetToFGAbCokernel is defined similarly to SetFromFGAbKernel)
        SetToFGAbCokernel(scoker, Bto_scoker);

        return scoker;
    end 
);

# FGAbPullback( f1, f2 )
InstallMethod( FGAbPullback, "for two FGAbHomomorphism", 
    [ IsFGAbHomomorphism and IsFGAbHomomorphismRep, IsFGAbHomomorphism and IsFGAbHomomorphismRep ], 
    function(f1, f2)
        local f1mat, f2mat, m, m_list, r, c1, c2, i, j, row,
              A1, A2, B, A1A2, f1f2, pb, frompbtoA1A2, fromA1A2toA1, fromA1A2toA2, p1, p2;

        if not IsIdenticalObj(Range(f1), Range(f2)) then
            Error("FGAbPullback: Range(f1) and Range(f2) must be the same.");
        fi;

        f1mat := UnderlyingMatrix(f1);
        f2mat := UnderlyingMatrix(f2);
        
        r  := DimensionsMat(f1mat)[1];
        c1 := DimensionsMat(f1mat)[2];
        c2 := DimensionsMat(f2mat)[2];

        # Construct the horizontal block matrix m = [ f1, -f2 ]
        if r = 0 or c1 + c2 = 0 then
            m := ZeroMatrix(Integers, r, c1 + c2);
        else
            m_list := [];
            for i in [1..r] do
                row := [];
                # Append elements from f1
                for j in [1..c1] do Add(row, f1mat[i,j]); od;
                # Append elements from -f2
                for j in [1..c2] do Add(row, -f2mat[i,j]); od;
                Add(m_list, row);
            od;
            m := Matrix(Integers, m_list);
        fi;

        A1 := Source(f1);
        A2 := Source(f2);
        B := Range(f1);
        A1A2 := FGAbDirectSum(A1, A2);
        
        # m = (f1, -f2): A1A2 -> B induced by f1 and f2.
        f1f2 := FGAbHomomorphismNC(A1A2, B, m);
        pb := FGAbKernel(f1f2);
        
        # pb -> A1A2 uses FromFGAbKernel(pb)
        frompbtoA1A2 := FromFGAbKernel(pb);
        # A1A2 -> Ai uses Projection(A1A2, i)
        fromA1A2toA1 := Projection(A1A2, 1);
        fromA1A2toA2 := Projection(A1A2, 2);
        
        # p1: pb -> A1
        p1 := CompositionMapping(fromA1A2toA1, frompbtoA1A2);
        # p2: pb -> A2
        p2 := CompositionMapping(fromA1A2toA2, frompbtoA1A2);
        
        SetFGAbPullbackInfo(pb, rec(projection := [p1, p2]));

        return pb;
    end 
);

InstallMethod( FromFGAbPullback, "for a FGAb and a Int", [ IsFGAb and IsFGAbRep, IsInt ],
	function(A, n)
		return FGAbPullbackInfo(A).projection[n];
	end 
);

# FGAbPushout( f1, f2 )
InstallMethod( FGAbPushout, "for two FGAbHomomorphism", 
    [ IsFGAbHomomorphism and IsFGAbHomomorphismRep, IsFGAbHomomorphism and IsFGAbHomomorphismRep ], 
    function(f1, f2)
        local f1mat, f2mat, m, m_list, r1, r2, c, i, j, row,
              A, B1, B2, B1B2, f1f2, po, to_po, emb1, emb2, i1, i2;

        if not IsIdenticalObj(Source(f1), Source(f2)) then
            Error("FGAbPushout: Source(f1) and Source(f2) must be the same.");
        fi;

        f1mat := UnderlyingMatrix(f1);
        f2mat := UnderlyingMatrix(f2);
        
        r1 := DimensionsMat(f1mat)[1];
        c  := DimensionsMat(f1mat)[2];
        r2 := DimensionsMat(f2mat)[1];

        # Construct the vertical block matrix m = [  f1 ]
        #                                         [ -f2 ]
        if r1 + r2 = 0 or c = 0 then
            m := ZeroMatrix(Integers, r1 + r2, c);
        else
            m_list := [];
            # Add rows from f1
            for i in [1..r1] do
                row := [];
                for j in [1..c] do Add(row, f1mat[i,j]); od;
                Add(m_list, row);
            od;
            # Add rows from -f2
            for i in [1..r2] do
                row := [];
                for j in [1..c] do Add(row, -f2mat[i,j]); od;
                Add(m_list, row);
            od;
            m := Matrix(Integers, m_list);
        fi;

        A := Source(f1);
        B1 := Range(f1);
        B2 := Range(f2);
        B1B2 := FGAbDirectSum(B1, B2);
        
        # m = (f1, -f2)^T : A -> B1B2 induced by f1 and f2.
        f1f2 := FGAbHomomorphismNC(A, B1B2, m);
        po := FGAbCokernel(f1f2);
        
        # B1B2 -> po uses ToFGAbCokernel(po)
        to_po := ToFGAbCokernel(po);
        # Bi -> B1B2 uses Embedding(B1B2, i)
        emb1 := Embedding(B1B2, 1);
        emb2 := Embedding(B1B2, 2);
        
        # i1: B1 -> po
        i1 := CompositionMapping(to_po, emb1);
        # i2: B2 -> po
        i2 := CompositionMapping(to_po, emb2);
        
		SetFGAbPushoutInfo(po, rec(embedding := [i1, i2]));

        return po;
    end 
);

InstallMethod( ToFGAbPushout, "for a FGAb and a Int", [ IsFGAb and IsFGAbRep, IsInt ],
	function(A, n)
		return FGAbPushoutInfo(A).embedding[n];
	end 
);