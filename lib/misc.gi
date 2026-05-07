##################################################################
# extension of some functions to MatrixObj

InstallMethod( SolutionIntMat, "for a MatrixObj and a VectorObj", [ IsMatrixObj, IsVectorObj ],
    function(m, v)
        local dim, r, c, lv, mat_list, vec_list, sol;
        
        dim := DimensionsMat(m);
        r := dim[1];
        c := dim[2];
        lv := Length(v);
        
        # 1. Check dimension compatibility
        # For equation x * m = v, the number of columns in m must equal the length of v.
        if c <> lv then
            Error("SolutionIntMat: matrix columns and vector length must be equal.\n");
        fi;
        
        # 2. Handle degenerate case where the matrix has 0 rows (0 x c matrix)
        if r = 0 then
            # The length of the unknown vector x must be 0.
            # x * m will always yield a c-dimensional zero vector.
            if c = 0 then
                # c = 0, v is clearly a 0-length vector, exactly matched.
                return ZeroVector(Integers, 0); 
            else
                # c > 0. We must check if v is an all-zero vector.
                # Since c > 0, the vector v is not empty and can be safely unpacked.
                vec_list := Unpack(v); 
                if ForAll(vec_list, val -> val = 0) then
                    return ZeroVector(Integers, 0);
                else
                    return fail;
                fi;
            fi;
        fi;
        
        # 3. Handle degenerate case where the matrix has 0 columns and r > 0 (r x 0 matrix)
        if c = 0 then
            # The equation x * m = v becomes: (length r vector x) * (r x 0 matrix) = (empty vector v).
            # This equation holds for any vector x of length r.
            # We naturally return a length r all-zero vector as a particular solution.
            return ZeroVector(Integers, r);
        fi;
        
        # 4. Normal case: neither rows nor columns are 0, safe to use Unpack
        mat_list := Unpack(m);
        vec_list := Unpack(v);
        
        sol := SolutionIntMat(mat_list, vec_list);
        
        if sol = fail then
            return fail;
        else
            return Vector(Integers, sol);
        fi;
    end 
);

InstallMethod( SmithNormalFormIntegerMat, "for a MatrixObj", [ IsMatrixObj ],
    function(m)
        local dim, r, c, mat_list, snf_list;
        
        dim := DimensionsMat(m);
        r := dim[1];
        c := dim[2];
        
        # 1. Handle degenerate cases (r = 0 or c = 0)
        # The Smith Normal Form of an empty matrix is simply an empty matrix 
        # of the exact same dimensions. No unpacking is needed.
        if r = 0 or c = 0 then
            return ZeroMatrix(Integers, r, c);
        fi;
        
        # 2. Normal case: safe to use Unpack
        mat_list := Unpack(m);
        snf_list := SmithNormalFormIntegerMat(mat_list);
        
        return Matrix(Integers, snf_list);
    end 
);

InstallMethod( SmithNormalFormIntegerMatTransforms, "for a MatrixObj", [ IsMatrixObj ],
    function(m)
        local dim, r, c, mat_list, trans_rec, 
              S_obj, P_obj, Q_obj, P_inv_obj, Q_inv_obj,
              P_inv_list, Q_inv_list;
        
        dim := DimensionsMat(m);
        r := dim[1];
        c := dim[2];
        
        # 1. Handle degenerate cases (r = 0 or c = 0)
        if r = 0 or c = 0 then
            # The Smith Normal Form S is an r x c zero matrix
            S_obj := ZeroMatrix(Integers, r, c);
            
            # The row transform matrix P must be an r x r identity matrix.
            # Its inverse is also an r x r identity matrix.
            if r = 0 then
                P_obj := ZeroMatrix(Integers, 0, 0);
                P_inv_obj := ZeroMatrix(Integers, 0, 0);
            else
                P_obj := Matrix(Integers, IdentityMat(r));
                P_inv_obj := Matrix(Integers, IdentityMat(r));
            fi;
            
            # The column transform matrix Q must be a c x c identity matrix.
            # Its inverse is also a c x c identity matrix.
            if c = 0 then
                Q_obj := ZeroMatrix(Integers, 0, 0);
                Q_inv_obj := ZeroMatrix(Integers, 0, 0);
            else
                Q_obj := Matrix(Integers, IdentityMat(c));
                Q_inv_obj := Matrix(Integers, IdentityMat(c));
            fi;
            
            return rec( normal := S_obj, 
                        rowtrans := P_obj, 
                        coltrans := Q_obj,
                        invrowtrans := P_inv_obj, 
                        invcoltrans := Q_inv_obj );
        fi;
        
        # 2. Normal case: safe to use Unpack
        mat_list := Unpack(m);
        trans_rec := SmithNormalFormIntegerMatTransforms(mat_list);
        
        # Calculate the inverses of the transformation matrices.
        # Since P and Q are unimodular integer matrices, their inverses 
        # are guaranteed to be integer matrices.
        P_inv_list := trans_rec.rowtrans^-1;
        Q_inv_list := trans_rec.coltrans^-1;
        
        # Wrap the resulting matrices from the record back into MatrixObjs
        S_obj := Matrix(Integers, trans_rec.normal);
        P_obj := Matrix(Integers, trans_rec.rowtrans);
        Q_obj := Matrix(Integers, trans_rec.coltrans);
        P_inv_obj := Matrix(Integers, P_inv_list);
        Q_inv_obj := Matrix(Integers, Q_inv_list);
        
        return rec( normal := S_obj, 
                    rowtrans := P_obj, 
                    coltrans := Q_obj,
                    invrowtrans := P_inv_obj, 
                    invcoltrans := Q_inv_obj );
    end 
);

#####################################################################
# GAP bug or extension of functions.

# GAP doesn't have left action of a MatrixObj on a VectorObj.
InstallMethod( \*,
  [ "IsPlistMatrixRep", "IsPlistVectorRep" ],
  function( M, v )
    local i, res, s;
    
    # 1. 维度校验：矩阵的列数 (RLPOS 即 Row Length) 必须等于向量的长度
    if ValueOption( "check" ) <> false and
       ( not IsIdenticalObj( M![BDPOS], v![BDPOS] ) or
         M![RLPOS] <> Length( v ) ) then
      Error( "<M> and <v> are not compatible" );
    fi;
    
    # 2. 初始化结果向量：长度应为矩阵的行数
    res := ZeroVector( NumberRows(M), M![EMPOS] );
    
    # 3. 核心计算：逐行计算点积 (M 的第 i 行 * 向量 v)
    for i in [1..NumberRows(M)] do
        s := M![ROWSPOS][i] * v; 
        if not IsZero(s) then
            res[i] := s; # 直接赋值非零的分量
        fi;
    od;
    
    # 4. 保持不可变性
    if not IsMutable(M) and not IsMutable(v) then
        MakeImmutable(res);
    fi;
    
    return res;
  end );

# Bug: nrgens := GeneratorsOfAdditiveGroup( A );
InstallMethod( ViewObj,
    "for an add. magma-with-zero with generators",
    [ IsAdditiveMagmaWithZero and HasGeneratorsOfAdditiveMagmaWithZero ], 1,
    function( A )
    local nrgens;
	# Here is the bug.
    nrgens := Length(GeneratorsOfAdditiveMagmaWithZero( A ));
    if nrgens = 0 then
      Print( "<trivial additive magma-with-zero>" );
    else
      Print( "<additive magma-with-zero with ",
             Pluralize( nrgens, "generator" ), ">" );
    fi;
    end );

InstallMethod( ViewObj,
    "for an add. magma-with-inverses with generators",
    [ IsAdditiveGroup and HasGeneratorsOfAdditiveGroup ], 1,
    function( A )
    local nrgens;
    nrgens := Length(GeneratorsOfAdditiveGroup( A ));
    if nrgens = 0 then
      Print( "<trivial additive magma-with-inverses>" );
    else
      Print( "<additive magma-with-inverses with ",
             Pluralize( nrgens, "generator" ), ">" );
    fi;
    end );