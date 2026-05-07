#################################################################
# elements

InstallMethod( ZeroOp, "for a FGAbElement", [ IsFGAbElement and IsFGAbElementRep ],
	function(v)
		local fam, ev, ew;

		fam := FamilyObj(v);
		ev := ExtRepOfObj(v);
		ew := ZeroVector(Integers, Length(ev));

		return ObjByExtRep(fam, ew);
	end 
);

InstallMethod( \+, "for a FGAbElement and a FGAbElement", [ IsFGAbElement and IsFGAbElementRep, IsFGAbElement and IsFGAbElementRep ],
	function(v1, v2)
		local fam1, fam2, ev1, ev2, ev;

		fam1 := FamilyObj(v1);
		fam2 := FamilyObj(v2);
		if not IsIdenticalObj(fam1, fam2) then
			Error("\\+: v1 and v2 must in the same FGAb.");
		fi;
		ev1 := ExtRepOfObj(v1);
		ev2 := ExtRepOfObj(v2);
		ev := ev1 + ev2;
		
		return ObjByExtRep(fam1, ev);
	end 
);

InstallMethod( AdditiveInverseOp, "for a FGAbElement", [ IsFGAbElement and IsFGAbElementRep ],
	function(v)
		local fam, ev;

		fam := FamilyObj(v);
		ev := ExtRepOfObj(v);
		
		return ObjByExtRep(fam, -ev);
	end 
);

InstallMethod( \*, "for an integer and a FGAbElement", [ IsInt, IsFGAbElement and IsFGAbElementRep ],
	function(n, v)
		local fam, ev;

		fam := FamilyObj(v);
		ev := ExtRepOfObj(v);
		
		return ObjByExtRep(fam, n * ev);
	end 
);

InstallMethod( \*, "for a FGAbElement and an integer", [ IsFGAbElement and IsFGAbElementRep, IsInt ],
	function(v, n)
		local fam, ev;

		fam := FamilyObj(v);
		ev := ExtRepOfObj(v);
		
		return ObjByExtRep(fam, ev * n);
	end 
);

InstallMethod( \=, "for a FGAbElement and a FGAbElement", [ IsFGAbElement and IsFGAbElementRep, IsFGAbElement and IsFGAbElementRep ],
	function(v1, v2)
		local fam1, fam2, ev1, ev2;

		fam1 := FamilyObj(v1);
		fam2 := FamilyObj(v2);
		if not IsIdenticalObj(fam1, fam2) then
			return false;
		fi;
		ev1 := ExtRepOfObj(v1);
		ev2 := ExtRepOfObj(v2);
		
		return ev1 = ev2;
	end 
);

InstallMethod( \<, "for a FGAbElement and a FGAbElement", [ IsFGAbElement and IsFGAbElementRep, IsFGAbElement and IsFGAbElementRep ],
	function(v1, v2)
		local fam1, fam2, ev1, ev2;

		fam1 := FamilyObj(v1);
		fam2 := FamilyObj(v2);
		if not IsIdenticalObj(fam1, fam2) then
			Error("\\<: Cannot be compared.");
		fi;
		ev1 := ExtRepOfObj(v1);
		ev2 := ExtRepOfObj(v2);
		
		return ev1 < ev2;
	end 
);

InstallMethod( \in, "for a FGAbElement and a FGAb", [ IsFGAbElement and IsFGAbElementRep, IsFGAb and IsFGAbRep ],
	function(v, A)
		return IsIdenticalObj(FamilyObj(v), FamilyObj(Representative(A)));
	end 
);

InstallMethod( ObjByExtRep, "for a FGAbElement", [ CategoryFamily( IsFGAbElement ), IsVector ],
	function(fam, v)
		local snf, D, f, finv, i;

		v := Vector(Integers, v);
		snf := fam!.snf;
		D := snf.normal;
		f := snf.rowtrans;
		finv := snf.invrowtrans;
		v := f * v;
		for i in [1..DimensionsMat(D)[1]] do
			if i <= DimensionsMat(D)[2] and D[i, i] <> 0 then
				v[i] := v[i] mod D[i, i];
			fi;
		od;
		v := finv * v;

		return Objectify(fam!.type, rec(vector := v));
	end 
);

InstallMethod( ExtRepOfObj, "for a FGAbElement", [ IsFGAbElement and IsFGAbElementRep ],
	function(v)
		return v!.vector;
	end 
);

InstallMethod( PrintObj, "for a FGAbElement", [ IsFGAbElement and IsFGAbElementRep ], 
	function(v)
		Print("< FGAbElement ", ExtRepOfObj(v),">");
	end 
);

InstallMethod( ViewObj, "for a FGAbElement", [ IsFGAbElement and IsFGAbElementRep ], 
	function(v)
		Print("< FGAbElement ", ExtRepOfObj(v),">");
	end 
);
