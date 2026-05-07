#################################################################
# homomorphisms

#! @ChapterInfo `FGAbHomomorphism`, `FGAbHomomorphism`
#! @Arguments m
DeclareCategory( "IsFGAbHomomorphism", IsAdditiveGroupGeneralMapping );

DeclareRepresentation( "IsFGAbHomomorphismRep", IsAttributeStoringRep );

# FGAbHomomorphism( A, B, m )
# It checks if LiftingMatrix exists.
#! @BeginGroup `FGAbHomomorphism`
#! @ChapterInfo `FGAbHomomorphism`, `FGAbHomomorphism`
#! @Arguments A, B, m
#! @Returns a `FGAbHomomorphism`
#! @Description `A` and `B` are `FGAb` as domain and codomain. Let `A` be $A_1\to A_0$, and `B` be $B_1\to B_0$, `m` is a matrix representing a map $A_0\to B_0$. `n` is a matrix representing a map $A_1\to B_1$. `n` is actually determined by `m` up to a chain homotopy. Therefore, we don't need to give `n` explicitly.
DeclareOperation( "FGAbHomomorphism", [ IsFGAb, IsFGAb, IsMatrixOrMatrixObj ] );

# FGAbHomomorphismNC( A, B, m )
#! @ChapterInfo `FGAbHomomorphism`, `FGAbHomomorphism`
#! @Arguments A, B, m
DeclareOperation( "FGAbHomomorphismNC", [ IsFGAb, IsFGAb, IsMatrixOrMatrixObj ] );

# FGAbHomomorphism( A, B, m, n )
# It checks if LiftingMatrix exists.
# n is the matrix A1 -> B1.
#! @ChapterInfo `FGAbHomomorphism`, `FGAbHomomorphism`
#! @Arguments A, B, m, n
#! @Description 
DeclareOperation( "FGAbHomomorphism", [ IsFGAb, IsFGAb, IsMatrixOrMatrixObj, IsMatrixOrMatrixObj ] );

# FGAbHomomorphismNC( A, B, m, n )
#! @ChapterInfo `FGAbHomomorphism`, `FGAbHomomorphism`
#! @Arguments A, B, m, n
DeclareOperation( "FGAbHomomorphismNC", [ IsFGAb, IsFGAb, IsMatrixOrMatrixObj, IsMatrixOrMatrixObj ] );
#! @EndGroup

#! @ChapterInfo `FGAbHomomorphism`, `FGAbHomomorphism`
#! @Arguments A
#! @Returns the underlying matirx $A_0\to B_0$ of a `FGAbHomomorphism`
DeclareAttribute( "UnderlyingMatrix", IsFGAbHomomorphism and IsFGAbHomomorphismRep );

# It returns the lifting f1:A1 -> B1 of f0:A0 -> B0.
#! @ChapterInfo `FGAbHomomorphism`, `FGAbHomomorphism`
#! @Arguments A
#! @Returns the underlying matirx $A_1\to B_1$ of a `FGAbHomomorphism`
DeclareAttribute( "LiftingMatrix", IsFGAbHomomorphism and IsFGAbHomomorphismRep );

# FGAbKernel( f )
#! @BeginGroup `FGAbKernel`
#! @ChapterInfo `FGAbHomomorphism`, limits and colimits
#! @Arguments f
#! @Description `FGAbKernel(f)` is a `FGAb` equipped with `FromFGAbKernel`.
DeclareAttribute( "FGAbKernel", IsFGAbHomomorphism and IsFGAbHomomorphismRep );
#! @ChapterInfo `FGAbHomomorphism`, limits and colimits
#! @Arguments A
#! @Description The input of `FromFGAbKernel` must be a output of `FGAbKernel`.
DeclareAttribute( "FromFGAbKernel", IsFGAb and IsFGAbRep );
#! @ChapterInfo `FGAbHomomorphism`, limits and colimits
#! @Arguments A
DeclareAttribute( "UnderlyingFGAbKernel", IsAdditiveGroup );
#! @EndGroup

# FGAbCokernel( f )
#! @BeginGroup `FGAbCokernel`
#! @ChapterInfo `FGAbHomomorphism`, limits and colimits
#! @Arguments f
#! @Description `FGAbKernel(f)` is a `FGAb` equipped with `ToFGAbCokernel`.
DeclareAttribute( "FGAbCokernel", IsFGAbHomomorphism and IsFGAbHomomorphismRep );
#! @ChapterInfo `FGAbHomomorphism`, limits and colimits
#! @Arguments A
#! @Description The input of `FromFGAbCokernel` must be a output of `FGAbCokernel`.
DeclareAttribute( "ToFGAbCokernel", IsFGAb and IsFGAbRep );
#! @EndGroup

# FGAbPullback( f1, f2 )
#! @BeginGroup `FGAbPullback`
#! @ChapterInfo `FGAbHomomorphism`, limits and colimits
#! @Arguments f1, f2
#! @Description `FGAbPullback(f1, f2)` is a `FGAb` equipped with `FromFGAbPullback`.
DeclareOperation( "FGAbPullback", [ IsFGAbHomomorphism and IsFGAbHomomorphismRep, IsFGAbHomomorphism and IsFGAbHomomorphismRep ] );

DeclareAttribute( "FGAbPullbackInfo", IsFGAb and IsFGAbRep );
#! @ChapterInfo `FGAbHomomorphism`, limits and colimits
#! @Arguments A, n
#! @Description `FromFGAbPullback(A, n)` is the projection to `Source(f1)` or `Source(f2)`.
DeclareOperation( "FromFGAbPullback", [ IsFGAb and IsFGAbRep, IsInt ] );
#! @EndGroup

# FGAbPushout( f1, f2 )
#! @BeginGroup `FGAbPushout`
#! @ChapterInfo `FGAbHomomorphism`, limits and colimits
#! @Arguments f1, f2
#! @Description `FGAbPushout(f1, f2)` is a `FGAb` equipped with `ToFGAbPushout`.
DeclareOperation( "FGAbPushout", [ IsFGAbHomomorphism and IsFGAbHomomorphismRep, IsFGAbHomomorphism and IsFGAbHomomorphismRep ] );

DeclareAttribute( "FGAbPushoutInfo", IsFGAb and IsFGAbRep );
#! @ChapterInfo `FGAbHomomorphism`, limits and colimits
#! @Arguments A, n
#! @Description `ToFGAbPushout(A, n)` is the projection to `Range(f1)` or `Range(f2)`.
DeclareOperation( "ToFGAbPushout", [ IsFGAb and IsFGAbRep, IsInt ] );
#! @EndGroup