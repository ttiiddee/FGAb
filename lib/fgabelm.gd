#################################################################
# elements

# IsExtLElement and IsExtRElement are elements that have a \* operation by external elements.
#! @BeginGroup `IsFGAbElement`
#! @ChapterInfo `FGAb`, `FGAbElement`
#! @Arguments v
DeclareCategory( "IsFGAbElement", IsAdditiveElementWithInverse and IsAdditiveElementWithZero and IsExtLElement and IsExtRElement );

#! @ChapterInfo `FGAb`, `FGAbElement`
#! @Arguments x
DeclareCategoryCollections( "IsFGAbElement" );

#! @ChapterInfo `FGAb`, `FGAbElement`
#! @Arguments x
DeclareCategoryCollections( "IsFGAbElementCollection" );
#! @EndGroup

DeclareRepresentation( "IsFGAbElementRep", IsComponentObjectRep );
