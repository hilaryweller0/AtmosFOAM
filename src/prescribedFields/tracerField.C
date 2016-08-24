#include "tracerField.H"

defineRunTimeSelectionTable(tracerField, dict);

autoPtr<tracerField> tracerField::New(const dictionary& dict)
{
    const word tracerFieldType(dict.lookup("type"));

    dictConstructorTable::iterator cstrIter =
        dictConstructorTablePtr_->find(tracerFieldType);

    if (cstrIter == dictConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "tracerField::New(const dictionary&)"
        ) << "Unknown type "
          << tracerFieldType << nl
          << "Valid types: " << endl
          << dictConstructorTablePtr_->sortedToc()
          << exit(FatalError);
    }

    return autoPtr<tracerField>
    (
       cstrIter()(dict)
    );
}

void tracerField::applyTo(volScalarField& T) const
{
    applyToInternalField(T);
    
    forAll(T.boundaryField(), patchI)
    {
        if (!isA<emptyPolyPatch>(T.mesh().boundaryMesh()[patchI]))
        {
//            applyToBoundary(T, patchI);
        }
    }
}

void tracerField::applyToInternalField(volScalarField& T) const
{
    const fvMesh& mesh = T.mesh();
    forAll(T, cellI)
    {
        const point& p = mesh.C()[cellI];
        T[cellI] = tracerAt(p, T.time());
    }
}
