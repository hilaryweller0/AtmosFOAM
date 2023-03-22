#include "velocityField.H"
#include "wallPolyPatch.H"

defineRunTimeSelectionTable(velocityField, dict);

autoPtr<velocityField> velocityField::New(const dictionary& dict)
{
    const word velocityFieldType(dict.lookupOrDefault<word>("type", "zeroVelocity"));

    dictConstructorTable::iterator cstrIter =
        dictConstructorTablePtr_->find(velocityFieldType);

    if (cstrIter == dictConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "velocityField::New(const dictionary&)"
        ) << "Unknown type "
          << velocityFieldType << nl
          << "Valid types: " << endl
          << dictConstructorTablePtr_->sortedToc()
          << exit(FatalError);
    }

    return autoPtr<velocityField>
    (
       cstrIter()(dict)
    );
}

void velocityField::applyTo(surfaceScalarField& phi) const
{
    applyToInternalField(phi);
    
    forAll(phi.boundaryField(), patchI)
    {
        if
        (
            !isA<emptyPolyPatch>(phi.mesh().boundaryMesh()[patchI])
         && !isA<wallPolyPatch>(phi.mesh().boundaryMesh()[patchI])
        )
        {
            applyToBoundary(phi, patchI);
        }
    }
}

