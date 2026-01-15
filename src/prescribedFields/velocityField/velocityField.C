#include "velocityField.H"
#include "wallPolyPatch.H"

namespace Foam
{

defineTypeNameAndDebug(velocityField, 0);
defineRunTimeSelectionTable(velocityField, dict);

autoPtr<velocityField> velocityField::New
(
    const dictionary& dict,
    const word velocityType
)
{
    const word velocityFieldType = velocityType == "" ?
            (dict.lookupOrDefault<word>("type", "zeroVelocity")) : velocityType;

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

void velocityField::applyTo(surfaceScalarField& phi, scalar time) const
{
    applyToInternalField(phi, time);
    
    forAll(phi.boundaryField(), patchI)
    {
        if
        (
            !isA<emptyPolyPatch>(phi.mesh().boundaryMesh()[patchI])
         && !isA<wallPolyPatch>(phi.mesh().boundaryMesh()[patchI])
        )
        {
            applyToBoundary(phi, patchI, time);
        }
    }
    
    project(phi);
}

}
