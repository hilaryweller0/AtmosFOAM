#include "tracerField.H"

defineRunTimeSelectionTable(tracerField, dict);

autoPtr<tracerField> tracerField::New
(
    const dictionary& dict,
    const advectable& velocityField
)
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
       cstrIter()(dict, velocityField)
    );
}

tracerField::tracerField(const advectable& velocityField)
:
velocityField(velocityField)
{}

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

void tracerField::applyTo(surfaceScalarField& Tf) const
{
    const fvMesh& mesh = Tf.mesh();
    forAll(Tf, faceI)
    {
        const point& p = mesh.Cf()[faceI];
        Tf[faceI] = tracerAt(velocityField.initialPositionOf(p, Tf.time()), Tf.time());
    }
}

void tracerField::applyToInternalField(volScalarField& T) const
{
    const fvMesh& mesh = T.mesh();
    forAll(T, cellI)
    {
        const point& p = mesh.C()[cellI];
        T[cellI] = tracerAt(velocityField.initialPositionOf(p, T.time()), T.time());
    }
}
