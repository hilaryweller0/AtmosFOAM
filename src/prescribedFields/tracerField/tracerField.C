#include "tracerField.H"
#include "fixedValueFvPatchField.H"
#include "fixedGradientFvPatchFields.H"
#include "sphericalMeshData.H"

namespace Foam
{
defineTypeNameAndDebug(tracerField, 0);
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

tracerField::tracerField
(
    const advectable& velocityField,
    const bool spherical__,
    const scalar earthRadius__
)
:
    velocityField(velocityField),
    spherical_(spherical__),
    earthRadius_(earthRadius__)
{
    if(spherical() && earthRadius() < SMALL)
    {
        WarningIn("tracerField::tracerField")
            << " spherical geometry specifed but earthRadius = "
            << earthRadius() << endl;
    }
}

void tracerField::applyTo(volScalarField& T) const
{
    applyToInternalField(T);
    
    forAll(T.boundaryField(), patchI)
    {
        applyToBoundary(T, patchI);
    }
}

void tracerField::applyTo(surfaceScalarField& Tf) const
{
    const fvMesh& mesh = Tf.mesh();

    // Face centres to apply the field at
    const vector* Cp { &mesh.Cf()[0] };

    if (spherical_)
    {
        // Get reference to spherical geometry if needed
        const sphericalMeshData& sphericalData = sphericalMeshData::New
        (
            mesh, earthRadius_
        );
        Cp = &sphericalData.faceCentres()[0];
    }

    forAll(Tf, faceI)
    {
        const point& p = Cp[faceI];
        Tf[faceI] = tracerAt
        (
            velocityField.initialPositionOf(p, Tf.time()),
            Tf.time()
        );
    }
    forAll(Tf.boundaryField(), patchI)
    {
        applyToBoundary(Tf, patchI);
    }
}

void tracerField::applyToInternalField(volScalarField& T) const
{
    const fvMesh& mesh = T.mesh();
    
    // Cell centres to apply the field at
    const vector* Cp { &mesh.C()[0] };

    if (spherical_)
    {
        // Get reference to spherical geometry if needed
        const sphericalMeshData& sphericalData = sphericalMeshData::New
        (
            mesh, earthRadius_
        );
        Cp = &sphericalData.cellCentres()[0];
    }

    forAll(T, cellI)
    {
        const point& p = Cp[cellI];
        T[cellI] = tracerAt
        (
            velocityField.initialPositionOf(p, T.time()),
            T.time()
        );
    }
}

void tracerField::applyToBoundary(volScalarField& T, const label patchI) const
{
    fvPatchField<scalar>& patch = T.boundaryFieldRef()[patchI];
    const fvPatch& pat = T.mesh().boundary()[patchI];

    // Patch centres to apply the field at
    const vector* Cp { &(T.mesh().Cf().boundaryField()[patchI][0]) };

    if (spherical_)
    {
        // Get reference to spherical geometry if needed
        const sphericalMeshData& sphericalData = sphericalMeshData::New
        (
            T.mesh(), earthRadius_
        );
        Cp = &(sphericalData.faceCentres()[pat.start()]);
    }

    if (patch.type() == "fixedValue")
    {
        forAll(patch, faceI)
        {
            const point& p = Cp[faceI];
            patch[faceI] = tracerAt
            (
                p,
                T.time()
            );
        }
    }
    else if (patch.type() == "fixedGradient" && hasGradient())
    {
        fixedGradientFvPatchField<scalar>& gradPatch
             = dynamic_cast<fixedGradientFvPatchField<scalar>& >(patch);
        forAll(gradPatch.gradient(), faceI)
        {
            const vectorField nf = gradPatch.patch().nf();
            gradPatch.gradient()[faceI] = nf[faceI] & gradientAt
            (
                gradPatch.patch().Cf()[faceI],
                T.time()
            );
        }
    }
}

void tracerField::applyToBoundary
(
    surfaceScalarField& T,
    const label patchI
) const
{
    fvsPatchField<scalar>& patch = T.boundaryFieldRef()[patchI];
    const fvPatch& pat = T.mesh().boundary()[patchI];

    // Patch centres to apply the field at
    const vector* Cp { &(T.mesh().Cf().boundaryField()[patchI][0]) };

    if (spherical_)
    {
        // Get reference to spherical geometry if needed
        const sphericalMeshData& sphericalData = sphericalMeshData::New
        (
            T.mesh(), earthRadius_
        );
        Cp = &(sphericalData.faceCentres()[pat.start()]);
    }

    forAll(patch, faceI)
    {
        const point& p = Cp[faceI];
        patch[faceI] = tracerAt
        (
            p,
            T.time()
        );
    }
}

bool tracerField::hasGradient() const
{
    return false;
}

vector tracerField::gradientAt
(
    const point& p,
    const Time& t
) const
{
    notImplemented
    (
        "tracerField::gradientAt(const point& p, const Time& t)"
    );
    return vector(0,0,0);
}
}
