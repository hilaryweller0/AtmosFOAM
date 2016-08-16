#include "velocityField.H"

defineRunTimeSelectionTable(velocityField, dict);

autoPtr<velocityField> velocityField::New(const dictionary& dict)
{
    const word velocityFieldType(dict.lookup("type"));

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
        if (!isA<emptyPolyPatch>(phi.mesh().boundaryMesh()[patchI]))
        {
            applyToBoundary(phi, patchI);
        }
    }
}

void velocityField::applyToInternalField(surfaceScalarField& phi) const
{
    phi.ref() = dimensionedScalar("phi", phi.dimensions(), scalar(0));
    const fvMesh& mesh = phi.mesh();
    forAll(phi, faceI)
    {
        const face& f = mesh.faces()[faceI];
        phi[faceI] = faceFlux(f, mesh, phi.time());
    }
}

void velocityField::applyToBoundary(surfaceScalarField& phi, const label patchI) const
{
    const fvMesh& mesh = phi.mesh();
    scalarField& bf = phi.boundaryFieldRef()[patchI];
    forAll(bf, faceI)
    {
        const face& f = mesh.boundaryMesh()[patchI][faceI];
        phi[faceI] = faceFlux(f, mesh, phi.time());
    }
}

scalar velocityField::faceFlux(const face& f, const fvMesh& mesh, const Time& t) const
{
    point p0 = mesh.points()[f.last()];
    point p1 = p0;

    scalar flux = 0;

    forAll(f, ip)
    {
        p0 = p1;
        p1 = mesh.points()[f[ip]];
        point pmid = 0.5*(p0 + p1);
        flux += streamfunctionAt(pmid, t) & (p0 - p1);
    }

    return flux;
}
