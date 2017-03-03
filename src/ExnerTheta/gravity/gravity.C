#include "gravity.H"

Foam::gravity::gravity(const fvMesh& mesh)
:
    mesh_(mesh),
    dict_(
        IOobject
        (
            "environmentalProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    g_(dict_.lookup("g")),
    unitFaceNormal_(mesh.Sf() / mesh.magSf() & (g_ / mag(g_)))
{}

Foam::dimensionedVector Foam::gravity::operator()() const
{
    return g_;
}

Foam::gravity::~gravity() {}
