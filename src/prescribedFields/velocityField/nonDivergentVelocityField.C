#include "nonDivergentVelocityField.H"
#include "polarPoint.H"

void nonDivergentVelocityField::applyToInternalField(surfaceScalarField& phi) const
{
    phi.ref() = dimensionedScalar("phi", phi.dimensions(), scalar(0));
    const fvMesh& mesh = phi.mesh();
    forAll(phi, faceI)
    {
        const face& f = mesh.faces()[faceI];
        phi[faceI] = faceFlux(f, mesh, phi.time());
        Info << "V/dt " << mesh.V()[mesh.faceOwner()[faceI]]/mesh.time().deltaTValue();
    }
}

void nonDivergentVelocityField::applyToBoundary(surfaceScalarField& phi, const label patchI) const
{
    const fvMesh& mesh = phi.mesh();
    scalarField& bf = phi.boundaryFieldRef()[patchI];
    forAll(bf, faceI)
    {
        const face& f = mesh.boundaryMesh()[patchI][faceI];
        bf[faceI] = faceFlux(f, mesh, phi.time());
    }
}

scalar nonDivergentVelocityField::faceFlux(const face& f, const fvMesh& mesh, const Time& t) const
{
    point p0 = mesh.points()[f.last()];
    point p1 = p0;

    scalar flux = 0;

    scalar minR = -1;
    scalar maxR = -1;
    forAll(f, ip)
    {
        p0 = p1;
        p1 = mesh.points()[f[ip]];
        point pmid = 0.5*(p0 + p1);
        //if (abs(pmid/mag(pmid) & (p0 - p1)) > 1e-8)
        //{
        
        if (minR < 0 || mag(p0) < minR) minR = mag(p0);
        if (mag(p0) > maxR) maxR = mag(p0);

        flux += streamfunctionAt(pmid, t) & ((p0 - p1)/mag(p0-p1));

        polarPoint p0p = convertToPolar(p0);
        polarPoint p1p = convertToPolar(p1);
        polarPoint pmidp = convertToPolar(pmid);
        Info << " mag(edge) " << mag(p0 - p1) << " flux " << (streamfunctionAt(pmid, t) & (p0 - p1));
        Info << " (lat,lon,r) (" << p0p.lat() << "," << p0p.lon() << "," << p0p.r() << ")--";
        Info << "(" << p1p.lat() << "," << p1p.lon() << "," << p1p.r() << ") mid ";
        Info << "(" << pmidp.lat() << "," << pmidp.lon() << "," << pmidp.r() << ")" << endl;
        
        //}
    }

    flux *= (sqr(maxR) - sqr(minR))/2;

    Info << "flux for " << f << " " << flux << " " << endl;

    return flux;
}
