offCentre = max
(
    ocAlpha,
    1 - rhof*mesh.magSf()
        /(runTime.deltaT()*mesh.deltaCoeffs()
            *max(mag(phi), dimensionedScalar(phi.dimensions(), SMALL)))
);
aphi = offCentre*phi;

Info << "ocAlpha = " << ocAlpha << " offCentre-ocAlpha goes from "
     << max(offCentre).value()-ocAlpha
     << " to " << min(offCentre).value()-ocAlpha << endl;
