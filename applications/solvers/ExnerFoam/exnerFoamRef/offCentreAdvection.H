CoNum = CourantNo(phi/rhof, runTime.deltaT());
surfaceScalarField Cof = maxInterp.interpolate(CoNum);

//offCentre = max(ocAlpha, 1 - 1/max(Cof, SMALL));
offCentre = max(ocAlpha, 1 - 1/(Cof + 0.2));
/*offCentre = max
(
    ocAlpha,
    1 - rhof*mesh.magSf()
        /(runTime.deltaT()*mesh.deltaCoeffs()
            *max(mag(phi), dimensionedScalar(phi.dimensions(), SMALL)))
);*/

aphi = offCentre*phi;

Info << "ocAlpha = " << ocAlpha << " offCentre-ocAlpha goes from "
     << max(offCentre).value()-ocAlpha
     << " to " << min(offCentre).value()-ocAlpha << endl;

ImEx = 0.5*(sign(maxInterp.interpolate(CoNum+SMALL) - CoLimitExp) + 1);

// Modify Imex for strong stratification
if (SIgravityWaves)
{
    volScalarField Ndt = runTime.deltaT()*Foam::sqrt(max
    (
        -(g & fvc::grad(thetaf))/theta,
        dimensionedScalar("", dimensionSet(0,0,-2,0,0), scalar(0))
    ));

    ImEx += (1 - ImEx)*0.5*(sign(maxInterp.interpolate(Ndt) - 1.8)+1);
}

