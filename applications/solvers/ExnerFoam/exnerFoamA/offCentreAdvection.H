CoNum = CourantNo(phi/rhof, runTime.deltaT());
surfaceScalarField Cof = runTime.deltaT()*mag(phi)*mesh.deltaCoeffs()
                            /(rhof*mesh.magSf());

offCentre = max(ocAlpha, 1 - 1/((max(Cof, maxInterp.interpolate(CoNum))))+0.1);
aphi = offCentre*phi;

Info << "offCentre goes from " << min(offCentre).value()
     << " to " << max(offCentre).value() << endl;

ImEx = 0.5*(sign(Cof - CoLimitExp)+1);
Info << "ImEx goes from " << min(ImEx).value()
         << " to " << max(ImEx).value() << endl;

Ndt = runTime.deltaT()*Foam::sqrt
(
    max
    (
        maxInterp.interpolate
        (
            -(g & fvc::grad(thetaf))/theta
        ),
        max
        (
            -gSfhat*fvc::snGrad(theta)/thetaf,
            dimensionedScalar("", dimensionSet(0,0,-2,0,0), scalar(0))
        )
    )
);

if (SIgravityWaves)
{
    offCentreGW = max(ocAlpha, 1 - 1/max(Ndt, SMALL));
    offCentrePG = max(ocAlpha, offCentreGW);

    Info << "offCentreGW goes from " << min(offCentreGW).value()
         << " to " << max(offCentreGW).value() << endl;
    Info << "offCentrePG goes from " << min(offCentrePG).value()
         << " to " << max(offCentrePG).value() << endl;
}

