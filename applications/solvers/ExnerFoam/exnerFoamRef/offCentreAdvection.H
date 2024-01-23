CoNum = CourantNo(phi/rhof, runTime.deltaT());
//CoNum = fvc::localMax(maxInterp.interpolate(CoNum));
surfaceScalarField Cof = max
(
    runTime.deltaT()*mag(phi)*mesh.deltaCoeffs()/(rhof*mesh.magSf()),
    maxInterp.interpolate(CoNum)
);

if (SIgravityWaves)
{
    offCentre = max(ocAlpha, 1 - 1/max(Cof + Ndt + 0.2, SMALL));
}
else
{
    offCentre = max(ocAlpha, 1 - 1/max(Cof + 0.2, SMALL));
}

aphi = offCentre*phi;

Info << "offCentre goes from " << min(offCentre).value()
     << " to " << max(offCentre).value() << endl;

ImEx = 1 - 0.5*(sign(CoLimitExp - Cof)+1)*0.5*(sign(CoLimitExp - Ndt)+1);
Info << "ImEx goes from " << min(ImEx).value()
         << " to " << max(ImEx).value() << endl;

