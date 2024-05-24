CoNum = CourantNo(phi/rhof, runTime.deltaT());
//CoNum = fvc::localMax(maxInterp.interpolate(CoNum));
surfaceScalarField Cof = max
(
    runTime.deltaT()*mag(phi)*mesh.deltaCoeffs()/(rhof*mesh.magSf()),
    //mag(gUnitNormal)*maxInterp.interpolate(CoNum)
    maxInterp.interpolate(CoNum)
);

if (SIgravityWaves)
{
    offCentre = max(ocAlpha, 1 - 1/max(Cof*(1 + Ndt) + offCentreSafety, SMALL));
}
else
{
    offCentre = max(ocAlpha, 1 - 1/max(Cof + offCentreSafety, SMALL));
}

aphi = offCentre*phi;

Info << "offCentre goes from " << min(offCentre).value()
     << " to " << max(offCentre).value() << endl;

//ImEx = 1 - 0.5*(sign(CoLimitExp - Cof)+1)*0.5*(sign(CoLimitExp - Ndt)+1);
ImEx = 0.5*(sign(Cof - CoLimitExp) + 1);
//ImEx = fvc::interpolate(fvc::localMax(ImEx));
Info << "ImEx goes from " << min(ImEx).value()
         << " to " << max(ImEx).value() << endl;

