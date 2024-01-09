CoNum = CourantNo(phi/rhof, runTime.deltaT());
surfaceScalarField Cof = max
(
    runTime.deltaT()*mag(phi)*mesh.deltaCoeffs()/(rhof*mesh.magSf()),
    maxInterp.interpolate(CoNum)
);

offCentre = max(0.5, 1 - 1/max(Cof, SMALL));
aphi = offCentre*phi;

Info << "offCentre goes from " << min(offCentre).value()
     << " to " << max(offCentre).value() << endl;

ImEx = 1 - 0.5*(sign(CoLimitExp - Cof)+1)*0.5*(sign(CoLimitExp - Ndtf)+1);
Info << "ImEx goes from " << min(ImEx).value()
         << " to " << max(ImEx).value() << endl;

if (SIgravityWaves)
{
    offCentrePG = max(scalar(0), 1 - 1/(2*Ndt));
    offCentrePG = max(ocAlpha, max(offCentre, offCentrePG));

//    Info << "offCentreGW goes from " << min(offCentreGW).value()
//         << " to " << max(offCentreGW).value() << endl;

    Info << "offCentrePG goes from " << min(offCentrePG).value()
         << " to " << max(offCentrePG).value() << endl;
}

