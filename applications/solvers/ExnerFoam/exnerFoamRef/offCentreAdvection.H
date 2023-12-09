CoNum = CourantNo(phi/rhof, runTime.deltaT());
surfaceScalarField Cof = runTime.deltaT()*mag(phi)*mesh.deltaCoeffs()
                            /(rhof*mesh.magSf());
                        //maxInterp.interpolate(CoNum);

offCentre = max(0.5, 1 - 1/(Cof + 0.2));

aphi = offCentre*phi;

Info << "offCentre goes from " << min(offCentre).value()
     << " to " << max(offCentre).value() << endl;

ImEx = //0.5*(sign(maxInterp.interpolate(CoNum+SMALL) - CoLimitExp) + 1);
       0.5*(sign(Cof - CoLimitExp)+1);

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
    offCentrePG = max(ocAlpha, 1 - 1/max(Ndt, SMALL));
    //offCentreGW = max(scalar(0), 1 - 1/max(Ndt, SMALL));
    Info << "offCentreGW goes from " << min(offCentreGW).value()
         << " to " << max(offCentreGW).value() << endl;

    // Modify Imex for strong stratification
    //ImEx += (1 - ImEx)*0.5*(sign(Ndt - 1.8) + 1);
}

