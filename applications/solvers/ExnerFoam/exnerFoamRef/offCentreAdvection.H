CoNum = CourantNo(phi/rhof, runTime.deltaT());
surfaceScalarField Cof = maxInterp.interpolate(CoNum);

offCentre = max(ocAlpha, 1 - 1/(Cof + 0.2));

aphi = offCentre*phi;

Info << "ocAlpha = " << ocAlpha << " offCentre-ocAlpha goes from "
     << max(offCentre).value()-ocAlpha
     << " to " << min(offCentre).value()-ocAlpha << endl;

ImEx = 0.5*(sign(maxInterp.interpolate(CoNum+SMALL) - CoLimitExp) + 1);

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

// Modify Imex for strong stratification
if (SIgravityWaves)
{
    //offCentrePG = max(ocAlpha, 1 - 1/max(2*Ndt, SMALL));
    offCentreGW = max(scalar(0), 1 - 1/max(Ndt, SMALL));
    Info << "offCentreGW goes from " << min(offCentreGW).value()
         << " to " << max(offCentreGW).value() << endl;

    ImEx += (1 - ImEx)*0.5*(sign(Ndt - 1.8) + 1);

/*    Info << "ocAlpha = " << ocAlpha << " offCentrePG-ocAlpha goes from "
         << max(offCentrePG).value()-ocAlpha
         << " to " << min(offCentrePG).value()-ocAlpha << endl;
*/
}

