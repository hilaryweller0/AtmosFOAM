CoNum = CourantNo(phi/rhof, runTime.deltaT());
surfaceScalarField Cof = maxInterp.interpolate(CoNum);

offCentre = max(ocAlpha, 1 - 1/max(Cof, SMALL));
aphi = offCentre*phi;

Info << "ocAlpha = " << ocAlpha << " offCentre-ocAlpha goes from "
     << max(offCentre).value()-ocAlpha
     << " to " << min(offCentre).value()-ocAlpha << endl;

ImEx = 0.5*(sign(maxInterp.interpolate(CoNum) - CoLimitExp) + 1);
