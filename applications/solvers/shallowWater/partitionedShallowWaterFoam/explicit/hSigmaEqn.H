// update h and sigma according to the continuity equation
h = h.oldTime() - dt*fvc::div(fluxSum);
hf = fvc::interpolate(h);

sumSigma == 0;
sumSigmaf == 0;
for(label ip = 0; ip < nParts; ip++)
{
    sigma[ip] = 
    (
        h.oldTime()*sigma[ip].oldTime()
      - dt*fvc::div(flux[ip], sigma[ip])
    )/h;
    sigmaf[ip] = fvc::interpolate(sigma[ip]);
    sumSigma += sigma[ip];
    sumSigmaf += sigmaf[ip];
}

Info << "1-sumSigma goes from " << min(1-sumSigma).value() << " to "
     << max(1-sumSigma).value() << endl;
Info << "1-sumSigmaf goes from " << min(1-sumSigmaf).value() << " to "
     << max(1-sumSigmaf).value() << endl;
Info << "sigma[0] goes from " << min(sigma[0]).value() << " to "
     << max(sigma[0]).value() << endl;
Info << "sigma[1] goes from " << min(sigma[1]).value() << " to "
     << max(sigma[1]).value() << endl;
Info << "sigmaf[0] goes from " << min(sigmaf[0]).value() << " to "
     << max(sigmaf[0]).value() << endl;

