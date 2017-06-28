// update sigma according to the continuity equation
sumSigma == 0;
for(label ip = 0; ip < nParts; ip++)
{
    sigma[ip] = (h.oldTime()*sigma[ip].oldTime() - dt*fvc::div(flux[ip]))/h;
    sumSigma += sigma[ip];
}
Info << "1-sumSigma goes from " << min(1-sumSigma).value() << " to "
     << max(1-sumSigma).value() << endl;
Info << "sigma[0] goes from " << min(sigma[0]).value() << " to "
     << max(sigma[0]).value() << endl;
Info << "sigma[1] goes from " << min(sigma[1]).value() << " to "
     << max(sigma[1]).value() << endl;

h.write();
u[0].write();
u[1].write();
sigma[0].write();
sigma[1].write();
FatalErrorIn("") << exit(FatalError);

