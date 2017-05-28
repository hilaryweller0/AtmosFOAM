// Update sigma for all partitions

// Update total density rate of change
atmosParts.dRhodt() = -fvc::div(atmosParts.flux());

// Check sumSigma
volScalarField sumSigma
(
    IOobject("sumSigma", Exner.mesh().time().timeName(), Exner.mesh()),
    Exner.mesh(), dimensionedScalar("sum", dimless, scalar(0))
);

// Loop over all partitions
for(label ip = 0; ip < atmosParts.size(); ip++)
{
    partition& parti = atmosParts[ip];
    
    // Advect density for this partition (without partition transfer terms)
    parti.dSigmaRhodt() = -fvc::div(parti.flux());

    volScalarField rho = parti.sigmaRho()/parti.sigma();

    parti.sigma() = (parti.sigmaRho().oldTime() + dt*
    (
        (1-offCentre)*parti.dSigmaRhodt().oldTime()
      + offCentre*parti.dSigmaRhodt()
    ))/rho;

    sumSigma += parti.sigma();
}

Info << "1-sumSigma goes from " << 1-max(sumSigma).value() << " to "
     << 1-min(sumSigma).value() << endl;
     
// Scale so that sigmas sum to 1
for(label ip = 0; ip < atmosParts.size(); ip++)
{
    atmosParts[ip].sigma() /= sumSigma;
}

