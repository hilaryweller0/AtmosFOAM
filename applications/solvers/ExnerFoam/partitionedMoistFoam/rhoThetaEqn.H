// flux form rho theta equation for each partition

for(label ip = 0; ip < atmosParts.size(); ip++)
{
    partition& parti = atmosParts[ip];
    const volScalarField& sigmaRhoi = parti.sigmaRho();
    volScalarField& thetai = parti.theta();
    const surfaceScalarField& fluxi = parti.flux();

    // Advect theta
    parti.dSigmaRhoThetadt() = -fvc::div(fluxi, thetai);
    thetai = 
    (
        sigmaRhoi.oldTime()*thetai.oldTime() + dt*
        (
            (1-offCentre)*parti.dSigmaRhoThetadt().oldTime()
          + offCentre*parti.dSigmaRhoThetadt()
        )
    )/sigmaRhoi;

    // Calculate and apply source term
    parti.T() = thetai*Exner;
    volScalarField heatSource = sigmaRhoi*thetai*parti.thetaSource();

    parti.dSigmaRhoThetadt() += heatSource;
    thetai += dt*offCentre*heatSource/sigmaRhoi;

    // Update temperature
    parti.T() = thetai*Exner;
}
