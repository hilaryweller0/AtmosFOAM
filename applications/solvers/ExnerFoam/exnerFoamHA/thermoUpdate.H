// Update the pressure and temperature
thermo.p() = pRef*pow(Exner, 1/kappa);
thermo.T() = theta*Exner;
thermo.he() == thermo.he(thermo.p(),thermo.T());
thermo.correct();

