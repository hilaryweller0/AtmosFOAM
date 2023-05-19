// Update the pressure and temperature
thermo.p() = pRef*pow(Exnerg +ExnerOffset, 1/kappa);
thermo.T() = theta*(Exnerg + ExnerOffset);
thermo.he() == thermo.he(thermo.p(),thermo.T());
thermo.correct();

