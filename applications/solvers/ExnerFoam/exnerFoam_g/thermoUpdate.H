// Update the pressure and temperature
thermo.p() = pRef*pow(Exnerg + gh/(Cp*theta), 1/kappa);
thermo.T() = theta*(Exnerg + gh/(Cp*theta));
thermo.he() == thermo.he(thermo.p(),thermo.T());
thermo.correct();

