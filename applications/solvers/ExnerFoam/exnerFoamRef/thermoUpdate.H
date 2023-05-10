// Update the pressure and temperature
thermo.p() = pRef*pow(Exnera + Exnerp, 1/kappa);
thermo.T() = (thetaa + thetap)*(Exnera + Exnerp);
thermo.he() == thermo.he(thermo.p(),thermo.T());
thermo.correct();

