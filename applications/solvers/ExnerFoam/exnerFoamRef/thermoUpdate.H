// Update the pressure and temperature
//Exnerp = pow(rho*R*theta/pRef, RbyCv) - Exnera;
thermo.p() == pRef*pow(Exnera + Exnerp, 1/kappa);
thermo.T() == theta*(Exnera + Exnerp);
//thermo.p() = rho*R*thermo.T();
thermo.he() == thermo.he(thermo.p(),thermo.T());
thermo.correct();

