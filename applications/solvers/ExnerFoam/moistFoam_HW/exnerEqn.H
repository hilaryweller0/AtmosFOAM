Psi == rho/Exner;
rhof = fvc::interpolate(rho);

// Explicit part of rate of change of flux
gradPcoeff = fvc::interpolate
(
    air.Cp()*theta*atmos.rhoR()
   /(atmos.sumDensity()*air.R()*atmos.volGas())
);
surfaceScalarField gradPcoeff2
(
    "gradPcoeff2", gradPcoeff*G*offCentre*rhof*dt
);
dVdt = -(H.delta() &
(
    fvc::interpolate(fvc::div(U, u))
 + (2*rhof*(Omega^(Uf-Ug)))
));

// Vu is the flux without g or dpdz (to get the BCs right)
Vu = G*
(
    V.oldTime()
  + (1-offCentre)*dt*dVdt.oldTime()
  + offCentre*dt*dVdt
);

// add the gravity terms
//V = Vu + G*offCentre*dt*rhof*gd;
// (orthogonal versin for consistency with partitionedMoistFoam)
V = Vu + G*offCentre*dt*rhof*gSf*H.magd()/mesh.magSf();

// Flux in the normal direction (including off diagonal part of pressure
// gradient)

for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
{
    //U = H.ddirToFlux(V)
    //  - H.ddirToFluxCorr(gradPcoeff2*fvc::snGrad(Exner)*H.magd());
    U = V*mesh.magSf()/H.magd();

    fvScalarMatrix ExnerEqn
    (
        fvm::ddt(Psi, Exner)
      - (1-offCentre)*dRhodt.oldTime()
      + offCentre*fvc::div(U)
      - fvm::laplacian(offCentre*gradPcoeff2, Exner)
    );

    if(ucorr == nOuterCorr-1 && nonOrth == nNonOrthCorr && corr == nCorr-1)
    {
        ExnerEqn.solve((Exner.name() + "Final"));
    }
    else
    {
        ExnerEqn.solve((Exner.name()));
    }

    if (nonOrth == nNonOrthCorr)
    {
        //U += ExnerEqn.flux()/offCentre;
        U -= gradPcoeff2*fvc::snGrad(Exner)*mesh.magSf();
        V -= gradPcoeff2*fvc::snGrad(Exner)*H.magd();
        //u = H.reconstructd(V/rhof);
        u = fvc::reconstruct(U/rhof);
        divu = fvc::div(U/rhof);
        Uf = fvc::interpolate(u);
        dVdt += rhof*gd
             - H.magd()*rhof*gradPcoeff*fvc::snGrad(Exner)
             - muSponge*V;
    }
}

