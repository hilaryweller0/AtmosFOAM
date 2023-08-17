// Explicit solution of the Momentum equation

// Exner from eqn of state
Exnerp = pow(rho*R*theta/pRef, RbyCv) - Exnera;

drhoUdt = -fvc::div(phi, U)
          //- rho*g*thetap/thetaa
          //- rho*Cp*theta*fvc::grad(Exnerp);
        - rho*fvc::weightedReconstruct
        (
            gSf*thetapf/thetaaf + Cp*thetaf*fvc::snGrad(Exnerp)*mesh.magSf(),
            0.25
        );


dPhidt = -fvc::flux(fvc::div(phi, U))
       - rhof*
       (
           gSf*thetapf/thetaaf + Cp*thetaf*fvc::snGrad(Exnerp)*mesh.magSf()
       );

phi = fvc::flux(rho.oldTime()*U.oldTime()) + runTime.deltaT()*
    (
        ocAlpha*dPhidt + (1-ocAlpha)*dPhidt.oldTime()
    );

fvVectorMatrix UEqn
(
    EulerDdtv.fvmDdt(rho, U)
  - ocAlpha*drhoUdt
  - (1-ocAlpha)*drhoUdt.oldTime()
);

UEqn.solve();

//phi = fvc::flux(rho*U);

