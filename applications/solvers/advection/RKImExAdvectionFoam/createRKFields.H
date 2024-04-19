Info << "Butcher Tableau\n" << endl;
const butcherTableau Bt
(
    mesh.schemes().subDict("ddtSchemes").lookup
    (
        "butcherTableau"
    )
);

/*// Fields to be stored at each Runge-Kutta sup-step
RKfield<scalar, fvPatchField, volMesh> dTdt
(
    Bt.nSteps(),
    0*T/dt
);*/

// Fields to be stored at each Runge-Kutta sup-step
RKfield<scalar, fvPatchField, volMesh> Ts
(
    Bt.nSteps(),
    0*T
);

RKfield<scalar, fvsPatchField, surfaceMesh> phis
(
    Bt.nSteps(),
    phi
);

volScalarField divSum("divSum", 0*T/dt);

