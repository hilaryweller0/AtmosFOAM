Info << "Butcher Tableau\n" << endl;
const butcherTableau Bt
(
    mesh.schemes().subDict("ddtSchemes").lookup
    (
        "butcherTableau"
    )
);

// Fields to be stored at each Runge-Kutta sub-step
PtrList<volScalarField> div0(T.size());
PtrList<PtrList<volScalarField> > div(T.size());
for(label iT = 0; iT < T.size(); iT++)
{
    div0.set(iT, divSum[iT]);

    div.set(iT, new PtrList<volScalarField>(Bt.nSteps()-1));

    for (int iRK = 0; iRK < Bt.nSteps()-1; iRK++)
    {
        div[iT].set(iRK, divSum[iT]);
    }
}

