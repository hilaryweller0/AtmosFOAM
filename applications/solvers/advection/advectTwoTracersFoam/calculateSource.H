IOdictionary dict
(
    IOobject
    (
        "advectionDict",
        mesh.time().system(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    )
);

const int sourcePresent = dict.lookupOrDefault<int>("sourcePresent", int(1));
const dimensionedScalar sourceLengthScale(dict.lookup("sourceLengthScale"));

// The vertical velocity used for calculating the source terms
const volScalarField w = U.component(vector::Z);

const volScalarField T1damping
(
    "T1damping",
    sourcePresent*0.5/sourceLengthScale*(mag(w)-w)
);
const volScalarField T2damping
(
    "T2damping",
    sourcePresent*0.5/sourceLengthScale*(mag(w)+w)
);

