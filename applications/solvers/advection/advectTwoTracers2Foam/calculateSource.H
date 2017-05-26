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
const dimensionedScalar sourceTemperatureScale(dict.lookup("sourceTemperatureScale"));
const dimensionedScalar velocityScale(dict.lookup("velocityScale"));

// The vertical velocity used for calculating the source terms
const volScalarField w = U.component(vector::Z);

const volScalarField rho1damping
(
    "rho1damping",
    sourcePresent*0.5*T/sourceTemperatureScale * (mag(w)-w)/(mag(w)+velocityScale)
);
const volScalarField rho2damping
(
    "rho2damping",
    sourcePresent*0.5*T/sourceTemperatureScale * (mag(w)+w)/(mag(w)+velocityScale)
);

