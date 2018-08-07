// Read in prognostic variables
Info<< "Reading field Exner\n" << endl;
volScalarField Exner
(
    IOobject
    (
        "Exner",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
Exner.oldTime();

// The partition names
const label nParts = 2;
wordList partNames(nParts);
partNames[0] = "stable.";
partNames[1] = "buoyant.";

// Set pointer lists for variables in all partitions. Read in later
PtrList<volScalarField> sigma(nParts);
PtrList<volScalarField> theta(nParts);
PtrList<surfaceVectorField> Uf(nParts);

PtrList<volScalarField> rho(nParts);
PtrList<volScalarField> sigmaRho(nParts);
PtrList<volVectorField> u(nParts);
PtrList<surfaceScalarField> flux(nParts);
//PtrList<surfaceScalarField> gradPcoeff(nParts);
PtrList<surfaceScalarField> dFluxdt(nParts);
PtrList<surfaceScalarField> drag(nParts);
List<PtrList<volScalarField> > massTransfer(nParts);

// Read in or set variables in partitions
for (label ip = 0; ip < nParts; ip++)
{
    sigma.set
    (
        ip,
        new volScalarField
        (
            IOobject
            (
                partNames[ip]+"sigma",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );
    sigma[ip].oldTime();

    theta.set
    (
        ip,
        new volScalarField
        (
            IOobject
            (
                partNames[ip]+"theta",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );
    theta[ip].oldTime();
    
    Uf.set
    (
        ip,
        new surfaceVectorField
        (
            IOobject
            (
                partNames[ip]+"Uf",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );
    
    rho.set
    (
        ip,
        new volScalarField
        (
            IOobject
            (
                partNames[ip]+"rho",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            pRef/(R*theta[ip])*pow(Exner, (1-kappa)/kappa)
        )
    );
    rho[ip].write();
    
    sigmaRho.set
    (
        ip,
        new volScalarField
        (
            IOobject
            (
                partNames[ip]+"sigmaRho",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            sigma[ip]*rho[ip],
            "zeroGradient"
        )
    );
    sigmaRho[ip].oldTime();
    sigmaRho[ip].write();
    
    u.set
    (
        ip,
        new volVectorField
        (
            IOobject
            (
                partNames[ip]+"u",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            H.reconstruct(Uf[ip] & mesh.Sf()),
            Uf[ip].boundaryField().types()
        )
    );
    
    flux.set
    (
        ip,
        new surfaceScalarField
        (
            IOobject
            (
                partNames[ip]+"flux",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            linearInterpolate(sigmaRho[ip])*(Uf[ip] & mesh.Sf())
        )
    );
    flux[ip].oldTime();
    
//    gradPcoeff.set
//    (
//        ip,
//        new surfaceScalarField
//        (
//            IOobject
//            (
//                partNames[ip]+"gradPcoeff",
//                runTime.timeName(),
//                mesh
//            ),
//            Cp*linearInterpolate(sigmaRho[ip])*linearInterpolate(theta[ip])
//        )
//    );
    
    dFluxdt.set
    (
        ip,
        new surfaceScalarField
        (
            IOobject
            (
                partNames[ip]+"dFluxdt",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar("dFluxdt", dimensionSet(1,0,-2,0,0), scalar(0)),
            "fixedValue"
        )
    );
    dFluxdt[ip].oldTime();
    
    drag.set
    (
        ip,
        new surfaceScalarField
        (
            IOobject
            (
                partNames[ip]+"drag",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar("drag", dimensionSet(1,0,-2,0,0), scalar(0))
        )
    );
    
    massTransfer[ip].setSize(nParts);

    for(label jp = 0; jp < nParts; jp++)
    {
        if (ip != jp)
        {
            massTransfer[ip].set
            (
                jp,
                new volScalarField
                (
                    IOobject
                    (
                        partNames[ip]+partNames[jp]+"massTransfer",
                        runTime.timeName(),
                        mesh
                    ),
                    mesh,
                    dimensionedScalar("", dimensionSet(1,-3,-1,0,0), scalar(0))
                )
            );
        }
    }
}

// Summed variables
volScalarField rhoSum
(
    IOobject
    (
        "rhoSum",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    sigmaRho[0]
);
for(label ip = 1; ip < nParts; ip++){ rhoSum += sigmaRho[ip]; }
rhoSum.write();

volScalarField rhoTheta("rhoSum", sigmaRho[0]*theta[0]);
for(label ip = 1; ip < nParts; ip++)
{
    rhoTheta += sigmaRho[ip]*theta[ip];
};

volScalarField Psi("Psi", rhoSum/Exner);
Psi.oldTime();

surfaceScalarField fluxSum
(
    IOobject("fluxSum", runTime.timeName(), mesh),
    flux[0]
);
for(label ip = 1; ip < nParts; ip++){ fluxSum += flux[ip]; }

surfaceScalarField gradPcoeffSum
(
    IOobject("gradPcoeff", runTime.timeName(), mesh),
    Cp*fvc::interpolate(rhoTheta)
);
//for(label ip = 1; ip < nParts; ip++)
//{
//    gradPcoeffSum += gradPcoeff[ip];
//}

#include "calculateDrag.H"
