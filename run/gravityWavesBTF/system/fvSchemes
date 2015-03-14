/*---------------------------------------------------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  1.4                                   |
|   \\  /    A nd           | Web:      http://www.openfoam.org               |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/

FoamFile
{
    version         2.0;
    format          ascii;

    root            "";
    case            "";
    instance        "";
    local           "";

    class           dictionary;
    object          fvSchemes;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

ddtSchemes
{
    default Euler;
}
offCentre 0.5;
SIgravityWaves true;

gradSchemes
{
    default         none;
    grad(theta)     Gauss midPoint;
}

divSchemes
{
    default             none;
    div(U,theta)    Gauss cubicUpwindCPCFit 3;
    div(U,u)        Gauss cubicUpwindCPCFit 3;
}

laplacianSchemes
{
    default         Gauss linear uncorrected;
}

interpolationSchemes
{
    default        none; 
    interpolate(theta)    midPoint;
    interpolate(rho)    midPoint;
    interpolate(convection(U,u)) midPoint;
    interpolate(u) midPoint;
    interpolate(grad(theta)) midPoint;
    interpolate(Exner)       midPoint;
    H                        midPoint;
}

snGradSchemes
{
    default         none;
    snGrad(Exner)   uncorrected;
    snGrad(theta)   uncorrected;
}

fluxRequired
{
    default         no;
    Exner;
}


// ************************************************************************* //
