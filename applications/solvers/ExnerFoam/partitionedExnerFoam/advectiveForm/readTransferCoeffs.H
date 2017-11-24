Info << "\nReading environmentalProperties" << endl;

IOdictionary transferProps
(
    IOobject
    (
        "transferProperties",
        runTime.system(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

const dimensionedScalar cloudRadiusMax(transferProps.lookup("cloudRadiusMax"));
const dimensionedScalar cloudRadiusMin(transferProps.lookup("cloudRadiusMin"));
const scalar dragCoeff(readScalar(transferProps.lookup("dragCoeff")));
const scalar heatTransferCoeff
(
    readScalar(transferProps.lookup("heatTransferCoeff"))
);
//const dimensionedScalar thermalConductivity
//(
//    envProperties.lookup("thermalConductivity")
//);
//const dimensionedScalar nu(envProperties.lookup("nu"));
//const dimensionedScalar heatDist(envProperties.lookup("heatDist"));

//const dimensionedScalar Ksigma(transferProps.lookup("Ksigma"));

