const dimensionedScalar cloudRadius(envProperties.lookup("cloudRadius"));
const scalar dragCoeff(readScalar(envProperties.lookup("dragCoeff")));
const scalar heatTransferCoeff
(
    readScalar(envProperties.lookup("heatTransferCoeff"))
);
//const dimensionedScalar thermalConductivity
//(
//    envProperties.lookup("thermalConductivity")
//);
//const dimensionedScalar nu(envProperties.lookup("nu"));
//const dimensionedScalar heatDist(envProperties.lookup("heatDist"));

const dimensionedScalar Ksigma(envProperties.lookup("Ksigma"));

