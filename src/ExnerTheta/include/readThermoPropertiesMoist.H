Info << "Reading additional thermophysical properties for moist air"
     << endl;

// Read constant for the Clausius-Clapyron relation
const scalar Pcc = readScalar(thermoDict.lookup("pClausiusClapyron"));

// Perfect gas properties of water vapour
const constTransport<hConstThermo<perfectGas<specie> > > waterVapour
(
    thermoDict.lookup("waterVapour")
);

// Latent heat release at reference temperature
const dimensionedScalar Lv0(thermoDict.lookup("latentHeatRef"));

// Specific heat of water vapour at constant pressure

const dimensionedScalar Cpv
(
    "Cpv",
    dimGasConstant,
    waterVapour.cp(pRef.value(), T0.value())/waterVapour.W()
);

// Specific heat of liquid water at constant pressure
const dimensionedScalar Cpl(thermoDict.lookup("specificHeatLiquidWater"));

// Gas constant of water vapour
const dimensionedScalar Rv("Rv", dimGasConstant, waterVapour.R());

// Ratio between molecular weight of water vapour and dry air
const scalar epsilon = (R/Rv).value();

// Speicific heat of air at constant volume
const dimensionedScalar Cvv("Cvv", Cpv - Rv);

Info << "Moist thermodynamic constants Cpv = " << Cpv
     << " epsilon = " << epsilon << " Rv = " << Rv << " Cvv = " << Cvv << endl;

