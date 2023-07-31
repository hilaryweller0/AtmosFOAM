Info<< "Reading thermophysical properties\n" << endl;

autoPtr<fluidThermo> pThermo
(
    fluidThermo::New(mesh)
);
fluidThermo& thermo = pThermo();

thermo.validate(args.executable(), "h", "e");

// Use of Exner needs reference T and p and a single value of R, Cp and Cv
const dimensionedScalar Tref
(
    "Tref", dimTemperature, readScalar(thermo.properties().lookup("Tref"))
);
const dimensionedScalar pRef
(
    "pRef", dimPressure, readScalar(thermo.properties().lookup("pRef"))
);

const constTransport<hConstThermo<perfectGas<specie> > > air
(
    thermo.properties().subDict("mixture")
);

const dimensionedScalar R("R", dimGasConstant, air.R());
const dimensionedScalar Cp
(
    "Cp", dimGasConstant, air.Cp(pRef.value(),Tref.value())
);
const dimensionedScalar kappa("kappa", R/Cp);
const scalar gamma = (1/(1-kappa)).value();
const dimensionedScalar Cv("Cv", Cp - R);
const dimensionedScalar speedSound("speedSound", sqrt(gamma*R*Tref));
const dimensionedScalar RbyCv("RbyCv", R/Cv);

Info << "Cp = " << Cp.value() << " speedSound = " << speedSound << " kappa = "
     << kappa << " R = " << air.R() << " RR = " << RR << " Cv = " << Cv 
     << " RbyCv = " << RbyCv << endl;

Info << "Maximum acoustic Courant number = "
     << max(mesh.deltaCoeffs())*runTime.deltaT()*speedSound << endl;

