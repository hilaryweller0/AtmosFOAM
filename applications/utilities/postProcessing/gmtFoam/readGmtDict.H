Info << "\nReading "
     << runTime.constant()/"gmtDicts"/args.additionalArgs()[0] << endl;

IOdictionary gmtDict
(
    IOobject
    (
        args.additionalArgs()[0],
        runTime.constant()/"gmtDicts",
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    )
);

// How to plot
// Compulsory inputs
const string projection(gmtDict.lookup("projection"));
const word projection2D(gmtDict.lookup("projection2D"));

// Fields to plot
const List<FieldToPlot> fieldsToPlot(gmtDict.lookup("FieldsToPlot"));

// Optional
const bool plotBoundaryDots = gmtDict.lookupOrDefault<bool>
(
    "plotBoundaryDots", true
);
const bool plotZeroValuePolys = gmtDict.lookupOrDefault<bool>("plotZeroValuePolys", true);
const word patchName = gmtDict.found("patch") ?
                       word(gmtDict.lookup("patch")) : "";
const word boundaryPatchName = gmtDict.lookupOrDefault<word>("boundaryPatch", "");

const string region = gmtDict.found("region") ?
                      string(gmtDict.lookup("region")) :
                      setRegion(mesh, projection2D);
const scalar rMin = gmtDict.lookupOrDefault<scalar>("rMin", scalar(0));
const scalar rMax = gmtDict.lookupOrDefault<scalar>("rMax", scalar(0));
const scalar lonMin = gmtDict.lookupOrDefault<scalar>("lonMin", scalar(0));
const scalar lonMax = gmtDict.lookupOrDefault<scalar>("lonMax", scalar(0));
const scalar latMin = gmtDict.lookupOrDefault<scalar>("latMin", scalar(0));
const scalar latMax = gmtDict.lookupOrDefault<scalar>("latMax", scalar(0));
const string boundaryMarks = gmtDict.found("boundaryMarks") ?
                      string(gmtDict.lookup("boundaryMarks")) : "0/0";
const bool separateLegend = gmtDict.found("separateLegend") ?
                         readBool(gmtDict.lookup("separateLegend")) : false;

const bool overlay = gmtDict.found("overlay") ?
                     readBool(gmtDict.lookup("overlay")) : false;
const bool lastPlot = gmtDict.found("lastPlot") ?
                     readBool(gmtDict.lookup("lastPlot")) : true;

IStringStream annotation = gmtDict.found("annotate") ?
                           IStringStream(gmtDict.lookup("annotate")) : "";

const word annotate = (annotation.str() == "") ? "" : word(annotation);
const scalar annotateX = (annotation.str() == "") ? 0 : readScalar(annotation);
const scalar annotateY = (annotation.str() == "") ? 0 : readScalar(annotation);

IStringStream vectorMax = gmtDict.found("writeVectorMax") ?
                          IStringStream(gmtDict.lookup("writeVectorMax")) : "";
const bool writeVectorMax = (vectorMax.str() == "") ? false : true;
const scalar vectorMaxX = writeVectorMax ? readScalar(vectorMax) : 0;
const scalar vectorMaxY = writeVectorMax ? readScalar(vectorMax) : 0;

const string xlabel = gmtDict.lookupOrDefault<string>("xlabel", "");
const string ylabel = gmtDict.lookupOrDefault<string>("ylabel", "");
const Pair<scalar> xlabelLoc = gmtDict.lookupOrDefault<Pair<scalar> >
        ("xlabelLoc", Pair<scalar>(scalar(0), scalar(0)));
const Pair<scalar> ylabelLoc = gmtDict.lookupOrDefault<Pair<scalar> >
        ("ylabelLoc", Pair<scalar>(scalar(0), scalar(0)));

Info << "done" << endl;

// output file name:
fileName epsFileName = fileName
(
    gmtDict.lookupOrDefault<fileName>
    (
        "overlayOnto",
        args.additionalArgs()[0]
    )
);

if (fieldsToPlot.size() == 0)
{
    FatalErrorIn("readGmtDict") << " nothing to plot";
    exit(FatalError);
}
else
{
    Info << "Plotting " << fieldsToPlot << " into "
         << epsFileName << ".eps for each time" << endl;
}
