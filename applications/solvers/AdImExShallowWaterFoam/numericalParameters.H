class numericalParameters
{
// Public data
public:
    const int nOuterCorrs;
    const int nPressureCorrs;
    const int nNonOrthogCorrs;
    const scalar alpha;
    const Switch opSplit;
    
// Construct from dictionary
    numericalParameters(dictionary dict)
    :
    nOuterCorrs(readLabel(dict.lookup("nOuterCorrs"))),
    nPressureCorrs(readLabel(dict.lookup("nPressureCorrs"))),
    nNonOrthogCorrs(readLabel(dict.lookup("nNonOrthogCorrs"))),
    alpha(readScalar(dict.lookup("timeOffCentre"))),
    opSplit(dict.lookup("operatorSplit"))
    {};
};

numericalParameters num(mesh.solution());
