class numericalParameters
{
// Public data
public:
    const int nNonLinIters;
    const int nCoriIters;
    const Switch FB;
    
// Construct from dictionary
    numericalParameters(dictionary dict)
    :
    nNonLinIters(readLabel(dict.lookup("nNonLinIters"))),
    nCoriIters(readLabel(dict.lookup("nCoriIters"))),
    FB(dict.lookup("forwardBackwardGravity"))
    {};
};

numericalParameters num(mesh.solution());
