class numericalParameters
{
// Public data
public:
    const int nItsU1;
    const int nItsU2;
    const int nItsh;
    
// Construct from dictionary
    numericalParameters(dictionary dict)
    :
    nItsU1(readLabel(dict.lookup("nItsAdv1"))),
    nItsU2(readLabel(dict.lookup("nItsAdv2"))),
    nItsh(readLabel(dict.lookup("nItsGrav")))
    {};
};

numericalParameters num(mesh.solution());
