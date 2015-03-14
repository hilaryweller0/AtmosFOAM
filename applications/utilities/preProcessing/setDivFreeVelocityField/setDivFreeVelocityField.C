#include "fvCFD.H"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"        

    // If p does not have any fixed boundaries, a reference value will be needed
    bool setRefp = true;
    for(label iBC = 0; iBC < mesh.boundary().size() && setRefp; iBC++)
    {
        if (p.boundaryField()[iBC].type() == "fixedValue")
        {
            setRefp = false;
        }
    }

    bool converged = false;
    const int ncorrMax = 100;
    for(int icorr = 0; icorr < ncorrMax && !converged; icorr++)
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(p) + fvc::div(phi)
        );
        if (setRefp) pEqn.setReference(0,0);
        const int nIters = pEqn.solve().nIterations();
        converged = nIters == 0;

        if (nIters == 0 || icorr == ncorrMax-1) phi += pEqn.flux();
    }

    phi.write();
}
