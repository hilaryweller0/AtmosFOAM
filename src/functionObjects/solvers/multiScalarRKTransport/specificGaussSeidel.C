else // Specific symmetric Gauss-Seidel for not density weighted
{
    volScalarField RHS = s_[is].oldTime() - dt*fvc::div(totalFlux[is]);
    s_[is] = s_[is].oldTime();
    volScalarField denom = 1 + dt*fvc::surfaceIntegrateOut(phi);
    const int nSweeps = 2;
    const labelList& own = mesh_.owner();
    const labelList& nei = mesh_.neighbour();
    for (int isw = 0; isw < nSweeps; isw++)
    {
        // Forward sweeps
        for(label icell = 0; icell < s_[is].size(); icell++)
        {
            scalar sumIn = 0;
            scalar V = mesh_.V()[icell];
            const cell& c = mesh_.cells()[icell];
            for(label ifac = 0; ifac < c.size(); ifac++)
            {
                label iF = c[ifac];
                label N = (own[iF] == icell) ? nei[iF] : own[iF];
                if
                (
                    ((phi[iF] < 0 && own[iF] == icell)
                    || (phi[iF] > 0 && nei[iF] == icell))
                    && N < mesh_.nCells()
                )
                {
                    sumIn += mag(phi[iF])*s_[is][N];
                }
            }

            s_[is][icell] = (RHS[icell] + dt.value()/V*sumIn)
                            /denom[icell];
        }
        // Backward sweeps
        for(label icell = s_[is].size()-1; icell >= 0 ; icell--)
        {
            scalar sumIn = 0;
            scalar V = mesh_.V()[icell];
            const cell& c = mesh_.cells()[icell];
            for(label ifac = 0; ifac < c.size(); ifac++)
            {
                label iF = c[ifac];
                label N = (own[iF] == icell) ? nei[iF] : own[iF];
                if
                (
                    ((phi[iF] < 0 && own[iF] == icell)
                    || (phi[iF] > 0 && nei[iF] == icell))
                    && N < mesh_.nCells()
                )
                {
                    sumIn += mag(phi[iF])*s_[is][N];
                }
            }

            s_[is][icell] = (RHS[icell] + dt.value()/V*sumIn)
                            /denom[icell];
        }
    }
}
