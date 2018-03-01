// Calculating the drag force between two partitions

if (mag(dragCoeff) > SMALL)
{
    // This code only works for nParts == 2
    if (nParts != 2)
    {
        FatalErrorIn("calculateDrag.H")
            << " drag is only formulated for nParts = 2, not nParts = "
            << nParts << exit(FatalError);
    }

    // Calculate the cloud radius and hence the drag
    cloudRadius = max(cloudRadiusMin, cloudRadiusMax*sigmaf[0]*sigmaf[1]);
    dragCommon = dragCoeff*sigmaRhof.sum()/cloudRadius * mag(Uf[0] - Uf[1]);
}
