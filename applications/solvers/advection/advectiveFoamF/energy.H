{
    scalar totalA = 0;
    scalar totalSum = 0;

    forAll(Tf, facei)
    {
        const scalar value = Tf[facei];
        const scalar area = mesh.delta()()[facei] & mesh.Sf()[facei];
        totalA += area;
        totalSum += value*area;
    }

    totalSum /= totalA;

    const dimensionedScalar diffNorm = (totalSum - totalInitialSum)/totalInitialSum;

    es << runTime.timeName() << " "
       << diffNorm.value() << endl;
}
