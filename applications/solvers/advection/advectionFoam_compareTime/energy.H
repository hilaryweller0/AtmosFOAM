{
    const dimensionedScalar totalSum = fvc::domainIntegrate(T)/totalV;
    const dimensionedScalar diffNorm = (totalSum - totalInitialSum)/totalInitialSum;
    es << runTime.timeName() << " "
       << diffNorm.value() << endl;
}
