// read height at which gravity is swithced off
const scalar gravityOffHeight = envProperties.lookupOrDefault<scalar>
(
    "gravityOffHeight", scalar(0)
);
if (gravityOffHeight > SMALL)
{
    forAll(gd, faceI)
    {
        if (mesh.Cf()[faceI].z() > gravityOffHeight)
        {
            gd[faceI] = 0;
        }
    }
    forAll(gd.boundaryField(), patchi)
    {
        forAll(gd.boundaryField()[patchi], faceI)
        {
            if (mesh.Cf().boundaryField()[patchi][faceI].z() > gravityOffHeight)
            {
                gd.boundaryField()[patchi][faceI] = 0;
            }
        }
    }
}

