IStringStream sGradType("interpGrad");
tmp<fv::sGradScheme<scalar> > tsGrad = fv::sGradScheme<scalar>::New(mesh, sGradType);
const fv::sGradScheme<scalar>& sGrad = tsGrad();

