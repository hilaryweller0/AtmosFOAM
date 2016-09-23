#ifndef analyticSurfaceInterpolationScheme_H
#define analyticSurfaceInterpolationScheme_H

#include "surfaceInterpolationScheme.H"
#include "velocityField.H"
#include "tracerField.H"

namespace Foam
{

template<class Type>
class analyticSurfaceInterpolationScheme
:
    public surfaceInterpolationScheme<Type>
{
private:
    const surfaceScalarField& faceFlux_;
    const IOdictionary tracerDict_;
    autoPtr<velocityField> velocityField_;
    autoPtr<tracerField> tracerField_;

public:
    TypeName("analyticSurfaceInterpolationScheme");

    declareRunTimeSelectionTable
    (
        tmp,
        analyticSurfaceInterpolationScheme,
        Mesh,
        (
            const fvMesh& mesh,
            Istream& schemeData
        ),
        (mesh, schemeData)
    );

	declareRunTimeSelectionTable
	(
		tmp,
		analyticSurfaceInterpolationScheme,
		MeshFlux,
		(   
			const fvMesh& mesh,
			const surfaceScalarField& faceFlux,
			Istream& schemeData
		),
		(mesh, faceFlux, schemeData)
	);

	analyticSurfaceInterpolationScheme
	(
		const fvMesh& mesh,
		const surfaceScalarField& faceFlux,
		Istream&
	)
	:
		analyticSurfaceInterpolationScheme<Type>(mesh, faceFlux)
	{}

	analyticSurfaceInterpolationScheme
	(
		const fvMesh& mesh,
		Istream& is
	)
	:
		analyticSurfaceInterpolationScheme
        (
			mesh,
			mesh.lookupObject<surfaceScalarField>
			(
				word(is)
			)
         )
	{}

	analyticSurfaceInterpolationScheme
	(
		const fvMesh& mesh,
		const surfaceScalarField& faceFlux
	)
	:
		surfaceInterpolationScheme<Type>(mesh),
		faceFlux_(faceFlux),
        tracerDict_
        (
			IOobject
			(
				"tracerFieldDict",
				mesh.time().constant(),
				mesh,
				IOobject::READ_IF_PRESENT,
				IOobject::NO_WRITE
			)
        ),
		velocityField_
		(
			velocityField::New
			(
				IOdictionary
				(
					IOobject
					(
						"velocityFieldDict",
						mesh.time().constant(),
						mesh,
						IOobject::READ_IF_PRESENT,
						IOobject::NO_WRITE
					)
				)
			)
		)
	{
		tracerField_ = tracerField::New
		(
			tracerDict_,
			dynamic_cast<advectable&>(velocityField_())
		);
	}

    virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
    interpolate(const GeometricField<Type, fvPatchField, volMesh>& vf) const
    {
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tsf
        (
            new GeometricField<Type, fvsPatchField, surfaceMesh>
            (
                IOobject("interpolate(" + vf.name() + ")", surfaceInterpolationScheme<Type>::mesh().time().timeName(), surfaceInterpolationScheme<Type>::mesh()),
                surfaceInterpolationScheme<Type>::mesh(),
                vf.dimensions()
            )
        );

        GeometricField<Type, fvsPatchField, surfaceMesh>& sf = tsf.ref();

        sf = pTraits<Type>::zero;
        tracerField_->applyTo(sf);

        return tsf;
    }

    virtual tmp<surfaceScalarField> weights
    (
        const GeometricField<Type, fvPatchField, volMesh>& vf
    ) const
    {

        return pos(this->faceFlux_);
    }
};

}

#endif
