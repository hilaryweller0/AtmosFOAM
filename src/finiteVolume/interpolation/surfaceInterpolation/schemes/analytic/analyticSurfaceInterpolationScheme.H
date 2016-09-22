#ifndef analyticSurfaceInterpolationScheme_H
#define analyticSurfaceInterpolationScheme_H

#include "surfaceInterpolationScheme.H"

namespace Foam
{

template<class Type>
class analyticSurfaceInterpolationScheme
:
    public surfaceInterpolationScheme<Type>
{
private:
    const surfaceScalarField& faceFlux_;

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
		surfaceInterpolationScheme<Type>(mesh),
		faceFlux_
		(
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
		faceFlux_(faceFlux)
	{}

    virtual tmp<surfaceScalarField> weights
    (
        const GeometricField<Type, fvPatchField, volMesh>&
    ) const
    {
        return pos(this->faceFlux_);
    }
};

}

#endif
