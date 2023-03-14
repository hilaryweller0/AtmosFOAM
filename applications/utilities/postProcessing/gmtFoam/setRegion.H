#include "boundBox.H"

inline const string convertToString(const scalar x)
{
    OStringStream os;
    os << x;
    return os.str();
}

inline const string setRegion(const fvMesh& mesh, const word& projection2D)
{
    const point& min = mesh.bounds().min();
    const point& max = mesh.bounds().max();

    Info << "Mesh bounding box = " << min << "  " << max << endl;

    string r = "";

    switch (projection2D[0])
    {
        case 'X':
        {
            r = convertToString(min[0])/convertToString(max[0]);
        } break;
        case 'Y':
        {
            r = convertToString(min[1])/convertToString(max[1]);
        } break;
        case 'Z':
        {
            r = convertToString(min[2])/convertToString(max[2]);
        } break;
        default:
        {
            r = 'g';
        }
    }
    switch (projection2D[1])
    {
        case 'X':
        {
            r = r/convertToString(min[0])/convertToString(max[0]);
        } break;
        case 'Y':
        {
            r = r/convertToString(min[1])/convertToString(max[1]);
        } break;
        case 'Z':
        {
            r = r/convertToString(min[2])/convertToString(max[2]);
        } break;
    }

    Info << "Region = " << r << endl;

    return r;
}

inline boundBox setRegion
(
    const fvMesh& mesh, const word& projection2D, const string& region,
    const scalar earthRadius
)
{
    boundBox bb = mesh.bounds();
    vector bbExpand = 0.01*bb.span();
    bb.min() -= bbExpand;
    bb.max() += bbExpand;

    if (region.size() > 4)
    {
        string regionStr = region;
        regionStr.replaceAll("/", " ");
        IStringStream bbStr(regionStr);
        scalar minx, miny, maxx, maxy = 0;
        bbStr >> minx >> maxx >> miny >> maxy;

        if (maxx > minx && maxy > miny)
        {
            switch (projection2D[0])
            {
                case 'X':
                    {
                        bb.min().x() = minx;
                        bb.max().x() = maxx;
                    } break;
                case 'Y':
                    {
                        bb.min().y() = minx;
                        bb.max().y() = maxx;
                    } break;
                case 'Z':
                    {
                        bb.min().z() = minx;
                        bb.max().z() = maxx;
                    } break;
            }

            switch (projection2D[1])
            {
                case 'X':
                    {
                        bb.min().x() = miny;
                        bb.max().x() = maxy;
                    } break;
                case 'Y':
                    {
                        bb.min().y() = miny;
                        bb.max().y() = maxy;
                    } break;
                case 'Z':
                    {
                        bb.min().z() = miny;
                        bb.max().z() = maxy;
                    } break;
            }

            if (projection2D == "latLon")
            {
                const scalar radToDeg = 180./constant::mathematical::pi;
                bb.min().z() = earthRadius*Foam::sin(miny/radToDeg);
                bb.max().z() = earthRadius*Foam::sin(maxy/radToDeg);
            }
            else if (projection2D == "cylinder")
            {
                bb.min().x() = minx;
                bb.max().x() = maxx;
                bb.min().z() = miny;
                bb.max().z() = maxy;
            }
        }
    }

    return bb;
}
