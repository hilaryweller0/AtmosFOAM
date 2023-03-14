#include "OStringStream.H"
#include "scalar.H"

namespace Foam
{

    inline const string stringScalar(const scalar x, const signed nsfs = 0)
    {
        OStringStream os;
        if (nsfs) os.precision(nsfs);
        os << x;
        return os.str();
    }

}
