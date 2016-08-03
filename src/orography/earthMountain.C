#include "earthMountain.H"
#include "fvCFD.H"
#include "addToRunTimeSelectionTable.H"
#include "gdal_priv.h"
#include "cpl_conv.h"

defineTypeNameAndDebug(earthMountain, 0);
addToRunTimeSelectionTable(Mountain, earthMountain, dict);

earthMountain::earthMountain(const IOdictionary& dict) :
    xResolution(readScalar(dict.lookup("xResolution")))
{
    GDALAllRegister();
    const string filename(dict.lookup("filename"));
    terrain = reinterpret_cast<GDALDataset*>(GDALOpen(filename.c_str(), GA_ReadOnly));
    if (terrain == NULL)
    {
        FatalErrorIn("earthMountain::<init>") << exit(FatalError);
    }

    GDALRasterBand* band = terrain->GetRasterBand(1);
    xSize = band->GetXSize();
    ySize = band->GetYSize();
    scanlines = reinterpret_cast<float*>(CPLMalloc(sizeof(float)*xSize));
    band->RasterIO(GF_Read, 0, 3000, xSize, 1, scanlines, xSize, 1, GDT_Float32, 0, 0);
}

scalar earthMountain::heightAt(const point& p) const
{
    label xi = floor(p.x() / xResolution);
    if (xi < 0 || xi >= xSize)
    {
        FatalErrorIn("earthMountain::heightAt") << "mesh domain larger than terrain data (xi " << xi << " outside bounds)" << exit(FatalError);
    }
    return scanlines[xi];
}

scalar earthMountain::gradientAt(const scalar x) const
{
    return 0;
}

scalar earthMountain::timeToCross(const scalar u0, const scalar H) const
{
    return 0;
}

earthMountain::~earthMountain() {
    delete terrain;
    free(scanlines);
}

scalar earthMountain::halfWidth() const
{
    return 0;
}
