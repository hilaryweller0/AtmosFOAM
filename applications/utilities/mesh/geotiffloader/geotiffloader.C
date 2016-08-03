#include <stdlib.h>
#include "fvCFD.H"
#include "gdal_priv.h"
#include "cpl_conv.h"

int main(int argc, char *argv[])
{
    GDALAllRegister();
    autoPtr<GDALDataset> terrain(reinterpret_cast<GDALDataset*>(GDALOpen("/home/jshaw/terrain/030363012518617/ASTGTM2_N47E010_dem.tif", GA_ReadOnly)));
    if (terrain.empty())
    {
        FatalErrorIn("geotiffloader::main") << exit(FatalError);
    }

    GDALRasterBand* band = terrain->GetRasterBand(1);
    int xSize = band->GetXSize();
    float *scanline = (float*) CPLMalloc(sizeof(float)*xSize);
    band->RasterIO(GF_Read, 0, 3000, xSize, 1, scanline, xSize, 1, GDT_Float32, 0, 0);

    for (label i = 0; i < xSize; i++)
    {
        Info << i << " " << scanline[i] << endl;
    }

    free(scanline);
    
    return EXIT_SUCCESS;
}

