
#include "tntn/RasterIO.h"
#include "fmt/format.h"
#include "tntn/logging.h"
#include "tntn/println.h"
#include "tntn/util.h"
#include "tntn/raster_tools.h"

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <iomanip>

#include <ogr_spatialref.h>
#include <gdal_priv.h>

#include "tntn/gdal_init.h"

namespace tntn {


// --------------------------------------------------------------------------------
// Import raster from GDAL

typedef std::unique_ptr<GDALDataset, void (*)(GDALDataset*)> GDALDataset_ptr;

struct TransformationMatrix
{
    union
    {
        double matrix[6];
        struct
        {
            double origin_x;
            double scale_x;
            double padding_0;
            double origin_y;
            double padding_1;
            double scale_y;
        };
    };
};

static void GDALClose_wrapper(GDALDataset* p)
{
    if(p != nullptr)
    {
        GDALClose(p);
    }
}

static bool get_transformation_matrix(GDALDataset* dataset, TransformationMatrix& gt)
{
    if(dataset->GetGeoTransform(gt.matrix) != CE_None)
    {
        TNTN_LOG_ERROR("Input raster is missing geotransformation matrix");
        return false;
    }

    if(fabs(gt.scale_x) != fabs(gt.scale_y))
    {
        TNTN_LOG_ERROR("Can not process rasters with non square pixels ({}x{})",
                       fabs(gt.scale_x),
                       fabs(gt.scale_y));
        return false;
    }

    return true;
}

static bool is_valid_projection(GDALDataset* dataset)
{
    // returned pointer should not be altered, freed or expected to last for long.
    const char* projection_wkt = dataset->GetProjectionRef();

    if(projection_wkt == NULL)
    {
        TNTN_LOG_ERROR("Input raster file does not provide spatial reference information");
        return false;
    }

    OGRSpatialReference raster_srs;

#if GDAL_VERSION_MAJOR == 2 && GDAL_VERSION_MINOR < 3
    raster_srs.importFromWkt(const_cast<char**>(&projection_wkt));
#else
    raster_srs.importFromWkt(projection_wkt);
#endif

    int matches_number;

    OGRSpatialReference web_mercator;
    web_mercator.importFromEPSG(3857);

    bool matched = false;

#if GDAL_VERSION_MAJOR == 2 && GDAL_VERSION_MINOR < 3
    OGRErr match_projection_error = raster_srs.AutoIdentifyEPSG();

    if(match_projection_error != 0)
    {
        TNTN_LOG_ERROR("Can not match projection to EPSG:3857");
        return false;
    }

    if(web_mercator.IsSame(&raster_srs))
    {
        matched = true;
    }
#else
    // Matched projections must be freed with OSRFreeSRSArray()
    OGRSpatialReferenceH* matched_projections =
            raster_srs.FindMatches(NULL, &matches_number, NULL);

    if(matched_projections == 0)
    {
        TNTN_LOG_ERROR("Can not match projection to EPSG:3857");
        return false;
    }

    for(int i = 0; i < matches_number; i++)
    {
        if(web_mercator.IsSame(static_cast<OGRSpatialReference*>(matched_projections[i])))
        {
            matched = true;
            break;
        }
    }

    if(matched_projections)
    {
        OSRFreeSRSArray(matched_projections);
    }
#endif

    if(!matched)
    {
        TNTN_LOG_ERROR("Can not match projection to EPSG:3857");
    }

    return matched;
}

bool load_raster_file(const std::string& file_name,
                      RasterDouble& target_raster,
                      bool validate_projection)
{
    initialize_gdal_once();

    TNTN_LOG_INFO("Opening raster file {} with GDAL...", file_name);

    GDALDataset_ptr dataset(static_cast<GDALDataset*>(GDALOpen(file_name.c_str(), GA_ReadOnly)),
                            &GDALClose_wrapper);

    if(dataset == nullptr)
    {
        TNTN_LOG_ERROR("Can't open input raster {}: ", file_name);
        return false;
    }

    TransformationMatrix gt;

    if(!get_transformation_matrix(dataset.get(), gt))
    {
        return false;
    }

    if(validate_projection && !is_valid_projection(dataset.get()))
    {
        println("input raster file must be in EPSG:3857 (Web Mercator) format");
        println("you can reproject raster terrain using GDAL");
        println("as follows: 'gdalwarp -t_srs EPSG:3857 input.tif output.tif'");

        return false;
    }

    int bands_count = dataset->GetRasterCount();

    if(bands_count == 0)
    {
        TNTN_LOG_ERROR("Can't process a raster file witout raster bands");
        return false;
    }
    else if(bands_count > 1)
    {
        TNTN_LOG_WARN("File {} has {} raster bands, processing with a raster band #1",
                      file_name,
                      bands_count);
    }

    // TODO: Perhaps make raster band number a parameter
    GDALRasterBand* raster_band = dataset->GetRasterBand(1);

    int raster_width = raster_band->GetXSize();
    int raster_height = raster_band->GetYSize();

    target_raster.set_cell_size(fabs(gt.scale_x));
    target_raster.allocate(raster_width, raster_height);
    target_raster.set_no_data_value(raster_band->GetNoDataValue());

    TNTN_LOG_INFO("reading raster data...");
    if(raster_band->RasterIO(GF_Read,
                             0,
                             0,
                             raster_width,
                             raster_height,
                             target_raster.get_ptr(),
                             raster_width,
                             raster_height,
                             GDT_Float64,
                             0,
                             0) != CE_None)
    {
        TNTN_LOG_ERROR("Can not read raster data");
        return false;
    }

    double x1 = gt.origin_x;
    double y1 = gt.origin_y;
    double x2 = gt.origin_x + raster_width * gt.scale_x;
    double y2 = gt.origin_y + raster_height * gt.scale_y;

    // Ensure raster's origin is exactly at the lower left corner
    target_raster.set_pos_x(std::min(x1, x2));
    target_raster.set_pos_y(std::min(y1, y2));

    if(gt.scale_x < 0)
    {
        raster_tools::flip_data_x(target_raster);
    }

    if(gt.scale_y > 0)
    {
        raster_tools::flip_data_y(target_raster);
    }

    return true;
}

} // namespace tntn
