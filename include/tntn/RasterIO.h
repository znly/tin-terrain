#pragma once

#include "Raster.h"
#include "tntn/raster_tools.h"
#include <string>

namespace tntn {

bool load_raster_file(const std::string& filename,
                      RasterDouble& raster,
                      bool validate_projection = true);
} // namespace tntn
