#pragma once

#include "../logger.hpp"
#include "triangulation.hpp"
#include "module_base.hpp"

#include "TPSpline.hpp"

#include <cstdlib>
#include <string>

#include <cmath>

#include <math.h>

/**
* \addtogroup modules
* @{
* \class tair_llra_lookup
* \brief Constant monthly linear lapse rate adjustment
*
* Constant monthly linear lapse rate adjustment for air temperature.
*
* Month | Lapse rate (degC/m)
* ------|----------------
* Jan | 0.0044
* Feb | 0.0059
* Mar | 0.0071
* Apr | 0.0078
* May | 0.0081
* Jun | 0.0082
* Jul | 0.0081
* Aug | 0.0081
* Sep | 0.0077
* Oct | 0.0068
* Nov | 0.0055
* Dec | 0.0047
*
* Depends:
* - None
*
* Provides:
* - Air temperature "t" [degC]
* - Air temperatue "tair_llra_lookup" [degC]
*
* Reference:
* > Liston, Glen E., and Kelly Elder. 2006. A meteorological distribution system for high-resolution terrestrial modeling (MicroMet). Journal of hydrometeorology 7: 217-234

*/
class tair_llra_lookup : public module_base
{
public:
    tair_llra_lookup(std::string ID);
    ~tair_llra_lookup();
    virtual void run(mesh_elem& elem, boost::shared_ptr<global> global_param);
};

/**
@}
*/