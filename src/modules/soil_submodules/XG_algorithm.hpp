#include "I_freeze_thaw_depths.hpp"
#include <cmath>

class XG_algorithm : public I_freeze_thaw_depths
{
public:
    ~XG_algorithm() {};
    XG_algorithm() {};

    virtual void run() override;

};
