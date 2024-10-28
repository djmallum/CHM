#include "I_freeze_thaw_depths.hpp"

class XG_algorithm
{
public:
    ~XG_algorithm() {};
    XG_algorithm(two_layer_DTO& _DTO) : DTO(_DTO) {};

    void run() override;

private:
    two_layer_DTO& DTO;

};
