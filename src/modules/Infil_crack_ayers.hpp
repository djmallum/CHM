#include "Ayers_infiltration.hpp"
#include "Crack.hpp"

template<typename Data>
class Infil_crack_ayers
{
public:
    Crack::Model<Data> crack;
    Ayers::Model<Data> ayers;


private:
    Status_Checker status_checker;
}
