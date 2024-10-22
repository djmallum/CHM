#include "I_K_estimate.hpp"

class K_estimate : public I_K_estimate
{
public:
    K_estimate(K_DTO& _DTO);
    ~K_estimate();

    void run(void) override;

private:

    K_DTO& DTO;
    struct params;
    struct Darcy_Vels;
    
};
