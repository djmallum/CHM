// k_estimator_mock.h
#pragma once
#include "submodule_base.hpp"
#include <gmock/gmock.h>

class KEstimatorMock : public I_K_estimate {
public:
    virtual ~KEstimatorMock() = default;
    
    MOCK_METHOD(void, run, (), (override));
    
    // Mock any other pure virtual functions from submodule_base if needed
    // For example:
    // MOCK_METHOD(void, initialize, (), (override));
};
