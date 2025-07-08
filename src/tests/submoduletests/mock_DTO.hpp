#include <gmock/gmock.h>
#include "soil_DTO.hpp"

class MockTwoLayerDTO : public two_layer_DTO
{
public:
    MOCK_METHOD(int,get_dt,(), (override));
    MOCK_METHOD(bool,get_new_day,(), (override));


	void set_new_day_true()
	{
		ON_CALL(*this, get_new_day).WillByDefault(::testing::Return(true));
	};

	void set_dt_3600()
	{
		ON_CALL(*this, get_dt).WillByDefault(::testing::Return(3600));
	};

};

class MockETDTO : public soil_ET_DTO
{
public:

};
