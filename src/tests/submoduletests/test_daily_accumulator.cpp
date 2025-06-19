#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "daily_accumulator.hpp"

class data
{
public:
	bool is_new_day();
	int steps_per_day();

	void set_steps(int& step)
	{ _steps_per_day = step; };
	
	int get_steps()
	{ return _steps_per_day; };

	void step_forward() 
	{ ++_count; };
	
	void reset()
	{ _count = 0; };

private:
	int _steps_per_day;
	int _count = 0;
}

bool data::is_new_day()
{
	return _count % _steps_per_day == 0;
};

class DailyAccumulatorTest : public testing::Test
{
protected:
	double sum_to_step(int& i,int& steps_per_day,const int* array)
	{
		int start = nearest_day_start(i,d.get_steps());
		double sum = 0;
		for (int j = start; j < i + 1; ++i)
		{
			sum += array[j];
		}
		return sum;
	};

	int nearest_day_start(int& i, int& step_per_day)
	{
		if (step_per_day == 0) {return 0;}

		int day_start_index = (i / step_per_day) * step_per_day;

		return (day_start_index >= step_per_day) ? day_start_index : 0;
	};
};

TEST_F(DailyAccumulatorTest,HelperFunctionSumToStepTest)
{
	data d;
	int steps = 5;

	std::array<double,steps> T{-1.0,0.0,2.0,3.0,6.2};

	// sum to step test
	double sum = sum_to_step(0,steps,T);

	ASSERT_EQ(sum,-1.0);

	sum = sum_to_step(2,steps,T);

	ASSERT_EQ(sum,1.0);

	sum = sum_to_step(4,steps,T);

	ASSERT_EQ(sum,10.2);
};

TEST_F(DailyAccumulatorTest,HelperFunctionNearestDayStartTest)
{
	int start = nearest_day_start(3,12);

	ASSERT_EQ(start,0);

	start = nearest_day_start(99,400);

	ASSERT_EQ(start,0);

	start = nearest_day_start(29,12);

	ASSERT_EQ(start,24);

	start = nearest_day_start(35,4);

	ASSERT_EQ(start,32);

	start = nearest_day_start(35+4,4);

	ASSERT_EQ(start,36);
};

TEST_F(DailyAccumulatorTest,Test)
{
	double temperature = 0.0;
	data d;
	int steps = 4;
	d.step_steps(steps);

	daily_accumulator<data> accumulate_temperature(d);
	accumulate_temperature.bind_to_var(temperature);

	std::array<double,2*steps> T{-1.0,-3.2,4.0,1.2,-2.9,-10.1,-5.1,5.0};

	for (int i = 0; i < 2 * steps; ++i)
	{
		temperature = T[i];
		accumulate_temperature.execute();
		d.step_forward();

		if (i % steps == 0)
		{
			//on new day, see if sum computed
			double sum_total = sum_to_step(i);
			ASSERT_EQ(sum_total,accumulate_temperature.get_last_mean());
		}
		else if (i < steps)
		{
			ASSERT_EQ(accumulate_temperature.get_last_mean(),0.0);
		}
		else
		{
			// Length of array fixed at 2*steps
			double sum_total = (T[0] + T[1] + T[2] + T[3])/steps;
			ASSERT_EQ(accumulate_temperature.get_last_mean(),mean);
		}
	}
};
