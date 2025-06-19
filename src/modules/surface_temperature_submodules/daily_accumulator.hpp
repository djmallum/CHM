#include "base_step.hpp"

template<typename T>
concept accumulator_data = requires(T& t) {
    // Must have is_new_day() const member function returning bool
    { t.is_new_day() } -> std::same_as<bool>;
    
    // Must have steps_per_day() const member function convertible to double
    { t.steps_per_day() } -> std::convertible_to<double>;
};

template<accumulator_data data>
class daily_accumulator : public base_step<data>
{
	const double* _target_var = nullptr;
	double _accumulator = 0.0;
	double mean_value = 0.0;
public:
	explicit daily_accumulator(data& _d) : base_step<data>(_d) {};
	explicit daily_accumulator(data& _d, double& value) : Base<data>(_d), mean_value(value) {};

	void bind_to_var(double& var)
	{
		_target_var = &var;
	};

	const double& get_last_mean() const
	{
		return mean_value;
	};

	void execute() override final
	{
		if (!_target_var) return;

		_accumulator += *_target_var;

		if (this->d.is_new_day())
		{
			mean_value = _accumulator / this->d.steps_per_day();
			_accumulator = 0.0;
		};
	};
};
