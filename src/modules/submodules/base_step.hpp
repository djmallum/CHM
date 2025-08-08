#pragma once
#include <stdexcept>
#include <concepts>

template<class data>
class base_step
{
public:
	explicit base_step() {}; 
	virtual ~base_step() = default;

	virtual void execute(data& d) = 0;
};


