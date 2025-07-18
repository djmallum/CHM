#pragma once
#include <stdexcept>
#include <iostream>
#include <concepts>

#define THROW_NULL_POINTER_EXCEPTION() \
	throw std::runtime_error( \
			std::string("Null pointer at ") + __File__ + ":" + std::to_string(__Line__) \
			)


template<class data>
class base_step
{
public:
	explicit base_step() {}; 
	virtual ~base_step() = default;

	virtual void execute(data& d) = 0;
};


