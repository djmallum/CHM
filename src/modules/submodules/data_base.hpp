#pragma once

#include "global.hpp"
#include "triangulation.hpp"
#include <optional>
#include <boost/shared_ptr.hpp>
#include <boost/property_tree/ptree.hpp>
#include <cstddef>
#include <cmath>
#include <stdexcept>
#include <cassert>
#include <cstdint>

/*
 * Base class for data objects defined in each module to cache and acess triangle specific data.
 *
 * data_base is the base class. It expects a single template parameter CacheType. CacheType is always 
 * derived from cache_base and this is enforced by the CacheRules concept.
 *
 * CacheType should be unique to each module, just the data class is unique to each object. The 
 * CacheType class should only store raw variables and initialized as NaN. This NaN initialization 
 * is done to allow member-by-member checking if the member has been cached and whether it should be 
 * accessed from face during each time step.
 *
 * data_base is written specifically to cater to the base_step class. Derived classes of base_step
 * have a single template parameter to a class/struct to access inputs and parameters. It is intended 
 * that the data class is used in this template. 
 *
 * data_base provides tools to write functions to that the base_step template expects to pass inputs
 * and parameters set at runtime. 
 *
 * Basic usage:
 *
 * // In the module header file define the Cache and data objects within the module definition
 * class Cache : public cache_base
 * {
 * 	   // Has to be NaN does not have to be this version of NaN
 * 	   double input_variable = std::numeric_limits::quiet_NaN();
 * };
 * 
 * //
 * class data : public data_base<Cache>
 * {
 * public:
 *     // Example of an input/provides
 * 	   double& input_variable();
 *     // Example of domain wide parameter set in the config
 * 	   double& config_parameter();
 *     // Example of parameter access from global instance
 * 	   int& get_dt();
 *     // Example of output function
 * 	   void output_variable(const double& in);
 *	   // Example of a parameter that varies on each face
 *	   double& spatial_parameter();
 * private:
 *
 *	   double spatial_parameter_ = face->veg_attribute("spatial_parameter"_s);
 * };
 * 
 * double& data::input_variable() const
 * {
 * 	   update_field(cache_->input_variable,
 * 	   		[this]() { return (*face)["input_variable"_s];} );
 * 
 * 	   return cache_->input_variable;
 * };
 * 
 * double& data::config_parameter() const
 * {
 * 	   static const double result = cfg_.get("config_parameter",1.0);
 * 
 * 	   return result;
 * };
 *
 * double& data::spatial_parameter() const
 * {
 *	   return spatial_parameter_;
 * };       
 * 
 * int& data::get_dt() const
 * {
 *     // Don't actually need to do 
 * 	   static const double dt = global_param->dt();
 * 
 * 	   return dt;
 * };
 * 
 * void data::output_variables(const double& out)
 * {
 * 	   set_output(cache_->output_variable,out);
 * };
 * 
 * void module_name::run(mesh_elem& face)
 * {
 * 	   // run whatever calculations here, using the data object 
 * 
 * 	   // manually reset the cache
 * 	   d.reset_cache();
 * };
 *
 */ 

struct cache_base 
{
    int64_t last_timestep = -1;
    bool is_stale(int64_t tn) const { return last_timestep != tn; };
};

template<typename C>
concept CacheRules = std::derived_from<C,cache_base>;

namespace pt = boost::property_tree;

template<CacheRules CacheType>
class data_base {
    
    void init_cache() const;

protected:
    
    data_base(const mesh_elem& face_in, const boost::shared_ptr<global> param, 
            const pt::ptree& cfg, const bool istest = false);
    ~data_base() {};
    
    const mesh_elem face{nullptr};
    const boost::shared_ptr<global> global_param;
    const pt::ptree& cfg_;
    mutable std::optional<CacheType> cache_;

    template<typename Value,typename Fetch>
    void update_field(Value& value, Fetch&& fetch) const;

    template<typename T>
    void set_output(T& output,const T& t) const;

public:
    void reset_cache() { cache_.reset(); };
    const std::optional<CacheType>& get_cache() { return cache_; }; 
};

template<CacheRules CacheType>
void data_base<CacheType>::init_cache() const {
    if (!cache_ || cache_->is_stale(global_param->timestep_counter)) {
        cache_.emplace();
        cache_->last_timestep = global_param->timestep_counter;
    }
}

template<CacheRules CacheType>
data_base<CacheType>::data_base(const mesh_elem& face_in, const boost::shared_ptr<global> param, 
        const pt::ptree& cfg, const bool istest) : face(face_in), global_param(param), cfg_(cfg)
{
	// Optional istest parameter only exists to skip these tests during tests of this class where we aren't testing whether the face object has been set correctly.
	// This means tests show that the underlying functions work as intended
    if (!face->is_valid() && !istest)
        throw std::invalid_argument("Face handle points to an invalid face");

    if (!global_param && !istest)
        throw std::invalid_argument("global parameter holder is null");
};

template<CacheRules CacheType>
template<typename Value,typename Fetch>
void data_base<CacheType>::update_field(Value& value, Fetch&& fetch) const {
    init_cache();

    if (std::isnan(value)) {
        value = fetch();
    }
};

template<CacheRules CacheType>
template<typename T>
void data_base<CacheType>::set_output(T& output,const T& t) const
{
    init_cache();

    output = t;
};

