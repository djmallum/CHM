#pragma once

#include "global.hpp"
#include "triangulation.hpp"
#include <optional>
#include <boost/shared_ptr.hpp>
#include <cstddef>
#include <cmath>
#include <stdexcept>

struct cache_base 
{
    size_t last_timestep = -1;
    bool is_stale(int tn) const { return last_timestep != tn; };
};

template<typename CacheType>
requires std::derived_from<CacheType,cache_base>
class data_base {
    
    void init_cache() const {
        if (!cache_ || cache_->is_stale(global_param->timestep_counter)) {
            cache_.emplace();
            cache_->last_timestep = global_param->timestep_counter;
        }
    }
    
protected:
    mutable std::optional<CacheType> cache_;
    mesh_elem face{nullptr};
    boost::shared_ptr<global> global_param;

    template<typename Value,typename Fetch>
    void update_field(Value& value, Fetch&& fetch) const {
        init_cache();

        if (std::isnan(value)) {
            value = fetch();
        }
    };

    template<typename T>
    void set_output(T& output,const T& t)
    {
        init_cache();

        output = t;
    };

public:
    void set_face(mesh_elem& face_in) { face = face_in; }
    void set_global(boost::shared_ptr<global> param) { global_param = param; };
    void reset_cache() { cache_.reset(); }
    const CacheType& get_cache() const 
    {
        if (!cache_)
           throw std::runtime_error("Cache accessed before it is initialized"); 
        return *cache_; 
    };
};
