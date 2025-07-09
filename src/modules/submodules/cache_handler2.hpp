#pragma once
#include <memory>
#include <type_traits>

template <typename CacheType, typename FieldEnum>
class cache_handler
{
    std::unique_ptr<CacheType> cache_;
    std::unique_ptr<bool[]> field_initialized_; // Tracks initialized fields
    static inline bool initializing_ = false;

public:
    class initializer_guard
    {
    public:
        initializer_guard() { initializing_ = true; }
        ~initializer_guard() { initializing_ = false; }
    };

    // Initialize cache and tracking array on first use
    void ensure_cache_ready(size_t total_fields)
    {
        if (!cache_)
        {
            cache_ = std::make_unique<CacheType>();
            field_initialized_ = std::make_unique<bool[]>(total_fields);
        }
    }

    // Main get_field: Uses enum to identify fields
    template <typename Fetch, typename T>
    T& get_field(Fetch&& fetch, FieldEnum field)
    {
        if (initializing_)
            return fetch();

        const size_t total_fields = static_cast<size_t>(FieldEnum::COUNT);
        ensure_cache_ready(total_fields);

        const size_t index = static_cast<size_t>(field);
        if (!field_initialized_[index])
        {
            cache_->*get_member_ptr(field) = fetch();
            field_initialized_[index] = true;
        }

        return cache_->*get_member_ptr(field);
    }

    // Reset a single field
    void reset_field(FieldEnum field)
    {
        if (field_initialized_)
            field_initialized_[static_cast<size_t>(field)] = false;
    }

    // Reset entire cache
    void reset()
    {
        cache_.reset();
        field_initialized_.reset();
    }

    CacheType& get_cache() { return *cache_; }

    static auto scoped_init() { return initializer_guard(); }

private:
    // Helper: Maps enum -> member pointer (must be specialized for each CacheType)
    auto get_member_ptr(FieldEnum field) -> decltype(&CacheType::air_temperature);
};

//     // Example CacheType and FieldEnum
//     struct MyCache
//     {
//         double air_temperature;
//         int humidity;
//         long pressure;
//     };
//     
//     enum class MyCacheFields
//     {
//         AIR_TEMPERATURE,
//         HUMIDITY,
//         PRESSURE,
//         COUNT // Must be last (used for array size)
//     };
//     
//     // Specialize get_member_ptr for MyCache
//     template <>
//     auto cache_handler<MyCache, MyCacheFields>::get_member_ptr(MyCacheFields field)
//     {
//         switch (field)
//         {
//             case MyCacheFields::AIR_TEMPERATURE: return &MyCache::air_temperature;
//             case MyCacheFields::HUMIDITY: return &MyCache::humidity;
//             case MyCacheFields::PRESSURE: return &MyCache::pressure;
//             default: throw std::invalid_argument("Unknown field");
//         }
//     }i
//     
//     // Fetch air_temperature (double)
//     double& temp = Cache.get_field(
//         [this]() -> double { return (*face)["air_temperature"_s]; },
//         MyCacheFields::AIR_TEMPERATURE
//     );
