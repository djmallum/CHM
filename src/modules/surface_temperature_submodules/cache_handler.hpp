#pragma once
#include <memory>

template <typename CacheType>
class cache_handler
{
    std::unique_ptr<CacheType> cache_;
    static bool initializing_ = false;
public:
    template<typename Fetch, typename Access>
    auto& get_field(Fetch&& fetch, Access&& access)
    {
        if (initializing_)
            return fetch();

        if (!cache)
        {
            cache_ = std::make_unique<CacheType>();
            access(cache_) = fetch();
        }

        return access(cache_);

    }

    void reset() { cache.reset() };
    static void start_init() { initialization = true; };
    static void end_init() { initialization = false; };

};
