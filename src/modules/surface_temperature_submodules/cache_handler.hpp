#pragma once
#include <memory>

template <typename CacheType>
class cache_handler
{
	std::unique_ptr<CacheType> cache_;
    static inline bool initializing_ = false;
public:

	class initializer_guard
	{
	public:
		initializer_guard() { initializing_ = true; };
		~initializer_guard() { initializing_ = false; };
	};

    template<typename Fetch, typename Access>
    auto& get_field(Fetch&& fetch, Access&& access)
    {
        if (initializing_)
            return fetch();

        if (!cache_)
        {
			cache_ = std::make_unique<CacheType>();
			access(*cache_) = fetch();
        }

        return access(*cache_);

    }

	template<typename Access>
	void set_field(Access&& access, double& value)
	{
		if (!cache_)
		{
			cache_ = std::make_unique<CacheType>();
		};
	
		access(*cache_) = value;
	};

    void reset() { cache_.reset(); };

	CacheType& get_cache()
	{ return *cache_; };

	statis auto scoped_init() { return initializer_guard(); };
};
