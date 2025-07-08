#include <gtest/gtest.h>
#include "cache_handler.hpp"

// Simple test cache structure
struct TestCache {
    double value;
    int count;
};

// Fixture for tests
class CacheHandlerTest : public ::testing::Test {
protected:
    
    void SetUp() override
    {
        cache_handler<TestCache>::end_init();
    };

    void TearDown() override
    {
        cache_handler<TestCache>::end_init();
    };

    cache_handler<TestCache> cache_;
    int fetch_count_ = 0;

    // Mock fetch function that tracks calls
    double& mock_fetch() {
        fetch_count_++;
        static double value = 3.14;
        return value;
    }
};

TEST_F(CacheHandlerTest, CachesValueAfterFirstAccess) {
    auto& result1 = cache_.get(
        [this]() -> double& { return mock_fetch(); },
        [](TestCache& c) -> double& { return c.value; }
    );
    
    auto& result2 = cache_.get(
        []() -> double& { throw std::runtime_error("Shouldn't be called!"); },
        [](TestCache& c) -> double& { return c.value; }
    );

    EXPECT_EQ(result1, 3.14);
    EXPECT_EQ(result2, 3.14);
    EXPECT_EQ(fetch_count_, 1);  // Fetch only called once
}

TEST_F(CacheHandlerTest, ResetForcesRefetch) {
    // First access
    cache_.get(
        [this]() -> double& { return mock_fetch(); },
        [](TestCache& c) -> double& { return c.value; }
    );

    cache_.reset();

    // Second access after reset
    cache_.get(
        [this]() -> double& { return mock_fetch(); },
        [](TestCache& c) -> double& { return c.value; }
    );

    EXPECT_EQ(fetch_count_, 2);  // Fetch called twice
}

TEST_F(CacheHandlerTest, MultipleCacheFields) {
    cache_.get(
        []() -> double& { static double v = 1.0; return v; },
        [](TestCache& c) -> double& { return c.value; }
    );

    cache_.get(
        []() -> int& { static int v = 42; return v; },
        [](TestCache& c) -> int& { return c.count; }
    );

    EXPECT_EQ(cache_.get(
        []() -> double& { throw "Shouldn't execute!"; },
        [](TestCache& c) -> double& { return c.value; }),
        1.0);

    EXPECT_EQ(cache_.get(
        []() -> int& { throw "Shouldn't execute!"; },
        [](TestCache& c) -> int& { return c.count; }),
        42);
}

TEST_F(CacheHandlerTest, InitializationByPassCacheTest)
{
    cache_handler<TestCache>::start_init();

    double* out;

    *out = cache_.get(
            []() -> double& { return mock_fetch(); },
            [](TestCache& c) -> double& { throw "Shouldn't call!"}
            );

    cache_.reset();

    EXPECT_EQ(mock_fetch(),*out);
    
};
    
