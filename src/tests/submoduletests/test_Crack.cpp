#include <concepts>
#include <functional>
#include <gtest/gtest.h>
#include <stdexcept>
#include "Crack.hpp"

using namespace Crack;

class CrackTestBase : public ::testing::Test
{
    size_t NONEcount = 0;
    size_t FIRST_MAJORcount = 0;
    size_t PRIOR_INFILTRATIONcount = 0;
    size_t RESETcount = 0;
    size_t NORMALcount = 0;

    typedef std::function<double()> func;

protected:

    struct output
    {
        State s;
        bool cont = true;
        func swe;
        LimitedPhase phase;
    };

    output prepareNONE()
    {
        State s;
        bool cont = true;
        switch (NONEcount)
        {
            case 0:
                // ice lens
                s.major_melt_count = 1;
                State::lenstemp = -5.0;
                s.daily_max_temp = State::lenstemp - 2.0;
                s.soil_saturation_at_freeze = 25.0;
        break;
            case 1:
                // Another ice lens
                s.major_melt_count = 50;
                State::lenstemp = -35.0;
                s.daily_max_temp = -100.0;
                s.soil_saturation_at_freeze = 500.0;
                break;
            case 2:
                // missed everything
                s.major_melt_count = 50;
                State::infDays = 6;
                State::lenstemp = -10.0;
                s.daily_max_temp = 5.0;
                State::major_melt_threshold = 5.0;
                cont = false;
                break;
        }
        
        ++NONEcount;
        return output{s,cont,[] () {return -99999.0;}};
    };
    
    output prepareFIRST_MAJOR()
    {
        State s;
        bool cont = true;
        std::function<double()> swe;
        switch (FIRST_MAJORcount)
        {
            case 0:
                State::lenstemp = -10.0;
     s.daily_max_temp = 0.0;
                s.major_melt_count = 0;
                // hack so I don't have to set dailt_melt_total
                State::major_melt_threshold = -5.0;
                State::infDays = 6;
                swe = []() { return 100.0;};
                break;
            case 1:
                // same same but different
                State::lenstemp = -5.0;
                s.daily_max_temp = 0.0;
                s.major_melt_count = 0;
                // hack so I don't have to set dailt_melt_total
                State::major_melt_threshold = -30.0;
                State::infDays = 6;
                swe = []() { return 22.0;};
                cont = false;
                break;
        }

        ++FIRST_MAJORcount;
        return output{s,cont,swe};
    };

    output preparePRIOR_INFILTRATION()
    {
        output o;
        switch (PRIOR_INFILTRATIONcount)
        {
            case 0:
                FIRST_MAJORcount = 0;
                o = prepareFIRST_MAJOR();
                State::major_melt_threshold = 5.0;
                o.swe = []() {return 0.0;};
                o.cont = false;
                break;
        }

        ++PRIOR_INFILTRATIONcount;
        return o;
    };

    output prepareRESET()
    {
        State s;
        bool cont = true;
        output o;
        switch (RESETcount)
        {
            case 0:
                FIRST_MAJORcount = 0;
                o = prepareFIRST_MAJOR();
                o.s.major_melt_count = 1;
                State::infDays = 6;
                State::major_melt_threshold = -5.0;
                o.s.init_SWE = 0.0;
                o.swe = []() {return 1.0;};
                o.s.daily_max_temp = 100.0;
                break;
            case 1:
                //same same but different
                FIRST_MAJORcount = 1;
                o = prepareFIRST_MAJOR();
                o.s.major_melt_count = 99;
                State::infDays = o.s.major_melt_count;
                State::major_melt_threshold = -5.0;
                o.s.init_SWE = 50.0;
                o.swe = []() {return 75.0;};
                cont = false;
                o.s.daily_max_temp = 100.0;
                break;
        };

        ++RESETcount;
        o.cont = cont;
        return o;                
    };
    
    output prepareNORMAL()
    {
        bool cont = true;
        output o;
        switch (NORMALcount)
        {
            case 0:
                RESETcount = 0;
                o = prepareRESET();
                State::major_melt_threshold = 5.0;
                break;
            case 1:
                RESETcount = 1;
                o = prepareRESET();
                State::major_melt_threshold = -5.0;
                o.s.init_SWE = 1e6;
                cont = false;
                break;
        }

        ++NORMALcount;
        o.cont = cont;
        return o;
    };

    template<LimitedPhase phase>
    output prepare()
    {

        output o;
        if constexpr(phase == LimitedPhase::NONE)
            o = prepareNONE();
        else if constexpr (phase == LimitedPhase::FIRST_MAJOR)
            o = prepareFIRST_MAJOR();
        else if constexpr (phase == LimitedPhase::PRIOR_INFILTRATION)
            o = preparePRIOR_INFILTRATION();
        else if constexpr (phase == LimitedPhase::RESET)
            o = prepareRESET();
        else if constexpr (phase == LimitedPhase::NORMAL)
            o = prepareNORMAL();

        return o;

    };
            
        
    std::string to_string(LimitedPhase phase)
    {
        switch (phase) {
            case LimitedPhase::NONE: 
                return "NONE";
            case LimitedPhase::FIRST_MAJOR:
                return "FIRST_MAJOR";
            case LimitedPhase::NORMAL:
                return "NORMAL";
            case LimitedPhase::RESET:
                return "RESET";
            case LimitedPhase::PRIOR_INFILTRATION:
                return "PRIOR_INFILTRATION";
        }
    }

    std::function<void()> build_phase_verify(const LimitedPhase& phase,const output& o)
    {
        auto phase_verify = [&phase,&o,this]() {
        EXPECT_EQ(o.phase,phase) << "Test of: " << to_string(phase);
        };
        return phase_verify;
    };
                                     
};

// ============================================================================
//                               phase_checker
// ============================================================================

class PhaseCheckerTest : public CrackTestBase
{
protected:
    template<LimitedPhase phase, typename ...Statements>
    void phase_test_loop(output& o,Statements... statements)
    {
        bool keep_going = true;
        phase_checker phase_check;
        while (keep_going)
        {
            o = prepare<phase>();
            o.phase = phase_check.get(o.s,o.swe());

            (statements(), ...);
            keep_going = o.cont;
        };
    };

    std::function<void()> build_check_count(const LimitedPhase& phase, const output& o)
    {
        auto check_count = [&o,&phase, this]() {
        EXPECT_GT(o.s.major_melt_count, State::infDays) << "Test of :" << to_string(phase);
        };

        return check_count;
    };

    phase_checker phase;
    output o;
    LimitedPhase tested_phase;

    template<LimitedPhase phase>
    void generic_phase_checker(const LimitedPhase& tested_phase, output& o)
    {
        auto phase_verify = build_phase_verify(tested_phase,o);

        phase_test_loop<phase>(o,phase_verify);
    };

};

TEST_F(PhaseCheckerTest,TestNONE)
{
    tested_phase = LimitedPhase::NONE;
    auto phase_verify = build_phase_verify(tested_phase,o);
    auto check_count = build_check_count(tested_phase,o);

    phase_test_loop<LimitedPhase::NONE>(o,phase_verify,check_count);
};

TEST_F(PhaseCheckerTest,TestFIRST_MAJOR)
{
    tested_phase = LimitedPhase::FIRST_MAJOR;

    generic_phase_checker<LimitedPhase::FIRST_MAJOR>(tested_phase,o);
};

TEST_F(PhaseCheckerTest,TestNORMAL)
{
    tested_phase = LimitedPhase::NORMAL;
    
    generic_phase_checker<LimitedPhase::NORMAL>(tested_phase,o);
};

TEST_F(PhaseCheckerTest,TestRESET)
{
    tested_phase = LimitedPhase::RESET;

    generic_phase_checker<LimitedPhase::RESET>(tested_phase,o);
};

TEST_F(PhaseCheckerTest,TestPRIOR_INFILTRATION)
{
    tested_phase = LimitedPhase::PRIOR_INFILTRATION;

    generic_phase_checker<LimitedPhase::PRIOR_INFILTRATION>(tested_phase,o);
};

// ============================================================================
//                                limited_inf
// ============================================================================

class LimitedInfTest : public CrackTestBase
{
protected:
    void SetUp() override
    {
        test_assertions.clear();
    };
    std::vector<std::function<void(double,output)>> test_assertions;
    
    std::function<double()> swe = []() {return 100.0;};
    std::function<void(State&,const double&)> init_inputs = [](State& s,const double& melt) {
        s.daily_melt_total.bind_target(melt);
        s.daily_melt_total.accumulate(false);
        s.daily_melt_total.accumulate(true);
    };

    output o;

    limited_inf inf_calc;

    template<LimitedPhase phase, std::invocable<double,output> ...Statements>
    void limited_test_loop(output& o,std::function<double()> swe, Statements... statements)
    {
        bool keep_going = true;
        limited_inf limited;
        while (keep_going)
        {
            o = prepare<phase>();

            auto melt = 0.0;

            if constexpr (phase == LimitedPhase::FIRST_MAJOR
                    || phase == LimitedPhase::RESET)
                melt = std::fabs(State::major_melt_threshold) * 2.0;
            else if constexpr (phase == LimitedPhase::NORMAL)
                if (State::major_melt_threshold > 0.0)
                    melt = State::major_melt_threshold / 2.0;
                else
                    melt = std::fabs(State::major_melt_threshold) * 2.0;
            else if constexpr (phase == LimitedPhase::PRIOR_INFILTRATION)
                melt = State::major_melt_threshold - fabs(State::major_melt_threshold)*0.5;
            else if constexpr (phase == LimitedPhase::NONE)
                melt = 1e6;

            init_inputs(o.s,melt);
            auto inf = limited.get(o.s,swe()); 

            (statements(inf,o), ...);
            keep_going = o.cont;
        };
    };

    
};

TEST_F(LimitedInfTest,FIRST_MAJOR_AND_RESET)
{
    test_assertions.push_back([this](double inf,output o) { 
            EXPECT_GT(o.s.index / swe(), 0.0) << "Test of: " << to_string(o.phase);
            });
    test_assertions.push_back([this](double inf,output o) { 
            EXPECT_GT(o.s.max_major_per_melt, 0.0) << "Test of: " << to_string(o.phase);
            EXPECT_LT(o.s.max_major_per_melt, swe()) << "Test of: " << to_string(o.phase);
            });
    test_assertions.push_back([this](double inf,output o) { 
            EXPECT_EQ(o.s.init_SWE,swe()) << "Test of: " << to_string(o.phase); 
            });
    test_assertions.push_back([this](double inf,output o) { 
            EXPECT_EQ(o.s.major_melt_count,1) << "Test of: " << to_string(o.phase);
            });

    test_assertions.push_back( [this](double inf, output o) {
            EXPECT_GT(inf,0.0) << "Test of: " << to_string(o.phase);
            } );

    limited_test_loop<LimitedPhase::FIRST_MAJOR>(o,swe,
            test_assertions.at(0),test_assertions.at(1),test_assertions.at(2),test_assertions.at(3),test_assertions.at(4));

    test_assertions.push_back( [this](double inf, output o) {
            EXPECT_GT(o.s.major_melt_count,1) << "Test of: " << to_string(o.phase);
            });
    // RESET is the same as FIRST_MAJOR, so change nothing and check again
    limited_test_loop<LimitedPhase::RESET>(o,swe,
            test_assertions.at(0),test_assertions.at(1),test_assertions.at(2),test_assertions.at(4));
};

TEST_F(LimitedInfTest,NORMAL)
{
    struct Holder
    {
        size_t melt_count = 0;
        double init_SWE = 0.0;
        double index = 0.0;
        double max_major = 0.0;
    } holder;

    test_assertions.push_back( [this,&holder](double inf, output o) {
            EXPECT_EQ(o.s.init_SWE,holder.init_SWE) << "Test of: " << to_string(o.phase);
            EXPECT_GT(inf,0.0) << "Test of: " << to_string(o.phase);
            EXPECT_GT(o.s.max_major_per_melt,0.0) << "Test of: " << to_string(o.phase);
            EXPECT_EQ(o.s.major_melt_count,holder.melt_count+1) << "Test of: " << to_string(o.phase);
            EXPECT_EQ(o.s.index,holder.index) << "Test of: " << to_string(o.phase);
            EXPECT_EQ(o.s.max_major_per_melt,holder.max_major) << "Test of: " << to_string(o.phase);
            });

    init_inputs = [&holder](State& s,const double melt) {
        s.daily_melt_total.bind_target(melt);
        s.daily_melt_total.accumulate(false);
        s.daily_melt_total.accumulate(true);

        holder = Holder{s.major_melt_count,s.init_SWE,
            s.index,s.max_major_per_melt};
    };

    limited_test_loop<LimitedPhase::NORMAL>(o,swe);

};

TEST_F(LimitedInfTest,PRIOR_INFILTRATION)
{
    test_assertions.clear();

    struct Holder
    {
        size_t melt_count = 0;
        double init_SWE = 0.0;
        double index = 0.0;
        double max_major = 0.0;
        double melt = 0.0;
    } holder;

    init_inputs = [&holder](State& s,const double melt) {
        s.daily_melt_total.bind_target(melt);
        s.daily_melt_total.accumulate(false);
        s.daily_melt_total.accumulate(true);

        holder = Holder{s.major_melt_count,s.init_SWE,
            s.index,s.max_major_per_melt,melt};
    };
    test_assertions.push_back( [&holder](double inf, output o) {
            EXPECT_EQ(o.s.major_melt_count,0);
            EXPECT_EQ(holder.melt,inf);
            });

    limited_test_loop<LimitedPhase::PRIOR_INFILTRATION>(o,swe,test_assertions[0]);
};

TEST_F(LimitedInfTest,NONE)
{
    test_assertions.clear();

    test_assertions.push_back( [](double inf, output o) {
            EXPECT_EQ(inf,0.0);
            });

    limited_test_loop<LimitedPhase::NONE>(o,swe,test_assertions[0]);
};

struct test_State : public State
{
    std::pair<bool,std::string> operator==(const test_State& other) const {
        bool output = true;
        std::string str;
        if (major_melt_count != other.major_melt_count)
            { output = false; str = "major_melt_count";}
        else if (index != other.index)
            { output = false; str = "index"; }
        else if (max_major_per_melt != other.max_major_per_melt)
            { output = false; str = "max_major_per_melt"; }
        else if (init_SWE != other.init_SWE)
            { output = false; str = "init_SWE"; }
        else if (daily_melt_total != other.daily_melt_total)
            { output = false; str = "daily_melt_total"; }
        else if (daily_rain_total != other.daily_rain_total)
            { output = false; str = "daily_rain_total"; }
        else if (daily_max_temp != other.daily_max_temp)
            { output = false; str = "daily_max_temp"; }
        else if (current.inf != other.current.inf)
            { output = false; str = "current.inf"; }
        else if (current.snow_inf != other.current.snow_inf)
            { output = false; str = "current.snow_inf"; }
        else if (current.runoff != other.current.runoff)
            { output = false; str = "current.runoff"; }
        else if (current.melt_runoff != other.current.melt_runoff)
            { output = false; str = "current.melt_runoff"; }
        else if (soil_saturation_at_freeze != other.soil_saturation_at_freeze)
            { output = false; str = "soil_saturation_at_freeze"; }

        return std::make_pair(output,str);
    };
};


//TEST_F(LimitedInfTest,OTHER)
//{
//    //test_assertions.push_back(
//    // Same outcome expected!
//    test_assertions.push_back([o_old = o](double inf, output o) { 
//            auto old = static_cast<test_State>(o_old.s);
//            auto current = static_cast<test_State>(o.s);
//            auto [eq,fault] = old == current;
//            EXPECT_TRUE(eq) << "Member at Fault: " << fault;
//            });
//    test_assertions.push_back([](double inf, output o) { EXPECT_DOUBLE_EQ(o.s.daily_melt_total.get_yesterday() * o.s.index,inf); });
//
//    limited_test_loop<LimitedPhase::RESET>(o,swe,test_assertions.at(0),test_assertions.at(1));
//    
//
//    // Again same outcome!
//    // 
//    test_assertions.push_back([](double inf, output o) { EXPECT_DOUBLE_EQ(o.s.max_major_per_melt,inf); });
//
//    limited_test_loop<LimitedPhase::NORMAL>(o,swe,test_assertions.at(0),test_assertions.at(1),test_assertions.at(2));
//
//    test_assertions.clear();
//
//    {
//        State s;
//        const double melt = 6.5;
//        s.daily_melt_total.bind_target(melt);
//        s.daily_melt_total.accumulate(false);
//        s.daily_melt_total.accumulate(true);
//        auto melt_check = [&s](double inf) { EXPECT_DOUBLE_EQ(s.daily_melt_total.get_yesterday(),inf);};
//        
//        output o = prepare<LimitedPhase::PRIOR_INFILTRATION>();
//        o.s = s;
//        auto inf = inf_calc.get(o.s,o.swe());
//
//        melt_check(inf);
//    }
//
//    test_assertions.push_back([](double inf,output o) { EXPECT_DOUBLE_EQ(inf,0.0); });
//
//    limited_test_loop<LimitedPhase::NONE>(o,swe,test_assertions[0]);
//
//};

// ============================================================================
//                               inf_calculator
// ============================================================================


class InfCalculator : public ::testing::Test
{
    // No test of the limited, public member function because its functionality 
    // is simply to call public members of other classes, its just a forwarding 
    // function and does not produce unique behaviour.
    // 
    // Testing that the function is called is enforcing a specific type of behaviour 
    // and it means that any changes to its private class would immediately break the test

protected:
    inf_calculator calc;
    State s;
    double melt = 2.2;
    double rain = 2.7;
    void SetUp() override {
        s = State::construct_test_state(melt,rain);
    };
};

TEST_F(InfCalculator,unlimited)
{
    auto inf = calc.unlimited(s);
    
    // No restriction on infiltration, soil is empty
    EXPECT_EQ(inf,melt);
    EXPECT_EQ(s.major_melt_count,1);
};

TEST_F(InfCalculator,restricted)
{
    auto inf = calc.restricted(s);

    // Fully stops infiltration, soil is full.
    EXPECT_EQ(inf,0.0);
    EXPECT_EQ(s.major_melt_count,0);
};



// ============================================================================
//                                 Processor
// ============================================================================



class ProcessorTest : public ::testing::Test
{
    static constexpr auto swe = 100.0;
protected:
    const std::vector<double> soil_saturation{0.0,45.0,100.0,150.99};
    double melt;
    double rain;

    State melt_no_rain() {
        melt = 1.6;
        rain = 0.0;

        return State::construct_test_state(melt, rain);
    };

    State melt_with_rain() {
        melt = 2.5;
        rain = 3.2;      

        return State::construct_test_state(melt, rain);
    };

    State rain_no_melt() { 
        melt = 0.0;
        rain = 1.3;

        return State::construct_test_state(melt, rain);
    };

    State no_rain_no_melt() { 
        melt = 0.0;
        rain = 0.0;

        return State::construct_test_state(melt,rain);
    }; 
    
};

TEST_F(ProcessorTest,NoMeltNoRainYesterdayProducesZeros)
{
    Processor process;
    auto s = no_rain_no_melt();

    process.new_day(s,100.0);

    EXPECT_DOUBLE_EQ(s.current.inf,0.0);
    EXPECT_DOUBLE_EQ(s.current.runoff,0.0);
    EXPECT_DOUBLE_EQ(s.current.snow_inf,0.0);
    EXPECT_DOUBLE_EQ(s.current.melt_runoff,0.0);
    EXPECT_DOUBLE_EQ(s.current.rain,0.0);
};

TEST_F(ProcessorTest, LimitedBranch)
{
    Processor process;
    inf_calculator inf;

    auto s = melt_with_rain();
    s.soil_saturation_at_freeze = 43.0;

    double swe = 55.0;

    // Force that it is not a has_lens case
    s.daily_max_temp = std::fabs(State::lenstemp) * 10;
    double expected_inf = inf.limited(s,swe);

    EXPECT_GT(expected_inf,0.0);
    
    process.new_day(s,swe);

    EXPECT_DOUBLE_EQ(s.current.snow_inf, expected_inf);
    EXPECT_DOUBLE_EQ(s.current.melt_runoff, melt - expected_inf);

    EXPECT_DOUBLE_EQ(s.current.runoff,s.current.melt_runoff);
    EXPECT_DOUBLE_EQ(s.current.inf,expected_inf + rain);

    // Pretend to force has_lens, however this will give the same answer as above
    // Requires small max temperature AND nonzero major_melt_count
    s.daily_max_temp = State::lenstemp - 5.0;

    expected_inf = inf.limited(s,swe);
    EXPECT_GT(expected_inf,0.0);

    process.new_day(s,swe);
    EXPECT_DOUBLE_EQ(s.current.snow_inf, expected_inf);
    EXPECT_DOUBLE_EQ(s.current.melt_runoff, melt - expected_inf);
    
    EXPECT_DOUBLE_EQ(s.current.runoff,s.current.melt_runoff);
    EXPECT_DOUBLE_EQ(s.current.inf,expected_inf+rain);

    // ACtually triggers lens case and restricts all infiltration

    s.major_melt_count = 1000;
    expected_inf = inf.limited(s,swe);
    EXPECT_DOUBLE_EQ(expected_inf,0.0);

    process.new_day(s,swe);
    EXPECT_DOUBLE_EQ(s.current.snow_inf, 0.0);
    EXPECT_DOUBLE_EQ(s.current.melt_runoff, melt);
    
    EXPECT_DOUBLE_EQ(s.current.runoff,melt+rain);
    EXPECT_DOUBLE_EQ(s.current.inf,0.0);
};

TEST_F(ProcessorTest,UnlimitedBranch)
{
    Processor process;
    inf_calculator inf;

    auto s = melt_with_rain();
    s.soil_saturation_at_freeze = 0.0;
    auto swe = 0.0;

    auto expected_inf = inf.unlimited(s);
    
    process.new_day(s,swe);

    EXPECT_DOUBLE_EQ(s.current.melt_runoff,0.0);
    EXPECT_DOUBLE_EQ(s.current.snow_inf,melt);
    EXPECT_DOUBLE_EQ(s.current.runoff,0.0);
    EXPECT_DOUBLE_EQ(s.current.inf,melt+rain);
    EXPECT_DOUBLE_EQ(s.current.rain,rain);

};

TEST_F(ProcessorTest,RestrictedBranch)
{
    Processor process;
    inf_calculator inf;

    auto s = melt_with_rain();
    s.soil_saturation_at_freeze = 100.0;
    auto swe = 0.0;

    auto expected_inf = inf.restricted(s);
    
    process.new_day(s,swe);

    EXPECT_DOUBLE_EQ(s.current.melt_runoff,melt);
    EXPECT_DOUBLE_EQ(s.current.snow_inf,0.0);
    EXPECT_DOUBLE_EQ(s.current.runoff,melt+rain);
    EXPECT_DOUBLE_EQ(s.current.inf,0.0);
    EXPECT_DOUBLE_EQ(s.current.rain,rain);

};

TEST_F(ProcessorTest,ThrowOnBadSoilSaturation)
{
    Processor process;
    auto s = melt_with_rain();
    auto swe = 0.0;
    s.soil_saturation_at_freeze = 1e6;

    EXPECT_THROW(process.new_day(s,swe),std::logic_error);

    s.soil_saturation_at_freeze = -22.0;

    EXPECT_THROW(process.new_day(s,swe),std::logic_error);
};

TEST_F(ProcessorTest,NoMeltWithRainHoldsRain)
{
    Processor process;
    auto s = rain_no_melt();
    double swe = 100.0;
    process.new_day(s,swe);

    EXPECT_DOUBLE_EQ(s.current.rain,rain);

    process.new_day(s,swe);

    EXPECT_DOUBLE_EQ(s.current.rain,rain * 2.0);

    process.new_day(s,swe);

    EXPECT_DOUBLE_EQ(s.current.rain,rain * 3.0);

    melt = 5.5;
    s.daily_melt_total.accumulate(false);
    s.daily_melt_total.accumulate(true);
    s.soil_saturation_at_freeze = 0.0;
    process.new_day(s,swe);

    EXPECT_DOUBLE_EQ(s.current.rain,rain);
    EXPECT_GT(s.current.inf,0.0);
};



// ============================================================================
//                                   Model
// ============================================================================


class Data
{
public:
    // ---- Required getters ----
    Crack::State& get_state() { return state_; }
    const bool is_newday() const { return is_newday_; }
    double snowmelt() const { return snowmelt_; }
    double rainfall() const { return rainfall_; }
    double swe() const { return swe_; }
    double air_temperature() const { return air_temperature_; }

    // ---- Required setters (mutators) ----
    void infiltrated(double v) { _inf = v; }
    void runoff(double v) { _runoff = v; }
    void snow_infiltrated(double v) { _snow_inf = v; }
    void melt_runoff(double v) { _melt_runoff = v; }
    void rain_on_snow(double v) { _rain = v; }

    // ---- Additional setters for the getters ----
    void set_is_newday(bool v) { is_newday_ = v; }
    void set_snowmelt(double v) { snowmelt_ = v; }
    void set_rainfall(double v) { rainfall_ = v; }
    void set_swe(double v) { swe_ = v; }
    void set_air_temperature(double v) { air_temperature_ = v; }

    // ---- Additional getters for the setters ----
    double get_infiltrated() const { return _inf; }
    double get_runoff() const { return _runoff; }
    double get_snow_infiltrated() const { return _snow_inf; }
    double get_melt_runoff() const { return _melt_runoff; }
    double get_rain_on_snow() const { return _rain; }

private:
    // Private members for all fields
    Crack::State state_{};

    bool is_newday_ = false;
    double snowmelt_ = 0.0;
    double rainfall_ = 0.0;
    double swe_ = 0.0;
    double air_temperature_ = 0.0;
    double _inf = 0.0;
    double _runoff = 0.0;
    double _snow_inf = 0.0;
    double _melt_runoff = 0.0;
    double _rain = 0.0;

};

class CrackModelTest : public ::testing::Test
{
protected:
    Model<Data> crack;
    static constexpr double melt = 7.5;
    static constexpr double rain = 2.0;
    Data d;

    void SetUp() override {
        auto& s = d.get_state();    
        s.daily_melt_total.bind_target(melt);
        s.daily_rain_total.bind_target(rain);
    };
    State::Current c{-9999.0,-1000.0,-50022.0,-55.0,-1.0};

    void current_equal_data(Data& d,State::Current& c,const std::string& str) {
        EXPECT_DOUBLE_EQ(d.get_rain_on_snow(),c.rain) << str;
        EXPECT_DOUBLE_EQ(d.get_runoff(),c.runoff) << str;
        EXPECT_DOUBLE_EQ(d.get_melt_runoff(),c.melt_runoff) << str;
        EXPECT_DOUBLE_EQ(d.get_infiltrated(),c.inf) << str;
        EXPECT_DOUBLE_EQ(d.get_snow_infiltrated(),c.snow_inf) << str;
    };
};

TEST_F(CrackModelTest,SameDay)
{
    d.set_is_newday(false);
    auto& s = d.get_state();
    auto t = 5.0; 
    s.daily_max_temp = 10.0;
    s.current = c;

    crack.execute(d);

    EXPECT_DOUBLE_EQ(s.current.rain,c.rain);
    EXPECT_DOUBLE_EQ(s.current.inf,c.inf);
    EXPECT_DOUBLE_EQ(s.current.snow_inf,c.snow_inf);
    EXPECT_DOUBLE_EQ(s.current.runoff,c.runoff);
    EXPECT_DOUBLE_EQ(s.current.melt_runoff,c.melt_runoff);

    current_equal_data(d,s.current,"Same Day");

    EXPECT_GT(s.daily_max_temp,t);
};

TEST_F(CrackModelTest,NewDay)
{
    d.set_is_newday(true);
    auto& s = d.get_state();

    auto t = 100.0;
    d.set_air_temperature(t);

    crack.execute(d);

    EXPECT_NE(s.current.rain,c.rain);
    EXPECT_NE(s.current.inf,c.inf);
    EXPECT_NE(s.current.snow_inf,c.snow_inf);
    EXPECT_NE(s.current.runoff,c.runoff);
    EXPECT_NE(s.current.melt_runoff,c.melt_runoff);

    current_equal_data(d,s.current,"New Day");

    EXPECT_DOUBLE_EQ(s.daily_max_temp,t);
};
