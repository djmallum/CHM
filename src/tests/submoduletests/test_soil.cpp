#include <gmock/gmock.h>
#include "soil_DTO.hpp"
#include "soil_classes.hpp"
#include <gtest/gtest.h>
#define diff 0.0001


class SoilComponentsTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Common setup code for all tests
        reset_DTO();
	    initializer init(DTO);    
		init.zero_single_step_vars();
    }
    
    void reset_DTO()
    {
        DTO.soil_storage_max = 750.0;
        DTO.soil_rechr_max = 250.0;
        DTO.soil_storage = 350.0;
        DTO.soil_rechr_storage = 125.0;
        DTO.excess_to_ssr = true;
        DTO.detention_snow_max = 0.0;
        DTO.detention_organic_max = 0.0;
        DTO.depression_max = 1.0;
        DTO.ground_water_max = 2.0;
        DTO.local_slope = 7.0 * 3.14159265/180.0;
        DTO.swe = 25.0;
        DTO.infil = 16.3497;
        
        DTO.pore_size_dist = 0.333;
        DTO.pore_size_dist_organic = 0.444;
        //const std::string soil_type = 
        DTO.porosity = 0.5;
        DTO.soil_index = 3.3;
        DTO.snow_grain_diameter = 3;
                                         
        DTO.Ksaturated_rechr = 0.0176;
        DTO.Ksaturated_lower = 6.9e-6;
        DTO.Ksaturated_ground_water = 6.9e-6;
        DTO.Ksaturated_organic = 6.9e-6;
        // Note: Ksaturated_snow is comput
        
        DTO.allow_runoff_from_infiltration = true;

        DTO.freeze_thaw_first_front = 0.0;
    };
    ::testing::NiceMock<MockTwoLayerDTO> DTO;
    void infil_setup()
    {
        DTO.thaw_fraction_rechr = 1.0;
        DTO.thaw_fraction_lower = 1.0;
    };

    void setup_depression(const double init_depression_storage,
		const double init_soil_excess);
    void test_depression(const double expected_excess,
		const double expected_storage,
		const double expected_runoff_to_storage,std::string errormsg);
};

TEST_F(SoilComponentsTest, SetLayerThawFractionTest) {
    DTO.set_new_day_true();
	DTO.set_dt_3600();
    
    initializer init(DTO);

    // Zero Porosity edge case
    DTO.porosity = 0.0;
    
    init.layer_thaw_fraction();
    
    EXPECT_DOUBLE_EQ(DTO.thaw_fraction_rechr, 0.0);
    EXPECT_DOUBLE_EQ(DTO.thaw_fraction_lower, 0.0);

    // Test allow_runoff_from_infiltration = true
    reset_DTO();
    DTO.allow_runoff_from_infiltration = false;
    init.layer_thaw_fraction();
    EXPECT_DOUBLE_EQ(DTO.thaw_fraction_rechr, 1.0);
    EXPECT_DOUBLE_EQ(DTO.thaw_fraction_lower, 1.0);

    // Test freeze_thaw_first_front = 0.0
    reset_DTO();
    DTO.allow_runoff_from_infiltration = true;
    DTO.freeze_thaw_first_front = 10.0;
    init.layer_thaw_fraction();
    EXPECT_DOUBLE_EQ(DTO.thaw_fraction_rechr, 0.0);
    EXPECT_DOUBLE_EQ(DTO.thaw_fraction_lower, 0.0);
    
    //reset
    // thaw front below all layers
    reset_DTO();
    DTO.thaw_front_depth = 10.0;  // Deeper than both layers
    init.layer_thaw_fraction(); 
    EXPECT_DOUBLE_EQ(DTO.thaw_fraction_rechr, 1.0);
    EXPECT_DOUBLE_EQ(DTO.thaw_fraction_lower, 1.0);

    // Thaw front in recharge layer
    reset_DTO();
    double rechr_depth = DTO.soil_rechr_max / DTO.porosity / 1000.0;
    DTO.thaw_front_depth = rechr_depth * 0.5;  // Halfway through recharge layer
    DTO.freeze_thaw_first_front = DTO.thaw_front_depth;
    
    init.layer_thaw_fraction(); 
    
    double expected_rechr = 0.5;
    EXPECT_TRUE(DTO.thaw_front_depth < rechr_depth); 
    EXPECT_DOUBLE_EQ(DTO.thaw_fraction_rechr, expected_rechr);
    EXPECT_DOUBLE_EQ(DTO.thaw_fraction_lower, 0.0) << "Expected to fail, possible Bug in CHRM.\n IGNORE";

    // Thaw front in lower layer
    reset_DTO();
    rechr_depth = DTO.soil_rechr_max / DTO.porosity / 1000.0;
    double soil_depth = DTO.soil_storage_max / DTO.porosity / 1000.0;
    DTO.thaw_front_depth = rechr_depth + (soil_depth - rechr_depth) * 0.5;  // Halfway through lower layer
    DTO.freeze_thaw_first_front = DTO.thaw_front_depth;
    
    init.layer_thaw_fraction(); 

    EXPECT_DOUBLE_EQ(DTO.thaw_fraction_rechr, 1.0);
    EXPECT_DOUBLE_EQ(DTO.thaw_fraction_lower, 0.5);

    // Edge case soil_storage max is 0
    reset_DTO();
DTO.soil_storage_max = 0.0;
    
    init.layer_thaw_fraction();

    EXPECT_DOUBLE_EQ(DTO.thaw_fraction_rechr, 0.0);
    EXPECT_DOUBLE_EQ(DTO.thaw_fraction_lower, 0.0);

    // off day
    reset_DTO(); 
    ON_CALL(DTO,get_new_day).WillByDefault(::testing::Return(false));
    double rechr =0.122334;
    double lower = 0.22454;
    DTO.thaw_fraction_rechr = rechr;
    DTO.thaw_fraction_lower = lower;

    init.layer_thaw_fraction();

    EXPECT_EQ(DTO.thaw_fraction_rechr,rechr);
    EXPECT_EQ(DTO.thaw_fraction_lower,lower);
} 

TEST_F(SoilComponentsTest, SetCondensationTest) {
    // nonzero SWE
    DTO.swe = 50.0;
    condensator C(DTO);
    double input = 153.4;
    DTO.actual_ET = input;

    C.set();

    EXPECT_EQ(DTO.actual_ET,input);
    EXPECT_EQ(DTO.condensation,0.0);

    // positive ET

    reset_DTO();
    DTO.swe = 0.0;
    DTO.actual_ET = input;

    C.set();

    EXPECT_EQ(DTO.actual_ET,input);
    EXPECT_EQ(DTO.condensation,0.0);

    // Negative ET
    
    reset_DTO();
    DTO.swe = 0.0;
    DTO.actual_ET = -input;

    C.set();

    EXPECT_EQ(DTO.actual_ET,0.0);
    EXPECT_EQ(DTO.condensation,input);
    
}



TEST_F(SoilComponentsTest, DistributeInfiltratedAmounts) {
    // Edge case soil_storage_max = 0
    DTO.soil_storage_max = 0.0;
    ::testing::NiceMock<MockTwoLayerDTO> copy_DTO;
    copy_DTO.soil_rechr_storage = DTO.soil_rechr_storage;
    copy_DTO.soil_storage = DTO.soil_storage;
    initializer init(DTO);

    infiltrator I(DTO);
    I.distribute();

    EXPECT_EQ(DTO.soil_storage,copy_DTO.soil_storage);
    EXPECT_EQ(DTO.soil_rechr_storage,copy_DTO.soil_rechr_storage);
    EXPECT_EQ(DTO.soil_excess_to_runoff,0.0);
    EXPECT_EQ(DTO.rechr_to_ssr,0.0);
    EXPECT_EQ(DTO.soil_to_ssr,0.0);
    EXPECT_EQ(DTO.excess,DTO.infil);

    reset_DTO();

    // allow_runoff_from_infiltration = true
    infil_setup();
	DTO.allow_runoff_from_infiltration = true;
    double potential = DTO.infil + DTO.condensation;
    double possible = potential * 0.9;
    DTO.soil_rechr_storage = DTO.soil_rechr_max - possible;
    init.zero_single_step_vars(); 
    I.distribute(); 

    EXPECT_NEAR(DTO.soil_excess_to_runoff,potential-possible,1e-4);
    EXPECT_EQ(DTO.soil_rechr_storage,DTO.soil_rechr_max);
    

	// possible > potential
    
    infil_setup();
    possible = potential * 1.2;
    DTO.soil_rechr_storage = DTO.soil_rechr_max - possible;
    DTO.allow_runoff_from_infiltration = true;
    init.zero_single_step_vars();
    I.distribute();

    EXPECT_EQ(DTO.soil_excess_to_runoff,0.0);
    EXPECT_EQ(DTO.soil_rechr_storage,DTO.soil_rechr_max - possible + potential);

	
	// possible < potential
    
    infil_setup();
    DTO.soil_rechr_storage = DTO.soil_rechr_max - 0.1;
    DTO.soil_storage = DTO.soil_storage_max - 0.2;
    DTO.K_soil_to_gw = 1e-4;
    DTO.allow_runoff_from_infiltration = false;
    DTO.excess_to_ssr = true;
    init.zero_single_step_vars();
    I.distribute();

    
    EXPECT_EQ(DTO.soil_excess_to_runoff,0.0);
    EXPECT_EQ(DTO.soil_rechr_storage,DTO.soil_rechr_max);
    EXPECT_EQ(DTO.soil_storage,DTO.soil_storage_max);
    EXPECT_EQ(DTO.soil_excess_to_gw,DTO.K_soil_to_gw);
    const double expected_excess = DTO.infil + DTO.condensation - 0.2-DTO.K_soil_to_gw; 
    EXPECT_NEAR(DTO.soil_to_ssr,expected_excess,1e-4);
	
    
    // nonzero swe
    infil_setup();
    DTO.soil_rechr_storage = DTO.soil_rechr_max - 0.1;
    DTO.soil_storage = DTO.soil_storage_max - 0.2;
    DTO.swe = 0.0;
    DTO.K_soil_to_gw = 1e-4;
    DTO.allow_runoff_from_infiltration = false;
    DTO.excess_to_ssr = false;
    DTO.K_rechr_to_ssr = 2e-2;
    init.zero_single_step_vars();
    I.distribute(); 

    EXPECT_NEAR(DTO.rechr_to_ssr,DTO.K_rechr_to_ssr * DTO.thaw_fraction_rechr,1e-4);
    EXPECT_NEAR(DTO.soil_rechr_storage,DTO.soil_rechr_max - DTO.rechr_to_ssr,1e-4);
    EXPECT_NEAR(DTO.soil_storage,DTO.soil_storage_max-DTO.rechr_to_ssr,1e-4);
    EXPECT_NEAR(DTO.soil_to_ssr,DTO.rechr_to_ssr,1e-4);

	// Test organize_soil_layers()
    // Add assertions here
}

TEST_F(SoilComponentsTest, ManageDetentionTest) {
    detention_layer detention(DTO);
    initializer init(DTO);
    
    DTO.runoff  = 5.0;
    DTO.excess = 20.0;
    DTO.swe = 0.0;
    DTO.detention_snow_max = 10.0;
    DTO.detention_organic_max = 20.0;
    DTO.detention_storage = 0.0;
    DTO.K_detention_to_runoff = 2.3e-3;
    detention.manage();
    
    EXPECT_EQ(DTO.detention_max,DTO.detention_organic_max);
    EXPECT_EQ(DTO.soil_excess_to_runoff,DTO.runoff+DTO.excess - DTO.detention_max+DTO.K_detention_to_runoff);
    EXPECT_EQ(DTO.detention_storage,DTO.detention_max-DTO.K_detention_to_runoff);
   
	// with snow
	
	init.zero_single_step_vars();
	DTO.detention_storage = 0.0;
	DTO.swe = 22.0;
	detention.manage();

	EXPECT_EQ(DTO.detention_max,DTO.detention_snow_max);
    EXPECT_EQ(DTO.soil_excess_to_runoff,DTO.runoff+DTO.excess - DTO.detention_max+DTO.K_detention_to_runoff);
    EXPECT_EQ(DTO.detention_storage,DTO.detention_max-DTO.K_detention_to_runoff);

}

void SoilComponentsTest::setup_depression(const double init_depression_storage,
		const double init_soil_excess)
{
	// depression_space < soil_excess_to_runoff
	// soil_excess_to_runoff > depression_max
	DTO.soil_excess_to_runoff = init_soil_excess;
	DTO.depression_storage = init_depression_storage;
	DTO.runoff_to_depression = 0.0;
};

void SoilComponentsTest::test_depression(const double expected_excess,
		const double expected_storage,
		const double expected_runoff_to_storage,std::string errormsg)
{
	EXPECT_EQ(DTO.soil_excess_to_runoff,expected_excess) << errormsg;
	EXPECT_EQ(DTO.depression_storage,expected_storage) << errormsg;
	EXPECT_EQ(DTO.runoff_to_depression,expected_runoff_to_storage) << errormsg;
};


TEST_F(SoilComponentsTest, ManageDepressionTest) {
	depression_layer depression(DTO);
	DTO.depression_max = 3.0;
	DTO.K_depression_to_gw = 0.0;
	struct initial_values
	{
		double init_Sd;
		double init_excess;
	};
	initial_values I;
	struct expected_values
	{
		double expected_excess;
		double expected_storage;
        double expected_runoff_to_storage;
	};
	expected_values E;
	std::string errormsg;
	
	errormsg = "soil_excess_to_runoff = 0.0";
	I.init_Sd = 1.234;
	I.init_excess = 0.0;
	E.expected_excess = I.init_excess;
	E.expected_storage = I.init_Sd;
	E.expected_runoff_to_storage = 0.0;
	setup_depression(I.init_Sd,I.init_excess);
	depression.manage();
	EXPECT_EQ(DTO.soil_excess_to_runoff,E.expected_excess) << errormsg;
	EXPECT_EQ(DTO.depression_storage,E.expected_storage) << errormsg;
	EXPECT_EQ(DTO.runoff_to_depression,E.expected_runoff_to_storage) << errormsg;

	errormsg = "depression_max = 0.0";
	DTO.depression_max = 0.0;
	I.init_Sd = 1.234;
	I.init_excess = 46.0;
	E.expected_excess = I.init_excess;
	E.expected_storage = I.init_Sd;
	E.expected_runoff_to_storage = 0.0;
	setup_depression(I.init_Sd,I.init_excess);
	depression.manage();
	EXPECT_EQ(DTO.soil_excess_to_runoff,E.expected_excess) << errormsg;
	EXPECT_EQ(DTO.depression_storage,E.expected_storage) << errormsg;
	EXPECT_EQ(DTO.runoff_to_depression,E.expected_runoff_to_storage) << errormsg;

	DTO.depression_max = 3.0;
	double depression_space;
	errormsg = "soil_excess_to_runoff / depression_max < 12.0";
	I.init_Sd = 1.234;
	I.init_excess = DTO.depression_max * 12.0 - 1.0;
	depression_space = (DTO.depression_max - I.init_Sd) * (1 - exp(-I.init_excess/DTO.depression_max));	
	E.expected_excess = I.init_excess - depression_space;
	E.expected_storage = I.init_Sd + depression_space;
	E.expected_runoff_to_storage = depression_space;
	setup_depression(I.init_Sd,I.init_excess);
	depression.manage();
	EXPECT_EQ(DTO.soil_excess_to_runoff,E.expected_excess) << errormsg;
	EXPECT_EQ(DTO.depression_storage,E.expected_storage) << errormsg;
	EXPECT_EQ(DTO.runoff_to_depression,E.expected_runoff_to_storage) << errormsg;

	// soil_excess_to_runoff / depression_max > 12.0
	errormsg = "soil_excess_to_runoff / depression_max < 12.0";
	I.init_Sd = 1.234;
	I.init_excess = DTO.depression_max * 12.0 + 1.0;
	depression_space = (DTO.depression_max - I.init_Sd) * (1 - exp(-12.0));	
	E.expected_excess = I.init_excess - depression_space;
	E.expected_storage = I.init_Sd + depression_space;
	E.expected_runoff_to_storage = depression_space;
	setup_depression(I.init_Sd,I.init_excess);
	depression.manage();
	EXPECT_EQ(DTO.soil_excess_to_runoff,E.expected_excess) << errormsg;
	EXPECT_EQ(DTO.depression_storage,E.expected_storage) << errormsg;
	EXPECT_EQ(DTO.runoff_to_depression,E.expected_runoff_to_storage) << errormsg;
	
	errormsg = "soil_excess_to_runoff < space";
	I.init_Sd = 1.234;
	depression_space = (DTO.depression_max - I.init_Sd) * (1 - exp(-12.0));	
	I.init_excess = depression_space * 0.8;
	E.expected_excess = 0.0;
	E.expected_storage = I.init_Sd + I.init_excess;
	E.expected_runoff_to_storage = I.init_excess;
	setup_depression(I.init_Sd,I.init_excess);
	depression.manage();
	EXPECT_EQ(DTO.soil_excess_to_runoff,E.expected_excess) << errormsg;
	EXPECT_EQ(DTO.depression_storage,E.expected_storage) << errormsg;
	EXPECT_EQ(DTO.runoff_to_depression,E.expected_runoff_to_storage) << errormsg;

	errormsg = "soil_excess_to_runoff = 0, but remaining depression storage";
	DTO.K_depression_to_gw = 1.2;
	I.init_Sd = 1.234;
	I.init_excess = 0.0;
	DTO.depression_to_gw = 0.0;
	E.expected_storage = I.init_Sd - DTO.K_depression_to_gw;
	double expected_Sd_to_gw = 
		DTO.K_depression_to_gw;
	setup_depression(I.init_Sd,0.0);
	depression.manage();
	EXPECT_EQ(DTO.soil_excess_to_runoff,E.expected_excess) << errormsg;
	EXPECT_EQ(DTO.depression_storage,E.expected_storage) << errormsg;
	EXPECT_EQ(DTO.depression_to_gw,expected_Sd_to_gw) << errormsg;
	
		
}

TEST_F(SoilComponentsTest, ManageGroundwaterTest) 
{
    initializer init(DTO);
	DTO.K_ground_water_out = 1.3e-2;
	DTO.ground_water_max = 321.234;
	double init_gw,init_from_depression,init_from_excess;
	
	std::string errormsg = "gw > gw_max > 0";
	init_gw = DTO.ground_water_max * 1.2;
	DTO.ground_water_storage = init_gw;
	init_from_depression = 123.0;
	init_from_excess = 148.9;
    DTO.soil_excess_to_gw = init_from_excess;
    DTO.depression_to_gw = init_from_depression;
    
	groundwater_layer groundwater(DTO);
	groundwater.manage();

	EXPECT_NEAR(DTO.ground_water_out,init_from_excess + init_from_depression + init_gw - DTO.ground_water_max + DTO.K_ground_water_out,diff) << errormsg;
	EXPECT_NEAR(DTO.ground_water_storage,DTO.ground_water_max - DTO.K_ground_water_out,diff) << errormsg;
	EXPECT_EQ(DTO.depression_to_gw,0.0) << errormsg;

	errormsg = "gw > gw_max = 0";
	init.zero_single_step_vars();
	DTO.ground_water_max = 0.0;
	init_from_depression = 99.3;
	init_gw = 323.4;
    DTO.ground_water_storage = init_gw;
	DTO.depression_to_gw = init_from_depression;
    DTO.soil_excess_to_gw = init_from_excess; 
	groundwater.manage();

	EXPECT_NEAR(DTO.ground_water_out,init_gw + init_from_excess + init_from_depression,diff) << errormsg;
	EXPECT_EQ(DTO.ground_water_storage,DTO.ground_water_max) << errormsg;
	EXPECT_EQ(DTO.depression_to_gw,0.0) << errormsg;

	errormsg = "gw < gw_max > 0";
	init.zero_single_step_vars();
	DTO.ground_water_max = 1050.0;
	init_gw = DTO.ground_water_max * 0.5;
	DTO.ground_water_storage = init_gw;
	DTO.depression_to_gw = 0.0; //eliminate its effect from the test
	groundwater.manage();

	EXPECT_NEAR(DTO.ground_water_out,(init_from_excess + init_from_depression + init_gw)/DTO.ground_water_max * DTO.K_ground_water_out,diff) << errormsg;
	EXPECT_NEAR(DTO.ground_water_storage,(init_from_excess + init_from_depression + init_gw) * (1 - DTO.K_ground_water_out/DTO.ground_water_max),diff) << errormsg;
}

TEST_F(SoilComponentsTest, ManageSubsurfaceRunoffTest) 
{
    initializer init(DTO);
	subsurface_runoff ssr(DTO);
	DTO.K_depression_to_ssr = 2.3e-3;
	double init_Sd;

    std::string errormsg = "Sd > K_Sd_to_ssr > 0.0";
	DTO.K_lower_to_ssr = 0.0; // zero so it is not computed
	init_Sd = DTO.K_depression_to_ssr * 1.1;
	DTO.depression_storage = init_Sd;
	ssr.manage();

	EXPECT_NEAR(DTO.depression_storage, init_Sd - DTO.K_depression_to_ssr,diff) << errormsg;
	EXPECT_NEAR(DTO.soil_to_ssr,DTO.K_depression_to_ssr,diff) << errormsg;
    
    errormsg = "K_Sd_to_ssr > Sd> 0.0";
	init.zero_single_step_vars();
    DTO.K_lower_to_ssr = 0.0; // zero so it is not computed
	init_Sd = DTO.K_depression_to_ssr * 0.9;
	DTO.depression_storage = init_Sd;
	ssr.manage();

	EXPECT_EQ(DTO.depression_storage,0.0) << errormsg;
	EXPECT_NEAR(DTO.soil_to_ssr,init_Sd,diff) << errormsg;

    errormsg = "K_Sd_to_ssr > Sd = 0.0";
	init.zero_single_step_vars();
	DTO.K_lower_to_ssr = 0.0; // zero so it is not computed
	init_Sd = 0.0;
	DTO.depression_storage = init_Sd;
	ssr.manage();

	EXPECT_EQ(DTO.depression_storage,0.0) << errormsg;
	EXPECT_EQ(DTO.soil_to_ssr,0.0) << errormsg;

    errormsg = "Sd > K_Sd_to_ssr = 0.0";
	init.zero_single_step_vars();
	DTO.K_lower_to_ssr = 0.0; // zero so it is not computed
	init_Sd = DTO.K_depression_to_ssr * 0.9;
	DTO.depression_storage = init_Sd;
	DTO.K_depression_to_ssr = 0.0;
	ssr.manage();

	EXPECT_EQ(DTO.depression_storage,0.0) << errormsg;
	EXPECT_NEAR(DTO.soil_to_ssr,init_Sd,diff) << errormsg;

    errormsg = "K_lower_to_ssr > 0 and < soil_lower_storage. thaw fraction is 1.0";
	init.zero_single_step_vars();
	double init_soil, init_rechr;
	init_rechr = 100.0;
	DTO.soil_rechr_storage = init_rechr;
	init_soil = init_rechr * 5.0;
	DTO.soil_storage = init_soil;
	DTO.thaw_fraction_lower = 1.0;
	DTO.K_lower_to_ssr = 1.994e-2; // zero so it is not computed
	DTO.K_depression_to_ssr = 0.0;
	ssr.manage();

	EXPECT_NEAR(DTO.soil_storage,init_soil - DTO.K_lower_to_ssr,diff*1e-2) << errormsg;
	EXPECT_NEAR(DTO.soil_to_ssr,DTO.K_lower_to_ssr,diff*1e-2) << errormsg;

    errormsg = "K_lower_to_ssr > 0 and > soil_lower_storage. thaw fraction is 1.0";
	init.zero_single_step_vars();
	init_rechr = 100.0;
	DTO.soil_rechr_storage = init_rechr;
	init_soil = init_rechr * 5.0;
	DTO.soil_storage = init_soil;
	DTO.thaw_fraction_lower = 1.0;
	DTO.K_lower_to_ssr = init_soil; // zero so it is not computed
	DTO.K_depression_to_ssr = 0.0;
	ssr.manage();

	EXPECT_NEAR(DTO.soil_storage,init_soil - (init_soil-init_rechr),diff*1e-2) << errormsg;
	EXPECT_NEAR(DTO.soil_to_ssr,init_soil-init_rechr,diff*1e-2) << errormsg;

    errormsg = "K_lower_to_ssr > 0 and < soil_lower_storage. thaw fraction is 0.4";
	init.zero_single_step_vars();
	init_rechr = 100.0;
	DTO.soil_rechr_storage = init_rechr;
	init_soil = init_rechr * 5.0;
	DTO.soil_storage = init_soil;
	DTO.thaw_fraction_lower = 0.4;
	DTO.K_lower_to_ssr = 4.6723e-1; // zero so it is not computed
	DTO.K_depression_to_ssr = 0.0;
	ssr.manage();

	EXPECT_NEAR(DTO.soil_storage,init_soil - (init_soil-init_rechr)*DTO.thaw_fraction_lower,diff*1e-2) << errormsg;
	EXPECT_NEAR(DTO.soil_to_ssr,(init_soil-init_rechr)*DTO.thaw_fraction_lower,diff*1e-2) << errormsg;
}

