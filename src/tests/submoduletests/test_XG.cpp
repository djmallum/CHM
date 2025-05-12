#include "XG_algorithm.hpp"
#include <vector>
#include <memory>
#include "gtest/gtest.h"

/*
 * XGStateTest: Wrapper class for tests
 * XGStateTest is effectively a mock of Infil_All module but done indirectly. Due to the complexity of the module classes, it was easier to write this.  
 * The member variables with the _ prefix are inputs that are supplied to the constructor of Crack.
 * Default values are given and used for most tests.
 * Other member variables are parameters that are also supplied to Crack unless it has the const specifier, then it is just useful for these tests.
 * Member functions are just tools to enable the tests.
 * Initialization of XGStateTest assumes that the frozen period has just begun. 
 * 
 */

class XGStateTest : public testing::Test
{
protected:
    XGStateTest()
    {
        init_vectors();
    };

    double surface_temperature = 0.0;
    std::unique_ptr<XG_algorithm::params> set_default_P();
    std::unique_ptr<XG_algorithm::state> set_default_S(const int& N);
    std::unique_ptr<XG_algorithm::params> P;
    std::unique_ptr<XG_algorithm::state> S;
    void distribute_moisture();

    // params
    int num_layers = 4;
    std::vector<double> depths;
    std::vector<double> porosity;
    std::vector<double> theta_default;
    std::vector<double> dry_soil_k;
    std::vector<double> saturated_frozen_soil_k;
    std::vector<double> saturated_unfrozen_soil_k;
    double Trigthrhld = 100.0;
    double theta_min = 0.001;
    double SWE_k = 0.35;
    double Zpf_init = 2.0;
    bool freeze_kw_ki_update = true;
    bool thaw_ki_kw_update = true;
    int k_update = 1;
    double soil_rechr_max = 350.0;
    double soil_moist_max = 625.0;
    int time_step_per_day = 24;
    double permafrost_depth = 2.0;
    bool calc_conductivity = false; 

    void init_vectors()
    {
        depths.resize(num_layers,0.5);
        porosity.resize(num_layers,0.5);
        theta_default.resize(num_layers,0.5);
        dry_soil_k.resize(num_layers,2.5);
        saturated_frozen_soil_k.resize(num_layers,1.98);
        saturated_frozen_soil_k[0] = 1.55;
        saturated_unfrozen_soil_k.resize(num_layers,1.67);
        saturated_unfrozen_soil_k[0] = 0.8;
    };    
};

// Use as first XGTest
//TEST_F(XGStateTest, NoRunDefaultOutput) {
//
//    XG_algorithm xg(surface_temperature,S,P);
//    
//    double thaw = xg.get_thaw_front_depth();
//
//    double freeze = xg.get_freeze_front_depth();
//
//    double first = xg.get_first_front_depth();
//
//    ASSERT_EQ(thaw,0.0);
//    ASSERT_EQ(freeze,0.0);
//    ASSERT_EQ(first,0.0);
//};
//

std::unique_ptr<XG_algorithm::state> XGStateTest::set_default_S(const int& N)
{
    std::unique_ptr<XG_algorithm::state> S = 
        std::make_unique<XG_algorithm::state>(N,*P);
    //S.set_layer_moisture_maximums(P)
    //    .set_thermal_conductivities(P)
    //    .set_freezethaw_ratios(P)
    //    .set_initial_freezethaw_depths(Zdf,Zdt);


    //S.nfront = 0;
    //S.Fz_low = 1;
    //S.Th_low = 1;
    //S.t_trend = 0.0;
    //S.Zd_front.assign(P.n,0.0);
    //S.XG_moist_d = 0.0;
    //S.XG_rechr_d = 0.0;
    //S.tc_composite.assign(P.n,0.0);
    //S.tc_composite2.assign(P.n,0.0);
    //S.XG_max.assign(P.n);
    //    
    //std::transform(P.por.begin(),P.por.end(),
    //        P.depths.begin(),
    //        S.XG_max.begin(), 
    //        [](double x,double y) {return x* y;}
    //        );
    //
    //S.theta = P.theta_default;
    //
    return S;

};

std::unique_ptr<XG_algorithm::params> XGStateTest::set_default_P(void)
{
    std::unique_ptr<XG_algorithm::params> P =
        std::make_unique<XG_algorithm::params>(depths,
                Trigthrhld,
                porosity,
                num_layers,
                theta_default,
                theta_min,
                dry_soil_k,
                saturated_frozen_soil_k,
                saturated_unfrozen_soil_k,
                SWE_k,
                Zpf_init,
                freeze_kw_ki_update,
                thaw_ki_kw_update,
                k_update,
                soil_rechr_max,
                soil_moist_max,
                time_step_per_day,
                calc_conductivity);

    return P;
};

TEST_F(XGStateTest,TestParams)
{
    P = set_default_P();

    EXPECT_EQ(P->por.size(),num_layers);
    EXPECT_EQ(P->getdepths().size(),num_layers);
    EXPECT_EQ(P->theta_default.size(),num_layers);
}

TEST_F(XGStateTest,MoistureMaxTest)
{
    P = set_default_P();
    S = set_default_S(P->N_Soil_layers);
    S->set_layer_moisture_maximums(*P);
    for (int i=0; i < P->N_Soil_layers; ++i)
    {
        EXPECT_EQ(S->XG_max[i],0.5*0.5*1000.0);
        EXPECT_EQ(S->theta[i],theta_default[i]);

        if (i == 0)
        {
            EXPECT_EQ(S->rechr_fract[i],1.0);
            EXPECT_EQ(S->moist_fract[i],0.0);
            EXPECT_EQ(S->default_fract[i],0.0);
        }
        else if (i == 1)
        {
            EXPECT_EQ(S->rechr_fract[i],
                    (soil_rechr_max - S->XG_max[0])/S->XG_max[1]);
            EXPECT_EQ(S->moist_fract[i],
                    1.0 - S->rechr_fract[i]);
            EXPECT_EQ(S->default_fract[i],0.0);
        }
        else if (i == 2)
        {
            EXPECT_EQ(S->rechr_fract[i],0.0);
            EXPECT_EQ(S->moist_fract[i],1.0);
            EXPECT_EQ(S->default_fract[i],0.0);
        }
        else if (i == 3)
        {
            EXPECT_EQ(S->rechr_fract[i],0.0);
            EXPECT_EQ(S->moist_fract[i],
                    (soil_moist_max - S->XG_max[1]*(1 - S->rechr_fract[1]) - S->XG_max[2])/S->XG_max[3]);
            EXPECT_EQ(S->default_fract[i],
                1 - S->moist_fract[i]);
        }
        else if (i == 4)
        {
            EXPECT_EQ(S->rechr_fract[i],0.0);
            EXPECT_EQ(S->moist_fract[i],0.0);
            EXPECT_EQ(S->default_fract[i],1.0);
        }
        

    }

    EXPECT_EQ(S->XG_rechr_d,
            0 + depths[0] + depths[1]*(soil_rechr_max - S->XG_max[0])/S->XG_max[1]);
    EXPECT_EQ(S->XG_moist_d,
            0 + depths[1] + depths[2] + depths[3]*(soil_moist_max - S->XG_max[1]*(1 - (soil_rechr_max - S->XG_max[0])/S->XG_max[1]) - S->XG_max[2])/S->XG_max[3]);

     

};

TEST_F(XGStateTest,MoistureMaxTest_Case2)
{
    num_layers = 1;
    soil_rechr_max = 100;
    soil_moist_max = 250;
    P = set_default_P();
    S = set_default_S(P->N_Soil_layers);
    S->set_layer_moisture_maximums(*P);
    
    EXPECT_EQ(S->XG_rechr_d,
            P->getdepths()[0]*soil_rechr_max/S->XG_max[0]);
    EXPECT_EQ(S->XG_moist_d,
            P->getdepths()[0]*soil_moist_max/S->XG_max[0]);
    EXPECT_EQ(S->rechr_fract[0],soil_rechr_max/S->XG_max[0]);
    EXPECT_EQ(S->moist_fract[0],(soil_moist_max - soil_rechr_max)/S->XG_max[0]);
    EXPECT_EQ(S->default_fract[0],1 - soil_moist_max/S->XG_max[0]);
};

TEST_F(XGStateTest,GetTtcAndGetFtcTest)
{
    //calc_conductivity off 
    P = set_default_P();
    S = set_default_S(P->N_Soil_layers);
    int i = 2;
    S->layer_h2o[i] = 23.0;

    double k_f = S->get_ftc(i,*P);

    double k_t = S->get_ttc(i,*P);

    double k_f_val = (1.0 - P->por[i])*P->soil_solid_km[i] + 
        S->layer_h2o[i]/1000.0*XG_algorithm::kw + (P->por[i] - S->layer_h2o[i]/1000.0)*XG_algorithm::ka;
    double k_t_val = (1.0 - P->por[i])*P->soil_solid_km[i] + 
        S->layer_h2o[i]/1000.0*XG_algorithm::ki + (P->por[i] - S->layer_h2o[i]/1000.0)*XG_algorithm::ka;

    EXPECT_EQ(k_f,k_f_val) << "off";
    EXPECT_EQ(k_t,k_t_val) << "off";
   
    // calc_conductivity on
    calc_conductivity = true;
    P = set_default_P();

    k_f = S->get_ftc(i,*P);
    k_t = S->get_ttc(i,*P);

    k_f_val = (P->soil_solid_km_kw[i] - P->soil_solid_km[i])*std::pow(S->layer_h2o[i]/(1000.0*P->por[i]),2) + P->soil_solid_km[i];
    k_t_val = P->soil_solid_km[i]*std::pow(P->soil_solid_km_ki[i]/P->soil_solid_km[i],S->layer_h2o[i]/(1000.0*P->por[i]));

    EXPECT_EQ(k_f,k_f_val) << "on";
    EXPECT_EQ(k_t,k_t_val) << "on";
   
};

TEST_F(XGStateTest,SetThermalConductivityTest)
{
	num_layers = 2;
	P = set_default_P();
    S = set_default_S(P->N_Soil_layers);
    S->ftc_contents[0] = 1; S->ttc_contents[0] = 1;
	S->ftc_contents[1] = 0; S->ttc_contents[1] = 0;
	const double soil_rechr = 100.0;
    const double soil_moist = 300.0;
    S->set_layer_moisture_maximums(*P)
        .set_thermal_conductivities(*P,soil_moist,soil_rechr);

    double XG_moist_compare = (S->rechr_fract[0]*soil_rechr/soil_rechr_max + S->moist_fract[0]*(soil_moist - soil_rechr)/(soil_moist_max - soil_rechr_max) + S->default_fract[0]*P->theta_default[0])*S->XG_max[0];
    EXPECT_DOUBLE_EQ(S->XG_moist[0],XG_moist_compare) << "Loop 0, XG moist";
	
	double check_XG_moist_compare = (S->rechr_fract[0]*soil_rechr/soil_rechr_max + S->moist_fract[0]*(soil_moist - soil_rechr)/(soil_moist_max - soil_rechr_max)) * S->XG_max[0];

	EXPECT_DOUBLE_EQ(S->theta[0],
			S->XG_moist[0] / S->XG_max[0]) << "Loop 0, theta";

	EXPECT_DOUBLE_EQ(S->layer_h2o[0],1000*S->theta[0]*P->por[0]) << "Loop 0, h2o";

	EXPECT_DOUBLE_EQ(S->ftc[0],S->get_ttc(0,*P)) << "Loop " << 0 << " ftc";

	EXPECT_DOUBLE_EQ(S->ttc[0],S->get_ftc(0,*P)) << "Loop 0, ttc";

	EXPECT_DOUBLE_EQ(S->pf[0],0.0);
	EXPECT_DOUBLE_EQ(S->pt[0],0.0);

    XG_moist_compare = (S->rechr_fract[1]*soil_rechr/soil_rechr_max + S->moist_fract[1]*(soil_moist - soil_rechr)/(soil_moist_max - soil_rechr_max) + S->default_fract[1]*P->theta_default[1])*S->XG_max[1];
    EXPECT_DOUBLE_EQ(S->XG_moist[1],XG_moist_compare) << "Loop 1, XG moist";
	
    check_XG_moist_compare += (S->rechr_fract[1]*soil_rechr/soil_rechr_max + S->moist_fract[1]*(soil_moist - soil_rechr)/(soil_moist_max - soil_rechr_max))*S->XG_max[1];

	EXPECT_DOUBLE_EQ(S->check_XG_moist,check_XG_moist_compare) << "check_XG_moist";
	
    EXPECT_DOUBLE_EQ(S->theta[1],
           S->XG_moist[1] / S->XG_max[1]) << "Loop 1, theta";

	EXPECT_DOUBLE_EQ(S->layer_h2o[1],1000*S->theta[1]*P->por[1]) << "Loop 1, h2o";

	EXPECT_DOUBLE_EQ(S->ftc[1],S->get_ftc(1,*P)) << "Loop " << 1 << " ftc";

	EXPECT_DOUBLE_EQ(S->ttc[1],S->get_ttc(1,*P)) << "Loop 1, ttc";

		

};

TEST_F(XGStateTest,SetThermalConductivityTest2)
{
	k_update = 0;
	num_layers = 1;
	P = set_default_P();
    S = set_default_S(P->N_Soil_layers);
    S->ftc_contents[0] = 1; S->ttc_contents[0] = 1;
	const double soil_rechr = 100.0;
    const double soil_moist = 300.0;
    S->set_layer_moisture_maximums(*P)
        .set_thermal_conductivities(*P,soil_moist,soil_rechr);

    double XG_moist_compare = (S->rechr_fract[0]*soil_rechr/soil_rechr_max + S->moist_fract[0]*(soil_moist - soil_rechr)/(soil_moist_max - soil_rechr_max) + S->default_fract[0]*P->theta_default[0])*S->XG_max[0];
    EXPECT_EQ(S->XG_moist[0],XG_moist_compare) << "Loop 0, XG moist";
	
	double check_XG_moist_compare = (S->rechr_fract[0]*soil_rechr/soil_rechr_max + S->moist_fract[0]*(soil_moist - soil_rechr)/(soil_moist_max - soil_rechr_max)) * S->XG_max[0];
	EXPECT_EQ(S->check_XG_moist,check_XG_moist_compare) << "Loop 0, check_XG_moist";

	EXPECT_EQ(S->theta[0],
			S->XG_moist[0] / S->XG_max[0]) << "Loop 0, theta";

	EXPECT_EQ(S->layer_h2o[0],1000*S->theta[0]*P->por[0]) << "Loop 0, h2o";

	EXPECT_EQ(S->ftc[0],S->get_ftc(0,*P)) << "Loop " << 0 << " ftc";

	EXPECT_EQ(S->ttc[0],S->get_ttc(0,*P)) << "Loop 0, ttc";
};

TEST_F(XGStateTest,SetThermaConductivityTest3)
{
	soil_moist_max = 0.0;
	num_layers = 1;
	P = set_default_P();
    S = set_default_S(P->N_Soil_layers);
    S->ftc_contents[0] = 1; S->ttc_contents[0] = 1;
	const double soil_rechr = 100.0;
	const double soil_moist = 300.0;
	S->set_layer_moisture_maximums(*P)
		.set_thermal_conductivities(*P,soil_moist,soil_rechr);

	EXPECT_EQ(S->XG_moist[0],S->default_fract[0]*P->theta_default[0]*S->XG_max[0]);

};

TEST_F(XGStateTest,SetFreezeThawRatiosTest)
{
	num_layers = 7;
    init_vectors();
	P = set_default_P();
    S = set_default_S(P->N_Soil_layers);
    for (int i = 0; i < P->N_Soil_layers; ++i)
	{
		S->ftc_contents[i] = 1; S->ttc_contents[i] = 1;
	}
	const double soil_rechr = 100.0;
    const double soil_moist = 300.0;
    
    S->set_layer_moisture_maximums(*P)
        .set_thermal_conductivities(*P,soil_moist,soil_rechr)
        .set_freezethaw_ratios(*P);

	for (int i = 0; i < P->N_Soil_layers; ++i)
	{
		if (i == 0)
		{	
			EXPECT_DOUBLE_EQ(S->pf[i], 0.0);
			EXPECT_DOUBLE_EQ(S->pt[i], 0.0);
		}
		else
		{
			EXPECT_DOUBLE_EQ(S->pf[i],std::sqrt(S->ftc[i-1]*S->layer_h2o[i]/(S->ftc[i]*S->layer_h2o[i-1])))
                << "Loop " << i;
			EXPECT_DOUBLE_EQ(S->pt[i],std::sqrt(S->ttc[i-1]*S->layer_h2o[i]/(S->ttc[i]*S->layer_h2o[i-1])))
                << "Loop " << i;
		}
	}

};

// TODO Mock state and param and move this to another file
class XGTest : public XGStateTest
{
protected:
    XGTest() {};
    const double soil_moist = 300;
    const double soil_rechr = 100;

    void do_my_set_up()
    {
        P = set_default_P();
        S = set_default_S(P->N_Soil_layers);
        
        S->set_layer_moisture_maximums(*P)
            .set_thermal_conductivities(*P,soil_moist,soil_rechr)
            .set_freezethaw_ratios(*P);
    }

};



TEST_F(XGTest,BasicOvertakenThawFrontTest)
{
	std::vector<double> TrigAcc{-Trigthrhld,Trigthrhld,
		Trigthrhld,-Trigthrhld,-Trigthrhld,-Trigthrhld};
	std::vector<double> t_trend(TrigAcc.size(),1.0);
	t_trend[0] *= -1.0;
	t_trend[3] *= -1.0;
	t_trend[4] *= -1.0;
	t_trend[5] *= -1.0;
	
	std::vector<double> Zd_front;

	std::vector<double> t_surface{-20.0,15.0,15.0,-10.0,-10.0,-8.5}; // Set this, essential for computing z_f and z_t, use it to control what kind of phase we are in (e.g., do fronts merge after a frozen/thaw front is overtaken by the surface thaw/freeze front. 

    do_my_set_up();

	for (int ii = 0; ii < TrigAcc.size(); ++ii)
	{
		
        S->is_newday = true;
		S->TrigAcc = TrigAcc.at(ii);
		S->t_trend = t_trend.at(ii);
		XG_algorithm XG(t_surface.at(ii),soil_moist,soil_rechr,*S,*P);
        XG.run();
		if (ii == 0 || ii == 1)
		{
            // Freeze at 0, Idle at 1
            //
            // Array at this stage:
            // [Z[0]] < 0
            // [0]
            // [0]
            // [0]
            // [0]
            
			EXPECT_TRUE(S->Zd_front.at(0) < 0.0) << "Loop: " << ii;
			for (int jj = 1; jj < S->Zd_front.size(); ++jj)
			{
				EXPECT_TRUE(S->Zd_front.at(jj) == 0.0) << "Loop: " << ii;
			}
            if (ii == 0)
			    Zd_front.push_back(S->Zd_front.at(0));
		}
		else if (ii == 2)
		{   
            // Thaw at 2
            //
            // Array at this stage:
            // [Z[1]] > 0
            // [Z[0]] < 0
            // [0]
            // [0]
            // [0]
			
            EXPECT_TRUE(S->Zd_front.at(0) > 0.0) << "Loop: " << ii;
			EXPECT_EQ(S->Zd_front.at(1),Zd_front.at(0)) << "Loop: " << ii;
			for (int jj = 2; jj < S->Zd_front.size(); ++jj)
			{
				EXPECT_TRUE(S->Zd_front.at(jj) == 0.0) << "Loop: " << ii;
			}
			Zd_front.push_back(S->Zd_front.at(0));
		}
		else if (ii == 3)
		{   
            // Idle at 3
            //
            // Array at this stage:
            // [0]
            // [Z[1]] > 0
            // [Z[0]] < 0
            // [0]
            // [0]
			
            EXPECT_EQ(S->Zd_front.at(0),0.0) << "Loop: " << ii;
			EXPECT_EQ(S->Zd_front.at(1),Zd_front.at(1)) << "Loop: " << ii;
			EXPECT_EQ(S->Zd_front.at(2),Zd_front.at(0)) << "Loop: " << ii;
			for (int jj = 3; jj < S->Zd_front.size(); ++jj)
			{
				EXPECT_TRUE(S->Zd_front.at(jj) == 0.0) << "Loop: " << ii;
			}
		}
		else if (ii == 4)
		{   
            // Freeze at 4
            //
            // Array at this stage:
            // [Z[2]] < 0
            // [Z[1]] > 0
            // [Z[0]] < 0
            // [0]
            // [0]
            EXPECT_TRUE(S->Zdf < S->Zdt) << "Loop: " << ii;
			EXPECT_TRUE(S->Zd_front.at(0) < 0.0) << "Loop: " << ii;
			EXPECT_EQ(S->Zd_front.at(1),Zd_front.at(1)) << "Loop: " << ii;
			EXPECT_EQ(S->Zd_front.at(2),Zd_front.at(0)) << "Loop: " << ii;
			for (int jj = 3; jj < S->Zd_front.size(); ++jj)
			{
				EXPECT_TRUE(S->Zd_front.at(jj) == 0.0) << "Loop: " << ii;
			}
			Zd_front.push_back(S->Zd_front.at(0));
		}
		else if (ii == 5)
        {
            // Freeze at 5
            //
            // Array at this stage:
            // [Z[3]] < 0
            // [0] 
            // [0]
            // [0]
            // [0]
            EXPECT_EQ(S->Zdf,-S->Zd_front[0]) << "Loop: " << ii;
		    EXPECT_TRUE(S->Zd_front[0] < 0) << "Loop: " << ii;
            EXPECT_EQ(S->Zdt,0.0) << "Loop: " << ii;
             
            for (int jj = 1; jj < S->Zd_front.size(); ++jj)
                EXPECT_TRUE(S->Zd_front.at(jj) == 0.0) << "Loop: " << ii;
        }
    }	
    //EXPECT_EQ(S->Zd_front[0],z_1);
    //for (int ii = 1; ii < S->Zd_front.size(); ++ii)
    //{
    //   EXPECT_TRUE(S->Zd_front[ii] == 0.0);
    //}

    //S->TrigAcc = P->Trigthrhld;
    //S->t_trend = 1.0;

    //EXPECT_TRUE(S->Zd_front[0] > 0.0);
    //std::cout << S->TrigState <<std::endl;
    //EXPECT_EQ(S->Zd_front[1],z_1);
    //for (int ii = 1; ii < S->Zd_front.size(); ++ii)
    //{
    //   EXPECT_TRUE(S->Zd_front[ii] == 0.0);
    //}      

    //double z_2 = S->Zd_front[0];

};

#define diff 0.000001
TEST_F(XGTest,FullZdFrontOrganizeTest)
{
    num_layers=6;
    do_my_set_up();
    std::vector<double> Z{-0.1,0.234,-0.758,0.9,-1.3,0.};
    S->Zd_front = Z;
    S->nfront=3;
    S->TrigState = -1;
    S->t_trend = -1;
    S->Zdt = S->Zd_front[1];
    S->Zdf = -S->Zd_front[0];
    double t_surface = -300.0;
    S->TrigAcc = -Trigthrhld;
    S->is_newday = true;

    XG_algorithm XG(t_surface,soil_moist,soil_rechr,*S,*P);
    XG.run();

    for (int i = 0; i < Z.size()-2; ++i)
    {
        EXPECT_NEAR(S->Zd_front.at(i),Z.at(i+2),diff) << "Loop: " << i;
    }
    for (int i = Z.size()-2; i < Z.size(); ++i)
        EXPECT_EQ(S->Zd_front.at(i),0.0) << "Loop: " << i;
    
    EXPECT_NEAR(S->Zdf,std::fabs(Z[2]),diff);
    EXPECT_NEAR(S->Zdt,Z[3],diff);

};


    



