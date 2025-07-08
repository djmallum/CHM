#include "gtest/gtest.h"
#include "CSVreader.hpp"
//reader comes from CSVReader.hpp as inline

TEST(CSVreaderTeat, CheckFunction) 
{
    double value = reader.getValue<double>("net_rain",1);

    ASSERT_DOUBLE_EQ(value,0.0);

};

TEST(CSVreaderTeat, CheckLaterValue) 
{   
    double value = reader.getValue<double>("net_rain",136768);

    ASSERT_DOUBLE_EQ(value,0.213297);

};
