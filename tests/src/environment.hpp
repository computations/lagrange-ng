#ifndef LAGRANGE_ENV_HPP
#define LAGRANGE_ENV_HPP

#include "gtest/gtest.h"
#include <fstream>
#include <iostream>

#define STRING(s) #s
#define STRINGIFY(s) STRING(s)

class LagrangeEnvironment: public ::testing::Environment{

    public:
    LagrangeEnvironment(){}

    std::ifstream get_datafile(){
        return std::ifstream(STRINGIFY(TREEPATH));
    }

};

extern LagrangeEnvironment * env;
#endif
