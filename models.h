/*
 * models.h
 *
 *  Created on: 12 марта 2014 г.
 *      Author: nlare
 */
#pragma once

#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <stack>
#include <omp.h>
#include "mersenne.cpp"

using namespace std;

#ifndef MODELS_H_
#define MODELS_H_

class Models    {
    public:
        Models(string, int);
        void show_choose_model(void);
        void set_model_name(string);
        void set_dimension(int);
        void set_grid_size(int);
        int N;
    private:
        string name;
        //int dimension;
        int grid_size;
//        double energy_per_spin;
//        double magnet_per_spin;
};

#endif /* MODELS_H_ */
