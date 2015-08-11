/*
 * =====================================================================================
 *
 *       Filename:  Models.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  20.12.2013 14:52:12
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */
#pragma once

#include "models.h"

using namespace std;

Models::Models(string _name, int _grid_size)   {
	  set_model_name(_name);
//    set_dimension(_dimension);
	  set_grid_size(_grid_size);
}
void Models::show_choose_model()  {
    cout << "You'r choose is " << name << setw(2) << "Grid = " << N << endl;
}
void Models::set_model_name(string _name)   {
    name = _name;
}
/*void Models::set_dimension(int _dimension)    {
    dimension = _dimension;
}*/
void Models::set_grid_size(int _grid_size)    {
    N = _grid_size;
}
