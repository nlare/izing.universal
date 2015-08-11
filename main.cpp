/*
 * =====================================================================================
 *
 *       Filename:  main.cpp
 *
 *    Description:
 *
 *        Version:  1.0
 *        Created:  07.12.2013 02:54:02
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  nlare
 *   Organization:  univer
 *
 * =====================================================================================
 */

// #include "izing-2d.cpp"
#include "izing-3d.cpp"
//#include "izing_data.pb.h"
//#include "izing_data.pb.cc"

using namespace std;

int main(int argc, char *argv[])    {
    long omp_time_start, omp_time_end, runtime_sec;
    stringstream command, time_filename;

    // Izing2D Model("Izing2D", 16, 1000, 100);
    //params - name, grid_size, mcs_max, statistic
    Izing3D Model("Izing3D", 128, 1000, 100, 1);

#ifndef _OPENMP
    cout << "OpenMP is not supported!" << endl;
    exit(0)
#endif

    // Model.setStreams(1);
    // Model.setTemperatureStatic(4.3, 4.51, 4.7, 5.0);
    Model.setTemperatureRange(1.8, 3.1);
    Model.setMethod("WS"); //

    cout << "this!" << endl;

    // omp_time_start = omp_get_wtime();

    Model.Start();
        // cerr << "Cannot start meashuring..." << endl;
    // }

    // if(Model.PlotGraph("E"))    {
    //     cout << "End of plot." << endl;
    // }   else    {
    //     cerr << "Error, graphs not gen. " << endl;
    // }
    // Model.PlotGraph("M");
    // Model.PlotGraph("X");
    // Model.PlotGraph("C_v");

    // omp_time_end = omp_get_wtime();

    // runtime_sec = omp_time_end - omp_time_start;

    // cout << "RUNTIME: [" << days << "'d " << hours << "'h " << minutes << "'m " << seconds << "'s ]" << std::endl;

    // getchar();

    return 0;
}
