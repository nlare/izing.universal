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
	// 1-й параметр - 3-х мерная (Izing3D) или 2-х мерная (Izing2D) модель Изинга будет использоваться
	// 2-й - размер решетки
	// 3-й - количество шагов Монте-Карло
	// 4-й - число потоков для обработки температур параллельно
	//
	//
    Izing3D Model("Izing3D", 32, 1000, 20, 1);

#ifndef _OPENMP
    cout << "OpenMP is not supported!" << endl;
    exit(0)
#endif

    // Model.setStreams(1);
    // Model.setTemperatureStatic(1.8, 4.2, 4.8, 6.5);
    Model.setTemperatureRange(1.8, 6.5, 4.2, 4.8);
    Model.setMethod("M"); //

//    cout << "this!" << endl;

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
