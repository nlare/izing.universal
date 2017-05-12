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

#include <list>
#include <algorithm>
#include <dirent.h>
#include "izing-3d.cpp"
// #include "izing-2d.cpp"
//#include "izing_data.pb.h"
//#include "izing_data.pb.cc"

using namespace std;

// namespace fs = std::filesystem;

int main(int argc, char *argv[])    {

    #ifndef _OPENMP
        cout << "OpenMP is not supported!" << endl;
        exit(0)
    #endif
    
    long omp_time_start, omp_time_end, runtime_sec;
    stringstream command, time_filename;

    int L, MCS, STAT, PP;
    double T_min, T_max, T_prec_min, T_prec_max;
    string simulation_method;

    bool uncorrect_parameters = false;

    // Izing2D Model("Izing2D", 16, 1000, 100);
    //params - name, grid_size, mcs_max, statistic
	// 1-й параметр - 3-х мерная (Izing3D) или 2-х мерная (Izing2D) модель Изинга будет использоваться
	// 2-й - размер решетки
	// 3-й - количество шагов Монте-Карло
	// 4-й - статистика
	// 5-й - число потоков для обработки температур параллельно
	
    DIR *path_to_files_with_parameters;
    struct dirent *rec_of_file;

    int count_of_input_files, count_of_parameters;
    std::string current_parameter;
    std::ifstream list_of_parameters;

    list<string> list_of_input_files;
    list<string>::iterator input_files_iterator;

    count_of_input_files = 0;
    count_of_parameters  = 0;

    stringstream buffer;
    string filename;

    if((path_to_files_with_parameters = opendir("input_parameters/")) != NULL)   {

        while((rec_of_file = readdir(path_to_files_with_parameters)) != NULL)   {

            count_of_input_files++;

            buffer.str("");
            buffer << rec_of_file->d_name;

            list_of_input_files.insert(list_of_input_files.begin(), buffer.str().c_str());

        }

    }

    bool increment_iterator = true;
    // Пройдем весь двусвязный список от начала до конца с помощью итератора
    for(auto it = list_of_input_files.begin(); it != list_of_input_files.end(); it++)  {
        // Делаем инкремент итератора +2 в начале для того чтобы исключить "." и ".." из списка имен файлов
        if(increment_iterator) {

            increment_iterator = false;
            // Нужно для перемещения итератора на определенную позицию списка
            advance(it,2);

        }

        count_of_input_files++;

        buffer.str("");
        buffer << "input_parameters/" << (string)*it;

        filename = buffer.str();
        replace(filename.begin(), filename.end(), '\n', ' ');

        // cout << filename << "\n";

        if(!(list_of_parameters.bad())) {

            list_of_parameters.close();

        }    

        list_of_parameters.open(filename);
    
        if(!(list_of_parameters.bad()))   {
            
            std::cout << "-------------------------------------- \n";
            std::cout << " ... READ FILE " << filename << "\n\n";
            // std::cout << "-------------------------------------- \n";

            count_of_parameters = 0;
    
            while(getline(list_of_parameters, current_parameter))   {
    
                count_of_parameters++;

                if(count_of_parameters == 1) {
                    
                    // replace(current_parameter.begin(), current_parameter.end(), '\n', ' ');
                    // current_parameter.erase(current_parameter.size() - 1);
                    simulation_method = current_parameter;
                    std::cout << count_of_parameters << " :: simulation_method = " << simulation_method << "\n";                
    
                }
    
                if(count_of_parameters == 2) { 
    
                    L   = std::stoi(current_parameter);
                    std::cout << count_of_parameters << " :: L = " << current_parameter << "\n";
    
                }
    
                if(count_of_parameters == 3) {
    
                    MCS = std::stoi(current_parameter);
                    std::cout << count_of_parameters << " :: MCS = " << current_parameter << "\n";
    
                }
    
                if(count_of_parameters == 4) {
    
                    STAT = std::stoi(current_parameter);
                    std::cout << count_of_parameters << " :: STAT = " << current_parameter << "\n";
    
                }
    
                if(count_of_parameters == 5) {
    
                    PP   = std::stoi(current_parameter);
                    std::cout << count_of_parameters << " :: PP = " << current_parameter << "\n";                
    
                }

                if(count_of_parameters == 6) {
    
                    T_min = std::stoi(current_parameter);
                    std::cout << count_of_parameters << " :: T_min = " << current_parameter << "\n";                
    
                }

                if(count_of_parameters == 7) {
    
                    T_max = std::stoi(current_parameter);
                    std::cout << count_of_parameters << " :: T_max = " << current_parameter << "\n";                
    
                }

                if(count_of_parameters == 8) {
    
                    T_prec_min = std::stoi(current_parameter);
                    std::cout << count_of_parameters << " :: T_prec_min = " << current_parameter << "\n";                
    
                }

                if(count_of_parameters == 9) {
    
                    T_prec_max = std::stoi(current_parameter);
                    std::cout << count_of_parameters << " :: T_prec_max = " << current_parameter << "\n";                
    
                }
    
            }
            
            cout << "\ncount_of_parameters = " << count_of_parameters << "\n\n";
            // cout << "-------------------------------------- \n";
            
            if(count_of_parameters < 7)    {

                cout << "One of the file have errors: \"" << filename << "\". Try to set correct count of parameters (7 OR 9!). ABORT.\n\n";
                cout << "EXAMPLE:\n";
                cout << "-------------- \n";
                cout << "M\n128\n1000\n20\n4\n1.5\n6.5\n4.4\n4.6\n";
                cout << "-------------- \n";
                cout << "DESCRIPTION:\n";
                cout << "-------------- \n";
                cout << "1 :: simulation_method   (Simulation Method e.g. Metropolis - [M], Wolf (StackBased - [WS]), Wolf (MassiveBased - [WM]), SwendsenWang - [SW])\n";
                cout << "2 :: L                   (Size of the lattice)\n";
                cout << "3 :: MCS                 (Count of Monte-Carlo steps)\n";
                cout << "4 :: STAT                (Statistic)\n";
                cout << "5 :: PP                  (Parallel threads (Work on STAT, no TEMP))\n";
                cout << "6 :: T_min               (Begin temperature of simulation)\n";
                cout << "7 :: T_max               (End temperature of simulation)\n";
                cout << "8 :: T_prec_min          (Begin of interval where steps for TEMP equal 0.01)\n";
                cout << "9 :: T_prec_max          (End of interval where steps for TEMP equal 0.01)\n";
                std::cout << "-------------- \n";
                exit(1);

            }

            if(count_of_parameters == 7) {

                uncorrect_parameters = false;

                if((simulation_method != "M") && (simulation_method != "WS") && (simulation_method != "WM") && (simulation_method != "SW"))   {

                    uncorrect_parameters = true;
                    cout << "1 :: simulation_method = " << simulation_method << ". Set correct simulation method parameter! Must be [M, WS, WM, SW]!\n";

                }

                if((L%2)!=0)    {
                    
                    uncorrect_parameters = true;
                    cout << "2 :: L = " << L << " .Set correct L! Must be devided on two!\n";
                    
                }

                if((STAT%PP)!=0)    {

                    uncorrect_parameters = true;
                    cout << "4,5 :: STAT = " << STAT << ", PP = " << PP << ". PP must to devide the whole STAT.";

                }

                if(T_min < 0 || T_max < 0)   {

                    uncorrect_parameters = true;
                    cout << "6 :: T_min < 0 OR 7 :: T_max < 0. Must be more than 0";

                }

                if((T_min > T_max)) {

                    uncorrect_parameters = true;
                    cout << "6,7 :: T_min cannot be bigger than T_max!";

                }

                if(uncorrect_parameters) exit(1);
                else cout << " ... CONFIG LOOKING FINE! ... \n";

            }

            if(count_of_parameters == 9) {

                uncorrect_parameters = false;

                if((simulation_method != "M") && (simulation_method != "WS") && (simulation_method != "WM") && (simulation_method != "SW"))   {

                    uncorrect_parameters = true;
                    cout << "1 :: simulation_method = " << simulation_method << ". Set correct simulation method parameter! Must be [M, WS, WM, SW]!\n";

                }

                if((L%2)!=0)    {
                    
                    uncorrect_parameters = true;
                    cout << "2 :: L = " << L << " .Set correct L! Must be devided on two!\n";
                    
                }

                if((STAT%PP)!=0)    {

                    uncorrect_parameters = true;
                    cout << "4,5 :: STAT = " << STAT << ", PP = " << PP << ". PP must to devide the whole STAT.";

                }

                if(T_min < 0 || T_max < 0)   {

                    uncorrect_parameters = true;
                    cout << "6 :: T_min < 0 OR 7 :: T_max < 0. Must be more than 0";

                }

                if((T_min > T_max)) {

                    uncorrect_parameters = true;
                    cout << "6,7 :: T_min cannot be bigger than T_max!";

                }

                if((T_prec_min < T_min) || (T_prec_min > T_max))    {

                    uncorrect_parameters = true;
                    cout << "8 :: T_prec_min must be inside common TEMP interval!";

                }

                if((T_prec_max < T_min) || (T_prec_max > T_max))    {

                    uncorrect_parameters = true;
                    cout << "9 :: T_prec_max must be inside common TEMP interval!";

                }

                if((T_prec_min > T_prec_max)) {

                    uncorrect_parameters = true;
                    cout << "8,9 :: T_prec_min cannot be bigger than T_prec_max!";

                }

                if(uncorrect_parameters) exit(1);
                else cout << " ... CONFIG LOOKING FINE! ... \n";
                
            }

            Izing3D Model("Izing3D", L, MCS, STAT, PP);

            Model.setMethod(simulation_method);

            Model.setTemperatureRange(T_min, T_max, T_prec_min, T_prec_max);

            Model.Start();

        }

    }

    exit(0);
	
    // Izing3D Model("Izing3D", 192, 1000, 20, 2);

    // Model.setStreams(1);
    // Model.setTemperatureStatic(1.8, 4.2, 4.8, 6.5);
    // Model.setTemperatureRange(2.1, 6.5, 4.4, 4.6);
    // Model.setMethod("M"); //

//    cout << "this!" << endl;

    // omp_time_start = omp_get_wtime();

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

    list_of_parameters.close();

    return 0;
}
