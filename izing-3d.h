/*
 * izing-3d.h
 *
 *  Created on: 22 апр. 2014 г.
 *      Author: nlare
 */

#include "models.cpp"

#ifndef IZING_3D_H_
#define IZING_3D_H_

class Izing3D: public Models  {
    public:
        Izing3D(string _name, int _grid_size, int _mks_max, int _statistic, int _streams);
        ~Izing3D(void);

        void setTemperatureRange(double, double, double, double);
        void setTemperatureRange(double, double);
        void setTemperatureStatic(double);
        void setTemperatureStatic(double, double, double, double);
        void setStreams(int);
        void setMethod(string);

        // void WriteRelax();
        void ShowCharacteristicOfModel();
        void FormSourceConf();
        void WriteTdCharToFile();
        // void WriteTimeToFile();

        // void Relaxation(int);
        void SaveSpinConf(int);
        // void ACF(int);
        // void StatisticProcessing();

        void Metropolis(double);
        void SwendsenWang(double);
        void WolfStackBased(double);
        void WolfMassiveBased(double);

        int Start();
        int PlotGraph(string);
        // void measure_only_relax();

    private:

        // int ***spin;

        double T_div[1000];
        double T, T_min, T_max, T_dr1, T_dr2;

        struct _data_t {

            double energy_per_spin;
            double magnetic_per_spin;

            double energy;
            double magnet;
            double C_v;
            double X;

            double relax;
            double correlation;

            int old_config;
            int new_config;

            double first_correlation_mcs;
            double first_correlation_magnet;

            double sum[5];

            _data_t()    {
                energy_per_spin = 0.0;
                magnetic_per_spin = 0.0;
            }

        } data[128];

        struct _statistic_t {

            double magnet;
            double energy;

            double relax;
            double correlation;

            double error_relax[2]; // 0 element - delta = (<mcs^2> - <mcs>^2)^1/2, 1 element - <mcs^2>, 2 element - <mcs>^2
            double error_correlation[2]; // 0 element - delta = (<mcs^2> - <mcs>^2)^1/2, 1 element - <mcs^2>, 2 element - <mcs>^2

        } statistic[128];

        struct trigger_t {
            bool correlation;
            bool relax_print;
            bool relax;

            trigger_t()   {
                correlation = true;
                relax_print = true;
                relax = true;
            }
        }   trigger[128];

        string algorithm;

        double first_mcs_magnet;

        int cluster_elements;

        int streams;
        int thread_num;

        int ***spin, ***spin_old, ***spin_thread;
        int statistic_max, mcs_max;
        int T_div_max_index;
        bool static_temp, static_procs_parallel;
        int static_procs_count;
        bool remove_relax, remove_appr, only_relax;
        char * name_relax, * name_appr, * name_temp, * name_U4;

        std::ofstream out_relax, out_appr, out_U4, out_time;

        int time_begin, time_end, runtime_sec;
        int days, hours, minutes, seconds;
        char * name_time_f;
        ofstream time_f;

        bool first_file_open;
           // out_appr - файл с аппроксимированными значениями

        void magnetic(int);
        void energy(int);

        int nearly_spins(int, int, int);
};

#endif /* IZING_3D_H_ */
