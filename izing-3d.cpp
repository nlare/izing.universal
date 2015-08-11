/*
 * izing.cpp
 *
 *  Created on: 03 марта 2014 г.
 *      Author: nlare
 */

#include "izing-3d.h"

using namespace std;

Izing3D::Izing3D(string _name, int _grid_size, int _mcs_max, int _statistic_max, int _streams): Models(_name, _grid_size)  {
    
    statistic_max = _statistic_max;
    mcs_max = _mcs_max;
    streams = _streams;

    stringstream ss;
    ss << "result/" << "TIME-3D-IZING-" << _grid_size << "-" << statistic_max << "-THR-" << streams << ".DAT";
    name_time_f = new char [strlen(ss.str().c_str())];
    strcpy(name_time_f, ss.str().c_str());

    time_f.open(name_time_f);

    time_begin = omp_get_wtime();

    omp_set_num_threads(streams);

    spin = new int **[_grid_size];
    for(int i=0;i<N;i++)    {
    	spin[i] = new int *[_grid_size];
    	for(int j=0;j<N;j++)	{
    		spin[i][j] = new int [_grid_size];
    	}
    }

    spin_old = new int **[_grid_size];
    for(int i=0;i<N;i++)    {
    	spin_old[i] = new int *[_grid_size];
    	for(int j=0;j<N;j++)	{
    		spin_old[i][j] = new int [_grid_size];
    	}
    }

    spin_thread = new int **[_grid_size];
    for(int i=0;i<N;i++)    {
        spin_thread[i] = new int *[_grid_size];
        for(int j=0;j<N;j++)    {
            spin_thread[i][j] = new int [_grid_size];
        }
    }

    T_min = 0;
    T_max = 0;
    T_dr1 = 0;
    T_dr2 = 0;

    for(int thread_num = 0; thread_num < streams; thread_num++)   {
        data[thread_num].energy_per_spin = 0.0;
        data[thread_num].magnetic_per_spin = 0.0;
        data[thread_num].new_config = 0;
        data[thread_num].old_config = 0;
    }

    static_procs_parallel = false;
    remove_relax = true;
    remove_appr = true;
    only_relax = false;
    first_file_open = true;

    // trig1 = true;
    // trig2 = true;
    // trig3 = true;

    cluster_elements = 0;
}

Izing3D::~Izing3D() {

    time_end = omp_get_wtime();

    runtime_sec = time_end - time_begin;

    // cout << "RUNTIME: [ " << runtime_sec << " ]\n"; 

    days = runtime_sec/60/60/24;
    hours = (runtime_sec/60/60)%24;
    minutes = (runtime_sec/60)%60;
    seconds = runtime_sec%60;

    time_f << std::fixed << "D\tH\tM\tS\n" << days << "\t" << hours << "\t" << minutes << "\t" << seconds << "\n";

    cout << "RUNTIME: [" << days << "'d " << hours << "'h " << minutes << "'m " << seconds << "'s ]" << std::endl;

    cout << endl << "Izing3D model measurement end" << endl;

    time_f.close();

    delete name_time_f;

    delete [] spin;
    delete [] spin_old;
    delete [] spin_thread;

}

void Izing3D::ShowCharacteristicOfModel()  {
    cout << "CHARACTERISTICS_OF_MODEL: " << endl;
    cout << "Grid = " << N << endl;
    cout << "Mcs = " << mcs_max << endl;
    cout << "statistic[thread_num] = " << statistic_max << endl;
    if(T_min!=0 && T_max!=0 && T_dr1!=0 && T_dr2!=0)	{
    	cout << "T_min = " << T_min << endl;
    	cout << "T_max = " << T_max << endl;
    	cout << "T_dr1 = " << T_dr1 << endl;
    	cout << "T_dr2 = " << T_dr2 << endl;
    }	else	cout << "T = " << T_div[0] << endl;
}

void Izing3D::FormSourceConf()   {

    for(int i = 0; i < N; i++)  {
        for(int j = 0; j < N; j++)  {
        	for(int k = 0; k < N; k++)  {
                spin[i][j][k] = 1; // ФОрмируем равновесную конфигурацию
        	}
        }
    }
}

int Izing3D::nearly_spins(int i, int j, int k)    {

    /*double sum, sum_i, sum_j;

    switch(i)   {
        case 0: sum_i = spin[N-1][j] + spin[i+1][j]; break;
        case NWO: sum_i = spin[N-2][j] + spin[0][j]; break;
        default: sum_i = spin[i+1][j] + spin[i-1][j]; break;
    }
    switch(j)   {
        case 0:sum_j = spin[i][N-1] + spin[i][0]; break;
        case NWO: sum_j = spin[i][N-2] + spin[i][0]; break;
        default: sum_j = spin[i][j+1] + spin[i][j-1]; break;
    }
    sum = sum_i + sum_j;
//    cout << "sum = " << sum << endl;
    return sum;*/

    int result;

    if(i==0)    {
        result=spin[N-1][j][k];
    }  else    {
        result=spin[i-1][j][k];
    }
    if(i==N-1)  {
        result=result+spin[0][j][k];
    }  else    {
        result=result+spin[i+1][j][k];
    }
    if(j==0)   {
        result=result+spin[i][N-1][k];
    }  else    {
        result=result+spin[i][j-1][k];
    }
    if(j==N-1) {
        result=result+spin[i][0][k];
    }  else    {
        result=result+spin[i][j+1][k];
    }
    if(k==0)   {
    	result=result+spin[i][j][N-1];
    }  else    {
    	result=result+spin[i][j][k-1];
    }
    if(k==N-1) {
    	result=result+spin[i][j][0];
    }  else    {
    	result=result+spin[i][j][k+1];
    }
    return result;
}

void Izing3D::energy(int thread_num)   {

    data[thread_num].energy_per_spin = 0.0;

    for(int i=0;i<N;i++)    {
        for(int j=0;j<N;j++)    {
        	for(int k=0;k<N;k++)	{
        		data[thread_num].energy_per_spin -= (spin[i][j][k] * nearly_spins(i,j,k));
        	}
        }
    }
//    cout << "::data[thread_num].energy_per_spin[]= " << data[thread_num].energy_per_spin << endl;
    data[thread_num].energy = data[thread_num].energy_per_spin/(6*N*N*N);
}

void Izing3D::magnetic(int thread_num)     {

    data[thread_num].magnetic_per_spin = 0.0;

    for(int i=0;i<N;i++)    {
        for(int j=0;j<N;j++)    {
        	for(int k=0;k<N;k++)	{
        		data[thread_num].magnetic_per_spin += spin[i][j][k];
        	}
        }
    }
//    cout << "::data[thread_num].magnetic_per_spin[]= " << data[thread_num].magnetic_per_spin << endl;
    data[thread_num].magnet = abs(data[thread_num].magnetic_per_spin/(N*N*N));
}

void Izing3D::setTemperatureRange(double _T_min, double _T_max, double _T_dr1, double _T_dr2) {

    struct step_t {
        float little;
        float big;
    } step;

    static_temp = false;

    step.little = 0.01;
    step.big = 0.1;

    T_min = _T_min;
    T_max = _T_max;
    T_dr1 = _T_dr1;
    T_dr2 = _T_dr2;

    int i = 0;
    if(T_max > T_min) {
        for(T=T_min; T<=T_dr1; T=T+step.big)   {
            T_div[i++] = T;
//            cout << i << ": " << T << endl;
        }
        for(T=T_dr1; T<=T_dr2;T=T+step.little)  {
            T_div[i++] = T;
//            cout << i << ": " << T << endl;
        }
        for(T=T_dr2; T<=T_max;T=T+step.big)  {
            T_div[i++] = T;
//            cout << i << ": " << T << endl;
        }
    }   else    {
        for(T=T_min; T>=T_dr2; T=T-step.big)   {
            T_div[i++] = T;
//            cout << i << ": " << T << endl;
        }
        for(T=T_dr2; T>=T_dr1;T=T-step.little)  {
            T_div[i++] = T;
//            cout << i << ": " << T << endl;
        }
        for(T=T_dr1; T>=T_max;T=T-step.big)  {
            T_div[i++] = T;
//            cout << i << ": " << T << endl;
        }
    }
    T_div_max_index = i;
}

void Izing3D::setTemperatureRange(double _T_min, double _T_max) {

    struct step_t {
        float big;
    } step;

    static_temp = false;

    step.big = 0.1;

    T_min = _T_min;
    T_max = _T_max;

    int i = 0;
    if(T_min < T_max)    {
        for(T=T_min; T<=T_max; T=T+step.big)   {
            T_div[i++] = T;
            cout << i << ": " << T << endl;
        }
    }   else  {
        for(T=T_min; T>=T_max; T=T-step.big)   {
            //cout << "T_min = " << T_min << endl;
            T_div[i++] = T;
            cout << i << ": " << T << endl;
        }
    }
    T_div_max_index = i;
}

void Izing3D::setTemperatureStatic(double _T)	{
	static_temp = true;
    T_div_max_index = 1;
	T_div[0] = _T;
}

void Izing3D::setTemperatureStatic(double _T1, double _T2, double _T3, double _T4)	{
	static_temp = true;
	static_procs_parallel = true;
	T_div[0] = _T1;
	T_div[1] = _T2;
	T_div[2] = _T3;
	T_div[3] = _T4;
	T_div_max_index = 4;
	static_procs_count = 4;
}

void Izing3D::setMethod(string _name)  {
    algorithm = _name;
}

void Izing3D::setStreams(int _streams)	{
	streams = _streams;
}

// void Izing2D::Relaxation(int _mcs)   {
//     double first_mcs_magnet;

//     // Запоминаем намагниченность на первом шаге
//     if(_mcs == 0)   {
//         first_mcs_magnet = data[thread_num].magnet;
//         //cout << "first_mcs_magnet = " << first_mcs_magnet << endl;
//     }
//     // Сравниваем намагниченность полученную на первом шаге, с намагниченностью полученной далее
//     // до момента, пока она не уменьшится на exp
//     if(data[thread_num].relax_mcs == 0 && data[thread_num].magnet <= first_mcs_magnet/exp(1))   {
//         data[thread_num].relax_mcs = _mcs;   // Время, с которого система начинает релаксировать
//         data[thread_num].relax_mcs_sum += data[thread_num].relax_mcs; // Статистика для усредения времени релаксации
//         /* Переменные, где считаем слагаемые для расчета погрешности как дисперсии
//          * ------------------------------ */
//         data[thread_num].error_relax_mcs[1] += _mcs; 
//         data[thread_num].error_relax_mcs[2] += _mcs*_mcs;
//         /* ------------------------------ */
//         data[thread_num].correlation_mcs = 0;
//         // cout << "Relaxation: relax_mcs = " << data[thread_num].relax_mcs << "\t" << endl;
//     }
// }

// void Izing2D::ACF(int _mcs) {

//     // int real_mcs = data[thread_num].relax_mcs; 

//     for(int i = 0; i < N; i++)  {
//         for(int j = 0; j < N; j++)  {
//              /* Считаем новое состояние системы исходя из старого, суммируя каждое новое полученное */
//             data[thread_num].new_config += spin_old[i][j]*spin[i][j];
//         }
//     }

//     if(data[thread_num].correlation_mcs == 0 && (double)data[thread_num].old_config <= (double)data[thread_num].new_config/2.71) {
//         // cout << "old_config: " << data[thread_num].old_config << ", data[thread_num].new_config: " << data[thread_num].new_config << endl;    
//         // trigger[thread_num].correlation = false;
        
//         data[thread_num].correlation_mcs = _mcs - data[thread_num].relax_mcs*10.0;

//         // cout << " Correlation: data[thread_num].correlation_mcs = " << data[thread_num].correlation_mcs << ", relax_mcs = " << data[thread_num].relax_mcs << endl;

//         // if(data[thread_num].correlation_mcs < data[thread_num].relax_mcs*10)    {
//             data[thread_num].correlation_mcs_sum += data[thread_num].correlation_mcs;
//             data[thread_num].error_correlation_mcs[1] += data[thread_num].correlation_mcs;
//             data[thread_num].error_correlation_mcs[2] += data[thread_num].correlation_mcs*data[thread_num].correlation_mcs;
//         // }
//     }
// }

void Izing3D::WriteTdCharToFile()   {

    double E = 0, M = 0, C_v = 0, X = 0, relax = 0, correlation = 0, error_relax = 0, error_correlation = 0, M2 = 0, M4 = 0, U4 = 0;

    struct _buffer  {

      double M;
      double E;
      double C_v;
      double X;

      double R;
      double C;
      double R_E;   //  relax error
      double R_E_SQUARE;
      double C_E;   //  correlation error
      double C_E_SQUARE;

      double M2;
      double M4;

    _buffer()   {
        M = 0.0;
        E = 0.0;
        C_v = 0.0;
        X = 0.0;

        R = 0.0;
        C = 0.0;
        R_E = 0.0;
        R_E_SQUARE = 0.0;
        C_E = 0.0;
        C_E_SQUARE = 0.0;

        M2 = 0.0;
        M4 = 0.0;
    } 

    } buffer;

    stringstream buff;

    buff.str("");
    buff << "result/" << algorithm << "-APPR-3D-" << N << "-" << statistic_max << "-THR-" << streams << ".DAT";
    cout << "Write to " << buff.str() << " ..." << endl;
    name_appr = new char[strlen(buff.str().c_str())];

    strcpy(name_appr, buff.str().c_str());

    buff.str("");
    buff << "result/U4-" << N << "-" << statistic_max << ".DAT";
    name_U4 = new char[strlen(buff.str().c_str())];

    strcpy(name_U4, buff.str().c_str());

    if(first_file_open && streams == 1)    {
        std::remove(name_appr);
        std::remove(name_U4);
        first_file_open = false;
    }

    std::cout << "Opening " << name_appr << "..." << std::endl;
    out_appr.open(name_appr, ios::app);
    std::cout << "Opening " << name_U4 << "..." << std::endl;
    out_U4.open(name_U4, ios::app);

    if(out_appr.bad() || out_U4.bad())   {
        cerr << "Cannot open files! Check permissions. " << endl;
    }
    out_appr.precision(4);

    for(int i = 0; i < streams; i++)   {

        buffer.M += data[i].sum[0];
        buffer.E += data[i].sum[1]; 
        buffer.C_v += data[i].sum[3];
        buffer.X += data[i].sum[2];

        buffer.R += statistic[i].relax;
        buffer.C += statistic[i].correlation;
        buffer.R_E += statistic[i].error_relax[0];
        buffer.R_E_SQUARE += statistic[i].error_relax[1];
        buffer.C_E += statistic[i].error_correlation[0];
        buffer.C_E_SQUARE += statistic[i].error_correlation[1];

        buffer.M2 += data[i].sum[2];
        buffer.M4 += data[i].sum[4];

    }
        
    M = buffer.M/(double)(mcs_max*(statistic_max));
    E = buffer.E/(double)(mcs_max*(statistic_max));
    C_v = (((buffer.C_v/(double)(mcs_max*(statistic_max))) - pow(buffer.E/(double)(mcs_max*(statistic_max)), 2.0))/T*T);
    X = (((buffer.X/(double)(mcs_max*(statistic_max))) - pow(buffer.M/(double)(mcs_max*(statistic_max)), 2.0))/T);

    M2 = buffer.M2/(double)(mcs_max*(statistic_max));
    M4 = buffer.M4/(double)(mcs_max*(statistic_max));

    std::cout << "\tM4 = " << M4 << "\tM2 = " << M2 << std::endl;

    // U4 = 1.0 - M4/(3*(M2*M2));
    U4 = 1.0 - M4/(3.0*M2*M2);
    // cout << "SR = " << statistic[thread_num].relax << ", SM = " <<  statistic_max << endl;

    relax = buffer.R/(double)(statistic_max);
    correlation = buffer.C/(double)(statistic_max);

    // cout << "statistic[thread_num].error_correlation[0] = " << statistic[thread_num].error_correlation[0] << \
    //         ", statistic[thread_num].error_correlation[1] = " << statistic[thread_num].error_correlation[1] << endl;

    error_relax = (sqrt(buffer.R_E_SQUARE/(double)(statistic_max) - buffer.R_E)/(double)(statistic_max));
    error_correlation = (sqrt(buffer.C_E_SQUARE - buffer.C_E))/(double)(statistic_max);
    
    cout << "( E = " << E << ", M = " << M << ", C_v = " << C_v << ", X = " << X << \
            ", R = " << relax << " ( " << error_relax << " ), Corr = " << correlation << \
            " ( " << error_correlation << " ),  U4 = " << U4 << ")" << endl;

    if(!out_appr.bad()) {

        out_appr << E << setw(10) << M << setw(10) << C_v << setw(10) << X << setw(10) << \
                    relax << setw(10) << error_relax << setw(10) << correlation << setw(10) << \
                    error_correlation << setw(10) << U4 << setw(10) << T << endl;

    }   else cerr << "Cannot open \"" << name_appr << "\", check permissions" << endl;

    if(!out_U4.bad())   {

        out_U4 << U4 << std::endl;

    }

    delete name_appr;

    out_appr.close();
    out_U4.close();
}

void Izing3D::SaveSpinConf(int thread_num) {

    for(int i = 0; i < N; i++)  {
        for(int j = 0; j < N; j++)  {
            for(int k = 0; k < N; k++)  {
            /* Запоминаем исходное состояние системы, когда система только начинает релаксировать. Нужно для расчета времени корелляции *
             * --------------------------------- */
                spin_old[i][j][k] = spin[i][j][k];
                data[thread_num].old_config += spin_old[i][j][k];    
            /* --------------------------------- */
            }
        } 
    }
}

void Izing3D::Metropolis(double _T)    {

    double dE, W, r;
    double exp_mcs;
        
    int spin_inverted;  // t_relax, M_begin, M_exp;
    int ci, cj, ck; // random var's
    long seed;
    // int thread_num;

    FormSourceConf();

    struct _time {
        float begin;
        float end;

        _time() {
            begin = 0.0;
            end = 0.0;
        }
    } times[streams];

    #pragma omp critical
    {
        cout << " Measurement... " << "( T = " << T << " )" << "[ N = " << N << ", Mcs = " << mcs_max << ", statistic = " << statistic_max << " ]" << endl;
    }

    for(int i = 0; i < streams; i++)    {

        data[i].sum[0] = 0.0;
        data[i].sum[1] = 0.0;
        data[i].sum[2] = 0.0;
        data[i].sum[3] = 0.0;

        statistic[i].relax = 0.0;
        statistic[i].correlation = 0.0;
        statistic[i].error_relax[0] = 0.0;
        statistic[i].error_relax[1] = 0.0;
        statistic[i].error_correlation[0] = 0.0;
        statistic[i].error_correlation[1] = 0.0;

    }

    omp_set_num_threads(streams);

    #pragma omp parallel for schedule(dynamic,1) private(ci,cj,ck,r,spin_inverted,dE) firstprivate(W)
    for(int start_conf = 0; start_conf < statistic_max; start_conf++)    { // собираем статистику, считаем конфигурации

        int thread_num = omp_get_thread_num();
        double first_mcs_magnet;

        #pragma omp critical
        {
            data[thread_num].relax = 0;
            data[thread_num].correlation = 0;

            data[thread_num].new_config = 0;
            data[thread_num].old_config = 0;

            trigger[thread_num].relax = true;
            trigger[thread_num].correlation = true;
            trigger[thread_num].relax_print = true;
        }

        times[thread_num].begin = omp_get_wtime();

        for(int mcs = 0; mcs < mcs_max; mcs++)    { // шагаем по Монте-Карло

            seed = time(0) + (start_conf%100) + (mcs%1000);

            CRandomMersenne Mersenne(seed);

            Mersenne.RandomInit(seed);  // используем "истинный" рандом

            for(int i=0;i<N;i++)    {
                for(int j=0;j<N;j++)    {
                	for(int k=0;k<N;k++)    {
                        ci = Mersenne.IRandom(0,N-1);
                        cj = Mersenne.IRandom(0,N-1);
                        ck = Mersenne.IRandom(0,N-1);

                        r = (Mersenne.IRandom(0,999))/1000.0;
                // Производим случайное пробное изменение в начальной конфигурации, переворачиваем спин
                        if(r < 0.5) { spin_inverted = spin[ci][cj][ck]; } else { spin_inverted = -spin[ci][cj][ck]; }
                // Считаем изменение энергии
                        dE = (spin[ci][cj][ck] - spin_inverted) * nearly_spins(ci, cj, ck);
                // W - вероятность перехода
                        if(dE<0) { spin[ci][cj][ck] = spin_inverted; } else { W=1.0/exp(dE/_T); }
                // принимаем новую конфигурацию если r <= W
                        if((r = Mersenne.IRandomX(0,999)/1000.0)<=W) { spin[ci][cj][ck] = spin_inverted; }
                	}
                }
            }
            
            magnetic(thread_num);
            // cout << "magnet = " << data[thread_num].magnet << endl;

            if(mcs == 0)    {
                first_mcs_magnet = data[thread_num].magnet;
            }

            if(trigger[thread_num].relax && data[thread_num].magnet < first_mcs_magnet/2.71)    {
                
                // cout << "mcs = " << mcs << endl;

                statistic[thread_num].relax += mcs;

                // cout << "statistic[thread_num].relax = " << statistic[thread_num].relax << endl;

                statistic[thread_num].error_relax[0] += mcs;
                statistic[thread_num].error_relax[1] += mcs*mcs;

                trigger[thread_num].relax = false; 

                SaveSpinConf(thread_num);
            }

            if(mcs > data[thread_num].relax*10)    {
            // ACF function
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    for(int k = 0; k < N; k++)  {
                        data[thread_num].new_config += spin_old[i][j][k]*spin[i][j][k];
                    }
                }
            }

            // cout << "old_config = " << data[thread_num].old_config << ", new config = " << data[thread_num].new_config << endl;

            if(trigger[thread_num].correlation && (double)data[thread_num].new_config > (double)data[thread_num].old_config/2.71)  {

                    // cout << "correlation = " << mcs << endl;

                    statistic[thread_num].correlation += mcs;

                    statistic[thread_num].error_correlation[0] += mcs;
                    statistic[thread_num].error_correlation[1] += mcs*mcs;

                    trigger[thread_num].correlation = false;
                }
            }


            data[thread_num].sum[0] += data[thread_num].magnet;
            data[thread_num].sum[2] += data[thread_num].magnet*data[thread_num].magnet;

            energy(thread_num);

            // cout << "energy = " << data[thread_num].energy << endl;

            data[thread_num].sum[1] += data[thread_num].energy;
            data[thread_num].sum[3] += data[thread_num].energy*data[thread_num].energy;
        }

        times[thread_num].end = omp_get_wtime();

        #pragma omp critical
        {
            cout << "Metropolis: stat №" << start_conf << " pass, T = " << _T << " thread = " << thread_num <<  \
                    ", iteration time = " << times[thread_num].end - times[thread_num].begin << "'s" << endl;
        }
    }

    WriteTdCharToFile();
    // WriteTimeToFile();

    // data[thread_num].relax_mcs = data[thread_num].relax_mcs_sum/statistic_max;
    // data[thread_num].correlation_mcs = data[thread_num].correlation_mcs_sum/statistic_max;
//    cout << "data[thread_num].relax_mcs_sum" << data[thread_num].relax_mcs_sum << "; approx data[thread_num].relax_mcs = " << data[thread_num].relax_mcs_sum/50.0 << endl;
    // data[thread_num].error_relax_mcs[0] = (sqrt(abs(data[thread_num].error_relax_mcs[2]/(double)statistic_max - data[thread_num].error_relax_mcs[1])/(double)statistic_max));
    // data[thread_num].error_correlation_mcs[0] = (sqrt(abs(data[thread_num].error_correlation_mcs[2]/(double)statistic_max - data[thread_num].error_correlation_mcs[1])/(double)statistic_max));

	// #pragma omp critical
 //    {
 //    	cout << " Relax = " << data[thread_num].relax_mcs << ", Err = " << data[thread_num].error_relax_mcs[0] << ", Corr. time = " << data[thread_num].correlation_mcs << ", Err = " << data[thread_num].error_correlation_mcs[0] << endl;
 //    	out_relax << fixed << T << " " << data[thread_num].relax_mcs << " " << data[thread_num].error_relax_mcs[0] << " " << data[thread_num].correlation_mcs << " " << data[thread_num].error_correlation_mcs[0] << endl;
 //    	out_appr << fixed << data[thread_num].sum[0]/static_cast<double>(mcs_max*statistic_max) << setw(10) << data[thread_num].sum[1]/static_cast<double>(mcs_max*statistic_max);
 //    	out_appr << fixed << setw(10) << ((data[thread_num].sum[2]/static_cast<double>(mcs_max*statistic_max)) - pow(data[thread_num].sum[0]/static_cast<double>(mcs_max*statistic_max), 2.0))/T*T;
 //    	out_appr << fixed << setw(10) << ((data[thread_num].sum[3]/static_cast<double>(mcs_max*statistic_max)) - pow(data[thread_num].sum[1]/static_cast<double>(mcs_max*statistic_max), 2.0))/T;
 //    	out_appr << fixed << setw(10) << T << endl;
 //    }

 //    out_appr.close();
 //    out_temp.close();
 //    out_relax.close();
}

void Izing3D::SwendsenWang(double _T)	{

	int ri, rj, rk;	//coordinates of random choose spin in grid
	int count_of_discovered;	// var, if equal N, then grid full of clusters
	double W,r;	// probability
	int cluster_num;
    int seed;
//	int count = 0;

	struct _checked_spins	{
//		int i;	// i coordinate
//		int j;	// j coordinate
		bool is_discovered;
		int cluster_num;

		_checked_spins()	{
			is_discovered = false;
			cluster_num = 0;
		}

		void clear()	{
			is_discovered = false;
			cluster_num = 0;
		}
	} checked_spins[N][N][N];

    struct _time {
        float begin;
        float end;

        _time() {
            begin = 0.0;
            end = 0.0;
        }
    } times;

    for(int i = 0; i < streams; i++)    {

        data[i].sum[0] = 0.0;
        data[i].sum[1] = 0.0;
        data[i].sum[2] = 0.0;
        data[i].sum[3] = 0.0;

        statistic[i].relax = 0.0;
        statistic[i].correlation = 0.0;
        statistic[i].error_relax[0] = 0.0;
        statistic[i].error_relax[1] = 0.0;
        statistic[i].error_correlation[0] = 0.0;
        statistic[i].error_correlation[1] = 0.0;
    }

	W = 1.0 - exp(-2.0/_T);
	cout << "W = " << W << endl;

    times.begin = 0;
    times.end = 0;

	FormSourceConf();

	for(int start_conf = 0; start_conf < statistic_max; start_conf++)	{

        times.begin = omp_get_wtime();

        srand(time(NULL));

        seed = 1 + rand() % 999;

        CRandomMersenne Mersenne(seed);

        Mersenne.RandomInit(seed);

        cout << "Swendsen_Wang: stat №" << start_conf << " pass, T = " << _T << ", iteration time = " << times.end - times.begin << "'s" << endl;
		// cout << "stat №" << start_conf << " pass" << endl;
		for(int mcs = 0; mcs < mcs_max; mcs++)	{

		cluster_num = 0;
		count_of_discovered = 0;

		for(int i = 0; i < N; i++)	{
			for(int j = 0; j < N; j++)	{
				for(int k = 0; k < N; k++)	{
					checked_spins[i][j][k].clear();
				}
			}
		}

		do	{

			ri = Mersenne.IRandomX(0,N-1);
			rj = Mersenne.IRandomX(0,N-1);
			rk = Mersenne.IRandomX(0,N-1);

		//check spin is marked? (used before check?)
			if(!checked_spins[ri][rj][rk].is_discovered)	{

				checked_spins[ri][rj][rk].is_discovered = true;
				cluster_num++;
				checked_spins[ri][rj][rk].cluster_num = cluster_num;
				// cout << "Cluster Num = " << checked_spins[ri][rj][rk].cluster_num << endl;

				count_of_discovered++;

				for(int i = ri-1; i <= ri+1; i++)	{
					for(int j = rj-1; j <= rj+1; j++)	{
						for(int k = rk-1; k <= rk+1; k++)	{
							if(i > 0 && i < N-1 && j > 0 && j < N-1 && k > 0 && k < N-1)	{
								if(checked_spins[i][j][k].is_discovered)	{
								// cout << "spin already discovered, cannot use " << endl;
								}	else	{
								// Form cluster and flip with probability 1/2
                                    r = (double)Mersenne.IRandomX(0,999)/1000.0;

									if(spin[i][j][k] == spin[ri][rj][rk] && r > W)	{

									   checked_spins[i][j][k].is_discovered = true;
									   checked_spins[i][j][k].cluster_num = checked_spins[ri][rj][rk].cluster_num;
									   count_of_discovered++;

									}
								}
							} else continue; 
						}
					}
				}

				if(Mersenne.IRandomX(0,999)/1000.0 > 0.5)	{

					spin[ri][rj][rk] = -spin[ri][rj][rk];

					for(int i = ri-1; i <= ri+1; i++)	{
						for(int j = rj-1; j <= rj+1; j++)	{
							for(int k = rk-1; k <= rk+1; k++)	{
								if(checked_spins[i][j][k].cluster_num == checked_spins[ri][rj][rk].cluster_num)
									spin[i][j][k] = -spin[i][j][k];
							}
						}
					}
				}

				if(count_of_discovered == N*N*N)	{
					// cout << "All elements in grid pass!" << endl;
					break;
//					goto m1;
				}
			}

			} 	while(1);
//m1:
			magnetic(thread_num);

            data[thread_num].sum[0] += data[thread_num].magnet;
            data[thread_num].sum[2] += data[thread_num].magnet*data[thread_num].magnet;

//			energy(thread_num;

//			data[thread_num].sum[1] += data[thread_num].energy;
//			cout << count++ << ":" << data[thread_num].energy << endl;
//			cout << "mcs №" << mcs << " pass" << endl;
		}

        times.end = omp_get_wtime();
	}
}

// void Izing3D::Recursion()   {

// }

void Izing3D::WolfStackBased(double _T)   {

    int central, state; // state - изначальное состояние спина, central - центральный спин, уже перевернутый (антипараллельный) изначальному
    int ci, cj, ck; // random coordinates
    // #pragma omp parallel default(none
    // int i, j, k; 

    int flip_max = 10;
    int iteration_count = 0;

    float W, r, seed;

    CRandomMersenne Mersenne(23567215);

    struct _time {

        float begin;
        float end;

        _time() {
            begin = 0.0;
            end = 0.0;
        }
    } times[streams];

    // int spin_pp[streams][N][N][N];

    int ****spin_pp;

    spin_pp = new int ***[N];
    for(int i=0;i<N;i++)    {
        spin_pp[i] = new int **[N];
        for(int j=0;j<N;j++)    {
            spin_pp[i][j] = new int *[N];
            for(int k=0;k<N;k++)    {
                spin_pp[i][j][k] = new int [N];
            }
        }
    }


    cout << "This!\n";

//     struct _checked_spin   {
// //      int i;  // i coordinate
// //      int j;  // j coordinate
//         bool with_defect; // Спин уже, добавленный в кластер

//         _checked_spin()    {
//             with_defect = false;
//         }

//         void clear()    {
//             with_defect = false;
//         }
//     } checked_spin[N][N][N];

    bool checked_spin[N][N][N];

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

    } data_pp[128];

    // struct _count   {

    //     int i;
    //     int j;
    //     int k;

    //     int number;

    // } cluster[N*N*N];

    // struct _spin_in_cluster  {
    //     int i;
    //     int j;
    //     int k;
    // } spin_in_cluster;

    stack <int> element_i;
    stack <int> element_j;
    stack <int> element_k;

    W = 1.0 - exp(-2.0/_T);

    // int pp_i = omp_get_thread_num();
    // int thread_num = omp_get_thread_num();

    #pragma omp critical
    {
        cout << "----------------------------------------------" << endl;
        cout << "Wolf: Measurement... " << "( T = " << _T << " ) [ N = " << N << ", mcs = " << mcs_max << ", statistic = " << statistic_max << " ]" << endl;
        cout << "W = " << W << endl;
    }

    // FormSourceConf();

    // add defects 
    for(int i=0;i<N;i++)    {
        for(int j=0;j<N;j++)    {
            for(int k=0;k<N;k++)    {
                if(Mersenne.IRandomX(0,999)/1000.0 < 0.4)   {
                    checked_spin[i][j][k]= true;
                }
            }
        }
    }

    for(int pp_i = 0; pp_i < streams; pp_i++)    {

        for(int i=0;i<N;i++)    {
            for(int j=0;j<N;j++)    {
                for(int k=0;k<N;k++)    {
                    // spin[i][j][k] = (Mersenne.IRandom(0,999)/1000.0 > 0.5) ? 1 : -1;
                    spin_pp[pp_i][i][j][k] = 1;
                }
            }
        }       

        data[pp_i].sum[0] = 0.0;
        data[pp_i].sum[1] = 0.0;
        data[pp_i].sum[2] = 0.0;
        data[pp_i].sum[3] = 0.0;
        data[pp_i].sum[4] = 0.0;

        // statistic[pp_i].relax = 0.0;
        // statistic[pp_i].correlation = 0.0;
        // statistic[pp_i].error_relax[0] = 0.0;
        // statistic[pp_i].error_relax[1] = 0.0;
        // statistic[pp_i].error_correlation[0] = 0.0;
        // statistic[pp_i].error_correlation[1] = 0.0;
    }

    // #pragma omp parallel for\
    // schedule(static,(statistic_max/streams))\
    // private(ci,cj,ck,r,seed,element_i,element_j,element_k,central,state)\
    // firstprivate(W)\
    // shared(spin_pp,flip_max,times,checked_spin)
    for(int start_conf = 0; start_conf < statistic_max; start_conf++)   {

        // #pragma omp critical
        // std::cout << "flip_max = " << flip_max << ", mcs_max = "\
                  // << mcs_max <<  ", statistic_max = " << statistic_max\
                  // << ", W = " << W << std::endl;

        // double first_mcs_magnet;

        thread_num = omp_get_thread_num();

        trigger[thread_num].relax = true;
        // trigger[thread_num].correlation = true;
        // data[thread_num].new_config = 0;
        // data[thread_num].old_config = 0;

        times[thread_num].begin = omp_get_wtime();

        // std::cout << "!\n";
        // for(int i = 0; i < N; i++)
        //     for(int j = 0; j < N; j++)  {
        //         std::cout << std::endl;
        //         for(int k = 0; k < N; k++)
        //             std::cout << spin[i][j][k] << "\t";
        //     }

        for(int mcs = 0; mcs < mcs_max; mcs++)  {

            for(int f = 0; f < flip_max; f++)   {

            // for(int i = 0; i < N; i++)  {
            //     for(int j = 0; j < N; j++)  {
            //         for(int k = 0; k < N; k++)  {
            //             checked_spins[i][j][k].clear();
            //         }
            //     }
            // }

            seed = 1 + rand() % 1000;

            CRandomMersenne Mersenne(seed);

            Mersenne.RandomInit(seed);

            ci = Mersenne.IRandom(0,N-1);
            cj = Mersenne.IRandom(0,N-1);
            ck = Mersenne.IRandom(0,N-1);

            // central = spin[ci][cj];

            element_i.push(ci);
            element_j.push(cj);
            element_k.push(ck);

            // cout << "element_i = " << element_i.top() << endl;

            spin_pp[thread_num][ci][cj][ck] = -spin_pp[thread_num][ci][cj][ck]; /* Инвертируем спин */

            do {

                if(element_i.empty() && element_j.empty() && element_k.empty())  {
                    break;
                };

                ci = element_i.top();
                cj = element_j.top();
                ck = element_k.top();

                element_i.pop();
                element_j.pop();
                element_k.pop();

                    // for X direction: 

                    if(ci!=0)    {

                        if(spin_pp[thread_num][ci-1][cj][ck]==-spin_pp[thread_num][ci][cj][ck] && checked_spin[ci-1][cj][ck] != true)   {

                            r = (double)Mersenne.IRandomX(0,999)/1000.0;

                            if(W > r)   {
                                spin_pp[thread_num][ci-1][cj][ck] = -spin_pp[thread_num][ci-1][cj][ck];

                                element_i.push(ci-1);
                                element_j.push(cj);
                                element_k.push(ck);
                            }
                        }   
                    }
                    else    {

                        if(spin_pp[thread_num][N-1][cj][ck]==-spin_pp[thread_num][ci][cj][ck] && checked_spin[ci-1][cj][ck]  != true)   {

                            r = (double)Mersenne.IRandomX(0,999)/1000.0;

                            if(W > r)   {

                                spin_pp[thread_num][N-1][cj][ck] = -spin_pp[thread_num][N-1][cj][ck];

                                element_i.push(N-1);
                                element_j.push(cj);
                                element_k.push(ck);
                            }
                        }
                    }

                    if(ci!=N-1) {

                        if(spin_pp[thread_num][ci+1][cj][ck]==-spin_pp[thread_num][ci][cj][ck] && checked_spin[ci-1][cj][ck]  != true) {

                            r = (double)Mersenne.IRandomX(0,1000)/1000.0;

                            if(W > r)   {

                                spin_pp[thread_num][ci+1][cj][ck]=-spin_pp[thread_num][ci+1][cj][ck];

                                element_i.push(ci+1);
                                element_j.push(cj);
                                element_k.push(ck);
                            }
                        }
                    }
                    else {

                        if(spin_pp[thread_num][0][cj][ck]==-spin_pp[thread_num][ci][cj][ck] && checked_spin[ci-1][cj][ck]  != true)  {

                            r = (double)Mersenne.IRandomX(0,1000)/1000.0;

                            if(W > r)   {

                                spin_pp[thread_num][0][cj][ck]=-spin_pp[thread_num][0][cj][ck];

                                element_i.push(0);
                                element_j.push(cj);
                                element_k.push(ck);
                            }
                        }
                    }

                    // for Y direction:

                    if(cj!=0)   {

                        if(spin_pp[thread_num][ci][cj-1][ck]==-spin_pp[thread_num][ci][cj][ck] && checked_spin[ci-1][cj][ck]  != true)   {

                            r = (double)Mersenne.IRandomX(0,1000)/1000.0;

                            if(W > r)   {

                                spin_pp[thread_num][ci][cj-1][ck]=-spin_pp[thread_num][ci][cj-1][ck];

                                element_i.push(ci);
                                element_j.push(cj-1);
                                element_k.push(ck);
                            }
                        }
                    }
                    else    {

                        if(spin_pp[thread_num][ci][N-1][ck]==-spin_pp[thread_num][ci][cj][ck] && checked_spin[ci-1][cj][ck]  != true)    {

                            r = (double)Mersenne.IRandomX(0,1000)/1000.0;

                            if(W > r)   {

                                spin_pp[thread_num][ci][N-1][ck]=-spin_pp[thread_num][ci][N-1][ck];

                                element_i.push(ci);
                                element_j.push(N-1);
                                element_k.push(ck);
                            }
                        }
                    }

                    if(cj!=N-1) {

                        if(spin_pp[thread_num][ci][cj+1][ck]==-spin_pp[thread_num][ci][cj][ck] && checked_spin[ci-1][cj][ck]  != true)   {

                            r = (double)Mersenne.IRandomX(0,1000)/1000.0;

                            if(W > r)   {

                                spin_pp[thread_num][ci][cj+1][ck]=-spin_pp[thread_num][ci][cj+1][ck];

                                element_i.push(ci);
                                element_j.push(cj+1);
                                element_k.push(ck);
                            }
                        }
                    }
                    else    {

                        if(spin_pp[thread_num][ci][0][ck]==-spin_pp[thread_num][ci][0][ck] && checked_spin[ci-1][cj][ck]  != true)   {

                            r = (double)Mersenne.IRandomX(0,1000)/1000.0;

                            if(W > r)   {

                                spin_pp[thread_num][ci][0][ck]=-spin_pp[thread_num][ci][0][ck];

                                element_i.push(ci);
                                element_j.push(0);
                                element_k.push(ck);
                            }
                        }
                    }

                    // for Z direction:

                    if(ck!=0)   {

                        if(spin_pp[thread_num][ci][cj][ck-1]==-spin_pp[thread_num][ci][cj][ck] && checked_spin[ci-1][cj][ck]  != true) {

                            r = (double)Mersenne.IRandomX(0,1000)/1000.0;

                            if(W > r)   {

                                spin_pp[thread_num][ci][cj][ck-1] = -spin_pp[thread_num][ci][cj][ck-1];

                                element_i.push(ci);
                                element_j.push(cj);
                                element_k.push(ck-1);
                            }
                        }
                    }
                    else    {

                        if(spin_pp[thread_num][ci][cj][N-1]==-spin_pp[thread_num][ci][cj][ck] && checked_spin[ci-1][cj][ck]  != true)   {

                            r = (double)Mersenne.IRandomX(0,1000)/1000.0;

                            if(W > r)   {

                                spin_pp[thread_num][ci][cj][N-1] = -spin_pp[thread_num][ci][cj][N-1];

                                element_i.push(ci);
                                element_j.push(cj);
                                element_k.push(N-1);
                            }
                        }
                    }

                    if(ck!=N-1) {

                        if(spin_pp[thread_num][ci][cj][ck+1]==-spin_pp[thread_num][ci][cj][ck] && checked_spin[ci-1][cj][ck]  != true)   {

                            r = (double)Mersenne.IRandomX(0,1000)/1000.0;

                            if(W > r)   {

                                spin_pp[thread_num][ci][cj][ck+1] = -spin_pp[thread_num][ci][cj][ck+1];

                                element_i.push(ci);
                                element_j.push(cj);
                                element_k.push(ck+1);
                            }
                        }
                    }
                    else    {

                        if(spin_pp[thread_num][ci][cj][0]==-spin_pp[thread_num][ci][cj][ck] && checked_spin[ci-1][cj][ck]  != true)  {

                            r = (double)Mersenne.IRandomX(0,1000)/1000.0;

                            if(W > r)   {

                                spin_pp[thread_num][ci][cj][0] = -spin_pp[thread_num][ci][cj][0];

                                element_i.push(ci);
                                element_j.push(cj);
                                element_k.push(0);
                            }
                        }
                    }

                    // cout << "cluster_size = " << element_i.size() << endl;

            } while(1);

        }

        // if(omp_get_thread_num()==0) {
        //     int count;
        //     std::cout << count++ << "!\n";
        // }

        // magnetic(thread_num);

        data[thread_num].magnetic_per_spin = 0.0;

        for(int i=0;i<N;i++)    {
            for(int j=0;j<N;j++)    {
                for(int k=0;k<N;k++)    {
                    data[thread_num].magnetic_per_spin += spin_pp[thread_num][i][j][k];
                }
            }
        }

        data[thread_num].magnet = abs(data[thread_num].magnetic_per_spin/(N*N*N));

        // #pragma omp critical
        // std::cout << "magnet = " << data[thread_num].magnet << std::endl;

        // if(trigger[thread_num].relax && mcs == 0)    {
        //     first_mcs_magnet = data[thread_num].magnet;
        //     cout << "first_mcs_magnet = " << first_mcs_magnet << std::endl;
        // }

        // if(trigger[thread_num].relax && (data[thread_num].magnet < (first_mcs_magnet/2.71)))    {
                
        //     // cout << "mcs = " << mcs << endl;

        //     statistic[thread_num].relax += mcs;

        //     cout << "current magnet = " << data[thread_num].magnet << ", relax = " << statistic[thread_num].relax << endl;

        //     statistic[thread_num].error_relax[0] += mcs;
        //     statistic[thread_num].error_relax[1] += mcs*mcs;

        //     trigger[thread_num].relax = false; 

        //     // SaveSpinConf(thread_num);;
        // }

        // if(mcs > data[thread_num].relax*10)    {
            // ACF function
            // for (int i = 0; i < N; i++) {
            //     for (int j = 0; j < N; j++) {
            //         for(int k = 0; k < N; k++)  {
            //             data[thread_num].new_config += spin_old[i][j][k]*spin[i][j][k];
            //         }
            //     }
            // }

            // if(trigger[thread_num].correlation && (double)data[thread_num].new_config > (double)data[thread_num].old_config/2.71)  {

            //     // cout << "correlation = " << mcs << endl;

            //     statistic[thread_num].correlation += mcs;

            //     statistic[thread_num].error_correlation[0] += mcs;
            //     statistic[thread_num].error_correlation[1] += mcs*mcs;

            //     trigger[thread_num].correlation = false;
            // }
        // }

        // if(mcs > 400)    {

        data[thread_num].sum[0] += data[thread_num].magnet;
        data[thread_num].sum[2] += data[thread_num].magnet*data[thread_num].magnet;
        data[thread_num].sum[4] += data[thread_num].magnet*data[thread_num].magnet*data[thread_num].magnet*data[thread_num].magnet;

        // if(omp_get_thread_num() == 0)   {
            // std::cout << data_pp[thread_num].sum[0] << std::endl;
            // cout << "end mcs " << mcs << endl;
        // }

        // energy(thread_num);

        data[thread_num].energy_per_spin = 0.0;

        for(int i=0;i<N;i++)    {
            for(int j=0;j<N;j++)    {
                for(int k=0;k<N;k++)    {

                    int neighbour_spins = 0;

                    if(i==0)   neighbour_spins= spin_pp[thread_num][N-1][j][k];
                    else       neighbour_spins= spin_pp[thread_num][i-1][j][k];
                    if(i==N-1) neighbour_spins+=spin_pp[thread_num][0][j][k];
                    else       neighbour_spins+=spin_pp[thread_num][i+1][j][k];
                    if(j==0)   neighbour_spins+=spin_pp[thread_num][i][N-1][k];
                    else       neighbour_spins+=spin_pp[thread_num][i][j-1][k];
                    if(j==N-1) neighbour_spins+=spin_pp[thread_num][i][0][k];
                    else       neighbour_spins+=spin_pp[thread_num][i][j+1][k];
                    if(k==0)   neighbour_spins+=spin_pp[thread_num][i][j][N-1];
                    else       neighbour_spins+=spin_pp[thread_num][i][j][k-1];
                    if(k==N-1) neighbour_spins+=spin_pp[thread_num][i][j][0];
                    else       neighbour_spins+=spin_pp[thread_num][i][j][k+1];

                    data[thread_num].energy_per_spin -= (spin_pp[thread_num][i][j][k] * neighbour_spins);
                }
            }
        }

        data[thread_num].energy = data[thread_num].energy_per_spin/(6*N*N*N);

        data[thread_num].sum[1] += data[thread_num].energy;
        data[thread_num].sum[3] += data[thread_num].energy*data[thread_num].energy;

        // }

        }

        times[thread_num].end = omp_get_wtime();

        #pragma omp critical
        { 
            cout << "Wolf: stat №" << iteration_count << " pass, T = " << _T << ", iteration time = " << times[thread_num].end - times[thread_num].begin << "'s" << endl;
            iteration_count++;
        }
    }

    // for(int pp_i = 0; pp_i < streams; pp_i++)   {

    //     data[pp_i].sum[0] = data_pp[pp_i].sum[0];
    //     data[pp_i].sum[1] = data_pp[pp_i].sum[1];
    //     data[pp_i].sum[2] = data_pp[pp_i].sum[2];
    //     data[pp_i].sum[3] = data_pp[pp_i].sum[3];
    //     data[pp_i].sum[4] = data_pp[pp_i].sum[4];

    // }

    #pragma omp ordered
    {
        WriteTdCharToFile();
    }

    delete [] spin_pp;
}

void Izing3D::WolfMassiveBased(double _T)   {
    int central, state; // state - изначальное состояние спина, central - центральный спин, уже перевернутый (антипараллельный изначальному)
    int ci, cj, ck; // random coordinates
    int i, j, k;    

    int flip_max = 1;
    int cluster_size[2] = {0,0};

    int cluster[N*N*N][3];
    int n_cl = 0;

    float W, r, seed;

    seed = 1 + rand() % 1000;

    CRandomMersenne Mersenne(seed);

    Mersenne.RandomInit(seed);

    struct _time {

        float begin;
        float end;

        _time() {

            begin = 0.0;
            end = 0.0;
        }
    } times[streams];

    struct _checked_spin   {
//      int i;  // i coordinate
//      int j;  // j coordinate
        bool with_defect; // Спин уже, добавленный в кластер

        _checked_spin()    {
            with_defect = false; 
        }

        void clear()    {
            with_defect = false;
        }
    } checked_spin[N][N][N];

    // struct _count   {

    //     int i;
    //     int j;
    //     int k;

    //     int number;

    // } cluster[N*N*N];

    struct _spin_in_cluster  {
        int i;
        int j;
        int k;
    } spin_in_cluster;

    stack <int> element_i;
    stack <int> element_j;
    stack <int> element_k;

    srand(time(NULL));

    cout << _T << endl;

    W = 1.0 - exp(-2.0/_T);

    cout << "----------------------------------------------" << endl;
    cout << "Wolf: Measurement... " << "( T = " << _T << " ) [ N = " << N << ", mcs = " << mcs_max << ", statistic = " << statistic_max << " ]" << endl;
    cout << "W = " << W << endl;

    FormSourceConf();

    // for(int i=0;i<N;i++)    {
    //     for(int j=0;j<N;j++)    {
    //         for(int k=0;k<N;k++)    {
    //             spin[i][j][k] = (Mersenne.IRandom(0,999)/1000.0 > 0.5) ? 1 : -1;
    //         }
    //     }
    // }

    // for(int i=0;i<N;i++)    {
    //     for(int j=0;j<N;j++)    {
    //         for(int k=0;k<N;k++)    {
    //             cout << spin[i][j][k] << setw(4);
    //         }
    //     }
    // }
    // adds defects 
    for(int i=0;i<N;i++)    {
        for(int j=0;j<N;j++)    {
            for(int k=0;k<N;k++)    {
                if(Mersenne.IRandomX(0,999)/1000.0 < 0.4)   {
                    checked_spin[i][j][k].with_defect = true;
                }
            }
        }
    }


    for(int i = 0; i < streams; i++)    {
        data[i].sum[0] = 0.0;
        data[i].sum[1] = 0.0;
        data[i].sum[2] = 0.0;
        data[i].sum[3] = 0.0;

        statistic[i].relax = 0.0;
        statistic[i].correlation = 0.0;
        statistic[i].error_relax[0] = 0.0;
        statistic[i].error_relax[1] = 0.0;
        statistic[i].error_correlation[0] = 0.0;
        statistic[i].error_correlation[1] = 0.0;
    }

    omp_set_num_threads(streams);

    #pragma omp parallel for schedule(static,10) private(ci,cj,ck,r,cluster_size,cluster) firstprivate(W,flip_max,n_cl)
    for(int start_conf = 0; start_conf < statistic_max; start_conf++)   {

        std::cout << "flip_max = " << flip_max << ", statistic_max = " << statistic_max << std::endl;

        int thread_num = omp_get_thread_num();

        trigger[thread_num].relax = true;
        trigger[thread_num].correlation = true;
        data[thread_num].new_config = 0;
        data[thread_num].old_config = 0;

        times[thread_num].begin = omp_get_wtime();

        for(int mcs = 0; mcs < mcs_max; mcs++)  {

            for(int f = 0; f < flip_max; f++)   {

            // for(int i = 0; i < N; i++)  {
            //     for(int j = 0; j < N; j++)  {
            //         for(int k = 0; k < N; k++)  {
            //             checked_spins[i][j][k].clear();
            //         }
            //     }
            // }


            cluster_size[0] = 0;
            cluster_size[1] = 0;

            ci = Mersenne.IRandom(0,N-1);
            cj = Mersenne.IRandom(0,N-1);
            ck = Mersenne.IRandom(0,N-1);  

            cluster[0][0] = ci;
            cluster[0][1] = cj;
            cluster[0][2] = ck;

            // central = spin[ci][cj];

            // element_i.push(ci);
            // element_j.push(cj);
            // element_k.push(ck);

            n_cl = 1;

            spin[ci][cj][ck] = -spin[ci][cj][ck]; /* Инвертируем спин */


            // setStreams(1);

            // cout << "streams=" << streams << endl;

            omp_set_num_threads(streams);

            do {

                if(n_cl == 0 || n_cl > N*N-2)  {
                    break;
                };

                // #pragma omp parallel if (n_cl > streams-1) shared(n_cl)
                // {
                    // cout << "ci = " << ci << ", cj = " << cj << ", ck = " << ck << endl;                    

                    // int size = element_i.size();
                    // #if (size > 1)
                    // #pragma omp sections
                    // {
                        // cout << "omp_get_num_threads = " << omp_get_num_threads() << endl;
                    // #pragma omp master
                    {
                        n_cl--;

                        ci = cluster[n_cl][0];
                        cj = cluster[n_cl][1];
                        ck = cluster[n_cl][2];
                    }
                            // ci = element_i.top();
                            // cj = element_j.top();
                            // ck = element_k.top();

                            // element_i.pop();
                            // element_j.pop();
                            // element_k.pop();

                            // cout << "omp_get_num_threads = " << omp_get_num_threads() << \
                                    // ", element_i = " << element_i.size() << endl;
                        // }

                        // N = N/streams;
                    // }
                    // #endif

                    // cout << "omp_get_num_threads = " << omp_get_num_threads() << endl;

                    // for X direction: 

                    if(ci!=0)    {

                        if(spin[ci-1][cj][ck]==-spin[ci][cj][ck] && checked_spin[ci-1][cj][ck].with_defect != true)   {

                            r = (double)Mersenne.IRandomX(0,1000)/1000.0;
                            // cout << r << endl;
                            if(W > r)   {

                                spin[ci-1][cj][ck] = -spin[ci-1][cj][ck];
                                // #pragma omp critical
                                // {

                                cluster[n_cl][0] = ci-1;
                                cluster[n_cl][1] = cj;
                                cluster[n_cl][2] = ck;

                                n_cl++;
                                // }
                            }
                        }   
                    }
                    else    {

                        if(spin[N-1][cj][ck]==-spin[ci][cj][ck] && checked_spin[N-1][cj][ck].with_defect != true)   {

                        r = (double)Mersenne.IRandomX(0,1000)/1000.0;

                            if(W > r)   {

                                spin[N-1][cj][ck] = -spin[N-1][cj][ck];

                                cluster[n_cl][0] = N-1;
                                cluster[n_cl][1] = cj;
                                cluster[n_cl][2] = ck;

                                n_cl++;
                            }
                        }
                    }

                    if(ci!=N-1) {

                        if(spin[ci+1][cj][ck]==-spin[ci][cj][ck] && checked_spin[ci+1][cj][ck].with_defect != true) {

                            r = (double)Mersenne.IRandomX(0,1000)/1000.0;

                            if(W > r)   {

                                spin[ci+1][cj][ck]=-spin[ci+1][cj][ck];

                                cluster[n_cl][0] = ci+1;
                                cluster[n_cl][1] = cj;
                                cluster[n_cl][2] = ck;

                                n_cl++;
                            }
                        }
                    }
                    else {

                        if(spin[0][cj][ck]==-spin[ci][cj][ck] && checked_spin[0][cj][ck].with_defect != true)  {

                            r = (double)Mersenne.IRandomX(0,1000)/1000.0;

                            if(W > r)   {

                                spin[0][cj][ck]=-spin[0][cj][ck];

                                cluster[n_cl][0] = 0;
                                cluster[n_cl][1] = cj;
                                cluster[n_cl][2] = ck;

                                n_cl++;
                            }
                        }
                    }

                    // for Y direction

                    if(cj!=0)   {

                        if(spin[ci][cj-1][ck]==-spin[ci][cj][ck] && checked_spin[ci][cj-1][ck].with_defect != true)   {

                            r = (double)Mersenne.IRandomX(0,1000)/1000.0;

                            if(W > r)   {

                                spin[ci][cj-1][ck]=-spin[ci][cj-1][ck];

                                cluster[n_cl][0] = ci;
                                cluster[n_cl][1] = cj-1;
                                cluster[n_cl][2] = ck;

                                n_cl++;
                            }
                        }
                    }
                    else    {

                        if(spin[ci][N-1][ck]==-spin[ci][cj][ck] && checked_spin[ci][N-1][ck].with_defect != true)    {

                            r = (double)Mersenne.IRandomX(0,1000)/1000.0;

                            if(W > r)   {

                                spin[ci][N-1][ck]=-spin[ci][N-1][ck];

                                cluster[n_cl][0] = ci;
                                cluster[n_cl][1] = N-1;
                                cluster[n_cl][2] = ck;

                                n_cl++;
                            }
                        }
                    }

                    if(cj!=N-1) {

                        if(spin[ci][cj+1][ck]==-spin[ci][cj][ck] && checked_spin[ci][cj+1][ck].with_defect != true)   {

                            r = (double)Mersenne.IRandomX(0,1000)/1000.0;

                            if(W > r)   {

                                spin[ci][cj+1][ck]=-spin[ci][cj+1][ck];

                                cluster[n_cl][0] = ci;
                                cluster[n_cl][1] = cj+1;
                                cluster[n_cl][2] = ck;

                                n_cl++;
                            }
                        }
                    }
                    else    {

                        if(spin[ci][0][ck]==-spin[ci][0][ck] && checked_spin[ci][0][ck].with_defect != true)   {

                            r = (double)Mersenne.IRandomX(0,1000)/1000.0;

                            if(W > r)   {

                                spin[ci][0][ck]=-spin[ci][0][ck];

                                cluster[n_cl][0] = ci;
                                cluster[n_cl][1] = 0;
                                cluster[n_cl][2] = ck;

                                n_cl++;
                            }
                        }
                    }

                    // for Z direction

                    if(ck!=0)   {

                        if(spin[ci][cj][ck-1]==-spin[ci][cj][ck] && checked_spin[ci][cj][ck-1].with_defect != true) {

                            r = (double)Mersenne.IRandomX(0,1000)/1000.0;

                            if(W > r)   {

                                spin[ci][cj][ck-1] = -spin[ci][cj][ck-1];

                                cluster[n_cl][0] = ci;
                                cluster[n_cl][1] = cj;
                                cluster[n_cl][2] = ck-1;

                                n_cl++;
                            }
                        }
                    }
                    else    {

                        if(spin[ci][cj][N-1]==-spin[ci][cj][ck] && checked_spin[ci][cj][N-1].with_defect != true)   {

                            r = (double)Mersenne.IRandomX(0,1000)/1000.0;

                            if(W > r)   {

                                spin[ci][cj][N-1] = -spin[ci][cj][N-1];

                                cluster[n_cl][0] = ci;
                                cluster[n_cl][1] = cj;
                                cluster[n_cl][2] = N-1;

                                n_cl++;
                            }
                        }
                    }

                    if(ck!=N-1) {

                        if(spin[ci][cj][ck+1]==-spin[ci][cj][ck] && checked_spin[ci][cj][ck+1].with_defect != true)   {

                            r = (double)Mersenne.IRandomX(0,1000)/1000.0;

                            if(W > r)   {

                                spin[ci][cj][ck+1] = -spin[ci][cj][ck+1];

                                cluster[n_cl][0] = ci;
                                cluster[n_cl][1] = cj;
                                cluster[n_cl][2] = ck+1;

                                n_cl++;
                            }
                        }
                    }
                    else    {

                        if(spin[ci][cj][0]==-spin[ci][cj][ck] && checked_spin[ci][cj][0].with_defect != true)  {

                            r = (double)Mersenne.IRandomX(0,1000)/1000.0;

                            if(W > r)   {

                                spin[ci][cj][0] = -spin[ci][cj][0];

                                cluster[n_cl][0] = ci;
                                cluster[n_cl][1] = cj;
                                cluster[n_cl][2] = 0;
                                // cout <<  << n_cl << endl;
                                n_cl++;
                            }
                        }
                    }

                    // cluster_size[0] = element_i.size();
                    // #pragma omp critical
                    // {
                        // cout << "cluster_size = " << element_i.size() << endl;
                    // }

            // } // omp parallel       
                    // cout << "cluster_size = " << n_cl << endl; 

            } while(1);

        }

        // cout << "while end" << endl;

        magnetic(thread_num);

        if(mcs == 0)    {
            first_mcs_magnet = data[thread_num].magnet;
        }

        if(trigger[thread_num].relax && data[thread_num].magnet < first_mcs_magnet/2.71)    {
                
            // cout << "mcs = " << mcs << endl;

            statistic[thread_num].relax += mcs;

            // cout << "statistic[thread_num].relax = " << statistic[thread_num].relax << endl;

            statistic[thread_num].error_relax[0] += mcs;
            statistic[thread_num].error_relax[1] += mcs*mcs;

            trigger[thread_num].relax = false; 

            SaveSpinConf(thread_num);
        }

        if(mcs > data[thread_num].relax*10)    {
            // ACF function
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    for(int k = 0; k < N; k++)  {
                        data[thread_num].new_config += spin_old[i][j][k]*spin[i][j][k];
                    }
                }
            }

            // cout << "old_config = " << data[thread_num].old_config << ", new config = " << data[thread_num].new_config << endl;

            if(trigger[thread_num].correlation && (double)data[thread_num].new_config > (double)data[thread_num].old_config/2.71)  {

                // cout << "correlation = " << mcs << endl;

                statistic[thread_num].correlation += mcs;

                statistic[thread_num].error_correlation[0] += mcs;
                statistic[thread_num].error_correlation[1] += mcs*mcs;

                trigger[thread_num].correlation = false;
            }
        }

        data[thread_num].sum[0] += data[thread_num].magnet;
        data[thread_num].sum[2] += data[thread_num].magnet*data[thread_num].magnet;
            // cout << "end mcs " << mcs << endl;

        energy(thread_num);

        data[thread_num].sum[1] += data[thread_num].energy;
        data[thread_num].sum[3] += data[thread_num].magnet*data[thread_num].magnet;

        }

        times[thread_num].end = omp_get_wtime();

        cout << "Wolf: stat №" << start_conf << " pass, T = " << _T << ", iteration time = " << times[thread_num].end - times[thread_num].begin << "'s" << endl;
    }

    WriteTdCharToFile();   
}

int Izing3D::PlotGraph(string _td_char)   {

    ofstream plot_file;
    ofstream datafile;
    
    char * name;
    char * datafile_name;
    
    string command_to_plot;

    stringstream buff;

    buff.str("");
    buff << "result/" << algorithm << "-APPR-3D-" << N << "-" << statistic_max << ".DAT";

    datafile_name = new char[strlen(buff.str().c_str())];

    strcpy(datafile_name, buff.str().c_str());

    cout << "datafile name = " << datafile_name << endl;

    datafile.open(datafile_name, ios::app);

    if(!datafile.bad())   {
        
        if(strcmp(_td_char.c_str(), "E") == 0)    {

            buff.str("");
            buff << "plot/" << algorithm << "-3D-" << _td_char << N << "-" << statistic_max << ".plot";

            name = new char [strlen(buff.str().c_str())];

            strcpy(name, buff.str().c_str());

            cout << "plot_file = " << name << endl;

            plot_file.open(name);

            if(plot_file.bad()) {
                cerr << "Cannot open \"" << name_appr << "\", check permissions" << endl;
            }

            plot_file << "#!/usr/bin/gnuplot -persist\n" << \
                            "set terminal jpeg font arial 12 size 800,600\n" << \
                            "set output \"" << name << ".jpg\"\n" << \
                            "set grid x y\n" << \
                            "set xlabel \"T\"\n" << \
                            "set ylabel \"mcs/s\"\n" << \
                            "plot \"" << datafile_name << "\" using 9:2 title \"izing3d-" << N << "\" with lines lt rgb \"red\"";   


            plot_file.close();
        }

        if(strcmp(_td_char.c_str(), "M") == 0)    {

            buff.str("");
            buff << "plot/" << algorithm << "-3D-" << _td_char << N << "-" << statistic_max << ".plot";

            name = new char [strlen(buff.str().c_str())];

            strcpy(name, buff.str().c_str());

            cout << "plot_file = " << name << endl;

            plot_file.open(name);

            if(plot_file.bad()) {
                cerr << "Cannot open \"" << name_appr << "\", check permissions" << endl;
            }

            plot_file << "#!/usr/bin/gnuplot -persist\n" << \
                            "set terminal jpeg font arial 12 size 800,600\n" << \
                            "set output \"" << name << ".jpg\"\n" << \
                            "set grid x y\n" << \
                            "set xlabel \"T\"\n" << \
                            "set ylabel \"mcs/s\"\n" << \
                            "plot \"" << datafile_name << "\" using 9:1 title \"izing3d-" << N << "\" with lines lt rgb \"red\"";   

            plot_file.close();
        }

        if(strcmp(_td_char.c_str(), "C_v") == 0)    {

            buff.str("");
            buff << "plot/" << algorithm << "-3D-" << _td_char << N << "-" << statistic_max << ".plot";

            name = new char [strlen(buff.str().c_str())];

            strcpy(name, buff.str().c_str());

            cout << "plot_file = " << name << endl;

            plot_file.open(name);

            if(plot_file.bad()) {
                cerr << "Cannot open \"" << name_appr << "\", check permissions" << endl;
            }

            plot_file << "#!/usr/bin/gnuplot -persist\n" << \
                            "set terminal jpeg font arial 12 size 800,600\n" << \
                            "set output \"" << name << ".jpg\"\n" << \
                            "set grid x y\n" << \
                            "set xlabel \"T\"\n" << \
                            "set ylabel \"mcs/s\"\n" << \
                            "plot \"" << datafile_name << "\" using 9:3 title \"izing3d-" << N << "\" with lines lt rgb \"red\"";   

            plot_file.close();
        }

        if(strcmp(_td_char.c_str(), "X") == 0)    {

            buff.str("");
            buff << "plot/" << algorithm << "-3D-" << _td_char << N << "-" << statistic_max << ".plot";

            name = new char [strlen(buff.str().c_str())];

            strcpy(name, buff.str().c_str());

            cout << "plot_file = " << name << endl;

            plot_file.open(name);

            if(plot_file.bad()) {
                cerr << "Cannot open \"" << name_appr << "\", check permissions" << endl;
            }

            plot_file << "#!/usr/bin/gnuplot -persist\n" << \
                            "set terminal jpeg font arial 12 size 800,600\n" << \
                            "set output \"" << name << ".jpg\"\n" << \
                            "set grid x y\n" << \
                            "set xlabel \"T\"\n" << \
                            "set ylabel \"mcs/s\"\n" << \
                            "plot \"" << datafile_name << "\" using 9:4 title \"izing3d-" << N << "\" with lines lt rgb \"red\"";   

            plot_file.close();
        }

        buff.str("");
        buff << "gnuplot " << name;

        command_to_plot = buff.str();

        cout << "command_to_plot = " << command_to_plot << endl;

        system(command_to_plot.c_str());

    }   else cerr << "File with TD data[thread_num] not exist!" << endl;

    buff.str("");

    delete name;
    delete datafile_name;

    datafile.close();
    plot_file.close();

    return 0;
}


int Izing3D::Start()    {

#ifndef _OPENMP
    cout << "Run on one thread!" << endl;
#else
    // omp_set_num_threads(streams);
    cout << "number of threads = " << streams << endl;
#endif

    if(strcmp("SW", algorithm.c_str()) == 0)   {
        for(int i = 0; i < T_div_max_index; i++)   {
            // cout << "T_div[0] = " << T_div[0] << endl;
            T = T_div[i];
            SwendsenWang(T_div[i]);
        }
    }

    if(strcmp("M", algorithm.c_str()) == 0)   {
        for(int i = 0; i < T_div_max_index; i++)   {
            // cout << "T_div[0] = " << T_div[0] << endl;
            T = T_div[i];
            Metropolis(T_div[i]);
        }
    }

    if(strcmp("WS", algorithm.c_str()) == 0)   {

    // #pragma omp parallel for\
    // schedule(static,(statistic_max/streams))\
    // private(ci,cj,ck,r,seed,element_i,element_j,element_k,central,state)\
    // firstprivate(W)\
    // shared(spin_pp,flip_max,times,checked_spin)

        // #pragma omp parallel for num_threads(streams) schedule(static,(statistic_max/streams))
        for(int i = 0; i < T_div_max_index; i++)   {
            T = T_div[i];
            WolfStackBased(T_div[i]);
        }
    }

    if(strcmp("WM", algorithm.c_str()) == 0)   {
        // setStreams(2);
        // #pragma omp parallel for num_threads(streams)
        for(int i = 0; i < T_div_max_index; i++)   {
            // #pragma omp critical    
            // {
            //    cout << "T_div[i] = " << T_div[i] << endl;    
            // }
            T = T_div[i];
            WolfMassiveBased(T_div[i]);
        }
    }
    return 0;
}