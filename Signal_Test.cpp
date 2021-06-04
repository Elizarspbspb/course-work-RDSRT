#include <cmath>
#include <float.h>
#include "Signal_Test.h"
#include <iostream>
#include <time.h>

#define PI 3.14159265 

void Signal(float Frequency, float Amplitude, float Phase, float Signalsin[])
{
    float x = 0.0;
    float M = 1080;

    float FD = 0;
    float DT = 0;
    float T0 = 0;

    FD = M;
    T0 = 1 / FD;
    DT = T0 / M;
    float t = 0;
    for (int i = 0; i < M; i++)
    {
        x = sin(2*PI*T0*i*Frequency+Phase);
        x *= Amplitude;
        Signalsin[i] = x;
    }
}

void Noise_Test(double Amplitude, double Noise[], int Count)
{
    int a;
    a = 1.0 * Amplitude;
    double high = 0;
    double low = 0;
    double randNum = 1.0;
    for (int i0 = 0; i0 < Count; i0++)
    {
        do {
            high = 2.0 * (double(rand()) / double(RAND_MAX)) - 1.0; // диапазон от -1 до +1
            low = 2.0 * (double(rand()) / double(RAND_MAX)) - 1.0;
            randNum = high * high + low * low;
        } while (randNum >= 1.0 || randNum == 0.0);

        Noise[i0] = high;

        //Noise[i0] = (float(rand())) / RAND_MAX - (float(rand())) / RAND_MAX;
        //Noise[i0] = 2 * ((rand() / ((double)RAND_MAX)) - 0.5);
    }
}



/*void Noise_Test(double Amplitude, double Noise[], int Count)
{
    int a;
    a = 1.0 * Amplitude;
    double high = 0;
    double low = 0;
    double randNum = 1.0;
    for (int i0 = 0; i0 < Count; i0++)
    {
        do {
            high = 2.0 * (double(rand()) / double(RAND_MAX)) - 1.0; // диапазон от -1 до +1
            low = 2.0 * (double(rand()) / double(RAND_MAX)) - 1.0;
            //std::cout << "high = " << high << std::endl;
            //std::cout << "low = " << low << std::endl;
            randNum = high * high + low * low;
        } while (randNum >= 1.0 || randNum == 0.0);

        Noise[i0] = high;





        //Noise[i0] = (float(rand())) / RAND_MAX - (float(rand())) / RAND_MAX;
        //Noise[i0] = sqrt((1/2) * Noise[i0]);



        //Noise[i0] = (float(rand()))/ 32765 - (float(rand())) / 32765;
    }


    //Дисперсия шума
    //for (i0 = 0; i0 < 360 * Frequency; i0++)
    //for (i0 = 0; i0 < 360; i0++)
    //{
        //Noise[i0] = (rand() % Frequency - rand() % Frequency)*a; // допустимый
    //    Noise[i0] = (rand() % 4 - rand() % 4) * a; // отличный
        //Noise[i0] = (rand() % (Frequency - rand() % Frequency) - rand() % (Frequency - rand() % Frequency)) * a; // плохой шум, все утонуло
    //}
}*/

















#define PI 3.1415926536

float random() {
    srand(rand() * cos((double)rand()) * (int)time(NULL));
    return (float)rand() / RAND_MAX;
}

//double AWGN_generator() /* Генерация аддитивного белого гауссовского шума с нулевым средним и стандартным отклонением, равным 1. */
//void AWGNoise(double Amplitude, double Noise[], int Count, double my)  
void AWGNoise1(double Amplitude, double Noise[], int Count)
{
    int EbNo_dB = -5;		// SNR

    double polar = 1.0, rsquared = 1.0, var1 = 1.0, var2 = 1.0;
    double mu = 0.0; 
    double sigma = 1.0;
    static bool deviateAvailable = false;	//	flag
    static float storedDeviate;			//	deviate from previous calculation
    double dist = 1.0, angle = 1.0;
    for (int j = 0; j < Count; j++)
    {
        /*static bool b_cached = false;
        float randNum = 1.0, a = 1.0, b=1.0, tmp = 1.0, res;
        float sigmaI = (float)sqrt(2 * 2 * (float)pow(10.0, (EbNo_dB / 10)));
        do {
            //randNum = random();
            //a = randNum * 2 - 1;
            //randNum = random();
            //b = randNum * 2 - 1;
            //randNum = a * a + b * b;
            a = 2.0 * (double(rand()) / double(RAND_MAX)) - 1.0; // диапазон от -1 до +1
            b = 2.0 * (double(rand()) / double(RAND_MAX)) - 1.0;
            randNum = a * a + b * b;
        } while (randNum >= 1.0 || randNum == 0);
        tmp = sqrt((-2 * log(randNum)) / randNum);
        if (b_cached) {
            res = b * tmp;
        }
        else {
            res = a * tmp;
        }
        b_cached = !b_cached;
        Noise[j] = res / sigmaI;
        //return res / sigmaI;*/


        if (!deviateAvailable) {
            /*do {
                var1 = 2.0 * (double(rand()) / double(RAND_MAX)) - 1.0; // диапазон от -1 до +1
                var2 = 2.0 * (double(rand()) / double(RAND_MAX)) - 1.0;

                rsquared = var1 * var1 + var2 * var2;
                std::cout << "VAR1 = " << var1 << std::endl;
                std::cout << "VAR2 = " << var2 << std::endl;
            } while (rsquared >= 1.0 || rsquared == 0.0);

            //polar = sqrt(-2.0 * log(rsquared) / rsquared);

            //	store first deviate and set flag
            storedDeviate = var1 * polar;
            deviateAvailable = true;
            double my;
            //	return second deviate
            Noise[j] = var2 * polar * sigma + mu;
            //return my;
            //return(var2 * polar * sigma + mu);
            rsquared = 0;*/

            //dist = sqrt(-2.0 * log(double(rand()) / double(RAND_MAX)));
            angle = 2.0 * PI * (double(rand()) / double(RAND_MAX));
            //	calculate and store first deviate and set flag
            storedDeviate = dist * cos(angle);
            deviateAvailable = true;
            //	calcaulate return second deviate
            Noise[j] =  (dist * sin(angle) * sigma + mu);
        }
        else {
            deviateAvailable = false;
            Noise[j] = storedDeviate * sigma + mu;
            //return storedDeviate * sigma + mu;
        }
    }








    /*double temp1;
    double temp2;
    int p = 1;

    for (int i = 0; i < Count; i++)
    {
        while (p > 0)
        {
            temp2 = (float(rand()) / ((double)RAND_MAX)); // функция rand() генерирует целое число между 0 и  RAND_MAX,
            std::cout << temp2 << std::endl;
            if (temp2 == 0) // temp2 >= (RAND_MAX / 2)
                p = 1;
            else            // temp2 < (RAND_MAX / 2)
                p = -1;
        }
        temp1 = cos((2.0 * (double)PI*i) * float(rand()) / ((double)RAND_MAX));
        //temp1 = cos((2.0 * (double)PI) * float(rand()) / ((double)RAND_MAX));
        Noise[i] = (sqrt(-2.0 * log(temp2)) * temp1);
        //Noise[i] = temp1;
        p = 1;
    }*/
}

/*float gauss_rand(float mean, float stdev)
{
    int i;
    const int ORDER = 2 * 12; /* 12,24,36 etc. due to del^2/12 
    const double dev_norm = 1.4142136; /* sqrt(ORDER/12) 
    double rndno;

    rndno = -(ORDER >> 1);
    for (i = 0; i < ORDER; i++) 
    {
        rndno += (double)(rand() / (RAND_MAX + 1.0));
    }

    rndno *= stdev / dev_norm;
    rndno += mean;
    return((float)rndno);
}

void add_gaussian_noise(float** orig, int Ni, int Nj, float** noisy, float mean, float stdev)
{
    int i, j;
    static int kilroy = 0;
    unsigned int seed;

    if (!kilroy) {
        kilroy = 1;

        seed = (unsigned)time(NULL);

        // uncomment for the same noise process
        //  seed=0;
        srand(seed);
    }

    for (i = 0; i < Ni; i++)
        for (j = 0; j < Nj; j++)
            noisy[i][j] = orig[i][j] + gauss_rand(mean, stdev);
}*/



/*#define PI 3.1415926536
double AWGN_generator()
{// Генерация аддитивного белого гауссовского шума с нулевым средним и стандартным отклонением, равным 1. 

    double temp1;
    double temp2;
    double result;
    int p;

    p = 1;

    while (p > 0)
    {
        temp2 = (rand() / ((double)RAND_MAX)); /* функция rand() генерирует
                                                       целое число между 0 и  RAND_MAX,
                                                       которое определено в stdlib.h.
                                                   
        std::cout << "  Temp2 = " << temp2 << std::endl;
        if (temp2 == 0)
        {// temp2 >= (RAND_MAX / 2)
            std::cout << "      1       " << std::endl;
            p = 1;
        }// конец if
        else
        {// temp2 < (RAND_MAX / 2)
            std::cout << "      B       " << std::endl;
            p = -1;
        }// конец else

    }// конец while()

    temp1 = cos((2.0 * (double)PI) * rand() / ((double)RAND_MAX));
    std::cout << "  Temp1 = " << temp1 << std::endl;
    result = sqrt(-2.0 * log(temp2)) * temp1;
    std::cout << "  Result = " << result << std::endl << std::endl;

    return result;        // возвращаем сгенерированный сэмпл

}// конец AWGN_generator()*/
