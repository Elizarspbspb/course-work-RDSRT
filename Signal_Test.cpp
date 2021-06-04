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