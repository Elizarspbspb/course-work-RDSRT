#include <stddef.h>
#include <stdlib.h>

extern void Signal(float Frequency, float Amplitude, float Phase, float Signalsin[]);

extern void Noise_Test(double Amplitude, double Noise[], int Frequency);

extern void AWGNoise(double Amplitude, double Noise[], int Frequency);

extern double AWGN_generator();

extern void AWGNoise1(double Amplitude, double Noise[], int Frequency);