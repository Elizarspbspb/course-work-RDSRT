#include "Signal_Test.h"
#include "Matlab_all.h"
#include <iostream>
#include <cmath>                                  // для функции exp
#include <fstream>

#include <iostream>
#include <iomanip>
#include <cstring>
#include <bitset>

#define PI 3.14159265 
using namespace std;

int printInBinary(const unsigned char val)
{
	int mas[32];
	for (int i = 7; i >= 0; --i)
		if (val & (1 << i))
		{
			cout << "1";
			mas[i] = 1;
		}
		else
		{
			cout << "0";
			mas[i] = 0;
		}
	return *mas;
}


int main(int, const char* const [])
{
	std::ofstream sinus1("SIN1.txt", std::ios_base::out);
	std::ofstream sinus2("SIN2.txt", std::ios_base::out);
	std::ofstream sinus3("SIN3.txt", std::ios_base::out);
	std::ofstream sinus4("SIN4.txt", std::ios_base::out);
	int Count = 1080;

	int Amplitude1 = 10;
	int Frequency1 = 10;

	int Amplitude2 = 10;
	int Frequency2 = 50;

	int Amplitude3 = 10;
	int Frequency3 = 56;

	float* Signalsin1 = new float[Count];
	Signal(Frequency1, Amplitude1, 0, Signalsin1);
	for (int i = 0; i < Count; i++)
		sinus1 << Signalsin1[i] << "		" << i << "\n";

	float* Signalsin2 = new float[Count];
	Signal(Frequency2, Amplitude2, 0, Signalsin2);
	for (int i = 0; i < Count; i++)
		sinus2 << Signalsin2[i] << "		" << i << "\n";

	float* Signalsin3 = new float[Count];
	Signal(Frequency3, Amplitude3, 0, Signalsin3);
	for (int i = 0; i < Count; i++)
		sinus3 << Signalsin3[i] << "		" << i << "\n";

	////////////////////////////////////////////////////////////////////////////

		// Шум
	std::ofstream noise("NOISE.txt", std::ios_base::out);
	double* Noise = new double[Count];
	int Amplitude4 = (Amplitude1 + Amplitude2 + Amplitude3) / 3;
	std::cout << std::endl << " Noise " << std::endl;
	Noise_Test(Amplitude4, Noise, Count);
	double my = 0.0;
	for (int i = 0; i < Count; i++)
	{
		noise << Noise[i] << "		" << i << "\n";
	}


	////////////////////////////////////////////////////////////////////////////

	std::ofstream signal("SIGNAL.txt", std::ios_base::out);
	float* Signal = new float[Count];
	std::cout << std::endl << " Signal " << std::endl;
	for (int i = 0; i < Count; i++)
	{
		//Signal[i] = Noise[i] + Signalsin1[i] + Signalsin2[i] + Signalsin3[i];
		Signal[i] = Signalsin1[i] + Signalsin2[i] + Signalsin3[i];
		signal << Signal[i] << "		" << i << "\n";
	}

	const int CountTest = 1080;

	float mas_int[CountTest];
	float mas_work[CountTest];
	float power = 0.0;
	std::ofstream dftsn("DFTSignal.txt", std::ios_base::out);
	std::ofstream FCHt("FCH.txt", std::ios_base::out);
	std::ofstream ACHt("ACH.txt", std::ios_base::out);
	std::ofstream FCHtlg("FCHLG.txt", std::ios_base::out);
	std::ofstream ACHtlg("ACHLG.txt", std::ios_base::out);
	std::cout << std::endl << " DFTSignal " << std::endl;
	
	float N_SN = CountTest;
	float* DFT_SN = new float[Count];

	float real_SN = 0.0;
	float mnim_SN = 0.0;
	float Sum1_SN = 0.0;
	float Sum2_SN = 0.0;
	float param_SN = 0.0;
	float st = 0;
	int sr = 0;
	float* FCH = new float[Count];
	float* ACH = new float[Count];
	float* FCHlg = new float[Count];
	float* ACHlg = new float[Count];

		for (int j = 0; j< CountTest; j++) // Прием данных в буфер
			mas_int[j] = Signal[j];
		for (int i = 0; i < N_SN; i++) // Обработка данных из буфера
		{
			Sum1_SN = 0;
			Sum2_SN = 0;
			for (int j = 0; j < N_SN; j++)
			{
				param_SN = -2.0 * PI * j * i / N_SN;
				real_SN = mas_int[j] * cos(param_SN);
				mnim_SN = mas_int[j] * sin(param_SN);
				Sum1_SN = Sum1_SN + real_SN;
				Sum2_SN = Sum2_SN + mnim_SN;
			}
			param_SN = 2.0 * PI  * i / N_SN;
			real_SN = Signal[i] * cos(param_SN);
			mnim_SN = Signal[i] * sin(param_SN);

			//FCH[i] = atan(Sum2_SN / Sum1_SN); // ФЧХ
			FCH[i] = atan(mnim_SN / real_SN); // ФЧХ
			FCHt << FCH[i] << "		" << i << "\n";
			//FCHt << FCH[i] << "		" << i << "\n";

			ACH[i] = (sqrt(real_SN * real_SN + mnim_SN * mnim_SN)); // АЧХ
			ACHt << ACH[i] << "		" << i << "\n";
			//ACHt << ACH[i] << "		" << i << "\n";

			mas_work[i] = ((Sum1_SN * Sum1_SN + Sum2_SN * Sum2_SN)); // Спектр
			dftsn << mas_work[i] << "		" << i << "\n";

			ACHlg[i] = 20 * log10(mas_int[i]);
			ACHtlg << ACHlg[i] << "		" << i << "\n";

			power = power + ((Sum1_SN * Sum1_SN + Sum2_SN * Sum2_SN) / (CountTest)); // мощность
			//power = power + ((real_SN * real_SN + mnim_SN * mnim_SN) / (CountTest)); // мощность
		}
	std::cout << std::endl << " Power signal = " << power << std::endl;




	float mas_intN[CountTest];
	float mas_workN[CountTest];
	std::ofstream dftnoise("DFT_Noise.txt", std::ios_base::out);
	std::cout << std::endl << " DFT_Noise " << std::endl;

	float N_SNN = CountTest;
	float* DFT_SNN = new float[Count];

	float real_SNN = 0.0;
	float mnim_SNN = 0.0;
	float Sum1_SNN = 0.0;
	float Sum2_SNN = 0.0;
	float param_SNN = 0.0;
	float stN = 0;
	int srN = 0;

	for (int j = 0; j < CountTest; j++) // Прием данных в буфер
		mas_intN[j] = Noise[j];
	for (int i = 0; i < N_SNN; i++) // Обработка данных из буфера
	{
		Sum1_SNN = 0;
		Sum2_SNN = 0;
		for (int j = 0; j < N_SNN; j++)
		{
			param_SNN = -2.0 * PI * j * i / N_SNN;
			real_SNN = mas_intN[j] * cos(param_SNN);
			mnim_SNN = mas_intN[j] * sin(param_SNN);
			Sum1_SNN = Sum1_SNN + real_SNN;
			Sum2_SNN = Sum2_SNN - mnim_SNN;
		}
		param_SN = 2.0 * PI * i / N_SN;
		real_SN = Noise[i] * cos(param_SN);
		mnim_SN = Noise[i] * sin(param_SN);

		//mas_workN[i] = 20 * log10(sqrt(Sum1_SNN * Sum1_SNN + Sum2_SNN * Sum2_SNN)) / (CountTest / 2); // Работает для шума
		//mas_workN[i] = 20 * log10(sqrt(mas_intN[i])) / CountTest; // Работает для шума
		//mas_workN[i] = 20 * log10(sqrt((Sum1_SNN * Sum1_SNN + Sum2_SNN * Sum2_SNN))) / CountTest; // Работает для шума


		mas_workN[i] = ((Sum1_SNN * Sum1_SNN + Sum2_SNN * Sum2_SNN));
		mas_workN[i] = 20 * log10(sqrt(mas_workN[i])) / CountTest;

		//mas_workN[i] = ((real_SN * real_SN + mnim_SN * mnim_SN));
		//mas_workN[i] = 20 * log10(sqrt(mas_workN[i])) / CountTest;

		//mas_workN[i] = 20 * log10(sqrt(Noise[i]));

		//dftnoise << mas_workN[i] << "		" << log10(i) << "\n";
		dftnoise << mas_workN[i] << "		" << i << "\n";
	}



	////////////////////////////// ФИЛЬТР ////////////////////
	float z1 = 0, z2 = 0, z3 = 0, z4 = 0;
	float sum1 = 0, sum2 = 0;
	float* filtr = new float[Count + 2];
	float* exit = new float[Count + 2];
	int flag1 = 0, flag2 = 0, flag3 = 0, flag4 = 0;

	std::cout << std::endl;
	for (int i = 0; i < Count + 2; i++)
	{
		filtr[i] = 0.99 * Signal[i];
		sum1 = filtr[i] + z1 * (-1.902) + z2 * (1);
		z2 = z1;
		z1 = filtr[i];
		exit[i] = sum1;
		flag1 = 1;
		//std::cout << "Si" << "\t\tZ1" << "\t\tZ2" << "\t\tSum1" << "\t\tSum2" << "\t\tZ3" << "\t\tZ4" << endl << Signal[i] << "\t\t" << z1 << "\t" << z2 << "\t" << sum1 << "\t" << z3 + sum1 + sum2 << "\t" << z3 + sum1 + sum2 << "\t" << sum2 << endl << endl;
	}
	std::cout << std::endl;


	float zn1 = 0, zn2 = 0, zn3 = 0, zn4 = 0;
	float sumn1 = 0, sumn2 = 0;
	float* filtrn = new float[Count + 2];
	float* exitn = new float[Count + 2];
	int flagn1 = 0, flagn2 = 0, flagn3 = 0, flagn4 = 0;

	for (int i = 0; i < Count + 2; i++)
	{
		filtrn[i] = 0.99 * Noise[i];
		sumn1 = filtrn[i] + zn1 * (-1.902) + zn2 * 1;
		zn2 = zn1;
		zn1 = filtrn[i];
		exitn[i] = sumn1;
		flagn1 = 1;
		//std::cout << "Signal[i]" << "\tZ1" << "\tZ2" << "\tSum1" << "\tSum2" << "\tZ3" << "\tZ4" << endl << Signal[i] << "\t" << z1 << "\t" << z2 << "\t" << sum1 << "\t" << z3 + sum1 + sum2 << "\t" << z3 + sum1 + sum2 << "\t" << sum2 << endl << endl;
	}
	std::cout << std::endl;



	//////////////////////////////////  DFT _FILTR   /////////////////////////////////////////
	std::ofstream FCHtFiltr("FCHFiltr.txt", std::ios_base::out);
	std::ofstream ACHtFiltr("ACHFiltr.txt", std::ios_base::out);
	std::ofstream FCHtlgFiltr("FCHLGFiltr.txt", std::ios_base::out);
	std::ofstream ACHtlgFiltr("ACHLGFiltr.txt", std::ios_base::out);
	std::ofstream filtr1("Filtr.txt", std::ios_base::out);
	float* FCHFiltr = new float[Count];
	float* ACHFiltr = new float[Count];
	float* FCHlgFiltr = new float[Count];
	float* ACHlgFiltr = new float[Count];

	float mas_int_F[CountTest];
	float mas_work_F[CountTest];

	float N_SN_F = CountTest;
	float* DFT_SN_F = new float[Count];

	float real_SN_F = 0.0;
	float mnim_SN_F = 0.0;
	float Sum1_SN_F = 0.0;
	float Sum2_SN_F = 0.0;
	float param_SN_F = 0.0;
	int st_F = 0;
	int sr_F = 0;
	power = 0.0;
	for (int i = 0; i < N_SN_F; i++) // Обработка данных из буфера
	{
		Sum1_SN_F = 0;
		Sum2_SN_F = 0;
		for (int j = 0; j < N_SN_F; j++)
		{

			param_SN_F = -2.0 * PI * j * i / N_SN_F;
			real_SN_F = exit[j] * cos(param_SN_F);
			mnim_SN_F = exit[j] * sin(param_SN_F);
			Sum1_SN_F = Sum1_SN_F + real_SN_F;
			Sum2_SN_F = Sum2_SN_F + mnim_SN_F;
		}
		param_SN_F = 2.0 * PI * i / N_SN_F;
		real_SN_F = exit[i] * cos(param_SN_F);
		mnim_SN_F = exit[i] * sin(param_SN_F);

		FCHFiltr[i] = atan(mnim_SN_F / real_SN_F); // ФЧХ
		//FCHtFiltr << FCHFiltr[i] << "		" << i << "\n";
		FCHtFiltr << FCHFiltr[i] << "		" << i << "\n";

		ACHFiltr[i] = (sqrt(real_SN_F * real_SN_F + mnim_SN_F * mnim_SN_F)); // АЧХ
		//ACHtFiltr << ACHFiltr[i] << "		" << i << "\n";
		ACHtFiltr << ACHFiltr[i] << "		" << i << "\n";

		mas_work_F[i] = ((Sum1_SN_F * Sum1_SN_F + Sum2_SN_F * Sum2_SN_F)); // Спектр
		filtr1 << mas_work_F[i] << "		" << i << "\n";

		ACHlgFiltr[i] = 20 * log10(exit[i]);
		ACHtlgFiltr << ACHlgFiltr[i] << "		" << i << "\n";

		power = power + ((Sum1_SN_F * Sum1_SN_F + Sum2_SN_F * Sum2_SN_F) / (CountTest)); // мощность
		//power = power + ((real_SN_F * real_SN_F + mnim_SN_F * mnim_SN_F) / (CountTest)); // мощность
	}
	std::cout << std::endl << " Power signal after filtr111 = " << power << std::endl;




	std::ofstream filtr_noise("Filtr_Noise.txt", std::ios_base::out);
	float mas_int_FN[CountTest];
	float mas_work_FN[CountTest];

	float N_SN_FN = CountTest;
	float* DFT_SN_FN = new float[Count];

	float real_SN_FN = 0.0;
	float mnim_SN_FN = 0.0;
	float Sum1_SN_FN = 0.0;
	float Sum2_SN_FN = 0.0;
	float param_SN_FN = 0.0;
	int st_FN = 0;
	int sr_FN = 0;
	for (int i = 0; i < N_SN_FN; i++) // Обработка данных из буфера
	{
		Sum1_SN_FN = 0;
		Sum2_SN_FN = 0;
		for (int j = 0; j < N_SN_FN; j++)
		{

			param_SN_FN = -2.0 * PI * j * i / N_SN_FN;
			real_SN_FN = exitn[j] * cos(param_SN_FN);
			mnim_SN_FN = exitn[j] * sin(param_SN_FN);
			Sum1_SN_FN = Sum1_SN_FN + real_SN_FN;
			Sum2_SN_FN = Sum2_SN_FN - mnim_SN_FN;
		}
		//mas_work_FN[i] = 20*log10((Sum1_SN_FN * Sum1_SN_FN + Sum2_SN_FN * Sum2_SN_FN)) / (CountTest / 2); // Это хорошо для шума
		mas_work_FN[i] = ((Sum1_SN_FN * Sum1_SN_FN + Sum2_SN_FN * Sum2_SN_FN));
		mas_work_FN[i] = 20 * log10(sqrt(mas_work_FN[i])) / CountTest;

		filtr_noise << mas_work_FN[i] << "		" << i << "\n";
	}








	//////////////////////////////////     16  &  8  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

	float index = 0.0;
	for (int i = 0; i < Count; i++)
	{
		if (index < Signal[i])
		{
			index = Signal[i];
		}
	}
	//std::cout << "Max Amplitude = " << index << " ";

	std::ofstream signal16("SIGNAL16.txt", std::ios_base::out);
	float* Signal16 = new float[Count];
	std::cout << std::endl << " Signal16 " << std::endl;
	for (int i = 0; i < Count; i++)
	{
		Signal16[i] = Signal[i]/ index;

			Signal16[i] = Signal16[i] * (2^16) - 1;
			Signal16[i] = Signal16[i] / (2 ^ 16);
			Signal16[i] = Signal16[i] * (index);

		signal16 << Signal16[i] << "		" << i << "\n";
	}


	std::ofstream signal8("SIGNAL8.txt", std::ios_base::out);
	float* Signal8 = new float[Count];
	std::cout << std::endl << " Signal8 " << std::endl;
	for (int i = 0; i < Count; i++)
	{
		Signal8[i] = Signal[i] / index;

			Signal8[i] = Signal8[i] * (2 ^ 8) - 1;
			Signal8[i] = Signal8[i] / (2 ^ 8);
			Signal8[i] = Signal8[i] * (index);

		signal8 << Signal8[i] << "		" << i << "\n";
	}



	/*std::ofstream signal8("SIGNAL8.txt", std::ios_base::out);
	float* Signal8 = new float[Count];
	std::cout << std::endl << " Signal8 " << std::endl;
	for (int i = 0; i < Count; i++)
	{
		//Signal[i] = Noise[i] + Signalsin1[i] + Signalsin2[i] + Signalsin3[i];
		Signal8[i] = Signal[i] / index;
		//Signal8[i] = Signal[i];
		if (Signal8[i] >= 0)
		{
			Signal8[i] = Signal8[i] * (2 ^ 8) - 1;
			Signal8[i] = Signal8[i] / (2 ^ 8);
			Signal8[i] = Signal8[i] * (index);
		}
		else if (Signal8[i] < 0)
		{
			Signal8[i] = Signal8[i] * (2 ^ 8) - 1;
			Signal8[i] = Signal8[i] / (2 ^ 8);
			Signal8[i] = Signal8[i] * (index);
		}
		signal8 << Signal8[i] << "		" << i << "\n";
	}*/

	/*std::ofstream signal24("SIGNAL24.txt", std::ios_base::out);
	float* Signal24 = new float[Count];
	std::cout << std::endl << " Signal24 " << std::endl;
	for (int i = 0; i < Count; i++)
	{
		//Signal[i] = Noise[i] + Signalsin1[i] + Signalsin2[i] + Signalsin3[i];
		Signal24[i] = Signal[i] / index;
		//Signal24[i] = Signal[i];
		if (Signal24[i] >= 0)
		{
			Signal24[i] = Signal24[i] * (2 ^ 32) - 1;
			Signal24[i] = Signal24[i] / (2 ^ 32);
			Signal24[i] = Signal24[i] * (index);
		}
		else if (Signal24[i] < 0.0)
		{
			Signal24[i] = Signal24[i] * (2 ^ 32) - 1;
			Signal24[i] = Signal24[i] / (2 ^ 32);
			Signal24[i] = Signal24[i] * (index);
		}
		signal24 << Signal24[i] << "		" << i << "\n";
	}*/

	//std::cout << x<<16 << '\n';
	/*float f = 1.0;
	unsigned char* ucp =
		reinterpret_cast<unsigned char*>(&f);

	for (int i = sizeof(float) - 1; i >= 0; --i)
	{
		printInBinary(ucp[i]);
	}*/



	/*std::cout << '\n';
	float f1 = 1.0;
	int* ucp1 =
		reinterpret_cast<int*>(&f1);

	for (int i = sizeof(float) - 1; i >= 0; --i)
	{
		mas[i] = printInBinary(ucp1[i]);
		//mas[i] = ucp1[i];
		std::cout << ucp1[i] << " ";
	}
	std::cout << '\n';
	for (int i = 16; i >= 0; i--)
	{
		//std::cout << mas[i] << " ";
	}
	std::cout << '\n';

	// inspect memory from c[0] to c[sizeof f - 1]

	/*for (int i = sizeof(x) - 1; i >= 0; i--)
	{
		std::cout << std::bitset<8>(reinterpret_cast<char*>(&x)[i]);
		mas[i] = (reinterpret_cast<int*>(&x)[i]);
	}
	std::cout << '\n';
	for (int i = 31; i >= 0; i--)
	{
		std::cout << mas[i] << " ";
	}
	std::cout << '\n';*/







	///////////////////////////////    DFT Signal 8      \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

	std::ofstream Signal8_d("DFT_Signal8.txt", std::ios_base::out);
	std::cout << std::endl << " DFT_Signal8 " << std::endl;

	std::ofstream FCHt8("FCH8.txt", std::ios_base::out);
	std::ofstream ACHt8("ACH8.txt", std::ios_base::out);
	std::ofstream FCHtlg8("FCHLG8.txt", std::ios_base::out);
	std::ofstream ACHtlg8("ACHLG8.txt", std::ios_base::out);


	N_SN = CountTest;
	real_SN = 0.0;
	mnim_SN = 0.0;
	Sum1_SN = 0.0;
	Sum2_SN = 0.0;
	param_SN = 0.0;
	st = 0;
	sr = 0;
	power = 0.0;
	for (int j = 0; j < CountTest; j++) // Прием данных в буфер
	{
		mas_int[j] = Signal8[j];
		mas_work[j] = 0;
	}
		for (int i = 0; i < N_SN; i++) // Обработка данных из буфера
		{
			Sum1_SN = 0;
			Sum2_SN = 0;
			for (int j = 0; j < N_SN; j++)
			{
				param_SN = -2.0 * PI * j * i / N_SN;
				real_SN = mas_int[j] * cos(param_SN);
				mnim_SN = mas_int[j] * sin(param_SN);
				Sum1_SN = Sum1_SN + real_SN;
				Sum2_SN = Sum2_SN - mnim_SN;
			}

			/*//mas_work[i] = ((Sum1_SN * Sum1_SN) / (CountTest / 2)); // Вроде работает
			mas_work[i] = ((Sum1_SN * Sum1_SN + Sum2_SN * Sum2_SN) / (CountTest)); // Работает для синуосв
			//mas_work[i] = ((Sum1_SN * Sum1_SN + Sum2_SN * Sum2_SN) / (CountTest / 2)); // Работает для синуосв

			DFT_SN[sr] = mas_work[i];
			Signal8_d << DFT_SN[sr] << "		" << sr << "\n";
			sr++;*/

			param_SN = 2.0 * PI * i / N_SN;
			real_SN = Signal8[i] * cos(param_SN);
			mnim_SN = Signal8[i] * sin(param_SN);

			FCH[i] = atan(mnim_SN / real_SN); // ФЧХ
			FCHt8 << FCH[i] << "		" << i << "\n";

			ACH[i] = (sqrt(real_SN * real_SN + mnim_SN * mnim_SN)); // АЧХ
			ACHt8 << ACH[i] << "		" << i << "\n";

			mas_work[i] = ((Sum1_SN * Sum1_SN + Sum2_SN * Sum2_SN)); // Спектр
			Signal8_d << mas_work[i] << "		" << i << "\n";

			ACHlg[i] = 20 * log10(mas_int[i]);
			ACHtlg8 << ACHlg[i] << "		" << i << "\n";

			power = power + ((Sum1_SN * Sum1_SN + Sum2_SN * Sum2_SN) / (CountTest)); // мощность
		}
		std::cout << std::endl << " Power signal_8 = " << power << std::endl;



	////////////////////////////// ФИЛЬТР ////////////////////
	z1 = 0, z2 = 0, z3 = 0, z4 = 0;
	sum1 = 0, sum2 = 0;
	flag1 = 0, flag2 = 0, flag3 = 0, flag4 = 0;

	std::cout << std::endl;
	for (int i = 0; i < Count + 2; i++)
	{
		filtr[i] = 0.99 * Signal8[i];
		//z4 = z3;
		//z3 = sum2;
		//sum2 = sum1 + z3 * 1.883 + z4 * (-0.98);
		sum1 = filtr[i] + z1 * (-1.902) + z2 * 1;
		z2 = z1;
		z1 = filtr[i];

		//exit[i] = sum2;
		exit[i] = sum1;
		flag1 = 1;
		//std::cout << "Signal[i]" << "\tZ1" << "\tZ2" << "\tSum1" << "\tSum2" << "\tZ3" << "\tZ4" << endl << Signal[i] << "\t" << z1 << "\t" << z2 << "\t" << sum1 << "\t" << z3 + sum1 + sum2 << "\t" << z3 + sum1 + sum2 << "\t" << sum2 << endl << endl;
	}
	std::cout << std::endl;



	//////////////////////////////////  DFT _FILTR   /////////////////////////////////////////

	std::ofstream Signal8_dF("Filtr_SIgnal8.txt", std::ios_base::out);
	std::ofstream FCHtFiltr8("FCHFiltr8.txt", std::ios_base::out);
	std::ofstream ACHtFiltr8("ACHFiltr8.txt", std::ios_base::out);
	std::ofstream FCHtlgFiltr8("FCHLGFiltr8.txt", std::ios_base::out);
	std::ofstream ACHtlgFiltr8("ACHLGFiltr8.txt", std::ios_base::out);
	float* FCHFiltr8 = new float[Count];
	float* ACHFiltr8 = new float[Count];
	float* FCHlgFiltr8 = new float[Count];
	float* ACHlgFiltr8 = new float[Count];

	N_SN_F = CountTest;
	real_SN_F = 0.0;
	mnim_SN_F = 0.0;
	Sum1_SN_F = 0.0;
	Sum2_SN_F = 0.0;
	param_SN_F = 0.0;
	st_F = 0;
	sr_F = 0;
	power = 0.0;
	for (int i = 0; i < N_SN_F; i++) // Обработка данных из буфера
	{
		Sum1_SN_F = 0;
		Sum2_SN_F = 0;
		for (int j = 0; j < N_SN_F; j++)
		{

			param_SN_F = -2.0 * PI * j * i / N_SN_F;
			real_SN_F = exit[j] * cos(param_SN_F);
			mnim_SN_F = exit[j] * sin(param_SN_F);
			Sum1_SN_F = Sum1_SN_F + real_SN_F;
			Sum2_SN_F = Sum2_SN_F - mnim_SN_F;
		}

		/*//mas_work_F[i] = ((Sum1_SN_F * Sum1_SN_F) / (CountTest / 2)); // Вроде работает для синусов
		mas_work_F[i] = ((Sum1_SN_F * Sum1_SN_F + Sum2_SN_F * Sum2_SN_F) / (CountTest)); // Вроде работает для синусов
																	 
																	 //mas_work_F[i] = 20*log10(sqrt(Sum1_SN_F * Sum1_SN_F + Sum2_SN_F * Sum2_SN_F)) / (CountTest / 2); // Это хорошо для шума
		//mas_work_F[i] = 20 * log10(sqrt(Sum1_SN_F * Sum1_SN_F + Sum2_SN_F * Sum2_SN_F)) / (CountTest / 2); //
		//mas_work_F[i] = (sqrt(Sum1_SN_F * Sum1_SN_F + Sum2_SN_F * Sum2_SN_F) / (CountTest / 2)); // Старое

		DFT_SN_F[sr_F] = mas_work_F[i];
		Signal8_dF << DFT_SN_F[sr_F] << "		" << sr_F << "\n";
		sr_F++;*/

		FCHFiltr8[i] = atan(Sum2_SN_F / Sum1_SN_F); // ФЧХ
		FCHtFiltr8 << FCHFiltr8[i] << "		" << i << "\n";

		ACHFiltr8[i] = (sqrt(Sum1_SN_F * Sum1_SN_F + Sum2_SN_F * Sum2_SN_F)); // АЧХ
		ACHtFiltr8 << ACHFiltr8[i] << "		" << i << "\n";

		mas_work_F[i] = ((Sum1_SN_F * Sum1_SN_F + Sum2_SN_F * Sum2_SN_F)); // Спектр
		Signal8_dF << mas_work_F[i] << "		" << i << "\n";

		ACHlgFiltr8[i] = 20 * log10(exit[i]);
		ACHtlgFiltr8 << ACHlgFiltr8[i] << "		" << i << "\n";

		power = power + ((Sum1_SN_F * Sum1_SN_F + Sum2_SN_F * Sum2_SN_F) / (CountTest)); // мощность
	}
	std::cout << std::endl << " Power signal_8 after filtr = " << power << std::endl;




	float summ1 = 0, summ2 = 0;
	for (int i = 0; i < Count; i++)
	{
		summ1 = summ1 + abs(Signal[i]);
		summ2 = summ2 + abs(Signal16[i]);
	}
	summ1 = summ1 / Count;
	summ2 = summ2 / Count;
	std::cout << summ1 << " = ";
	std::cout << summ2 << " " << endl;
	summ2 = summ2 - summ1;
	float SNR = 0;
	SNR = 20 * log10(summ1/ summ2);
	std::cout << "SNR16 = " << SNR << " " << endl;

	summ2 = 0;
	for (int i = 0; i < Count; i++)
		summ2 = summ2 + abs(Signal8[i]);
	summ2 = summ2 / Count;
	std::cout << summ1 << " = ";
	std::cout << summ2 << " " << endl;
	summ2 = summ2 - summ1;
	SNR = 0;
	SNR = 20 * log10(summ1 / summ2);
	std::cout << "SNR8 = " << SNR << " " << endl;


	delete[] Signal;
	delete[] Noise;
	delete[] Signalsin1; 
	delete[] Signalsin2;
	delete[] Signalsin3;
	sinus1.close();
	sinus2.close();
	sinus3.close();
	sinus4.close();
	noise.close();
	signal.close();
	std::cout << std::endl << " Finish " << std::endl;
	
	return 0;
}