#pragma once

#include <cmath>
#include <stdexcept>

namespace AmpModelling {

	constexpr float PI = 3.141592653589793f;
	constexpr float SQRT2 = 1.414213562373095f;
	constexpr float DELTA = 1e-6f;


	enum FilterType {
		LowPass1P, LowPass, HighPass, HighPass1P, BandPass, BandStop, LowShelf, HighShelf
	};

	enum TriodeType {
		ECC83_EH, ECC83_JJ, ECC83_NK
	};

	//declaration of classes


	class ampStage {
	private:
		float a;
	public:
		ampStage(float gain);
		ampStage();

		float run(float input);

	};

	class ampStageL {
	private:
		float gain;
		float limitUp;
		float limitDown;
	public:
		ampStageL();
		ampStageL(float gain, float limit);
		ampStageL(float gain, float limitUp, float limitDown);
		float run(float input);
	};


	class BQFilter {
	public:
		float a0;
		float a1;
		float a2;
		float b1;
		float b2;

		float z1;
		float z2;

		float gain;

		BQFilter(float a0, float a1, float a2, float b1, float b2, float overallgain, FilterType type);
		BQFilter(FilterType type, float Fs, float Fc, float Q, float gain, float overallGain);
		float run(float input);
	};

	class triode {
	private:
		float diff(float Ia);
	public:
		//static parameters of the triode
		float KP;
		float mu;
		float KVB;
		float X;
		float KG1;
		//TODO add dynamic parameters
		//component values around the triode
		float Ra;
		float Rk;
		float Vb;

		float Ia;

		triode(TriodeType type, float Ra, float Rk, float Vb);

		float E1(float Ia);

		float IaApprox(float Ia);

	};

	//declaration of misc functions

	float diff(float (*f)(float), float x);

}
