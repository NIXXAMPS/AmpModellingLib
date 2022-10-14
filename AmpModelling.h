#pragma once

#include <cmath>
#include <stdexcept>
constexpr auto PI = 3.141592653589793f;
constexpr auto SQRT2 = 1.414213562373095f;


	enum FilterType {
		LowPass1P, LowPass, HighPass, HighPass1P, BandPass, BandStop, LowShelf, HighShelf
	};

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
	private:
		float a0;
		float a1;
		float a2;
		float b1;
		float b2;

		float z1;
		float z2;
		
		float gain;
	public:
		BQFilter(float a0, float a1, float a2, float b1, float b2);
		BQFilter(float a0, float a1, float a2, float b1, float b2, float overallgain, FilterType type);
		BQFilter(FilterType type, float Fs, float Fc, float Q, float gain, float overallGain);
		float run(float input);
	};


