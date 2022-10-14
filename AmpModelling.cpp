#include "AmpModelling.h"


//implementation of classes
	ampStage::ampStage(float gain) {
		this->a = gain;
	}

	ampStage::ampStage() {
		this->a = 1;
	}
	
	float ampStage::run(float input) {
		return a * input;
	}
	


	ampStageL::ampStageL() {
		this->gain = 1;
		this->limitDown = -1;
		this->limitUp = 1;
	}

	ampStageL::ampStageL(float gain, float limit) {
		this->gain = gain;
		this->limitUp = limit;
		this->limitDown = -limit;
	}

	ampStageL::ampStageL(float gain, float limitUp, float limitDown) {
		this->gain = gain;
		this->limitUp = limitUp;
		this->limitDown = limitDown;
	}


	float ampStageL::run(float input) {
		float temp = input * this->gain;
		if (temp > this->limitUp)return limitUp;
		if (temp < this->limitDown)return limitDown;
		return temp;
	}



	BQFilter::BQFilter(float a0, float a1, float a2, float b1, float b2, float overallgain, FilterType type) {
		this->a0 = a0;
		this->a1 = a1;
		this->a2 = a2;
		this->b1 = b1;
		this->b2 = b2;
		this->z1 = 0;
		this->z2 = 0;
		this->gain = overallgain;
	}

	BQFilter::BQFilter(FilterType type, float Fs, float Fc, float Q, float gain, float overallGain) {
		float norm;
		float K = std::tan(PI * Fc / Fs);
		this->gain = overallGain;
		
		switch (type)
		{
		case LowPass1P:
			this->b1 = std::exp(-2.0f * PI * (Fc / Fs));
			this->a0 = 1.0f - this->b1;
			this->b1 = -this->b1;
			this->a1 = this->a2 = this->b2 = 0;
			break;
		case LowPass:
			norm = 1 / (1 + K / Q + K * K);
			this->a0 = K * K * norm;
			this->a1 = 2 * this->a0;
			this->a2 = this->a0;
			this->b1 = 2 * (K * K - 1) * norm;
			this->b2 = (1 - K / Q + K * K) * norm;
			break;
		case HighPass:
			norm = 1 / (1 + K / Q + K * K);
			this->a0 = 1 * norm;
			this->a1 = -2 * this->a0;
			this->a2 = this->a0;
			this->b1 = 2 * (K * K - 1) * norm;
			this->b2 = (1 - K / Q + K * K) * norm;
			break;
		case HighPass1P:
			this->b1 = -std::exp(-2 * PI * (0.5f - Fc / Fs));
			this->a0 = 1 + this->b1;
			this->b1 = -this->b1;
			this->a1 = this->a2 = this->b2 = 0;
			break;
		case BandPass:
			norm = 1 / (1 + K / Q + K * K);
			this->a0 = K / Q * norm;
			this->a1 = 0;
			this->a2 = -this->a0;
			this->b1 = 2 * (K * K - 1) * norm;
			this->b2 = (1 - K / Q + K * K) * norm;
			break;
		case BandStop:
			norm = 1 / (1 + K / Q + K * K);
			this->a0 = (1 + K * K) * norm;
			this->a1 = 2 * (K * K - 1) * norm;
			this->a2 = this->a0;
			this->b1 = this->a1;
			this->b2 = (1 - K / Q + K * K) * norm;
			break;
		case LowShelf:
			if (gain >= 0) {
				norm = 1 / (1 + SQRT2 * K + K * K);
				this->a0 = (1 + std::sqrt(2 * gain) * K + gain * K * K) * norm;
				this->a1 = 2 * (gain * K * K - 1) * norm;
				this->a2 = (1 - std::sqrt(2 * gain) * K + gain * K * K) * norm;
				this->b1 = 2 * (K * K - 1) * norm;
				this->b2 = (1 - SQRT2 * K + K * K) * norm;
			}
			else {
				norm = 1 / (1 + std::sqrt(2 * gain) * K + gain * K * K);
				this->a0 = (1 + SQRT2 * K + K * K) * norm;
				this->a1 = 2 * (K * K - 1) * norm;
				this->a2 = (1 - SQRT2 * K + K * K) * norm;
				this->b1 = 2 * (gain * K * K - 1) * norm;
				this->b2 = (1 - std::sqrt(2 * gain) * K + gain * K * K) * norm;
			}
			break;
		case HighShelf:
			if (gain >= 0) {
				norm = 1 / (1 + SQRT2 * K + K * K);
				this->a0 = (gain + std::sqrt(2 * gain) * K + K * K) * norm;
				this->a1 = 2 * (K * K - gain) * norm;
				this->a2 = (gain - std::sqrt(2 * gain) * K + K * K) * norm;
				this->b1 = 2 * (K * K - 1) * norm;
				this->b2 = (1 - SQRT2 * K + K * K) * norm;
			}
			else {
				norm = 1 / (gain + std::sqrt(2 * gain) * K + K * K);
				this->a0 = (1 + SQRT2 * K + K * K) * norm;
				this->a1 = 2 * (K * K - 1) * norm;
				this->a2 = (1 - SQRT2 * K + K * K) * norm;
				this->b1 = 2 * (K * K - gain) * norm;
				this->b2 = (gain - std::sqrt(2 * gain) * K + K * K) * norm;
			}
			break;
		default:
			throw std::invalid_argument("filter type error");
			break;
		}
		this->z1 = 0;
		this->z2 = 0;
	}

	float BQFilter::run(float input) {
		float y = input * this->a0 + this->z1;
		this->z1 = input * this->a1 + this->z2 - this->b1 * y;
		this->z2 = input * this->a2 - this->b2 * y;
		return y;
	}


	//implementation of misc functions