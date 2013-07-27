#ifndef __UTILS_H__
#define __UTILS_H__

#include <cmath>
#include <windows.h>

class Timer
{
public:
	Timer() {
		this->init();
	}

	~Timer(){};

	bool start() {

		bool bSuccess = false;

		if(!m_bTimerRunning && m_bTimerSupported)
		{
			m_startCount = 0;
			m_stopCount = 0;
			m_interval = 0;

			if(QueryPerformanceCounter((LARGE_INTEGER*)&m_startCount))
			{
				m_bTimerRunning = true;
				bSuccess = true;
			}
		}

		return bSuccess;
	}

	
	bool stop() {

		bool bSuccess = false;

		if(m_bTimerRunning && m_bTimerSupported)
		{
			if(QueryPerformanceCounter((LARGE_INTEGER*)&m_stopCount))
			{
				m_bTimerRunning = false;
				bSuccess = true;
			}
		}

		return bSuccess;
	}

	void reset() {
		this->init();
	}

	double interval_S() {
		return ((double)(m_stopCount - m_startCount) - m_adjustCount) / (double)m_frequency;
	}

	double interval_mS() {
		return (((m_stopCount - m_startCount) - m_adjustCount) * 1000.0) / (double)m_frequency;
	}

	double interval_uS() {
		return (((m_stopCount - m_startCount) - m_adjustCount) * 1000000.0) / (double)m_frequency;
	}

	double resolution_S() {
		return 1.0 / (double)m_frequency;
	}

	double resolution_mS() {
		return 1000.0 / (double)m_frequency;
	}

	double resolution_uS() {
		return 1000000.0 / (double)m_frequency;
	}

	double correction_uS() {
		return (m_adjustCount * 1000000.0) / (double)m_frequency;
	}

	void init() {

		m_frequency = 0;
		m_adjustCount = 0;
		m_bTimerSupported = false;
		m_bTimerRunning = false;

		if(QueryPerformanceFrequency((LARGE_INTEGER*)&m_frequency))
		{
			m_bTimerSupported = true;

			// Measure the 'Stop' function call overhead
			const int iNumSamples = 10;
			__int64 samples[iNumSamples];
			__int64 countTot = 0;
			double dAvCount = 0.0;
			double dAvDeviance = 0.0;

			for(int i = 0; i < iNumSamples; i++)
			{
				this->start();
				this->stop();

				samples[i] = m_stopCount - m_startCount;
				countTot += samples[i];
			}

			dAvCount = (double)countTot / (double)iNumSamples;

			// Get the average deviance
			for(int i = 0; i < iNumSamples; i++)
			{
				dAvDeviance += fabs(((double)samples[i]) - dAvCount);
			}

			// Average deviance only required for debug
			dAvDeviance /= iNumSamples;
			m_adjustCount = (__int64)dAvCount;
		}
	}

	void set_start_time(__int64 & c) {
		c = 0;
		QueryPerformanceCounter((LARGE_INTEGER*) &c);
	}

	double elapsed_time(__int64 & c) {
		__int64 end = 0;
		QueryPerformanceCounter((LARGE_INTEGER*)&end);
		return (double)((end - c) - m_adjustCount) / (double)m_frequency;
	}

private:
	__int64 m_interval;
	__int64 m_frequency;
	__int64 m_startCount;
	__int64 m_stopCount;
	__int64 m_adjustCount;
	bool m_bTimerSupported;
	bool m_bTimerRunning;
};

#endif