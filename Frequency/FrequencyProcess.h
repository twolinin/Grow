///-----------------------------------------------
// 
// Written by Jyun-Hong Lin
// 
//-----------------------------------------------
//
// Frequency
//

#ifndef FrequencyProcess_H
#define FrequencyProcess_H

#include "SequenceWorkItem.h"
#include "BWTAlgorithms.h"
#include "BWTIndexSet.h"
#include "SGVisitors.h"

// Parameter object for the reassembler
struct FrequencyParameters
{    
    BWTIndexSet indices;
    
    int kmerLength;
	int numThreads;

};

class FrequencyResult
{
    public:
        FrequencyResult(){}
        //SeedByReadHashMap seedHash;
};

class ReadFrequencyBasicElements
{
	public:

		ReadFrequencyBasicElements(const FrequencyParameters params);
		~ReadFrequencyBasicElements();
		
		std::string backtrackRead(const BWTIndexSet indices, int64_t inputIdx);
		// return string frequency
		size_t getFrequency(const BWTIndexSet indices, std::string query);
		
		void showReadFrequency(size_t queryIndex, int kmerSize);
		void showReadFrequency(std::string sequence, int kmerSize);
		
	protected:
		
		FrequencyParameters m_params;
};

class FrequencyProcess : public ReadFrequencyBasicElements
{
	public:
		
		FrequencyProcess(const FrequencyParameters params);
		~FrequencyProcess();

		FrequencyResult PBFrequency(const SequenceWorkItem& workItem);
		
		FrequencyResult process(const SequenceWorkItem& workItem)
		{
			return PBFrequency(workItem);
		}
		
	protected:

};

class FrequencyPostProcess : public ReadFrequencyBasicElements
{
    public:
        
		FrequencyPostProcess(const FrequencyParameters params);
        ~FrequencyPostProcess();

        void process(const SequenceWorkItem& item, FrequencyResult& result);
    
    protected:
		
};



#endif
