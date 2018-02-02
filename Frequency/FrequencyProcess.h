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

//typedef std::tr1::hash<int> Int64Hasher;
typedef std::tr1::hash<size_t> SizetHasher;

typedef SparseHashMap<size_t, size_t, SizetHasher> KmerFreqHashMap;
typedef std::pair<size_t, size_t> sizetPair;
typedef std::vector<sizetPair> sizetVec;


// Parameter object for the reassembler
struct FrequencyParameters
{    
    BWTIndexSet indices;
    
    int kmerLength;
	int numThreads;
};

struct KmerInfo
{
	std::string kmerStr;
	
	bool kmerStrand;
	bool isRepeat;
	
	std::string readConnectContigID;
	int readConnectContigLength;
	bool readConnectContigSide;
	
	int readIndex;
	int position;
	int secondPosition;
	
	int queryPos;
	int sbjctPos;
	int pairNumber;
};

struct ReadScanAndOverlapPosition
{
	int scanQueryStart;
	int scanQueryEnd;
	int scanSbjctStart;
	int scanSbjctEnd;
	int cutContigToEndPointLength;
	int kmerCount;
	int sbjctLength;
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
		std::pair<size_t, size_t> BacktrackReadIdx(const BWTIndexSet indices, size_t idx, bool activeHash);
		
		size_t getSeqIdx(const BWTIndexSet indices, std::string query);
		
		// return string frequency
		size_t getFrequency(const BWTIndexSet indices, std::string query);
		//void findRepeatRegion(std::string sequence, int kmerSize, int repeatThreshold);
		
		void showReadFrequency(size_t queryIndex, int kmerSize);
		void showReadFrequency(std::string sequence, int kmerSize);
		
	protected:
		
		SparseHashMap<size_t, sizetPair, SizetHasher> recordBacktrackIndex;
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
		
		int repeatThreshold(int kmerSize);
		void showFrequency(std::string originSequence);
		bool checkChimera(size_t readIndex);
		bool checkChimera(std::string querySeq);
		
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
