///-----------------------------------------------
// 
// Written by Jyun-Hong Lin
// 
//-----------------------------------------------
//
// ReadLongReadReassemble
//

#ifndef ReassemblerProcess_H
#define ReassemblerProcess_H

#include "SequenceWorkItem.h"
#include "BWTAlgorithms.h"
#include "BWTIndexSet.h"
#include "SGVisitors.h"
#include "ContigGraph.h"

typedef std::tr1::hash<int> Int64Hasher;
typedef std::tr1::hash<size_t> SizetHasher;

// collect all first step seeds. store seeds, read and contig information 
struct SeedSequenceInfo;
typedef std::vector<SeedSequenceInfo> SeedSequenceInfoVec ;
// collect all thread result. used in filter error reads
struct OverlapSeedVecInfo;
typedef SparseHashMap<int, std::vector<OverlapSeedVecInfo>, Int64Hasher> SeedByReadHashMap;
// build relation between contigs after filter error read
struct OverlapRelation;
typedef SparseHashMap<std::string, std::vector<OverlapRelation>, StringHasher> ContigRelationHashMap;
// collect all possible read
typedef SparseHashMap<int, std::vector<OverlapRelation>, Int64Hasher> ReadConnectContigHashMap;
// read index and it's seeds(kmers)
struct KmerInfo;
typedef std::vector<KmerInfo> KmerInfoVec;
typedef std::pair<int, KmerInfoVec> ReadKmerPair;
typedef std::vector<ReadKmerPair> ReadKmerVec;
// all kmers infomation
typedef SparseHashMap<std::string, std::vector<KmerInfo>, StringHasher> KmerHashMap;
// speed up backtrack
typedef std::pair<size_t, size_t> sizetPair;

struct ReassemblerSeedFeature
{
    public:
        ReassemblerSeedFeature(size_t startPos, std::string str, bool repeat, size_t kmerSize, size_t repeatCutoff);
         
        size_t seedStartPos;
        size_t seedEndPos;
        size_t seedLength;
        std::string seedStr;
        bool isRepeat;
        
        // estimated by calling estimateBestKmerSize
        size_t startBestKmerSize;
        size_t endBestKmerSize;
        size_t startKmerFreq;
        size_t endKmerFreq;
        
    private:
        size_t freqUpperBound;
        size_t freqLowerBound;

};

struct SeedSequenceInfo
{
        void setContigSeq(std::string contigID, std::string contigPartialSeq, int contigLength, int contigStartPos, int contigEndPos, bool contigSide);
        void setReadSeq(std::string seedString, int frequency, int readIndex, int readLocation, bool strand);
        void printSeedInfo();

		std::string contigPartialSeq;
        std::string contigID;
        int contigLength;
        int contigStartPos;
        int contigEndPos;
        int seedLength;
        bool contigSide; 
        
        std::string seedString;
        int frequency;
        int readIndex;
        int readLocation;
        bool strand;

};

struct ResultCount
{
	ResultCount():totalCount(0),partialCount(0),passCount(0),zeroSeed(0),seedsThreshold(0),seedDifferentStrand(0),seedOrderDiscordant(0),abnormalDistance(0),repeatSeeds(0),repeatStatus(0),
	abnormalOverlapLength(0),insufficientLength(0),palindrome(0),chimera(0),sucsess(0){}
	
	int totalCount;
	int partialCount;
	int passCount;
	int zeroSeed;
    int seedsThreshold;
    int seedDifferentStrand;
    int seedOrderDiscordant;
    int abnormalDistance;
    int repeatSeeds;
	int repeatStatus;
    int abnormalOverlapLength;
    int insufficientLength;
	int palindrome;
	int chimera;
	int sucsess;

	int readSeedsThreshold;
	int readRepeatSeeds;
	int readSeedDifferentStrand;
	int readSeedOrderDiscordant;
	int readAbnormalDistance;
	int readAbnormalOverlapLength;
};

struct OverlapSeedVecInfo
{
    std::string contig;
	SeedSequenceInfoVec SeedInfoVec;
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

struct contigReadPosition
{
    int64_t firstContigStartPos;
    int64_t lastContigStartPos;
    int64_t firstReadStartPos;
    int64_t lastReadStartPos;
    int64_t minReadPos;
    int64_t maxReadPos;
    int64_t contigEndPointOnReadPos;
    int64_t offsetOverlapLength;
};

struct OverlapRelation
{
    std::string firstContig;
    std::string secondContig;
    bool firstSide;
    bool secondSide;
    bool firstStrand;
    bool secondStrand;
    // overlap parameter
    int firContigStart;
    int firContigEnd;
    int secContigStart;
    int secContigEnd;
        
    std::string readFragment;
    int readLength;
    int readIndex;
    
    int overlapLength;
    int firstContigLength;
    int secondContigLength;
    
    bool connect;// false : overlap, true : nonOverlap
    
	bool readDirection;
    bool isSameStrand;
	
	// two read step
	int firIdx;
	int secIdx;
	
	ReadScanAndOverlapPosition singleOverlap;
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

// Parameter object for the reassembler
struct ReassemblerParameters
{    
    BWTIndexSet indices;
    BWTIndexSet contigIndices;
    
    int numKmerRounds;
    int kmerLength;
	int numThreads;
	
    // tree search parameters
    int maxLeaves;
    int minOverlap;
    int maxOverlap;

    // Read
    int minKmerLength;
    int FMWKmerThreshold;
    int seedKmerThreshold;
    int collectedSeeds;

    size_t maxSeedInterval;
    int readAvgLen;
	int searchRange;

    // first, reassemble overlap
    // second, reassemble non overlap by one read
	// third, reassemble non overlap by two read
    bool isFirst;
    bool isSecond;
	bool isThird;
};

class ReassemblerResult
{
    public:
        ReassemblerResult(){}
        SeedByReadHashMap seedHash;
};

class ReadReassemblerBasicElements
{
	public:

		ReadReassemblerBasicElements(const ReassemblerParameters params);
		~ReadReassemblerBasicElements();
		
		std::vector<ReassemblerSeedFeature> seedingByDynamicKmer(const std::string readSeq);
		std::vector<ReassemblerSeedFeature> seedingByDynamicKmer(const std::string readSeq, size_t scanStart, size_t scanEnd);
		// third step reassembler, take kmer as seed
		std::vector<ReassemblerSeedFeature> collectKmer(const std::string readSeq, int scanStart, int scanEnd);
		// return kmer freq of beginning and ending kmers
		std::pair<size_t, size_t> refineRepeatSeed(const std::string readSeq, size_t& seedStartPos, size_t& seedEndPos);
		// ratio of A,T,C,G
		bool isLowComplexity (std::string seq, float threshold);
		// backtrack to '$' and return index
		std::pair<size_t, size_t> BacktrackReadIdx(const BWTIndexSet indices, size_t idx, bool activeHash);
		// backtrack origin read
		std::string backtrackRead(const BWTIndexSet indices, int64_t inputIdx);
		// return string frequency
		size_t getFrequency(const BWTIndexSet indices, std::string query);
		// if read is origin read, this function will return index
		size_t getSeqIdx(const BWTIndexSet indices, std::string query);
		
		void showReadFrequency(size_t queryIndex, int kmerSize);
		
	protected:
		
		SparseHashMap<size_t, sizetPair, SizetHasher> recordBacktrackIndex;
		ReassemblerParameters m_params;
};

class ReassemblerProcess : public ReadReassemblerBasicElements
{
	public:
		
		ReassemblerProcess(const ReassemblerParameters params);
		~ReassemblerProcess();

		ReassemblerResult PBReassembler(const SequenceWorkItem& workItem);
		
		ReassemblerResult process(const SequenceWorkItem& workItem)
		{
			return PBReassembler(workItem);
		}
		
	protected:

		void findSeedOnLongRead(const BWTIndexSet indices, std::string contigID, std::vector<ReassemblerSeedFeature> inputSeedVec, SeedSequenceInfoVec &result, int contigLength, std::string contigPartialSeq, bool contigSide);
};

class ReassemblerPostProcess : public ReadReassemblerBasicElements
{
    public:
        
		ReassemblerPostProcess(const ReassemblerParameters params, ContigGraph* pGraph);
        ~ReassemblerPostProcess();

        void process(const SequenceWorkItem& item, ReassemblerResult& result);
        void filterErrorRead();
        void buildGraphByOneRead();
        void buildGraphByTwoRead();
        
    protected:
		
		void collectSeeds();
		void filterErrorConnect(ReadKmerVec &inputReadVec, ContigRelationHashMap &inputContigPair, bool isKmer);
		void layoutGraph(ContigRelationHashMap &inputContigPair);
		
		///************************************************///
        ///******** filter error read function **********///
        ///************************************************///
		
		//
		void combineEveryCheckFunction(OverlapSeedVecInfo firVec, OverlapSeedVecInfo secVec, std::string firContig, std::string secContig, int readIndex, int section);
        // check two contigs side and strand status
        bool checkContigRelation(bool firstSide, bool firstStrand, bool secondSide, bool secondStrand);
        // check each seed status and length
        bool checkContigSeed(OverlapSeedVecInfo input, int lengthThreshold);
        // check order by position, increase or decrease
        bool checkSeedsOrientations(OverlapSeedVecInfo input);
        // the distance between two contig and read
        bool checkSeedsDistance(OverlapSeedVecInfo input);
        // using read index check repeat
        bool checkRepeatByRead(OverlapSeedVecInfo input, bool show);
        // using contig index check repeat
        bool checkRepeatByContig(OverlapSeedVecInfo input);
		//
		bool checkPalindrome(std::string readSeq);
		//
		bool checkChimera(size_t readIndex);
		bool checkChimera(std::string querySeq);
        // the distance between first and last seed
        bool checkOverlapLength(OverlapSeedVecInfo input);
		// using in two read connect partial. record non overlap position in pair structure
		bool checkNonOverlapPartialLength(OverlapSeedVecInfo input, ReadScanAndOverlapPosition &position);
		
		///************************************************///
        ///**************** tool function *****************///
        ///************************************************///
		
		
		void repeatThresholdCaculate();
		// show kmer distribution, each filter count
		void showDetailInformation();
		// first and second step reassembler, store two contigs connect relation
        void insertContigRelationHash(std::string firstContig, bool firstSide, bool firstStrand, std::string secondContig, bool secondSide, bool secondStrand, int firstContigLength, int secondContigLength, std::vector<OverlapSeedVecInfo> contigInfoVec);
		//
		void insertContigRelationHash(ReadKmerPair readVecIter, OverlapSeedVecInfo resultSeedVec, ContigRelationHashMap &inputContigPair);
		// third step reassembler, store a read overlap on contig information 
        void insertSingleRelationHash(OverlapSeedVecInfo input, ReadScanAndOverlapPosition scanRange);
        // found first and last seed position on contig, so and read 
        contigReadPosition OnReadRange(SeedSequenceInfoVec inputSeedInfoVec, int64_t readLength);
		//
		ReadScanAndOverlapPosition reScanContigOnReadRange(std::string partialContig, std::string originRead);
        // return star and end position on contig
        std::pair<int,int> overlapPosition(bool side, int contigLength, int overlapLength);
        // count all seeds second connect contig. A to B have most connect, B should have most connect to A.
        std::string mostConnectContig(std::vector<OverlapRelation> inputRelationVec);
		// return most seeds connect read 
		int mostConnectRead(OverlapSeedVecInfo input, bool show);
		// return most seeds strand 
		bool mostSeedStrand(OverlapSeedVecInfo input);
		// return most read connect contig 
		std::string mostReadConnectContig(OverlapSeedVecInfo input);
        // check two contigs are overlap or not
        int connectVote(size_t readIdx, std::vector<OverlapRelation> inputRelationVec);
		// check two contigs side are not connect, store this pair in hash if not connect
        bool checkConnect(std::string firstID, bool firstSide,std::string secondID, bool secondSide);
        //
		void filterErrorStrand(OverlapSeedVecInfo input, OverlapSeedVecInfo &result, bool mostStrand);
		// 
		void filterErrorSeeds(OverlapSeedVecInfo input, OverlapSeedVecInfo &result, int mostRead);
		// third step reassembler, construct kmer hash map and read seeds vector
		void buildKmerHashAndReadKmerVec( std::string readConnectContig, bool connectContigSide, int connectContigLength, int readIndex, std::vector<ReassemblerSeedFeature> inputSeedVec, KmerHashMap &allKmerHash, ReadKmerVec &readVec);
		//
		bool checkSameStrand(int startReadIndex, int targetReadIndex, SeedSequenceInfoVec brigdeSeedVec);
		// return overlap length
		int overlapLength(OverlapSeedVecInfo input);
		// second read strand will same of the first 
		
		
		std::string concordantNonOverlapFragment(OverlapSeedVecInfo input);
		
		///************************************************///
        ///******************* Variable *******************///
        ///************************************************///
		
        ContigGraph* m_pGraph;
        // collect full origin read connect 
        SeedByReadHashMap collectSeedHashMap;
        // use to know how many contig connect target contig, include side infomation
        ContigRelationHashMap connectContigHashMap;
        // use to check two read whether can connect
        ReadConnectContigHashMap singleContigHashMap;
        // check this contig pair connect 
        SparseHashMap<std::string,bool,StringHasher> connectStatus;
        
		///*** kmer ***///
		
		// third step reassemble, all kmer appear on the read information
		KmerHashMap collectAllKmerHashMap;
		// third step reassemble, read index and it's seeds vector	
		ReadKmerVec AllReadKmerVec;
		// use to know how many contig connect target contig, include side infomation
        ContigRelationHashMap contigPairByKmer;
		
		///*** gkmer ***///
		
		// third step reassemble, all gkmer appear on the read information
        KmerHashMap collectAllGKmerHashMap;
		// third step reassemble, read index and it's seeds vector	
		ReadKmerVec AllReadGKmerVec;
		// use to know how many contig connect target contig, include side infomation
        ContigRelationHashMap contigPairByGKmer;
        
		///************************************************///
        ///********* construct kmer distribution **********///
        ///************************************************/// 
 
		int seedCount_19;
        int seedCount_15;
		int seedCount_11;
		int kmerSizeFreqByRead_19[100000];
        int kmerSizeFreqByRead_15[100000];
		int kmerSizeFreqByRead_11[100000];
        void countTotalSeedFreq(std::string seed);
		
		int differenceCount_19;
        int differenceCount_15;
		int differenceCount_11;
		int differenceValue_19[100000];
        int differenceValue_15[100000];
		int differenceValue_11[100000];
		
		///************************************************///
        ///********* result statistics **********///
        ///************************************************/// 
		
		ResultCount onePB,twoPB;

};



#endif
