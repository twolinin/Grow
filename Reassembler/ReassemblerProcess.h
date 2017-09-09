///-----------------------------------------------
// 
// Written by Jyun-Hong Lin
// 
//-----------------------------------------------
//
// TGSLongReadReassemble
//

#ifndef ReassemblerProcess_H
#define ReassemblerProcess_H

#include "SequenceWorkItem.h"
#include "BWTAlgorithms.h"
#include "BWTIndexSet.h"
#include "SGVisitors.h"
#include "ContigGraph.h"

typedef std::tr1::hash<int> Int64Hasher;

// collect all first step seeds. store seeds, tgs and contig information 
struct SeedSequenceInfo;
typedef std::vector<SeedSequenceInfo> SeedSequenceInfoVec ;
// collect all thread result. used in filter error tgss
struct OverlapSeedVecInfo;
typedef SparseHashMap<int, std::vector<OverlapSeedVecInfo>, Int64Hasher> SeedByTGSHashMap;
// build relation between contigs after filter error tgs
struct OverlapRelation;
typedef SparseHashMap<std::string, std::vector<OverlapRelation>, StringHasher> ContigRelationHashMap;
// collect all possible tgs
typedef SparseHashMap<int, std::vector<OverlapRelation>, Int64Hasher> TGSConnectContigHashMap;
// tgs index and it's seeds(kmers)
struct KmerInfo;
typedef std::vector<KmerInfo> KmerInfoVec;
typedef std::pair<int, KmerInfoVec> TGSKmerPair;
typedef std::vector<TGSKmerPair> TGSKmerVec;
// all kmers infomation
typedef SparseHashMap<std::string, std::vector<KmerInfo>, StringHasher> KmerHashMap;

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
        void setContigSeq(std::string contigID, std::string contigPartialSeq, int contigLength, int contigStartPos, int contigEndPos, bool contigSide)
        {
            this->contigID         = contigID;
            this->contigPartialSeq = contigPartialSeq;
			this->contigLength     = contigLength;
            this->contigStartPos   = contigStartPos;
            this->contigEndPos     = contigEndPos;
            this->contigSide       = contigSide;
            seedLength = contigEndPos - contigStartPos + 1;
        }
        
        void setTGSSeq(std::string seedString, int frequency, int tgsIndex, int tgsLocation, bool strand)
        {
            this->seedString     = seedString;
            this->frequency      = frequency;
            this->tgsIndex    = tgsIndex;
            this->tgsLocation = tgsLocation;
            this->strand = strand;
        }    
        
        void printSeedInfo()
        {
            std::cout<< contigID       << "\t" 
                     << contigLength   << "\t"
                     << contigStartPos << "\t"
                     << contigEndPos   << "\t"
                     << seedLength     << "\t"
                     << (contigSide ? "Head" : "Tail") << "\t"
                     << tgsIndex    << "\t"
                     << tgsLocation << "\t"
                     << (strand ? "NonRC" : "RC") << "\t"
                     << frequency      << "\t"
                     //<< seedString     << "\t"
                     << "\n";
        }
        
		std::string contigPartialSeq;
        std::string contigID;
        int contigLength;
        int contigStartPos;
        int contigEndPos;
        int seedLength;
        bool contigSide; 
        
        std::string seedString;
        int frequency;
        int tgsIndex;
        int tgsLocation;
        bool strand;

};

struct ResultCount
{
	ResultCount():totalCount(0),partialCount(0),passCount(0),zeroSeed(0),seedsThreshold(0),seedDifferentStrand(0),seedOrderDiscordant(0),abnormalDistance(0),repeatSeeds(0),repeatStatus(0),
	abnormalOverlapLength(0),insufficientLength(0),palindrome(0){}
	
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
};

struct OverlapSeedVecInfo
{
    std::string contig;
	SeedSequenceInfoVec SeedInfoVec;
};

struct TGSScanAndOverlapPosition
{
	int scanTGSStart;
	int scanTGSEnd;
	int tgsFragmentStart;
	int tgsFragmentEnd;
	int cutContigToEndPointLength;
	int kmerCount;
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
        
    std::string tgsFragment;
    int tgsLength;
    int tgsIndex;
    
    int overlapLength;
    int firstContigLength;
    int secondContigLength;
    
    bool connect;// false : overlap, true : nonOverlap
    
	bool tgsDirection;
    bool isSameStrand;
	
	TGSScanAndOverlapPosition singleOverlap;
};

struct KmerInfo
{
	std::string kmerStr;
	
	bool kmerStrand;
	bool isRepeat;
	
	std::string tgsConnectContigID;
	int tgsConnectContigLength;
	bool tgsConnectContigSide;
	
	int tgsIndex;
	int position;
	int secondPosition;
};

struct contigTGSPosition
{
    int64_t firstContigStartPos;
    int64_t lastContigStartPos;
    int64_t firstTGSStartPos;
    int64_t lastTGSStartPos;
    int64_t minTGSPos;
    int64_t maxTGSPos;
    int64_t contigEndPointOnTGSPos;
    int64_t offsetOverlapLength;
};

// Parameter object for the reassembler
struct ReassemblerParameters
{    
    BWTIndexSet indices;
    BWTIndexSet contigIndices;
    
    int numKmerRounds;
    int kmerLength;

    // tree search parameters
    int maxLeaves;
    int minOverlap;
    int maxOverlap;

    // TGS
    int minKmerLength;
    int FMWKmerThreshold;
    int seedKmerThreshold;
    int collectedSeeds;

    size_t maxSeedInterval;
    int tgsAvgLen;
	int searchRange;

    // first, reassemble overlap
    // second, reassemble non overlap by one tgs
	// third, reassemble non overlap by two tgs
    bool isFirst;
    bool isSecond;
	bool isThird;
};

class ReassemblerResult
{
    public:
        ReassemblerResult(){}
        SeedByTGSHashMap seedHash;
};

class TGSReassemblerBasicElements
{
	public:

		TGSReassemblerBasicElements(const ReassemblerParameters params);
		~TGSReassemblerBasicElements();
		
		std::vector<ReassemblerSeedFeature> seedingByDynamicKmer(const std::string readSeq);
		// third step reassembler, take kmer as seed
		std::vector<ReassemblerSeedFeature> collectKmer(const std::string readSeq, int scanStart, int scanEnd);
		// return kmer freq of beginning and ending kmers
		std::pair<size_t, size_t> refineRepeatSeed(const std::string readSeq, size_t& seedStartPos, size_t& seedEndPos);
		// ratio of A,T,C,G
		bool isLowComplexity (std::string seq, float threshold);
		// backtrack to '$' and return index
		std::pair<size_t, size_t> BacktrackTGSIdx(const BWTIndexSet indices, size_t idx);
		// backtrack origin read
		std::string backtrackTGS(const BWTIndexSet indices, int64_t inputIdx);
		
	protected:
		
		ReassemblerParameters m_params;
};

class ReassemblerProcess : public TGSReassemblerBasicElements
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

class ReassemblerPostProcess : public TGSReassemblerBasicElements
{
    public:
        
		ReassemblerPostProcess(const ReassemblerParameters params, ContigGraph* pGraph);

        ~ReassemblerPostProcess();

        void process(const SequenceWorkItem& item, ReassemblerResult& result);
        
        void filterErrorTGS();
        
        void buildGraphByOneTGS();
        
        void buildGraphByTwoTGS();
        
    protected:

		///************************************************///
        ///******** filter error tgs function **********///
        ///************************************************///
		
        // check two contigs side and strand status
        bool checkContigRelation(bool firstSide, bool firstStrand, bool secondSide, bool secondStrand);
        // check each seed status and length
        bool checkContigSeed(OverlapSeedVecInfo input, int lengthThreshold);
        // check order by position, increase or decrease
        bool checkSeedsOrientations(OverlapSeedVecInfo input);
        // the distance between two contig and tgs
        bool checkSeedsDistance(OverlapSeedVecInfo input);
        // using tgs index check repeat
        bool checkRepeatByTGS(OverlapSeedVecInfo input);
        // using contig index check repeat
        bool checkRepeatByContig(OverlapSeedVecInfo input);
		//
		bool checkPalindrome(std::string readSeq);
        // the distance between first and last seed
        bool checkOverlapLength(OverlapSeedVecInfo input);
		// using in two tgs connect partial. record non overlap position in pair structure
		bool checkNonOverlapPartialLength(OverlapSeedVecInfo input, TGSScanAndOverlapPosition &position);
		
		///************************************************///
        ///**************** tool function *****************///
        ///************************************************///
		
		// first and second step reassembler, store two contigs connect relation
        void insertContigRelationHash(std::string firstContig, bool firstSide, bool firstStrand, std::string secondContig, bool secondSide, bool secondStrand, int firstContigLength, int secondContigLength, std::vector<OverlapSeedVecInfo> contigInfoVec);
		// third step reassembler, store a tgs overlap on contig information 
        void insertSingleRelationHash(OverlapSeedVecInfo input, TGSScanAndOverlapPosition scanRange);
        // found first and last seed position on contig, so and tgs 
        contigTGSPosition OnTGSRange(SeedSequenceInfoVec inputSeedInfoVec, int64_t tgsLength);
		//
		TGSScanAndOverlapPosition reScanContigOnTGSRange(std::string partialContig, std::string originTGS);
        // return star and end position on contig
        std::pair<int,int> overlapPosition(bool side, int contigLength, int overlapLength);
        // count all seeds second connect contig. A to B have most connect, B should have most connect to A.
        std::string mostConnectContig(std::vector<OverlapRelation> inputRelationVec);
		// return most seeds connect tgs 
		int mostConnectTGS(OverlapSeedVecInfo input);
		// return most seeds strand 
		bool mostSeedStrand(OverlapSeedVecInfo input);
		// return most tgs connect contig 
		std::string mostTGSConnectContig(OverlapSeedVecInfo input);
        // check two contigs are overlap or not
        int connectVote(size_t tgsIdx, std::vector<OverlapRelation> inputRelationVec);
		// check two contigs side are not connect, store this pair in hash if not connect
        bool checkConnect(std::string firstID, bool firstSide,std::string secondID, bool secondSide);
        //
		std::vector<ReassemblerSeedFeature> filterRepeatSeed(std::vector<ReassemblerSeedFeature> inputSeedVec);
		//
		void filterErrorStrand(OverlapSeedVecInfo input, OverlapSeedVecInfo &result, bool mostStrand);
		// 
		void filterErrorSeeds(OverlapSeedVecInfo input, OverlapSeedVecInfo &result, int mostTGS);
		//develop
		void filterErrorSeeds2(OverlapSeedVecInfo input, OverlapSeedVecInfo &result, std::string mostConnectContig);
		// third step reassembler, construct kmer hash map and tgs seeds vector
		void buildKmerHashAndTGSKmerVec( std::string tgsConnectContig, bool connectContigSide, int connectContigLength, int tgsIndex, std::vector<ReassemblerSeedFeature> inputSeedVec);
		//
		bool checkSameStrand(int startTGSIndex, int targetTGSIndex, SeedSequenceInfoVec brigdeSeedVec);
		// return overlap length
		int overlapLength(OverlapSeedVecInfo input);
		// second tgs strand will same of the first 
		
		
		std::string concordantNonOverlapFragment(OverlapSeedVecInfo input);
		
		///************************************************///
        ///******************* Variable *******************///
        ///************************************************///
		
        ContigGraph* m_pGraph;
        // collect full origin tgs connect 
        SeedByTGSHashMap collectSeedHashMap;
        // use to know how many contig connect target contig, include side infomation
        ContigRelationHashMap connectContigHashMap;
        // use to check two tgs whether can connect
        TGSConnectContigHashMap singleContigHashMap;
        // check this contig pair connect 
        SparseHashMap<std::string,bool,StringHasher> connectStatus;
        // third step reassemble, all kmer appear on the tgs information
        KmerHashMap collectAllKmerHashMap;
        // third step reassemble,	tgs index and it's seeds vector	
		TGSKmerVec AllTGSKmerVec;
 
		///************************************************///
        ///********* construct kmer distribution **********///
        ///************************************************/// 
 
        int seedCount_15;
		int seedCount_11;
        int kmerSizeFreqByTGS_15[100000];
		int kmerSizeFreqByTGS_11[100000];
		int kmerFreqByContig_15[100000];
        void countTotalSeedFreq(std::string seed);

};



#endif
