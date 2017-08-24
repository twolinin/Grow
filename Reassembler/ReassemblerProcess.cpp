///-----------------------------------------------
// 
// Written by Jyun-Hong Lin
// 
//-----------------------------------------------
//
// PacBioLongReadReassembler
//

#include "ReassemblerProcess.h"
#include "ContigVertex.h"
#include "ContigEdge.h"

		///************************************************///
        ///******* TGS Reassemble Basic Elements *********///
        ///************************************************///

TGSReassemblerBasicElements::TGSReassemblerBasicElements(const ReassemblerParameters params): m_params(params)
{
}
TGSReassemblerBasicElements::~TGSReassemblerBasicElements()
{
}

std::vector<ReassemblerSeedFeature> TGSReassemblerBasicElements::seedingByDynamicKmer(const std::string readSeq)
{
    std::vector<ReassemblerSeedFeature> seedVec;
    const size_t smallKmerSize = m_params.kmerLength;
    const size_t largeKmerSize = smallKmerSize+6;
    
    // set dynamic kmer as largest kmer initially, which will reduce size when no seeds can be found within a maximum interval
    size_t dynamicKmerSize = largeKmerSize;
    size_t kmerThreshold = m_params.seedKmerThreshold;
    
    // prevention of short reads
    if(readSeq.length() < dynamicKmerSize) return seedVec;

    // search for solid kmers and group consecutive solids kmers into one seed
    // reduce kmer size and kmerThreshold if no seeds can be found within maxSeedInterval
    for(size_t i = 0 ; i+dynamicKmerSize <= readSeq.length() ; i++)
    {
        std::string kmer = readSeq.substr(i, dynamicKmerSize);

        size_t fwdKmerFreqs = BWTAlgorithms::countSequenceOccurrencesSingleStrand(kmer, m_params.indices);
        size_t rvcKmerFreqs = BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseComplement(kmer), m_params.indices);
        size_t kmerFreqs = fwdKmerFreqs+rvcKmerFreqs;
        
         //std::cout << i << ": " << kmer << "\t" << kmerFreqs <<":" << fwdKmerFreqs << ":" << rvcKmerFreqs << "\n";
        
        if(kmerFreqs >= kmerThreshold && fwdKmerFreqs>=1 && rvcKmerFreqs>=1)
        {    
            // skip low-complexity seeds, e.g., AAAAAAAATAAAA

            if(isLowComplexity(kmer, 0.7)) continue;
            
            size_t seedStartPos = i, 
            seedLen = 0;
            
            // Group consecutive solid kmers into one seed if possible
            size_t maxKmerFreq = kmerFreqs;
            for(i++ ; i+dynamicKmerSize <= readSeq.length(); i++)
            {
                kmer = readSeq.substr(i, dynamicKmerSize);

                fwdKmerFreqs = BWTAlgorithms::countSequenceOccurrencesSingleStrand(kmer, m_params.indices);
                rvcKmerFreqs = BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseComplement(kmer), m_params.indices);
                kmerFreqs = fwdKmerFreqs + rvcKmerFreqs;
        
                // std::cout << i << ": " << kmer << "\t" << kmerFreqs <<":" << fwdKmerFreqs << ":" << rvcKmerFreqs << "\n";

                if(isLowComplexity(kmer, 0.7)) break;

                maxKmerFreq = std::max(maxKmerFreq,kmerFreqs);

                if(kmerFreqs >= (size_t) kmerThreshold && fwdKmerFreqs>=1 && rvcKmerFreqs>=1)
                    seedLen++;
                else 
                    break;
            }

            size_t seedEndPos = seedStartPos+seedLen+dynamicKmerSize-1;
                                                                    
            // refine repeat seed only if it's small size
            if(dynamicKmerSize < largeKmerSize && maxKmerFreq > 40)
            {
                refineRepeatSeed(readSeq, seedStartPos, seedEndPos);
            }
            
            if(maxKmerFreq > 80)
            {
                // push concatenated seeds into seed vector
                ReassemblerSeedFeature newSeed(seedStartPos, readSeq.substr(seedStartPos, seedEndPos-seedStartPos+1), true, dynamicKmerSize, kmerThreshold+10);
                // newSeed.setBestKmerSize(dynamicKmerSize);
                //newSeed.estimateBestKmerSize(m_params.indices.pBWT);
                seedVec.push_back(newSeed);
            }
            else
            {
                // push concatenated seeds into seed vector
                ReassemblerSeedFeature newSeed(seedStartPos, readSeq.substr(seedStartPos, seedEndPos-seedStartPos+1), false, dynamicKmerSize, kmerThreshold+10);
                //newSeed.setBestKmerSize(dynamicKmerSize);
                //newSeed.estimateBestKmerSize(m_params.indices.pBWT);
                seedVec.push_back(newSeed);
            }
            
            // std::cout << "Seed:\t " << seedVec.back().seedLength << "\t" << seedVec.back().seedStr << "\n";

            // jump to the index after new seed
            i = i + dynamicKmerSize - 2;
            
            // reset kmer size and threshold tail to strict criterion
            kmerThreshold = m_params.seedKmerThreshold;
            dynamicKmerSize = largeKmerSize;
        }// end of sufficient kmerThreshold

        /* reduce kmerThreshold or kmer size if no seeds can be found within maxSeedInterval
        It happened when kmer become larger in 2nd round
        The regions with ultra-high errors are difficult to contain seeds and were mostly not corrected in the first round
        reduce kmer size and threshold if no seeds can be found within maxSeedInterval*/
        size_t prevSeedEndPos = seedVec.empty()? 0 : seedVec.back().seedStartPos + seedVec.back().seedStr.length();
        int distToPrevSeed = i + 1 - prevSeedEndPos;
        
        if(distToPrevSeed >= (int)m_params.maxSeedInterval)
        {
            if(kmerThreshold > 10)
            {
                kmerThreshold = 10;
                i = prevSeedEndPos;
            }
    
            // previous seed is a correct repeat seed
            if( dynamicKmerSize > 17 && !seedVec.empty() && !seedVec.back().isRepeat)
            {
                dynamicKmerSize = 17;
                i = prevSeedEndPos;
            }
        }
            
    }// end of for
    
    return seedVec;
}

std::vector<ReassemblerSeedFeature> TGSReassemblerBasicElements::collectKmer(const std::string readSeq, int scanStart, int scanEnd)
{
    std::vector<ReassemblerSeedFeature> seedVec;
    const size_t kmerSize = m_params.kmerLength;

    for(size_t i = scanStart ; i+kmerSize <= (size_t)scanEnd ; i++)
    {
        std::string kmer = readSeq.substr(i, kmerSize);
        std::string reverseKmer = reverseComplement(kmer);
        
        size_t fwdKmerFreqsUsingTGS = BWTAlgorithms::countSequenceOccurrencesSingleStrand(kmer, m_params.indices);
        size_t rvcKmerFreqsUsingTGS = BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseKmer, m_params.indices);
        size_t KmerFreqsUsingTGS = fwdKmerFreqsUsingTGS + rvcKmerFreqsUsingTGS;
        
		if(KmerFreqsUsingTGS<2)continue;
		
		if(isLowComplexity(kmer, 0.7)) continue;
		
        ReassemblerSeedFeature newSeed(i, readSeq.substr(i, kmerSize), false, kmerSize, 0);
        seedVec.push_back(newSeed);    
    }
    return seedVec;
}

std::pair<size_t, size_t> TGSReassemblerBasicElements::refineRepeatSeed(const std::string readSeq, size_t& seedStartPos, size_t& seedEndPos)
{
    // initially set to max unisnged int value
    size_t newSeedStartPos = (size_t)-1;
    size_t newSeedEndPos = (size_t)-1;
    size_t startKmerFreq=0, endKmerFreq=0;
    
    const int minRepeatFreq = 40, minFreqDiff = 30;
    
    size_t kmerSize = m_params.kmerLength;
    
    std::string kmer = readSeq.substr(seedStartPos, kmerSize);
    int initKmerFreq = BWTAlgorithms::countSequenceOccurrences(kmer, m_params.indices);
    int prevKmerFreq = initKmerFreq;

    // first kmer is a repeat
    if(initKmerFreq > minRepeatFreq)
    {
        newSeedStartPos = seedStartPos;
        startKmerFreq = initKmerFreq;
    }
    
    
    // identify breakpoints of large freq difference between two kmers    
    for(size_t i=seedStartPos+1 ; i+kmerSize-1 <= seedEndPos; i++)
    {
        kmer = readSeq.substr(i, kmerSize);
        int currKmerFreq = BWTAlgorithms::countSequenceOccurrences(kmer, m_params.indices);

        // std::cout << i << ": " << kmer << "\t" << currKmerFreq << "\n";

        // error kmers within repeats often lead to large freq diff
        bool isLargeFreqDiff = currKmerFreq - prevKmerFreq > minFreqDiff;
        
        // PB36993_4517.fa, TTATGTAAGGAGTATTTGAT
        // some error kmers are with moderate frequency and no freq diff can be observed
        // pick up the first repeat kmer as starting point
        bool isRepeatKmer = (newSeedStartPos == (size_t)-1) && (currKmerFreq >= (int)minRepeatFreq);
        if(isLargeFreqDiff || isRepeatKmer)
        {
            // capture the highest repeat start pos
            bool isBetterRepeatKmer = (startKmerFreq!=0 && currKmerFreq > (int)startKmerFreq);
            if(newSeedStartPos == (size_t)-1 || isBetterRepeatKmer)
            {
                newSeedStartPos = i;
                startKmerFreq = currKmerFreq;
            }
        }
            
        // repeat end is reached
        // PB36993_4517.fa, AGGCTTGTCTGTAATCGGG
        if(prevKmerFreq - currKmerFreq > minFreqDiff /*|| currKmerFreq < minFreqDiff*/)
        {
            // do not enter unless start pos was found
            // if(newSeedStartPos != (size_t)-1)
            // {
                newSeedEndPos = i + kmerSize -2;
                endKmerFreq = prevKmerFreq;
                break;
            // }
        }
            
        prevKmerFreq = currKmerFreq;
    }
    
    if(newSeedStartPos == (size_t)-1)
    {
        newSeedStartPos = seedStartPos;
        startKmerFreq = initKmerFreq;
    }
        
    if(newSeedEndPos == (size_t)-1)
    {
        newSeedEndPos = seedEndPos;
        endKmerFreq = prevKmerFreq;
    }
    
    // std::cout << newSeedStartPos << "\t" << newSeedEndPos << "\n";

    seedStartPos = newSeedStartPos;
    seedEndPos = newSeedEndPos;
    return std::make_pair(startKmerFreq, endKmerFreq);
}

bool TGSReassemblerBasicElements::isLowComplexity (std::string seq, float threshold)
{
    size_t seqLen = seq.length();
    size_t countG =0 ;
    size_t countC =0 ;
    size_t countT =0 ;
    size_t countA =0 ;

    for (size_t i=0; i<seqLen; i++)
    {
        switch(seq[i]){
            case 'A': countA ++ ;break;
            case 'T': countT ++ ;break;
            case 'C': countC ++ ;break;
            case 'G': countG ++ ;break;
            default:  assert(false);
        }
    }

    if((float)(countG+countC)/seqLen>0.9) return true;
    
    if (  ((float) countA/seqLen >= threshold ) || ((float) countT/seqLen >= threshold)
       || ((float) countC/seqLen >= threshold ) || ((float) countG/seqLen >= threshold) )
        return true;

    return false;

}

std::pair<size_t, size_t> TGSReassemblerBasicElements::BacktrackTGSIdx(const BWTIndexSet indices, size_t InputIdx)
{
    size_t backtrackLength = 0;
    size_t idx = InputIdx;
    
    while(1)
    {
        char b = indices.pBWT->getChar(idx);
        size_t new_idx = indices.pBWT->getPC(b) + indices.pBWT->getOcc(b, idx - 1);
        if(b == '$')
        {
            idx = indices.pSSA->lookupLexoRank(new_idx);
            return std::make_pair(idx,backtrackLength);
        }
        else idx = new_idx;
        
        backtrackLength++;
        if( idx == InputIdx )
        {
            return std::make_pair(0,backtrackLength);
        }        
    }
}

std::string TGSReassemblerBasicElements::backtrackTGS(const BWTIndexSet indices, int64_t inputIdx)
{
    std::string backtrackTGS;
            
    while(1)
    {
        char b = indices.pRBWT->getChar(inputIdx);
        int64_t new_idx = indices.pRBWT->getPC(b) + indices.pRBWT->getOcc(b, inputIdx - 1);
        
        if(b == '$')break;
        else inputIdx = new_idx;
        
        backtrackTGS += b ;
    }
    return backtrackTGS;
}

		///************************************************///
        ///*********** TGS Reassembler Process ************///
        ///************************************************///

ReassemblerProcess::ReassemblerProcess(const ReassemblerParameters params) 
: TGSReassemblerBasicElements(params)
{
}
ReassemblerProcess::~ReassemblerProcess()
{
}

ReassemblerResult ReassemblerProcess::PBReassembler(const SequenceWorkItem& workItem)
{    
    //assume contig length bigger than 2n
    //we will use head and tail n length to find seeds
    ReassemblerResult result;
    SeedSequenceInfoVec seedVec;
    
    std::vector<ReassemblerSeedFeature> headSeedVec,tailSeedVec;
    // get origin contig sequence
    std::string contigSeq = workItem.read.seq.toString();
    
    // skip tgs, if tgs length less than 2n
    if( (int)contigSeq.length() < m_params.tgsAvgLen*2 ) return result;

    // cut head and tail n length
    std::string headContigSeq = contigSeq.substr(0,m_params.tgsAvgLen);
    std::string tailContigSeq = contigSeq.substr(contigSeq.length()-m_params.tgsAvgLen,m_params.tgsAvgLen);

	// using recheck overlap position
	std::string headPartialSeq = contigSeq.substr(0,400);
    std::string tailPartialSeq = contigSeq.substr(contigSeq.length()-400,400);
	
    // find seeds on head and tail n length
    headSeedVec = seedingByDynamicKmer(headContigSeq);
    tailSeedVec = seedingByDynamicKmer(tailContigSeq);

    // push head and tail seeds to seed vector
    findSeedOnLongRead( m_params.indices, workItem.read.id, headSeedVec, seedVec, contigSeq.length(), headPartialSeq, true);
    findSeedOnLongRead( m_params.indices, workItem.read.id, tailSeedVec, seedVec, contigSeq.length(), tailPartialSeq, false);
    
    // insert vector seeds to hash
    for(SeedSequenceInfoVec::iterator iter = seedVec.begin(); iter != seedVec.end(); ++iter )
    {
        SeedByTGSHashMap::iterator pbIterator = result.seedHash.find((*iter).tgsIndex);
        
        // filter short seed
        if( (*iter).seedLength<=15 )continue;

        //this tgs exist
        if( pbIterator!=result.seedHash.end() )
        {
            pbIterator->second[0].SeedInfoVec.push_back((*iter));
        }
        // this tgs not exist
        else
        {
            OverlapSeedVecInfo pushContig;
            std::vector<OverlapSeedVecInfo> contigVec;
            pushContig.contig = (*iter).contigID;
            pushContig.SeedInfoVec.push_back((*iter));
            contigVec.push_back(pushContig);
            result.seedHash.insert(std::make_pair((*iter).tgsIndex,contigVec));
        }
    }
    return result;
}

void ReassemblerProcess::findSeedOnLongRead(const BWTIndexSet indices, std::string contigID, std::vector<ReassemblerSeedFeature> inputSeedVec, std::vector<SeedSequenceInfo> &result, int contigLength, std::string contigPartialSeq, bool contigSide)
{
    for( std::vector<ReassemblerSeedFeature>::iterator iter = inputSeedVec.begin(); iter != inputSeedVec.end(); ++iter )
    {    
        SeedSequenceInfo currentSeedInfo;
        size_t contigStartPos = contigSide ? (*iter).seedStartPos : contigLength - m_params.tgsAvgLen + (*iter).seedStartPos;
        size_t contigEndPos   = contigSide ? (*iter).seedEndPos   : contigLength - m_params.tgsAvgLen + (*iter).seedEndPos;
        
        BWTInterval interval  = BWTAlgorithms::findInterval(indices,(*iter).seedStr) ;
        BWTInterval rinterval = BWTAlgorithms::findInterval(indices,reverseComplement((*iter).seedStr)) ;
        
        size_t fwdSeedFreqs = BWTAlgorithms::countSequenceOccurrencesSingleStrand((*iter).seedStr, m_params.indices);
        size_t rvcSeedFreqs = BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseComplement((*iter).seedStr), m_params.indices);
        size_t SeedFreqs = fwdSeedFreqs + rvcSeedFreqs;

        currentSeedInfo.setContigSeq(contigID, contigPartialSeq, contigLength, contigStartPos, contigEndPos, contigSide);                          

        for(size_t j = interval.lower ; j <= (size_t)interval.upper ; j++)
        {
            std::pair<size_t, size_t> currentSeed = BacktrackTGSIdx(indices,j);
            currentSeedInfo.setTGSSeq((*iter).seedStr, SeedFreqs, currentSeed.first, currentSeed.second, true);
            result.push_back(currentSeedInfo);
        }
            
        for(size_t j = rinterval.lower ; j <= (size_t)rinterval.upper ; j++)
        {
            std::pair<size_t, size_t> currentSeed = BacktrackTGSIdx(indices,j);
            currentSeedInfo.setTGSSeq(reverseComplement((*iter).seedStr), SeedFreqs, currentSeed.first, currentSeed.second, false);
            result.push_back(currentSeedInfo);
        }
    }
}

		///************************************************///
        ///********* TGS Reassembler Post Process *********///
        ///************************************************///

ReassemblerPostProcess::ReassemblerPostProcess(const ReassemblerParameters params, ContigGraph* pGraph) 
: TGSReassemblerBasicElements(params),m_pGraph(pGraph)
{
    collectSeedHashMap.set_deleted_key(-1);
    connectContigHashMap.set_deleted_key("");
    
    singleContigHashMap.set_deleted_key(-1);
    collectAllKmerHashMap.set_deleted_key("");
    
    seedCount_15 = 0;
	seedCount_11 = 0;
}
ReassemblerPostProcess::~ReassemblerPostProcess()
{
}
		
void ReassemblerPostProcess::process(const SequenceWorkItem& item, /*const*/ ReassemblerResult& result)
{    
    for(SeedByTGSHashMap::iterator currentPbIdx = result.seedHash.begin() ; currentPbIdx != result.seedHash.end() ; ++currentPbIdx)
    {
        for(SeedSequenceInfoVec::iterator iter = currentPbIdx->second[0].SeedInfoVec.begin(); iter !=currentPbIdx->second[0].SeedInfoVec.end() ; ++iter)
            if( seedCount_15<100000 && seedCount_11<100000) countTotalSeedFreq((*iter).seedString);
        
        SeedByTGSHashMap::iterator pbIterator = collectSeedHashMap.find( currentPbIdx->first );
        //this tgs already exist
        if( pbIterator != collectSeedHashMap.end() ) 
        {
            // push current contig seeds to hash
            pbIterator->second.push_back(currentPbIdx->second[0]);
        }
        //this tgs not exist
        else
        {
            collectSeedHashMap.insert(std::make_pair(currentPbIdx->first,currentPbIdx->second));
        }
    }
}

void ReassemblerPostProcess::filterErrorTGS()
{
    std::cout<< "\n------------filter error tgs------------\n";
	
	ResultCount onePB = ResultCount();
	ResultCount twoPB = ResultCount();

    int tmp = 0;
    for(SeedByTGSHashMap::iterator pbIterator = collectSeedHashMap.begin() ; pbIterator!=collectSeedHashMap.end() ; ++pbIterator )
    {    
        std::cout<< std::flush << '\r' << "total overlap between contig and tgs : (" << ++tmp << "/" << collectSeedHashMap.size() << ")";
        
        std::vector<OverlapSeedVecInfo>::iterator firstContigIter  ;
        std::vector<OverlapSeedVecInfo>::iterator secondContigIter ;

		//std::vector<OverlapSeedVecInfo>::iterator t = pbIterator->second.begin();
		//std::cout<< (*t).SeedInfoVec.size() << "\n";
		
        if( pbIterator->second.size()==2 )
        {   
            onePB.totalCount++;
            // this tgs overlap two contigs
        
            std::vector<OverlapSeedVecInfo>::iterator contigIter = pbIterator->second.begin();
            firstContigIter  = contigIter;
            secondContigIter = ++contigIter;
        }
        else if(pbIterator->second.size()>2)
        {   
			onePB.partialCount++;
            // multiple align,
            // this tgs overlap three or more contigs
            // select the most seeds and the second most seeds
            
            size_t firstSize  = 0;
            size_t secondSize = 0;
			
            for(std::vector<OverlapSeedVecInfo>::iterator contigIter = pbIterator->second.begin() ; contigIter!=pbIterator->second.end() ; ++contigIter)
            {
                if( firstSize < (*contigIter).SeedInfoVec.size())
                {
                    secondSize = firstSize;
                    secondContigIter = firstContigIter;

                    firstSize = (*contigIter).SeedInfoVec.size();
                    firstContigIter = contigIter;
                }
                else if( secondSize < (*contigIter).SeedInfoVec.size())
                {
                    secondSize = (*contigIter).SeedInfoVec.size();
                    secondContigIter = contigIter;
                }
            }
        }
        else 
        {
            if( !m_params.isThird ) continue;
            twoPB.totalCount++;

            TGSScanAndOverlapPosition singleOverlap;
            std::vector<OverlapSeedVecInfo>::iterator contigIter = pbIterator->second.begin();

            // each contig have three or more seeds on this tgs
            if((*contigIter).SeedInfoVec.size()<=3){ twoPB.seedsThreshold++;        continue; }
            // filter tgs if it concurrently exist different strand seed and all seeds length are rather than 19
			if(checkContigSeed(*contigIter,19))    { twoPB.seedDifferentStrand++;   continue; }
            // check repeat by contig index.
			if(checkRepeatByContig(*contigIter))   { twoPB.repeatSeeds++;           continue; }
            // check seeds orientations concordant or discordant 
			if(checkSeedsOrientations(*contigIter)){ twoPB.seedOrderDiscordant++;   continue; }
            // check distance between to seeds
			if(checkSeedsDistance(*contigIter))    { twoPB.abnormalDistance++;      continue; }
            // check overlap length 
			if(checkOverlapLength(*contigIter))    { twoPB.abnormalOverlapLength++; continue; }
            // check nonOverlap length
			if(checkNonOverlapPartialLength(*contigIter,singleOverlap)){twoPB.insufficientLength++; continue; }

            twoPB.passCount++;
            
            insertSingleRelationHash(*contigIter,singleOverlap);
        
            continue;
        }
        
        if( m_params.isThird ) continue;
        
        std::string firstContig  = (*firstContigIter).contig;
        std::string secondContig = (*secondContigIter).contig;
        bool firstStrand         = (*firstContigIter).SeedInfoVec[0].strand;
        bool secondStrand        = (*secondContigIter).SeedInfoVec[0].strand;
        bool firstContigSide     = (*firstContigIter).SeedInfoVec[0].contigSide;
        bool secondContigSide    = (*secondContigIter).SeedInfoVec[0].contigSide;
        int firstContigLength    = (*firstContigIter).SeedInfoVec[0].contigLength;
        int secondContigLength   = (*secondContigIter).SeedInfoVec[0].contigLength;
		
        //each contig have three or more seeds on this tgs
        if((*firstContigIter).SeedInfoVec.size()>=3 && (*secondContigIter).SeedInfoVec.size()>=3)
        {  
            // prevent large repeat, filter tgs if it concurrently exist same/different side and strand  
            if(checkContigRelation(firstContigSide,firstStrand,secondContigSide,secondStrand))           { onePB.repeatStatus++;          continue; }
            // check repeat by tgs index, old version. It can't detect every repeat seed.
            //if(checkRepeatByTGS((*firstContigIter)) || checkRepeatByTGS((*secondContigIter)))continue;
            // check repeat by contig index
            if(checkRepeatByContig((*firstContigIter)) || checkRepeatByContig((*secondContigIter)))      { onePB.repeatSeeds++;           continue; }
            // filter tgs if it concurrently exist different strand seed and all seeds length are rather than 19
            if(checkContigSeed((*firstContigIter),19) || checkContigSeed((*secondContigIter),19))        { onePB.seedDifferentStrand++;   continue; }
            //check seeds orientations and seeds distance between contig and tgs 
            if(checkSeedsOrientations((*firstContigIter)) || checkSeedsOrientations((*secondContigIter))){ onePB.seedOrderDiscordant++;   continue; }
            //check distance between to seeds
            if(checkSeedsDistance(*firstContigIter) || checkSeedsDistance(*secondContigIter))            { onePB.abnormalDistance++;      continue; }
            //check overlap length 
            if(checkOverlapLength((*firstContigIter)) || checkOverlapLength((*secondContigIter)))        { onePB.abnormalOverlapLength++; continue; }
            
            onePB.passCount++;
            
            std::vector<OverlapSeedVecInfo> tmp;
                        
            tmp.push_back((*firstContigIter));
            tmp.push_back((*secondContigIter));
                
            // build relation hash, use to check how many contig are connect the other contig
            insertContigRelationHash(firstContig,firstContigSide,firstStrand,secondContig,secondContigSide,secondStrand,firstContigLength,secondContigLength,tmp);
            insertContigRelationHash(secondContig,secondContigSide,secondStrand,firstContig,firstContigSide,firstStrand,secondContigLength,firstContigLength,tmp);
        }
        else onePB.seedsThreshold++;
    }
	
    // 15mer distribution using tgs
    std::sort(kmerSizeFreqByTGS_15,kmerSizeFreqByTGS_15+seedCount_15);
    std::sort(kmerFreqByContig_15,kmerFreqByContig_15+seedCount_15);
    std::cout << "\ntgs's 15mer freq by tgs indices :\n"  
              << kmerSizeFreqByTGS_15[seedCount_15/10*1] << "\t" << kmerSizeFreqByTGS_15[seedCount_15/10*2] << "\t"
              << kmerSizeFreqByTGS_15[seedCount_15/10*3] << "\t" << kmerSizeFreqByTGS_15[seedCount_15/10*4] << "\t"
              << kmerSizeFreqByTGS_15[seedCount_15/10*5] << "\t" << kmerSizeFreqByTGS_15[seedCount_15/10*6] << "\t"
              << kmerSizeFreqByTGS_15[seedCount_15/10*7] << "\t" << kmerSizeFreqByTGS_15[seedCount_15/10*8] << "\t"
              << kmerSizeFreqByTGS_15[seedCount_15/10*9] << "\t" << kmerSizeFreqByTGS_15[seedCount_15-1]    << "\n";
                                        
    std::cout << "tgs's 15mer freq by contig indices :\n" 
              << kmerFreqByContig_15[seedCount_15/10*1] << "\t" << kmerFreqByContig_15[seedCount_15/10*2] << "\t"
              << kmerFreqByContig_15[seedCount_15/10*3] << "\t" << kmerFreqByContig_15[seedCount_15/10*4] << "\t"
              << kmerFreqByContig_15[seedCount_15/10*5] << "\t" << kmerFreqByContig_15[seedCount_15/10*6] << "\t"
              << kmerFreqByContig_15[seedCount_15/10*7] << "\t" << kmerFreqByContig_15[seedCount_15/10*8] << "\t"
              << kmerFreqByContig_15[seedCount_15/10*9] << "\t" << kmerFreqByContig_15[seedCount_15-1]    << "\n";    
    
    std::cout<<"\n";
    if( m_params.isThird )
    {
        std::cout<< "tgs align one contig : "         << twoPB.totalCount            << "\n"
                 << "pass tgs number      : "         << twoPB.passCount             << "\n"
                 << "case 1. two flank are repeats   : " << twoPB.repeatSeeds           << "\n"
                 << "case 2. seeds num less than two : " << twoPB.seedsThreshold        << "\n"
                 << "case 3. seeds different strand  : " << twoPB.seedDifferentStrand   << "\n"
                 << "case 4. seeds order discordant  : " << twoPB.seedOrderDiscordant   << "\n"
                 << "case 5. abnormal overlap length : " << twoPB.abnormalOverlapLength << "\n"
                 << "case 6. abnormal distance between two seeds : " << twoPB.abnormalDistance   << "\n"
                 << "case 7. insufficient non overlap length     : " << twoPB.insufficientLength << "\n";
    }    
    if( m_params.isFirst || m_params.isSecond )
    {
        std::cout<< "tgs align two contig    : "      << onePB.totalCount            << "\n"
                 << "align more than two contig : "      << onePB.partialCount          << "\n"
                 << "pass tgs number      : "         << onePB.passCount             << "\n"
                 << "case 1. two reads are same side : " << onePB.repeatStatus          << "\n"
                 << "case 2. two flank are repeats   : " << onePB.repeatSeeds           << "\n"
                 << "case 3. seeds num less than two : " << onePB.seedsThreshold        << "\n"
                 << "case 4. seeds different strand  : " << onePB.seedDifferentStrand   << "\n"
                 << "case 5. seeds order discordant  : " << onePB.seedOrderDiscordant   << "\n"
                 << "case 6. abnormal overlap length : " << onePB.abnormalOverlapLength << "\n"
                 << "case 7. abnormal distance between two seeds : " << onePB.abnormalDistance << "\n"; 
    }
}

void ReassemblerPostProcess::buildGraphByOneTGS()
{
    std::cout<< "\n---------build graph by one tgs---------\n";
    std::cout<< "contig side number:\t" << connectContigHashMap.size() << "\n";
    
    for(ContigRelationHashMap::iterator contigSideIter = connectContigHashMap.begin() ; contigSideIter!=connectContigHashMap.end() ; ++contigSideIter )
    {
        // contigA have most tgs to contigB, so and contigB
        
        std::vector<OverlapRelation>::iterator iter ;
        std::string endContig = mostConnectContig(contigSideIter->second);

        if( endContig=="" ) continue;
        else
        {
            std::string startContig;
            int tmpOverlapLength = -1;

            int ec_same = 0;
            int ec_rev  = 0;
            bool ec_vote = true;
            
            for(std::vector<OverlapRelation>::iterator tmpIter = contigSideIter->second.begin(); tmpIter != contigSideIter->second.end(); ++tmpIter )
            {
                if((*tmpIter).secondContig == endContig)
                {
                    if((*tmpIter).isSameStrand) ec_same++;
                    else ec_rev++;
                }
            }
            
            if(ec_same == ec_rev) continue;
            if(ec_same < ec_rev) ec_vote = false;
            
            for(std::vector<OverlapRelation>::iterator tmpIter = contigSideIter->second.begin(); tmpIter != contigSideIter->second.end(); ++tmpIter )
            {
                if( (*tmpIter).secondContig == endContig )
                {    
                    if((*tmpIter).isSameStrand != ec_vote) continue;
                    
                    if( ((*tmpIter).overlapLength > tmpOverlapLength) && ((*tmpIter).firstContig == startContig || tmpOverlapLength == -1) )
                    {
                        iter = tmpIter;
                        
                        tmpOverlapLength = (*tmpIter).overlapLength;
                        
                        std::string insertKey = (*tmpIter).secondContig+((*iter).secondSide ? "_Head":"_Tail");
                    
                        ContigRelationHashMap::iterator contigIterator = connectContigHashMap.find(insertKey);
                        
                        if( contigIterator != connectContigHashMap.end() ) 
                            startContig = mostConnectContig(contigIterator->second);
                    }
                }
            }
            
            if( startContig != (*iter).firstContig ) continue;
        }
        
        std::string firstContig  = (*iter).firstContig;
        std::string secondContig = (*iter).secondContig;
        bool firstSide    = (*iter).firstSide;
        bool secondSide   = (*iter).secondSide;
        bool firstStrand  = (*iter).firstStrand;
        bool secondStrand = (*iter).secondStrand;
        
        int connectAndOverlap = connectVote((*iter).tgsIndex,contigSideIter->second);
        // 0:same number of overlap and connect 1:overlap 2:nonOverlap
        if( connectAndOverlap == 0) continue;
        // check two contigs side are not connect, store this pair in hash if not connect 
        if(checkConnect(firstContig,firstSide,secondContig,secondSide))continue;
        
        if( m_params.isFirst && connectAndOverlap!=1 )continue;
        
        if( m_params.isSecond && connectAndOverlap!=2 )continue;
        
            std::cout//<< contigSideIter->first << "\t"
                     //<< contigSideIter->second.size() << "\t"
                     << (connectAndOverlap-1 ? "->...<-" : "-><-   ") << "\t"
                     << (*iter).tgsIndex              << "\t"
                     //<< (*iter).tgsLength             << "\t"
                     //<< (*iter).overlapLength            << "\t"
                     << firstContig                      << "\t"
                     << (firstSide    ? "Head" : "Tail") << "\t"
                     << (firstStrand  ? "NonRC" : "RC")  << "\t"
                     //<< (*iter).firstContigLength        << "\t"
                     //<< (firstSide ? "anti" : "sense")   << "\t"
                     << secondContig                     << "\t"
                     << (secondSide   ? "Head" : "Tail") << "\t"
                     << (secondStrand ? "NonRC" : "RC")  << "\t"
                     //<< (*iter).secondContigLength       << "\t"
                     //<< (secondSide ? "anti" : "sense")  << "\t"
                     << ((*iter).isSameStrand ? "same" : "reverse")      << "\t"
                     << "\n";
            
            SeqCoord firstSeqCoord((*iter).firContigStart,(*iter).firContigEnd,(*iter).firstContigLength); 
            SeqCoord secondSeqCoord((*iter).secContigStart,(*iter).secContigEnd,(*iter).secondContigLength);
            
            ContigVertex* pVL = ((ContigVertex*)m_pGraph->getVertex(firstContig));
            ContigVertex* pVR = ((ContigVertex*)m_pGraph->getVertex(secondContig));
            
            ContigEdge* edgeLR = new ContigEdge(pVR,
                                                firstSide ? ED_ANTISENSE : ED_SENSE,
                                                (*iter).isSameStrand ? EC_SAME : EC_REVERSE,
                                                firstSeqCoord,
                                                firstStrand,
                                                secondStrand,
                                                (*iter).tgsFragment,
                                                (*iter).overlapLength);
            ContigEdge* edgeRL = new ContigEdge(pVL,
                                                secondSide ? ED_ANTISENSE : ED_SENSE,
                                                (*iter).isSameStrand ? EC_SAME : EC_REVERSE,
                                                secondSeqCoord,
                                                secondStrand,
                                                firstStrand,
                                                (*iter).tgsFragment,
                                                (*iter).overlapLength);
            
            edgeLR->setTwin(edgeRL);
            edgeRL->setTwin(edgeLR);
     
            m_pGraph->addEdge(pVL,edgeLR);
            m_pGraph->addEdge(pVR,edgeRL);    
    }
}

void ReassemblerPostProcess::buildGraphByTwoTGS()
{
    std::cout << "\n---------build graph by two tgs---------\n";
	
	std::sort(kmerSizeFreqByTGS_11,kmerSizeFreqByTGS_11+seedCount_11);
    std::cout << "tgs's 11mer freq by tgs indices :\n"  
              << kmerSizeFreqByTGS_11[seedCount_11/10*1] << "\t" << kmerSizeFreqByTGS_11[seedCount_11/10*2] << "\t"
              << kmerSizeFreqByTGS_11[seedCount_11/10*3] << "\t" << kmerSizeFreqByTGS_11[seedCount_11/10*4] << "\t"
              << kmerSizeFreqByTGS_11[seedCount_11/10*5] << "\t" << kmerSizeFreqByTGS_11[seedCount_11/10*6] << "\t"
              << kmerSizeFreqByTGS_11[seedCount_11/10*7] << "\t" << kmerSizeFreqByTGS_11[seedCount_11/10*8] << "\t"
              << kmerSizeFreqByTGS_11[seedCount_11/10*9] << "\t" << kmerSizeFreqByTGS_11[seedCount_11-1]    << "\n";
    
	
	ResultCount twoPB = ResultCount();
	
    int tmp = 0;

    // collect all seeds on tgs fragment
    for(TGSConnectContigHashMap::iterator tgsIter = singleContigHashMap.begin() ; tgsIter!=singleContigHashMap.end() ; ++tgsIter )
    {
		std::cout<< std::flush << '\r' << "collect seeds on tgs fragment : (" << ++tmp << "/" << singleContigHashMap.size() << ")";
        std::string originTGS  = backtrackTGS( m_params.indices, tgsIter->first );
        std::string connectContig = tgsIter->second[0].firstContig;
        bool connectContigSide    = tgsIter->second[0].firstSide;
        int connectContigLength   = tgsIter->second[0].firstContigLength;
		// check palindrome tgs
		if(checkPalindrome(originTGS)) 
		{ 
			//std::cout<< tgsIter->first << "\n";
			twoPB.palindrome++; 
			continue; 
		}
		// collect seeds on tgs
		std::vector<ReassemblerSeedFeature> tgsSeedVec = collectKmer(originTGS,tgsIter->second[0].singleOverlap.scanTGSStart, tgsIter->second[0].singleOverlap.scanTGSEnd);
		// filter repeat seed
		std::vector<ReassemblerSeedFeature> filterTGSSeedVec = filterRepeatSeed(tgsSeedVec);
        buildKmerHashAndTGSKmerVec( connectContig, connectContigSide, connectContigLength, tgsIter->first, filterTGSSeedVec);
		//buildKmerHashAndTGSKmerVec( connectContig, connectContigSide, connectContigLength, tgsIter->first, tgsSeedVec);

		twoPB.passCount++;		
    }
	std::cout<< "\n";
	std::cout<< "filter palindrome : " << twoPB.palindrome << "\n\n";

    // filter error relation connect between two tgss
    for(TGSKmerVec::iterator tgsVecIter = AllTGSKmerVec.begin() ; tgsVecIter!=AllTGSKmerVec.end() ; ++tgsVecIter)
    {
        // use to record already connect tgss
		SparseHashMap<int, int, Int64Hasher> connectTGSCount;
        // use to find current tgs connect contig
        TGSConnectContigHashMap::iterator selfTGSIter = singleContigHashMap.find((*tgsVecIter).first);
        // current tgs connect contig
        std::string selfContig = selfTGSIter->second[0].firstContig; 
		// use to collect useful seeds
		OverlapSeedVecInfo currentSeedVec;
		// after filter, select most seeds connect tgs
		OverlapSeedVecInfo resultSeedVec;
		// look current tgs seeds
        for(KmerInfoVec::iterator selfSeedIter = tgsVecIter->second.begin() ; selfSeedIter != tgsVecIter->second.end() ; ++selfSeedIter)
        {
			KmerHashMap::iterator seedInfo = collectAllKmerHashMap.find((*selfSeedIter).kmerStr);
            // confirm this seed exist in hash
            if( seedInfo != collectAllKmerHashMap.end() )
            {
				SeedSequenceInfo currentSeed;
				// use to prevent repeat kmer
				int connectOtherTGSCount = 0;
				// look for other tgs if exist current seed
				for(std::vector<KmerInfo>::iterator otherKmerIter = seedInfo->second.begin() ; otherKmerIter != seedInfo->second.end() ; ++otherKmerIter )
                {	
					// prevent two tgs connect same contig
					if( (*otherKmerIter).tgsConnectContigID == selfContig ) continue;
                    
					connectOtherTGSCount++;

					currentSeed.contigID       = (*otherKmerIter).tgsConnectContigID;
					currentSeed.contigLength   = (*otherKmerIter).tgsConnectContigLength;
					currentSeed.contigSide     = (*otherKmerIter).tgsConnectContigSide;
					currentSeed.seedLength     = (*otherKmerIter).kmerStr.length();
					currentSeed.seedString     = (*otherKmerIter).kmerStr;
					currentSeed.strand         = (*otherKmerIter).kmerStrand;
					currentSeed.tgsLocation = (*otherKmerIter).position;
					currentSeed.tgsIndex    = (*otherKmerIter).tgsIndex;
					// contigStartPos is self tgs position
					currentSeed.contigStartPos = (*selfSeedIter).position;
                }
				// filter repeat seed
				if( connectOtherTGSCount>1 || connectOtherTGSCount==0 )continue;
				// push valid seed to vector
				currentSeedVec.SeedInfoVec.push_back(currentSeed);
            }
        }
		
		if( currentSeedVec.SeedInfoVec.size() == 0 ){ twoPB.zeroSeed++;        continue; }
		
		filterErrorSeeds(currentSeedVec,resultSeedVec,mostConnectTGS(currentSeedVec));
		
		if( resultSeedVec.SeedInfoVec.size() < 10 ) { twoPB.seedsThreshold++;        continue; }
		// filter tgs if it concurrently exist different strand seed and all seeds length are rather than 14
        if(checkContigSeed(resultSeedVec,14))       { twoPB.seedDifferentStrand++;   continue; }
		// check repeat by tgs index, old version. It can't detect every repeat seed.
		if(checkRepeatByTGS(resultSeedVec))      { twoPB.repeatSeeds++;           continue; }
		// check seeds order concordant or discordant 
        if(checkSeedsOrientations(resultSeedVec))   { twoPB.seedOrderDiscordant++;   continue; }
		// check distance between to seeds
        if(checkSeedsDistance(resultSeedVec))       { twoPB.abnormalDistance++;      continue; }
		// check overlap length 
        if(checkOverlapLength(resultSeedVec))       { twoPB.abnormalOverlapLength++; continue; }
		
		TGSConnectContigHashMap::iterator targetTGSIter = singleContigHashMap.find(resultSeedVec.SeedInfoVec[0].tgsIndex);
		
		OverlapRelation connectResult;

		connectResult.firstContig        = selfTGSIter->second[0].firstContig;
		connectResult.firstContigLength  = selfTGSIter->second[0].firstContigLength;
		connectResult.firstSide          = selfTGSIter->second[0].firstSide;
		connectResult.firstStrand        = selfTGSIter->second[0].firstStrand;
		
		connectResult.secondContig       = targetTGSIter->second[0].firstContig;
		connectResult.secondContigLength = targetTGSIter->second[0].firstContigLength;
		connectResult.secondSide         = targetTGSIter->second[0].firstSide;
		connectResult.secondStrand       = (resultSeedVec.SeedInfoVec[0].strand ?  targetTGSIter->second[0].firstStrand : !targetTGSIter->second[0].firstStrand );

		connectResult.firContigStart = 0;
		connectResult.firContigEnd   = 0;
		connectResult.secContigStart = 0;
		connectResult.secContigEnd   = 0;
		connectResult.overlapLength  = 0;
		
		connectResult.connect = true;
		
		connectResult.isSameStrand = checkSameStrand( (*tgsVecIter).first, resultSeedVec.SeedInfoVec[0].tgsIndex, resultSeedVec.SeedInfoVec);

		// prevent large repeat, filter tgs if it concurrently exist same/different side and strand 
		if( checkContigRelation( connectResult.firstSide, connectResult.firstStrand, connectResult.secondSide, connectResult.secondStrand )) {twoPB.repeatStatus++;continue;}
		
		// check two contigs side are not connect, store this pair in hash if not connect 
		if( checkConnect(connectResult.firstContig,connectResult.firstSide,connectResult.secondContig,connectResult.secondSide))continue;
		
		std::string firstTGSFragment;
		std::string secondTGSFragment = concordantNonOverlapFragment(resultSeedVec);
		std::string originTGS  = backtrackTGS( m_params.indices, (*tgsVecIter).first );
		
		if( selfTGSIter->second[0].singleOverlap.scanTGSStart == 0 ) 
		{
			SeedSequenceInfoVec::iterator first = resultSeedVec.SeedInfoVec.begin();
			SeedSequenceInfoVec::iterator last = resultSeedVec.SeedInfoVec.end();
			last--;
			
			int start = (*first).contigStartPos > (*last).contigStartPos ? (*first).contigStartPos : (*last).contigStartPos;
			int end   = selfTGSIter->second[0].singleOverlap.scanTGSEnd;
			firstTGSFragment  = originTGS.substr( start, end - start );
			
			//firstTGSFragment  = originTGS.substr( selfTGSIter->second[0].singleOverlap.scanTGSStart, selfTGSIter->second[0].singleOverlap.scanTGSEnd );
			connectResult.tgsFragment = secondTGSFragment + firstTGSFragment;
			
			//std::cout<< "start1\n";
			//std::cout<< (*first).contigStartPos << "\t" << (*last).contigStartPos <<"\n";
			//std::cout<< start << "\t" << end << "\t idx " << (*tgsVecIter).first << "\n";
		}
		else 
		{
			SeedSequenceInfoVec::iterator first = resultSeedVec.SeedInfoVec.begin();
			SeedSequenceInfoVec::iterator last = resultSeedVec.SeedInfoVec.end();
			last--;
			
			int start = selfTGSIter->second[0].singleOverlap.scanTGSStart;
			int end   = (*first).contigStartPos < (*last).contigStartPos ? (*first).contigStartPos : (*last).contigStartPos;
			
			firstTGSFragment  = originTGS.substr( start, end - start );
			
			//firstTGSFragment  = originTGS.substr( selfTGSIter->second[0].singleOverlap.scanTGSStart, originTGS.length() - selfTGSIter->second[0].singleOverlap.scanTGSStart );
			connectResult.tgsFragment = firstTGSFragment + secondTGSFragment;
			
			//std::cout<< "end1\n";
			//std::cout<< (*first).contigStartPos << "\t" << (*last).contigStartPos <<"\n";
			//std::cout<< start << "\t" << end << "\t idx " << (*tgsVecIter).first << "\n";
		}
		
		std::cout    << connectResult.firstContig                      << "\t"
                     << (connectResult.firstSide    ? "Head" : "Tail") << "\t"
                     << (connectResult.firstStrand  ? "NonRC" : "RC")  << "\t"
                     << connectResult.firstContigLength        << "\t"
                     //<< (firstSide ? "anti" : "sense")   << "\t"
                     << connectResult.secondContig                     << "\t"
                     << (connectResult.secondSide   ? "Head" : "Tail") << "\t"
                     << (connectResult.secondStrand ? "NonRC" : "RC")  << "\t"
                     << connectResult.secondContigLength       << "\t"
                     //<< (secondSide ? "anti" : "sense")  << "\t"
                     << (connectResult.isSameStrand ? "same" : "reverse")      << "\t"
					 << (*tgsVecIter).first << "\t"
					 << resultSeedVec.SeedInfoVec[0].tgsIndex << "\t" 
					 << resultSeedVec.SeedInfoVec.size() << "\t"
                     << "\n";
		
		SeqCoord firstSeqCoord(connectResult.firContigStart,connectResult.firContigEnd,connectResult.firstContigLength); 
		SeqCoord secondSeqCoord(connectResult.secContigStart,connectResult.secContigEnd,connectResult.secondContigLength);
				
		ContigVertex* pVL = ((ContigVertex*)m_pGraph->getVertex(connectResult.firstContig));
		ContigVertex* pVR = ((ContigVertex*)m_pGraph->getVertex(connectResult.secondContig));
				
		ContigEdge* edgeLR = new ContigEdge(pVR,
											connectResult.firstSide ? ED_ANTISENSE : ED_SENSE,
											connectResult.isSameStrand ? EC_SAME : EC_REVERSE,
											firstSeqCoord,
											connectResult.firstStrand,
											connectResult.secondStrand,
											connectResult.tgsFragment,
											connectResult.overlapLength);
		ContigEdge* edgeRL = new ContigEdge(pVL,
											connectResult.secondSide ? ED_ANTISENSE : ED_SENSE,
											connectResult.isSameStrand ? EC_SAME : EC_REVERSE,
											secondSeqCoord,
											connectResult.secondStrand,
											connectResult.firstStrand,
											connectResult.tgsFragment,
											connectResult.overlapLength);
				
		edgeLR->setTwin(edgeRL);
		edgeRL->setTwin(edgeLR);
		 
		m_pGraph->addEdge(pVL,edgeLR);
		m_pGraph->addEdge(pVR,edgeRL);
	}
	
	std::cout<<"\n";
	std::cout<< "pass tgs number         : "      << twoPB.passCount             << "\n"
             << "case 1. repeat status      : "      << twoPB.repeatStatus          << "\n"
			 << "case 2. zero seed          : "      << twoPB.zeroSeed              << "\n"
             << "case 3. seeds num less than ten : " << twoPB.seedsThreshold        << "\n"
             << "case 4. seeds different strand  : " << twoPB.seedDifferentStrand   << "\n"
             << "case 5. seeds order discordant  : " << twoPB.seedOrderDiscordant   << "\n"
             << "case 6. abnormal overlap length : " << twoPB.abnormalOverlapLength << "\n"
             << "case 7. abnormal distance between two seeds : " << twoPB.abnormalDistance << "\n\n"; 

}

		///******** filter error tgs function **********///

bool ReassemblerPostProcess::checkContigRelation(bool firstSide, bool firstStrand, bool secondSide, bool secondStrand)
{
    // same side and same strand
    if(  firstSide &&  secondSide &&  firstStrand &&  secondStrand)    return true;
    if(  firstSide &&  secondSide && !firstStrand && !secondStrand)    return true;
    if( !firstSide && !secondSide &&  firstStrand &&  secondStrand)    return true;
    if( !firstSide && !secondSide && !firstStrand && !secondStrand)    return true;
    
    // different side and different strand
    if(  firstSide && !secondSide &&  firstStrand && !secondStrand)    return true;
    if(  firstSide && !secondSide && !firstStrand &&  secondStrand)    return true;
    if( !firstSide &&  secondSide &&  firstStrand && !secondStrand)    return true;
    if( !firstSide &&  secondSide && !firstStrand &&  secondStrand)    return true;
    
    return false;
}
		
bool ReassemblerPostProcess::checkContigSeed(OverlapSeedVecInfo input, int lengthThreshold)
{
    // filter tgs if it have different seed status
    int shortSeedCount = 0;
    int reverseSeed = 0;
    bool allSmall = true;
    
    for(SeedSequenceInfoVec::iterator iter = input.SeedInfoVec.begin(); iter!=input.SeedInfoVec.end() ; ++iter)
    {
        if((*iter).seedLength<=lengthThreshold)shortSeedCount++;
        else allSmall = false;
        
        if((*iter).contigSide!=input.SeedInfoVec[0].contigSide || (*iter).strand != input.SeedInfoVec[0].strand)
            reverseSeed++;
    }
    
    if(allSmall || reverseSeed >0)return true;

    return false;
}

bool ReassemblerPostProcess::checkSeedsOrientations(OverlapSeedVecInfo input)
{
    SeedSequenceInfoVec::iterator first = input.SeedInfoVec.begin();
    SeedSequenceInfoVec::iterator last  = first++;
    
    bool dir = ((*last).tgsLocation > (*first).tgsLocation);
    
    while(last!=input.SeedInfoVec.end()-1)
    {    
        if(dir != ((*last).tgsLocation > (*first).tgsLocation) )
            return true;            
                
        first++;
        last++;
    }
    return false;
}

bool ReassemblerPostProcess::checkSeedsDistance(OverlapSeedVecInfo input)
{
    SeedSequenceInfoVec::iterator first = input.SeedInfoVec.begin();
    SeedSequenceInfoVec::iterator last  = first++;
    
    while(last!=input.SeedInfoVec.end()-1)
    {    
        int contigDistanceBetweenTwoSeeds = std::abs((*last).contigStartPos - (*first).contigStartPos);
        int tgsDistanceBetweenTwoSeeds = std::abs((*last).tgsLocation - (*first).tgsLocation);

        if(std::abs(tgsDistanceBetweenTwoSeeds-contigDistanceBetweenTwoSeeds)>100 &&
           std::abs(tgsDistanceBetweenTwoSeeds-contigDistanceBetweenTwoSeeds)>contigDistanceBetweenTwoSeeds*0.2)
            return true;
        
        first++;
        last++;
    }
    return false;
}
		
bool ReassemblerPostProcess::checkRepeatByTGS(OverlapSeedVecInfo input)
{
    int slideWindow = 15;
    int repeatKmerCount_70 = 0;
    int totalKmerNumber = 0;
    
    for(SeedSequenceInfoVec::iterator iter = input.SeedInfoVec.begin(); iter!=input.SeedInfoVec.end() ; ++iter)
    {
        for( int i = 0 ; i< (int)((*iter).seedString.length()) - slideWindow + 1 ; i++ )
        {
            std::string kmer = (*iter).seedString.substr(i, slideWindow);

            size_t fwdKmerFreqs = BWTAlgorithms::countSequenceOccurrencesSingleStrand(kmer, m_params.indices);
            size_t rvcKmerFreqs = BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseComplement(kmer), m_params.indices);
            size_t kmerFreqs = fwdKmerFreqs+rvcKmerFreqs;
            
            totalKmerNumber++;
           
            // repeat kmer
            if((int)kmerFreqs>=kmerSizeFreqByTGS_15[seedCount_15/10*7])repeatKmerCount_70++;
        }
    }

    if( repeatKmerCount_70/totalKmerNumber>=0.8 )return true;

    return false;
}

bool ReassemblerPostProcess::checkRepeatByContig(OverlapSeedVecInfo input)
{
    SeedSequenceInfoVec::iterator test = input.SeedInfoVec.begin();
    
    int uniqueOrOverlap = 0;
    int unique = 0;
    int total  = 0;
    
    int start1 = 0;
    int start2 = 0;
    int last1  = 0;
    int last2  = 0;
    
    for(SeedSequenceInfoVec::iterator iter = input.SeedInfoVec.begin(); iter!=input.SeedInfoVec.end() ; ++iter)
    {
        std::string kmer = (*iter).seedString;

        size_t fwdKmerFreqs = BWTAlgorithms::countSequenceOccurrencesSingleStrand(kmer, m_params.contigIndices);
        size_t rvcKmerFreqs = BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseComplement(kmer), m_params.contigIndices);
        size_t kmerFreqs = fwdKmerFreqs+rvcKmerFreqs;
            
        //    contig1 ----------*-----------*--------*--------
        //                          --------*--------*-------------*-------contig2
        //        pacibo        -------------------------------------
        // contig index freq    1           2        2             1 
        
        if(kmerFreqs<=2)uniqueOrOverlap++;
        
        //    contig1 ----------*--*---*--           ----*----*--*-------contig2
        //        pacibo      -------------------------------------
        //    freq              1  1   1                 1    1  1 
        
        if(kmerFreqs==1)unique++;
        total++;
        
        last2 = last1;
        last1 = kmerFreqs;
        
        if(start1==0)
        {
            start1 = start2;
            start2 = kmerFreqs;
        }
    }

    if(m_params.isFirst  && uniqueOrOverlap==0 ) return true;
    
    if(m_params.isSecond && unique==0 ) return true;
    
    if( start1 > 1 || start2 > 1 ) return true;
    
    if( last1 > 1 || last2 > 1 ) return true;
    
    return false;
}

int Cmp(const void *lhs, const void *rhs) {
    return ((const int *)lhs)[1]-((const int *)rhs)[1];
}

bool ReassemblerPostProcess::checkPalindrome(std::string readSeq)
{
	SparseHashMap<std::string, KmerInfo, StringHasher> readKmerHashMap;
	size_t kmerSize = 11;
	
    for(size_t i = 0 ; i+kmerSize <= (size_t)readSeq.length() ; i++)
    {
        std::string forwardkmer = readSeq.substr(i, kmerSize);
		std::string reverseKmer = reverseComplement(forwardkmer);
		
		size_t fwdKmerFreqs = BWTAlgorithms::countSequenceOccurrencesSingleStrand(forwardkmer, m_params.indices);
        size_t rvcKmerFreqs = BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseKmer, m_params.indices);
        size_t kmerFreqs = fwdKmerFreqs+rvcKmerFreqs;
		
		if( ((int)kmerFreqs > kmerSizeFreqByTGS_11[seedCount_11/10*7])) continue;
		
		SparseHashMap<std::string, KmerInfo, StringHasher>::iterator findKmerIter = readKmerHashMap.find(forwardkmer);
		
		if( findKmerIter == readKmerHashMap.end() ) 
		{
			KmerInfo tmpKmer;
			tmpKmer.position = (int)i;
			tmpKmer.secondPosition = 0;
			tmpKmer.isRepeat = false;
			readKmerHashMap.insert(std::make_pair(forwardkmer,tmpKmer));
		}
		else
		{
			findKmerIter->second.isRepeat = true;
		}
		
		findKmerIter = readKmerHashMap.find(reverseKmer);
		
		if( findKmerIter!=readKmerHashMap.end() )
		{
			if(!findKmerIter->second.isRepeat && findKmerIter->second.secondPosition == 0 )
				findKmerIter->second.secondPosition = (int)i;
			else if(!findKmerIter->second.isRepeat)
			{
				findKmerIter->second.isRepeat = true;
			}
		}
	}

	int tmp = 0;
	int count = 0;

	int positionArray[readSeq.length()][2];
	
	for(SparseHashMap<std::string, KmerInfo, StringHasher>::iterator kmerIter = readKmerHashMap.begin() ; kmerIter != readKmerHashMap.end() ; ++kmerIter )
	{
		if(!kmerIter->second.isRepeat && kmerIter->second.secondPosition != 0)
		{
			positionArray[tmp][0] = kmerIter->second.position;
			positionArray[tmp][1] = kmerIter->second.secondPosition;
			tmp++;
		}
	}
	
	qsort(positionArray,tmp,sizeof(positionArray[0]),Cmp);
	
	for(int i = 1; i<tmp ; i++)
	{
		if(std::abs(positionArray[i][1] - positionArray[i-1][1]) < 40 && positionArray[i][0] < positionArray[i-1][0]) 
		{	
			count++;
		}
		else 
		{
			count = 0;
		}
		if(count>=6)return true;
	}
	
	return false;
}

bool ReassemblerPostProcess::checkOverlapLength(OverlapSeedVecInfo input)
{
    SeedSequenceInfoVec::iterator first = input.SeedInfoVec.begin();
    SeedSequenceInfoVec::iterator last  = input.SeedInfoVec.end();
	--last;
	
    if(std::abs((int)(*first).contigStartPos - (int)(*last).contigStartPos)<400) return true;
    if(std::abs((int)(*first).tgsLocation - (int)(*last).tgsLocation)<400) return true;
    
    return false;
}

bool ReassemblerPostProcess::checkNonOverlapPartialLength(OverlapSeedVecInfo input, TGSScanAndOverlapPosition &position)
{
    bool contigSide  = input.SeedInfoVec[0].contigSide;
    int tgsIndex  = input.SeedInfoVec[0].tgsIndex;
    int contigLength = input.SeedInfoVec[0].contigLength;
    std::string originTGS = backtrackTGS( m_params.indices, tgsIndex );
	
	TGSScanAndOverlapPosition reScanPosition = reScanContigOnTGSRange( input.SeedInfoVec[0].contigPartialSeq, originTGS );
	
    int thresholdLength = 400; 
    
    SeedSequenceInfoVec::iterator first = input.SeedInfoVec.begin();
    SeedSequenceInfoVec::iterator last  = input.SeedInfoVec.end();
    last--;

    // contig head side
    if(contigSide)
    {
        // tgs dir ---->
        if( (*first).tgsLocation < (*last).tgsLocation )
        {
            if( (int)(*first).tgsLocation - (int)(*first).contigStartPos < thresholdLength ) return true;
            position.scanTGSStart = 0;
			position.scanTGSEnd   = (*first).tgsLocation - (*first).contigStartPos;
			position.tgsFragmentStart = 0;
			position.tgsFragmentEnd   = (*first).tgsLocation;
			//position.cutContigToEndPointLength = (*first).contigStartPos;
			
			if(reScanPosition.kmerCount>10)
			{
				position.scanTGSEnd     = reScanPosition.scanTGSStart;
				position.tgsFragmentEnd = reScanPosition.scanTGSStart;
			}
        }
        // tgs dir <----
        else
        {
            if( (int)originTGS.length() - ( (*first).tgsLocation + (*first).contigStartPos ) < thresholdLength ) return true;
            position.scanTGSStart = (*first).tgsLocation + (*first).contigStartPos;
			position.scanTGSEnd   = originTGS.length();
			position.tgsFragmentStart = (*first).tgsLocation;
			position.tgsFragmentEnd   = originTGS.length();
			//position.cutContigToEndPointLength = (*first).contigStartPos;
			
			if(reScanPosition.kmerCount>10)
			{
				position.scanTGSEnd     = reScanPosition.scanTGSEnd;
				position.tgsFragmentEnd = reScanPosition.scanTGSEnd;
			}
        }
    }
    // contig tail side
    else
    {
        int toEndLength = contigLength - (*last).contigStartPos;
		// tgs dir ---->
		if( (*first).tgsLocation < (*last).tgsLocation )
        {
			if( ( (int)originTGS.length() - ( (*last).tgsLocation + toEndLength ) ) < thresholdLength ) return true;
            position.scanTGSStart = (*last).tgsLocation + toEndLength;
			position.scanTGSEnd   = originTGS.length();
			position.tgsFragmentStart = (*last).tgsLocation;
			position.tgsFragmentEnd   = originTGS.length();
			//position.cutContigToEndPointLength = toEndLength;
			
			if(reScanPosition.kmerCount>10)
			{
				position.scanTGSEnd     = reScanPosition.scanTGSEnd;
				position.tgsFragmentEnd = reScanPosition.scanTGSEnd;
			}
        }
		// tgs dir <----
        else
        {
            if( (*last).tgsLocation - toEndLength < thresholdLength ) return true;
			position.scanTGSStart = 0;
			position.scanTGSEnd   = (*first).tgsLocation - toEndLength;
			position.tgsFragmentStart = 0;
			position.tgsFragmentEnd   = (*last).tgsLocation;
			//position.cutContigToEndPointLength = toEndLength;
			
			if(reScanPosition.kmerCount>10)
			{
				position.scanTGSEnd     = reScanPosition.scanTGSStart;
				position.tgsFragmentEnd = reScanPosition.scanTGSStart;
			}
        }
    }
    
	//std::cout<<position.scanTGSStart << "\t" << position.scanTGSEnd <<"   range\n";
    return false;
}

		///**************** tool function *****************///
		
void ReassemblerPostProcess::insertContigRelationHash(std::string firstContig, bool firstSide, bool firstStrand, std::string secondContig, bool secondSide, bool secondStrand, int firstContigLength, int secondContigLength, std::vector<OverlapSeedVecInfo> contigInfoVec)
{
    OverlapRelation contigPair;
    std::vector<OverlapRelation> insertVec;
    std::string insertKey = firstContig+(firstSide ? "_Head":"_Tail");

    std::string originTGS = backtrackTGS( m_params.indices, (int)contigInfoVec[0].SeedInfoVec[0].tgsIndex );
    contigTGSPosition f   = OnTGSRange( contigInfoVec[0].SeedInfoVec , (int64_t)originTGS.length() );
    contigTGSPosition s   = OnTGSRange( contigInfoVec[1].SeedInfoVec , (int64_t)originTGS.length() );

    contigPair.firstContig   = firstContig;
    contigPair.secondContig  = secondContig;
    contigPair.firstSide     = firstSide;
    contigPair.secondSide    = secondSide;
    contigPair.firstStrand   = firstStrand;
    contigPair.secondStrand  = secondStrand;
    
    contigPair.tgsLength  = originTGS.length();
    contigPair.tgsIndex   = contigInfoVec[0].SeedInfoVec[0].tgsIndex;
    contigPair.connect       = (f.maxTGSPos < s.minTGSPos || s.maxTGSPos < f.minTGSPos);
    
    contigPair.firstContigLength  = firstContigLength;
    contigPair.secondContigLength = secondContigLength;

    if( (f.firstTGSStartPos < f.lastTGSStartPos && s.firstTGSStartPos < s.lastTGSStartPos) || 
        (f.firstTGSStartPos > f.lastTGSStartPos && s.firstTGSStartPos > s.lastTGSStartPos))
        contigPair.isSameStrand = true;
    else
        contigPair.isSameStrand = false;

    //non overlap
    if(contigPair.connect)
    {
        if(f.minTGSPos > s.maxTGSPos)
            contigPair.tgsFragment = originTGS.substr(s.maxTGSPos,f.minTGSPos - s.maxTGSPos);
        else
            contigPair.tgsFragment = originTGS.substr(f.maxTGSPos,s.minTGSPos - f.maxTGSPos); 
        
        contigPair.overlapLength  = 0;
        contigPair.firContigStart = 0;
        contigPair.firContigEnd   = 0;
        contigPair.secContigStart = 0;
        contigPair.secContigEnd   = 0;
    }
    // overlap
    else
    {  
        contigPair.overlapLength = std::abs(f.contigEndPointOnTGSPos-s.contigEndPointOnTGSPos) + f.offsetOverlapLength + s.offsetOverlapLength;
        
        std::pair<int,int> firOverlap = overlapPosition(firstSide,contigPair.firstContigLength,contigPair.overlapLength);
        std::pair<int,int> secOverlap = overlapPosition(secondSide,contigPair.secondContigLength,contigPair.overlapLength);
        
        contigPair.firContigStart = firOverlap.first;
        contigPair.firContigEnd   = firOverlap.second;
        contigPair.secContigStart = secOverlap.first;
        contigPair.secContigEnd   = secOverlap.second;
    }

    // prevent different tgs build same relation edge  
    ContigRelationHashMap::iterator contigIterator = connectContigHashMap.find(insertKey);
    //this contigSide already exist
    if( contigIterator != connectContigHashMap.end() ) 
    {
        // push this contigSide to hash
        contigIterator->second.push_back(contigPair);
    }
    //this contigSide not exist
    else
    {
        insertVec.push_back(contigPair);
        connectContigHashMap.insert(std::make_pair(insertKey,insertVec));
    }
}

void ReassemblerPostProcess::insertSingleRelationHash(OverlapSeedVecInfo input, TGSScanAndOverlapPosition scanRange)
{
    OverlapRelation contigPair;
    std::vector<OverlapRelation> insertVec;
    
    std::string originTGS = backtrackTGS( m_params.indices, (int)input.SeedInfoVec[0].tgsIndex );
    contigTGSPosition f   = OnTGSRange( input.SeedInfoVec, (int)originTGS.length() );
	
    contigPair.firstContig   = input.contig;
    contigPair.firstSide     = input.SeedInfoVec[0].contigSide;
    contigPair.firstStrand   = input.SeedInfoVec[0].strand;
    contigPair.tgsLength  = originTGS.length();
    contigPair.tgsIndex   = input.SeedInfoVec[0].tgsIndex;
    contigPair.firstContigLength  = input.SeedInfoVec[0].contigLength;
    //contigPair.scanStart = scanRange.first;
    //contigPair.scanEnd   = scanRange.second;
	contigPair.singleOverlap = scanRange;
	
	contigPair.tgsDirection = f.firstTGSStartPos < f.lastTGSStartPos;
	
    TGSConnectContigHashMap::iterator tgsIterator = singleContigHashMap.find(contigPair.tgsIndex);
    
    if( tgsIterator != singleContigHashMap.end() ) 
    {
        tgsIterator->second.push_back(contigPair);
    }
    else
    {
        insertVec.push_back(contigPair);
        singleContigHashMap.insert(std::make_pair(contigPair.tgsIndex,insertVec));
    }
}

contigTGSPosition ReassemblerPostProcess::OnTGSRange(SeedSequenceInfoVec inputSeedInfoVec, int64_t tgsLength)
{
    contigTGSPosition result;
    
    int64_t firstContigStartPos = inputSeedInfoVec[0].contigStartPos;
    int64_t lastContigStartPos  = inputSeedInfoVec[0].contigStartPos;
    int64_t firstTGSStartPos = inputSeedInfoVec[0].tgsLocation;
    int64_t lastTGSStartPos  = inputSeedInfoVec[0].tgsLocation;
    int64_t contigEndPointOnTGSPos;

    result.offsetOverlapLength = 0;
    
	std::string originTGS = backtrackTGS( m_params.indices, (int)inputSeedInfoVec[0].tgsIndex );
	TGSScanAndOverlapPosition reScanPosition = reScanContigOnTGSRange( inputSeedInfoVec[0].contigPartialSeq, originTGS );
	
	
    for(SeedSequenceInfoVec::iterator iter =inputSeedInfoVec.begin(); iter!=inputSeedInfoVec.end() ; ++iter)
    {
        if( (*iter).contigStartPos < firstContigStartPos )
        {
            firstContigStartPos = (*iter).contigStartPos;
            firstTGSStartPos = (*iter).tgsLocation;
        }
        if( (*iter).contigStartPos + (*iter).seedLength > lastContigStartPos )
        {
            lastContigStartPos = (*iter).contigStartPos + (*iter).seedLength;
            lastTGSStartPos = (*iter).tgsLocation + (*iter).seedLength;
        }
    }
    
    //head
    if(inputSeedInfoVec[0].contigSide)
    {
        //SENSE
        if( firstTGSStartPos < lastTGSStartPos )
        {
            contigEndPointOnTGSPos = std::max( firstTGSStartPos - firstContigStartPos, (int64_t)0 );
            
			if(reScanPosition.kmerCount>10)
				contigEndPointOnTGSPos  = reScanPosition.scanTGSStart;
			
            //if( firstContigStartPos > firstTGSStartPos )
            //    result.offsetOverlapLength = std::abs(firstTGSStartPos - firstContigStartPos);
        }
        //ANTISENSE
        else 
        {    
            contigEndPointOnTGSPos = std::min( firstTGSStartPos + firstContigStartPos, tgsLength );
            
			if(reScanPosition.kmerCount>=10)
				contigEndPointOnTGSPos  = reScanPosition.scanTGSEnd;
			
            //if( firstTGSStartPos - firstContigStartPos > tgsLength)
            //    result.offsetOverlapLength = std::abs(firstTGSStartPos + firstContigStartPos - tgsLength );
        }
    }
    //tail
    else
    {
        int64_t contigLength = inputSeedInfoVec[0].contigLength;
        //SENSE
        if( firstTGSStartPos < lastTGSStartPos ) 
        {
            contigEndPointOnTGSPos = std::min( lastTGSStartPos + (contigLength - lastContigStartPos), tgsLength );
            
			if(reScanPosition.kmerCount>=10)
				contigEndPointOnTGSPos  = reScanPosition.scanTGSEnd;
			
            //if( lastTGSStartPos + (contigLength - lastContigStartPos) > tgsLength)
            //    result.offsetOverlapLength = std::abs( lastTGSStartPos + (contigLength - lastContigStartPos) - tgsLength );
        }
        //ANTISENSE
        else     
        {
            contigEndPointOnTGSPos = std::max( lastTGSStartPos - (contigLength - lastContigStartPos) , (int64_t)0 );
            
			if(reScanPosition.kmerCount>10)
				contigEndPointOnTGSPos  = reScanPosition.scanTGSStart;
			
            //if( lastTGSStartPos - (contigLength - lastContigStartPos) < 0)
            //    result.offsetOverlapLength = std::abs(lastTGSStartPos - (contigLength - lastContigStartPos));
        }
    }

    result.minTGSPos = std::min(lastTGSStartPos,firstTGSStartPos);
    result.minTGSPos = std::min(result.minTGSPos,contigEndPointOnTGSPos);
    
    result.maxTGSPos = std::max(lastTGSStartPos,firstTGSStartPos);
    result.maxTGSPos = std::max(result.maxTGSPos,contigEndPointOnTGSPos);
    
    result.firstContigStartPos = firstContigStartPos;
    result.lastContigStartPos  = lastContigStartPos;
    result.firstTGSStartPos = firstTGSStartPos;
    result.lastTGSStartPos  = lastTGSStartPos;
    result.contigEndPointOnTGSPos = contigEndPointOnTGSPos;
	
    return result;
}

TGSScanAndOverlapPosition ReassemblerPostProcess::reScanContigOnTGSRange(std::string partialContig, std::string originTGS)
{
	TGSScanAndOverlapPosition result;
	
	SparseHashMap<std::string, KmerInfo, StringHasher> contigPartialKmerHashMap;
	std::vector<TGSScanAndOverlapPosition> rePositionVec;
	
	int minPosition = 0;
	
	size_t kmerSize = 11;
	
	result.kmerCount = 0;
	result.scanTGSStart = 0;
	result.scanTGSEnd   = 0;

    for(size_t i = 0 ; i+kmerSize <= (size_t)partialContig.length() ; i++)
    {
		std::string forwardkmer = partialContig.substr(i, kmerSize);

		SparseHashMap<std::string, KmerInfo, StringHasher>::iterator forwardKmerIter = contigPartialKmerHashMap.find(forwardkmer);
		
		if(forwardKmerIter!=contigPartialKmerHashMap.end())
		{
			forwardKmerIter->second.isRepeat = true;
		}
		else
		{
			KmerInfo tmp;
			tmp.isRepeat = false;
			tmp.position = -1;
			contigPartialKmerHashMap.insert(std::make_pair(forwardkmer,tmp));
		}
	}
	
	for(size_t i = 0 ; i+kmerSize <= originTGS.length() ; i++)
    {
		std::string forwardkmer = originTGS.substr(i, kmerSize);
		std::string reversekmer = reverseComplement(forwardkmer);
		
		SparseHashMap<std::string, KmerInfo, StringHasher>::iterator forwardKmerIter = contigPartialKmerHashMap.find(forwardkmer);
		if(forwardKmerIter!=contigPartialKmerHashMap.end())
		{
			if(!forwardKmerIter->second.isRepeat && forwardKmerIter->second.position == -1)
			{
				forwardKmerIter->second.position = i;
			}
			else if(forwardKmerIter->second.position!= -1)
			{
				forwardKmerIter->second.isRepeat = true;
				forwardKmerIter->second.position = -1;
			}
		}
		
		SparseHashMap<std::string, KmerInfo, StringHasher>::iterator reverseKmerIter = contigPartialKmerHashMap.find(reversekmer);
		if(reverseKmerIter!=contigPartialKmerHashMap.end())
		{
			if(!reverseKmerIter->second.isRepeat && reverseKmerIter->second.position == -1)
			{
				reverseKmerIter->second.position = i;
			}
			else if(reverseKmerIter->second.position!= -1)
			{
				reverseKmerIter->second.isRepeat = true;
				reverseKmerIter->second.position = -1;
			}
		}
	}
	
	int tmp = 0;
	int contigOnTGSPosition[partialContig.length()];
	
	for(SparseHashMap<std::string, KmerInfo, StringHasher>::iterator kmerIter = contigPartialKmerHashMap.begin() ; kmerIter!=contigPartialKmerHashMap.end() ; ++kmerIter)
	{
		if(!kmerIter->second.isRepeat && kmerIter->second.position != -1)
		{
			contigOnTGSPosition[tmp] = kmerIter->second.position;
			tmp++;
		}
	}
	
	std::sort(contigOnTGSPosition,contigOnTGSPosition+tmp);
	
	for(int i =1;i<tmp;i++)
	{
		if( (std::abs(contigOnTGSPosition[i] - contigOnTGSPosition[i-1])) < 40 )
		{
			if(minPosition==0)
			{
				minPosition = contigOnTGSPosition[i-1];
				TGSScanAndOverlapPosition tmpPosition;
				tmpPosition.scanTGSStart = contigOnTGSPosition[i-1];
				tmpPosition.scanTGSEnd   = contigOnTGSPosition[i];
				tmpPosition.kmerCount = 1;
				rePositionVec.push_back(tmpPosition);
			}
			else
			{
				std::vector<TGSScanAndOverlapPosition>::iterator tmpIter = rePositionVec.end();
				--tmpIter;
				(*tmpIter).scanTGSEnd = contigOnTGSPosition[i];
				(*tmpIter).kmerCount++;
			}
		}
		else minPosition = 0;	
	}
	
	for(std::vector<TGSScanAndOverlapPosition>::iterator tmpIter = rePositionVec.begin() ; tmpIter != rePositionVec.end() ; ++tmpIter )
	{
		if( (*tmpIter).kmerCount > result.kmerCount )
		{
			result.kmerCount = (*tmpIter).kmerCount;
			result.scanTGSStart = (*tmpIter).scanTGSStart;
			result.scanTGSEnd = (*tmpIter).scanTGSEnd;
		}
	}

	return result;
}

std::pair<int,int> ReassemblerPostProcess::overlapPosition(bool side, int contigLength, int overlapLength)
{
    if(side) return std::make_pair(0,(overlapLength-1));
    else     return std::make_pair((contigLength-overlapLength),(contigLength-1));
}

std::string ReassemblerPostProcess::mostConnectContig(std::vector<OverlapRelation> inputRelationVec)
{
    std::string result;
    SparseHashMap<std::string, int, StringHasher> numberOfConnect;
    int first = 0;
    int second = 0;
    
    for(std::vector<OverlapRelation>::iterator iter = inputRelationVec.begin(); iter != inputRelationVec.end(); ++iter )
    {
        SparseHashMap<std::string, int, StringHasher>::iterator contigIter  = numberOfConnect.find((*iter).secondContig);
        
        if( contigIter != numberOfConnect.end() ) contigIter->second++;
        else numberOfConnect.insert(std::make_pair((*iter).secondContig,1));
    }
    
    for(SparseHashMap<std::string, int, StringHasher>::iterator numIter = numberOfConnect.begin(); numIter != numberOfConnect.end(); ++numIter )
    {
        if( numIter->second > first )
        {
            second = first;
            first = numIter->second;
            result = numIter->first;
        }
        else if( numIter->second > second )
        {
            second = numIter->second;
        }
    }
    if(first==second) return "";
    return result;
}
  
int ReassemblerPostProcess::mostConnectTGS(OverlapSeedVecInfo input)
{
	SparseHashMap<int, int, Int64Hasher> tgsAndCount;
	int mostconnectTGSIndex = 0;
	int maxSeedCount = 0;
	
	for( SeedSequenceInfoVec::iterator tmpSeedIter = input.SeedInfoVec.begin() ; tmpSeedIter != input.SeedInfoVec.end() ; ++tmpSeedIter )
	{	
		SparseHashMap<int, int, Int64Hasher>::iterator findTGSIndex = tgsAndCount.find((*tmpSeedIter).tgsIndex);
		
		if( findTGSIndex!= tgsAndCount.end() )
			findTGSIndex->second++;
		else
			tgsAndCount.insert(std::make_pair((*tmpSeedIter).tgsIndex,1));
	}
	
	for( SparseHashMap<int, int, Int64Hasher>::iterator tmp = tgsAndCount.begin() ; tmp != tgsAndCount.end() ; tmp++ )
	{
		if( tmp->second > maxSeedCount )
		{
			mostconnectTGSIndex = tmp->first;
			maxSeedCount  = tmp->second;
		}
	}
	return mostconnectTGSIndex;
}

int ReassemblerPostProcess::connectVote(size_t tgsIdx, std::vector<OverlapRelation> inputRelationVec)
{
    // 0:same number of overlap and connect 
	// 1:overlap case 
	// 2:nonOverlap case
	
	size_t overlap = 0;
    size_t connect = 0;
    size_t tgsStatus = 0;
    
    for(std::vector<OverlapRelation>::iterator iter = inputRelationVec.begin(); iter != inputRelationVec.end(); ++iter )
    {
        if((*iter).tgsIndex == (int)tgsIdx)
        {
            if((*iter).connect) tgsStatus = 2;
            else tgsStatus = 1;
        }
        
        if((*iter).connect) connect++;
        else overlap++;
    }

    if( tgsStatus==1 && connect < overlap ) return 1;
    if( tgsStatus==2 && connect > overlap ) return 2;
    
    return 0;

}

bool ReassemblerPostProcess::checkConnect(std::string firstID, bool firstSide,std::string secondID, bool secondSide)
{
    std::string firstInsert  = firstID  + (firstSide ? "_Head":"_Tail");
    std::string secondInsert = secondID + (secondSide ? "_Head":"_Tail");
    
    SparseHashMap<std::string,bool,StringHasher>::iterator firstIter  = connectStatus.find(firstInsert);
    SparseHashMap<std::string,bool,StringHasher>::iterator secondIter = connectStatus.find(secondInsert);
    
    //one or two contigSide already connect
    if( firstIter != connectStatus.end() || secondIter != connectStatus.end()) 
    {
        return true;
    }
    else
    {
        connectStatus.insert(std::make_pair(firstInsert,true));
        connectStatus.insert(std::make_pair(secondInsert,true));
        return false;
    }
}

std::vector<ReassemblerSeedFeature> ReassemblerPostProcess::filterRepeatSeed(std::vector<ReassemblerSeedFeature> inputSeedVec)
{
    std::vector<ReassemblerSeedFeature> seedVec;
	
	int repeatCount = 0;
	int reCount = 0;
	
    for(std::vector<ReassemblerSeedFeature>::iterator iter = inputSeedVec.begin(); iter != inputSeedVec.end(); ++iter)
    {
		if(repeatCount>0){repeatCount--;continue;}
		
		size_t fwdSeedFreqByTGS = BWTAlgorithms::countSequenceOccurrencesSingleStrand((*iter).seedStr, m_params.indices);
        size_t rvcSeedFreqByTGS = BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseComplement((*iter).seedStr), m_params.indices);
        size_t SeedFreqByTGS = fwdSeedFreqByTGS + rvcSeedFreqByTGS;
		
		size_t fwdSeedFreqByContig = BWTAlgorithms::countSequenceOccurrencesSingleStrand((*iter).seedStr, m_params.contigIndices);
        size_t rvcSeedFreqByContig = BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseComplement((*iter).seedStr), m_params.contigIndices);
        size_t SeedFreqByContig = fwdSeedFreqByContig + rvcSeedFreqByContig;
		
		//if((int)SeedFreqByTGS ==2)std::cout<<". ";
		//else std::cout<<SeedFreqByTGS<<" ";
		
		// filter repeat using tgs index and contig index
        if( ((int)SeedFreqByTGS > 4 /*kmerSizeFreqByTGS_15[seedCount_15/10*7]*/) || ((int)SeedFreqByContig > 2) ) 
		{
			int localCount = 1;
			repeatCount = 1;
			
			while(!seedVec.empty() && localCount>0 && reCount>0)
			{
				localCount--;
				reCount--;
				seedVec.pop_back();	
			}
			
			reCount = 0;
			continue;
		}
		
		reCount++;
        seedVec.push_back((*iter));    
		
    }
    return seedVec;
}

void ReassemblerPostProcess::filterErrorSeeds(OverlapSeedVecInfo input, OverlapSeedVecInfo &result, int mostTGS)
{
	for( SeedSequenceInfoVec::iterator tmpSeedIter = input.SeedInfoVec.begin() ; tmpSeedIter != input.SeedInfoVec.end() ; ++tmpSeedIter )
	{
		if( (*tmpSeedIter).tgsIndex != mostTGS ) continue;
		result.SeedInfoVec.push_back((*tmpSeedIter));
	}
}

void ReassemblerPostProcess::buildKmerHashAndTGSKmerVec( std::string tgsConnectContig, bool connectContigSide, int connectContigLength, int tgsIndex, std::vector<ReassemblerSeedFeature> inputSeedVec)
{
    KmerInfoVec currentTGSKmerVec;
    
    for( std::vector<ReassemblerSeedFeature>::iterator iter = inputSeedVec.begin(); iter != inputSeedVec.end(); ++iter )
    {    
        KmerInfo currentKmer;
        
        currentKmer.tgsConnectContigSide   = connectContigSide;
        currentKmer.tgsConnectContigID     = tgsConnectContig;
		currentKmer.tgsConnectContigLength = connectContigLength;
        currentKmer.tgsIndex = tgsIndex;
        currentKmer.kmerStr     = (*iter).seedStr;
        currentKmer.position    = (*iter).seedStartPos;
        
        currentTGSKmerVec.push_back(currentKmer);
        
        KmerHashMap::iterator findforwardKmerIter = collectAllKmerHashMap.find((*iter).seedStr);
        
        if( findforwardKmerIter!=collectAllKmerHashMap.end() )
        {
            currentKmer.kmerStrand = true;
            findforwardKmerIter->second.push_back(currentKmer);
        }
        else
        {
            std::vector<KmerInfo> tmpVec;
            currentKmer.kmerStrand = true;
            tmpVec.push_back(currentKmer);
            collectAllKmerHashMap.insert(std::make_pair((*iter).seedStr,tmpVec));
        }
        
        KmerHashMap::iterator findreverseKmerIter = collectAllKmerHashMap.find(reverseComplement((*iter).seedStr));
        
        currentKmer.kmerStr = reverseComplement((*iter).seedStr);
        
        if( findreverseKmerIter!=collectAllKmerHashMap.end() )
        {
            currentKmer.kmerStrand = false;
            findreverseKmerIter->second.push_back(currentKmer);
        }
        else
        {
            std::vector<KmerInfo> tmpVec;
            currentKmer.kmerStrand = false;
            tmpVec.push_back(currentKmer);
            collectAllKmerHashMap.insert(std::make_pair(reverseComplement((*iter).seedStr),tmpVec));
        }
    }
    AllTGSKmerVec.push_back(std::make_pair(tgsIndex,currentTGSKmerVec));
}

bool ReassemblerPostProcess::checkSameStrand(int startTGSIndex, int targetTGSIndex, SeedSequenceInfoVec brigdeSeedVec)
{
	bool strand;
	
	TGSConnectContigHashMap::iterator startTGSIter  = singleContigHashMap.find(startTGSIndex);
	TGSConnectContigHashMap::iterator targetTGSIter = singleContigHashMap.find(targetTGSIndex);
	
	std::vector<OverlapRelation>::iterator startFirst = startTGSIter->second.begin();
	std::vector<OverlapRelation>::iterator targetFirst = targetTGSIter->second.begin();


	if( (*startFirst).tgsDirection == (*targetFirst).tgsDirection ) strand = true;
	else strand = false;
		
	if(!brigdeSeedVec[0].strand) strand = !strand;

	return strand;

}

int ReassemblerPostProcess::overlapLength(OverlapSeedVecInfo input)
{
    SeedSequenceInfoVec::iterator first = input.SeedInfoVec.begin();
    SeedSequenceInfoVec::iterator last  = input.SeedInfoVec.end();
	--last;

    return std::abs((int)(*first).contigStartPos - (int)(*last).contigStartPos);
}

std::string ReassemblerPostProcess::concordantNonOverlapFragment(OverlapSeedVecInfo input)
{
	TGSConnectContigHashMap::iterator targetTGSIter = singleContigHashMap.find(input.SeedInfoVec[0].tgsIndex);
	std::string secondTGSRead = backtrackTGS( m_params.indices, input.SeedInfoVec[0].tgsIndex );	
	std::string fragment;
	
	SeedSequenceInfoVec::iterator first = input.SeedInfoVec.begin();
	SeedSequenceInfoVec::iterator last = input.SeedInfoVec.end();
	last--;
	
	if( targetTGSIter->second[0].singleOverlap.scanTGSStart == 0 )
	{
		//int fragmentLength = targetTGSIter->second[0].singleOverlap.scanTGSEnd - overlapLength(input);
		//fragment = secondTGSRead.substr( overlapLength(input), fragmentLength );
		int start = (*first).tgsLocation < (*last).tgsLocation ? (*first).tgsLocation : (*last).tgsLocation;
		int end   = targetTGSIter->second[0].singleOverlap.scanTGSEnd;
		fragment = secondTGSRead.substr( start, end - start );
		
		//std::cout<< "start2\n";
		//std::cout<< (*first).tgsLocation << "\t" << (*last).tgsLocation <<"\n";
		//std::cout<< start << "\t" << end << "\t idx " << input.SeedInfoVec[0].tgsIndex << "\n";
	}
	else
	{
		//int fragmentLength = targetTGSIter->second[0].singleOverlap.scanTGSEnd - targetTGSIter->second[0].singleOverlap.scanTGSStart - overlapLength(input);
		//fragment = secondTGSRead.substr( targetTGSIter->second[0].singleOverlap.scanTGSStart, fragmentLength );
		
		int start = targetTGSIter->second[0].singleOverlap.scanTGSStart;
		int end   = (*first).tgsLocation > (*last).tgsLocation ? (*first).tgsLocation : (*last).tgsLocation;
		fragment = secondTGSRead.substr( start, end - start );
		
		//std::cout<< "end2\n";
		//std::cout<< (*first).tgsLocation << "\t" << (*last).tgsLocation <<"\n";
		//std::cout<< start << "\t" << end << "\t idx " << input.SeedInfoVec[0].tgsIndex << "\n";
	}
	
	if( !input.SeedInfoVec[0].strand ) return reverseComplement(fragment);
	return fragment;
}

void ReassemblerPostProcess::countTotalSeedFreq(std::string seed)
{
    int slideWindow = 15;

    for( int i = 0 ; i < (int)seed.length() - slideWindow + 1 ; i++ )
    {
        std::string kmer = seed.substr(i, slideWindow);

        size_t fwdKmerFreqsUsingTGS = BWTAlgorithms::countSequenceOccurrencesSingleStrand(kmer, m_params.indices);
        size_t rvcKmerFreqsUsingTGS = BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseComplement(kmer), m_params.indices);
        size_t kmerFreqsUsingTGS = fwdKmerFreqsUsingTGS + rvcKmerFreqsUsingTGS;
        
        size_t fwdKmerFreqsUsingContig = BWTAlgorithms::countSequenceOccurrencesSingleStrand(kmer, m_params.contigIndices);
        size_t rvcKmerFreqsUsingContig = BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseComplement(kmer), m_params.contigIndices);
        size_t kmerFreqsUsingContig = fwdKmerFreqsUsingContig + rvcKmerFreqsUsingContig;
        
        if(seedCount_15>=100000)break;
       
        kmerSizeFreqByTGS_15[seedCount_15] = kmerFreqsUsingTGS;
        kmerFreqByContig_15[seedCount_15] = kmerFreqsUsingContig;
        seedCount_15++;
    }
	
	slideWindow = 11;
	
	for( int i = 0 ; i < (int)seed.length() - slideWindow + 1 ; i++ )
    {
        std::string kmer = seed.substr(i, slideWindow);

        size_t fwdKmerFreqsUsingTGS = BWTAlgorithms::countSequenceOccurrencesSingleStrand(kmer, m_params.indices);
        size_t rvcKmerFreqsUsingTGS = BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseComplement(kmer), m_params.indices);
        size_t kmerFreqsUsingTGS = fwdKmerFreqsUsingTGS + rvcKmerFreqsUsingTGS;
        
        if(seedCount_11>=100000)break;
       
        kmerSizeFreqByTGS_11[seedCount_11] = kmerFreqsUsingTGS;
        seedCount_11++;
    }
}

		///********** Reassembler Seed Feature **********///

ReassemblerSeedFeature::ReassemblerSeedFeature(size_t startPos, std::string str, bool repeat, size_t kmerSize, size_t repeatCutoff)
:seedStartPos(startPos), seedStr(str), isRepeat(repeat), freqUpperBound(repeatCutoff), freqLowerBound(15)
{
    seedEndPos = seedStartPos + seedStr.length() -1;
    seedLength = seedStr.length();

    if(repeat) startBestKmerSize = endBestKmerSize = seedLength>31?31: seedLength;				
    else startBestKmerSize = endBestKmerSize = kmerSize;
}





