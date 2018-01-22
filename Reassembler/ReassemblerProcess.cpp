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

bool showReassemblerProcess = false;

        ///************************************************///
        ///*********** Read Reassembler Process ************///
        ///************************************************///

ReassemblerProcess::ReassemblerProcess(const ReassemblerParameters params) 
: ReadReassemblerBasicElements(params)
{
}
ReassemblerProcess::~ReassemblerProcess()
{
}

ReassemblerResult ReassemblerProcess::PBReassembler(const SequenceWorkItem& workItem)
{    
    // using growing kmer find out the overlap between contig and long read
	
	ReassemblerResult result;
    SeedSequenceInfoVec seedVec;
    
    std::vector<ReassemblerSeedFeature> headSeedVec,tailSeedVec;
    // get origin contig sequence
    std::string contigSeq = workItem.read.seq.toString();
    
    //assume contig length bigger than 2n
    //we will use head and tail n length to find seeds
    // skip short read/contig
    if( (int)contigSeq.length() < m_params.searchRange*2 ) return result;

    // cut head and tail n length
    std::string headContigSeq = contigSeq.substr(0,m_params.searchRange);
    std::string tailContigSeq = contigSeq.substr(contigSeq.length()-m_params.searchRange,m_params.searchRange);

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
        SeedByReadHashMap::iterator pbIterator = result.seedHash.find((*iter).readIndex);
        
        //std::cout << (*iter).seedLength << " ";
        
        // filter short seed
        if( (*iter).seedLength<=15 )continue;

        //this read exist
        if( pbIterator!=result.seedHash.end() )
        {
            pbIterator->second[0].SeedInfoVec.push_back((*iter));
        }
        // this read not exist
        else
        {
            OverlapSeedVecInfo pushContig;
            std::vector<OverlapSeedVecInfo> contigVec;
            pushContig.contig = (*iter).contigID;
            pushContig.SeedInfoVec.push_back((*iter));
            contigVec.push_back(pushContig);
            result.seedHash.insert(std::make_pair((*iter).readIndex,contigVec));
        }
    }
    return result;
}

void ReassemblerProcess::findSeedOnLongRead(const BWTIndexSet indices, std::string contigID, std::vector<ReassemblerSeedFeature> inputSeedVec, std::vector<SeedSequenceInfo> &result, int contigLength, std::string contigPartialSeq, bool contigSide)
{
    for( std::vector<ReassemblerSeedFeature>::iterator iter = inputSeedVec.begin(); iter != inputSeedVec.end(); ++iter )
    {    
        SeedSequenceInfo currentSeedInfo;
        size_t contigStartPos = contigSide ? (*iter).seedStartPos : contigLength - m_params.searchRange + (*iter).seedStartPos;
        size_t contigEndPos   = contigSide ? (*iter).seedEndPos   : contigLength - m_params.searchRange + (*iter).seedEndPos;
        
        BWTInterval interval  = BWTAlgorithms::findInterval(indices,(*iter).seedStr) ;
        BWTInterval rinterval = BWTAlgorithms::findInterval(indices,reverseComplement((*iter).seedStr)) ;
        
        size_t SeedFreqs = getFrequency(m_params.indices,(*iter).seedStr);

        currentSeedInfo.setContigSeq(contigID, contigPartialSeq, contigLength, contigStartPos, contigEndPos, contigSide);                          

        for(size_t j = interval.lower ; j <= (size_t)interval.upper ; j++)
        {
            std::pair<size_t, size_t> currentSeed = BacktrackReadIdx(indices,j, (m_params.numThreads == 1) );
            currentSeedInfo.setReadSeq((*iter).seedStr, SeedFreqs, currentSeed.first, currentSeed.second, true);
            result.push_back(currentSeedInfo);
        }
            
        for(size_t j = rinterval.lower ; j <= (size_t)rinterval.upper ; j++)
        {
            std::pair<size_t, size_t> currentSeed = BacktrackReadIdx(indices,j, (m_params.numThreads == 1) );
            currentSeedInfo.setReadSeq(reverseComplement((*iter).seedStr), SeedFreqs, currentSeed.first, currentSeed.second, false);
            result.push_back(currentSeedInfo);
        }
    }
}

        ///************************************************///
        ///********* Read Reassembler Post Process *********///
        ///************************************************///

ReassemblerPostProcess::ReassemblerPostProcess(const ReassemblerParameters params, ContigGraph* pGraph) 
: ReadReassemblerBasicElements(params),m_pGraph(pGraph)
{
    collectSeedHashMap.set_deleted_key(-1);
    connectContigHashMap.set_deleted_key("");
    
    singleContigHashMap.set_deleted_key(-1);
    collectAllKmerHashMap.set_deleted_key("");
    
    onePB = ResultCount();
    twoPB = ResultCount();
    
    seedCount_15 = 0;
    seedCount_11 = 0;
}
ReassemblerPostProcess::~ReassemblerPostProcess()
{
}
        
void ReassemblerPostProcess::process(const SequenceWorkItem& item, /*const*/ ReassemblerResult& result)
{    
    // because ReassemblerProcess may use multiple thread
	// so, this process is collected all overlaps between contig and long read
	
	for(SeedByReadHashMap::iterator currentPbIdx = result.seedHash.begin() ; currentPbIdx != result.seedHash.end() ; ++currentPbIdx)
    {
        SeedByReadHashMap::iterator pbIterator = collectSeedHashMap.find( currentPbIdx->first );
        
        if( pbIterator != collectSeedHashMap.end() ) 
        {
            // this read already exist
            // push current contig seeds to hash
            pbIterator->second.push_back(currentPbIdx->second[0]);
        }
        else
        {
            //this read not exist
            collectSeedHashMap.insert(std::make_pair(currentPbIdx->first,currentPbIdx->second));
        }
        // sample kmer freqs
        for(SeedSequenceInfoVec::iterator iter = currentPbIdx->second[0].SeedInfoVec.begin(); iter !=currentPbIdx->second[0].SeedInfoVec.end() ; ++iter)
        {
			std::string read = backtrackRead( m_params.indices, currentPbIdx->first );
			//if( readNum < 100 )countTotalSeedFreq(read);
			if( seedCount_15<100000 && seedCount_11<100000) countTotalSeedFreq(read);
			//if( seedCount_15<100000 && seedCount_11<100000) countTotalSeedFreq((*iter).seedString);
		}
	}
}

void ReassemblerPostProcess::filterErrorRead()
{
    // confirm the overlap between contig and long read
	
	std::cout<< "\n------------filter error connect between contig and read------------\n";

	clock_t p1,p2;
	p1 = clock();
	
    // 15mer distribution using read
	std::sort(kmerSizeFreqByRead_11,kmerSizeFreqByRead_11+seedCount_11);
    std::sort(kmerSizeFreqByRead_15,kmerSizeFreqByRead_15+seedCount_15);
    
	//for(int i = 0 ; i < seedCount_15; i += seedCount_15/1000) std::cout << kmerSizeFreqByRead_15[i] << "\n";
	//showDetailInformation();
	
	// show process number 
    int tmp = 0;
    // check the seeds in the read to determine if the read overlaps with contig
    for(SeedByReadHashMap::iterator pbIterator = collectSeedHashMap.begin() ; pbIterator!=collectSeedHashMap.end() ; ++pbIterator )
    {    
        if(showReassemblerProcess)std::cout<< std::flush << '\r' << "total overlap between contig and read : (" << ++tmp << "/" << collectSeedHashMap.size() << ")";
        
		//showReadFrequency(pbIterator->first,15);
		
        std::vector<OverlapSeedVecInfo>::iterator firstContigIter  ;
        std::vector<OverlapSeedVecInfo>::iterator secondContigIter ;
        
        if( pbIterator->second.size()==2 )
        {   
            // this read overlap two contigs
            onePB.totalCount++;
            std::vector<OverlapSeedVecInfo>::iterator contigIter = pbIterator->second.begin();
            firstContigIter  = contigIter;
            secondContigIter = ++contigIter;
        }
        else if(pbIterator->second.size()>2)
        {   
            // multiple align,
            // this read overlap three or more contigs
            // select the most seeds and the second most seeds
            onePB.partialCount++;
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
        // other case is one read(pacbio/nanopore) connect one contig
        else if( !m_params.isThird ) continue;
        
        // use two reads to determine if two contigs are connected
        if( m_params.isThird )
        {
            // one read align one contig
            if( pbIterator->second.size()==1 )
            {
                std::vector<OverlapSeedVecInfo>::iterator contigIter = pbIterator->second.begin();
                std::string nothing;
                
                combineEveryCheckFunction(*contigIter,*contigIter,nothing,nothing,pbIterator->first,2);
            }
            else
            {
                std::string nothing;
                
                combineEveryCheckFunction(*firstContigIter,*firstContigIter,nothing,nothing,pbIterator->first,2);
                combineEveryCheckFunction(*secondContigIter,*secondContigIter,nothing,nothing,pbIterator->first,2);
            }
        }
        // use a read to determine if two contigs are connected
        if( m_params.isFirst || m_params.isSecond )
        {
            OverlapSeedVecInfo firstResultVec;
            OverlapSeedVecInfo secondResultVec;
            
            filterErrorStrand((*firstContigIter),firstResultVec,mostSeedStrand((*firstContigIter)));
            filterErrorStrand((*secondContigIter),secondResultVec,mostSeedStrand((*secondContigIter)));
            
            std::string firContig = (*firstContigIter).contig;
            std::string secContig = (*secondContigIter).contig;
    
            combineEveryCheckFunction(firstResultVec,secondResultVec,firContig,secContig,pbIterator->first,1);
        }
    }
    
    showDetailInformation();
	
	p2 = clock();
    std::cout<<  (p2 - p1) / CLOCKS_PER_SEC <<" sec\n";
}

void ReassemblerPostProcess::buildGraphByOneRead()
{
    // a long read across two contigs
	
	std::cout<< "\n---------build graph by one read---------\n";
    std::cout<< "contig side number:\t" << connectContigHashMap.size() << "\n";
    
    for(ContigRelationHashMap::iterator contigSideIter = connectContigHashMap.begin() ; contigSideIter!=connectContigHashMap.end() ; ++contigSideIter )
    {
        std::vector<OverlapRelation>::iterator iter ;
        std::string endContig = mostConnectContig(contigSideIter->second);

        if( endContig=="" ) continue;
        else
        {
            std::string startContig;
            int tmpOverlapLength = -1;

            int ec_same = 0;
            int ec_rev  = 0;
            bool ec_vote = true;  // default direct is ec_same
            
            // find most connect direct
            for(std::vector<OverlapRelation>::iterator tmpIter = contigSideIter->second.begin(); tmpIter != contigSideIter->second.end(); ++tmpIter )
                if((*tmpIter).secondContig == endContig)
                    (*tmpIter).isSameStrand ? ec_same++ : ec_rev++;
            
            // it can't determine ec_same or ec_rev
            if(ec_same == ec_rev) continue;
            // change direct to ec_rev
            if(ec_same < ec_rev) ec_vote = false;
            // contigA have most read to contigB, so and contigB
            for(std::vector<OverlapRelation>::iterator tmpIter = contigSideIter->second.begin(); tmpIter != contigSideIter->second.end(); ++tmpIter )
            {
                // check this contig is the most connected target
                if( (*tmpIter).secondContig == endContig )
                {    
                    // check the direction of the contig connection belongs to the direction of the many contigs connected
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
        
        int connectAndOverlap = connectVote((*iter).readIndex,contigSideIter->second);
        // 0:same number of overlap and connect 1:overlap 2:nonOverlap
        if( connectAndOverlap == 0) continue;
        // check two contigs side are not connect, store this pair in hash if not connect 
        if(checkConnect(firstContig,firstSide,secondContig,secondSide))continue;
        
        if( m_params.isFirst && connectAndOverlap!=1 )continue;
        if( m_params.isSecond && connectAndOverlap!=2 )continue;
        
        onePB.sucsess++;
        
            std::cout//<< contigSideIter->first << "\t"
                     //<< contigSideIter->second.size() << "\t"
                     << (connectAndOverlap-1 ? "->...<-" : "-><-   ") << "\t"
                     << (*iter).readIndex              << "\t"
                     //<< (*iter).readLength             << "\t"
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
            
            ContigEdge* edgeLR = new ContigEdge(pVR, firstSide ? ED_ANTISENSE : ED_SENSE, (*iter).isSameStrand ? EC_SAME : EC_REVERSE ,firstSeqCoord ,firstStrand ,secondStrand ,(*iter).readFragment ,(*iter).overlapLength);
            ContigEdge* edgeRL = new ContigEdge(pVL, secondSide ? ED_ANTISENSE : ED_SENSE, (*iter).isSameStrand ? EC_SAME : EC_REVERSE, secondSeqCoord, secondStrand, firstStrand, (*iter).readFragment, (*iter).overlapLength);
            
            edgeLR->setTwin(edgeRL);
            edgeRL->setTwin(edgeLR);
     
            m_pGraph->addEdge(pVL,edgeLR);
            m_pGraph->addEdge(pVR,edgeRL);    
    }
}

void ReassemblerPostProcess::buildGraphByTwoRead()
{
    // each long read overlaps with a single contig,
	// so this step is to check the overlap between the two long reads
	
	clock_t p1,p2,p3,p4;
    
	p1 = clock(); 
    collectSeeds();
	p2 = clock();
    std::cout<<  (p2 - p1) / CLOCKS_PER_SEC <<" sec\n";
	std::cout<< "vector size :" << AllReadKmerVec.size() << "\n";
    
	filterErrorConnect(AllReadGKmerVec,contigPairByGKmer,false);
    filterErrorConnect(AllReadKmerVec,contigPairByKmer,true);
    p3 = clock();
    std::cout<<  (p3 - p2) / CLOCKS_PER_SEC <<" sec\n";
    
	std::cout<<"\nusing gkmer\n";
    layoutGraph(contigPairByGKmer);
	std::cout<<"\nusing fixed-size kmer\n";
	layoutGraph(contigPairByKmer);
    p4 = clock();
    std::cout<<  (p4 - p3) / CLOCKS_PER_SEC <<" sec\n";
    
    std::cout<<"\n";
    std::cout<< "pass read number         : "                    << twoPB.passCount                 << "\n"
             << "case 1. zero seed          : "                  << twoPB.zeroSeed                  << "\n"
             //<< "case 2. repeat status      : "                 << twoPB.repeatStatus              << "\n"
             << "case 2. overlap in repeat  : "                  << twoPB.readRepeatSeeds           << "\n"
             << "case 3. seeds different strand  : "             << twoPB.readSeedDifferentStrand   << "\n"
             << "case 4. seeds order discordant  : "             << twoPB.readSeedOrderDiscordant   << "\n"
             << "case 5. abnormal overlap length : "             << twoPB.readAbnormalOverlapLength << "\n"
             << "case 6. seeds num less than threshold : "       << twoPB.readSeedsThreshold        << "\n"
             << "case 7. abnormal distance between two seeds : " << twoPB.readAbnormalDistance      << "\n"        
             << "case 8. chimeric read : "                       << twoPB.chimera                   << "\n"
             << "reassemble number     : "                       << twoPB.sucsess                   << "\n";
    std::cout<<"\n";
}

void ReassemblerPostProcess::collectSeeds()
{
	// find overlap between two read
	
	std::cout << "\n---------collect seeds on read---------\n";

    ResultCount twoPB = ResultCount();
    int tmp = 0;
    // collect all seeds on read fragment
    for(ReadConnectContigHashMap::iterator readIter = singleContigHashMap.begin() ; readIter!=singleContigHashMap.end() ; ++readIter )
    {
        if(showReassemblerProcess)std::cout<< std::flush << '\r' << "collect seeds on read fragment : (" << ++tmp << "/" << singleContigHashMap.size() << ")";
        std::string originRead  = backtrackRead( m_params.indices, readIter->first );
        std::string connectContig = readIter->second[0].firstContig;
        bool connectContigSide    = readIter->second[0].firstSide;
        int connectContigLength   = readIter->second[0].firstContigLength;
        int startPos = readIter->second[0].singleOverlap.scanReadStart < readIter->second[0].singleOverlap.scanReadEnd ? readIter->second[0].singleOverlap.scanReadStart : readIter->second[0].singleOverlap.scanReadEnd;
        int endPos   = readIter->second[0].singleOverlap.scanReadStart > readIter->second[0].singleOverlap.scanReadEnd ? readIter->second[0].singleOverlap.scanReadStart : readIter->second[0].singleOverlap.scanReadEnd;
        // collect seeds on read
        std::vector<ReassemblerSeedFeature> readGKmerVec = seedingByDynamicKmer(originRead, startPos, endPos);
        std::vector<ReassemblerSeedFeature> readKmerVec  = collectKmer(originRead, startPos, endPos);
        // build AllReadKmerVec
        buildKmerHashAndReadKmerVec( connectContig, connectContigSide, connectContigLength, readIter->first, readGKmerVec, collectAllGKmerHashMap, AllReadGKmerVec);
		buildKmerHashAndReadKmerVec( connectContig, connectContigSide, connectContigLength, readIter->first, readKmerVec, collectAllKmerHashMap, AllReadKmerVec);
        twoPB.passCount++;        
    }
}

void ReassemblerPostProcess::filterErrorConnect(ReadKmerVec &inputReadVec, ContigRelationHashMap &inputContigPair, bool isKmer)
{
    std::cout << "\n---------filter error connect between two read---------\n\n";
    
    // filter error relation connect between two reads
    for(ReadKmerVec::iterator readVecIter = inputReadVec.begin() ; readVecIter != inputReadVec.end() ; ++readVecIter)
    {    
        // use to record already connect reads
        SparseHashMap<int, int, Int64Hasher> connectReadCount;
        // use iterator to find current read connect contig
        ReadConnectContigHashMap::iterator selfReadIter = singleContigHashMap.find((*readVecIter).first);
        // current read connect contig string
        std::string selfContig = selfReadIter->second[0].firstContig; 
        // use to collect useful seeds
        OverlapSeedVecInfo currentSeedVec;
        // after filter, select most seeds connect read
        OverlapSeedVecInfo resultSeedVec;
        // look current read seeds
        for(KmerInfoVec::iterator selfSeedIter = readVecIter->second.begin() ; selfSeedIter != readVecIter->second.end() ; ++selfSeedIter)
        {
            KmerHashMap::iterator seedInfo = collectAllKmerHashMap.find((*selfSeedIter).kmerStr);
            // confirm this seed exist in hash
            if( seedInfo != collectAllKmerHashMap.end() )
            {
                SeedSequenceInfo currentSeed;
                // use to prevent repeat kmer
                int connectOtherReadCount = 0;
                // look for other read if exist current seed
                for(std::vector<KmerInfo>::iterator otherKmerIter = seedInfo->second.begin() ; otherKmerIter != seedInfo->second.end() ; ++otherKmerIter )
                {    
                    // prevent two read connect same contig
                    if( (*otherKmerIter).readConnectContigID == selfContig ) continue;
                    connectOtherReadCount++;
                    currentSeed.contigID       = (*otherKmerIter).readConnectContigID;
                    currentSeed.contigLength   = (*otherKmerIter).readConnectContigLength;
                    currentSeed.contigSide     = (*otherKmerIter).readConnectContigSide;
                    currentSeed.seedLength     = (*otherKmerIter).kmerStr.length();
                    currentSeed.seedString     = (*otherKmerIter).kmerStr;
                    currentSeed.strand         = (*otherKmerIter).kmerStrand;
                    currentSeed.readLocation = (*otherKmerIter).position;
                    currentSeed.readIndex    = (*otherKmerIter).readIndex;
                    // contigStartPos is self read position
                    currentSeed.contigStartPos = (*selfSeedIter).position;
                    // push valid seed to vector
                    currentSeedVec.SeedInfoVec.push_back(currentSeed);
                }
            }
        }
        
        if( currentSeedVec.SeedInfoVec.size() == 0 ){ twoPB.zeroSeed++; continue; }
        
        filterErrorSeeds(currentSeedVec,resultSeedVec,mostConnectRead(currentSeedVec,false));

        ReadConnectContigHashMap::iterator t = singleContigHashMap.find(resultSeedVec.SeedInfoVec[0].readIndex);
         
        // filter read if it concurrently exist different strand seed and all seeds length are rather than 14
        if(checkContigSeed(resultSeedVec,14))       { twoPB.readSeedDifferentStrand++;   continue; }
        // check seeds order concordant or discordant 
        if(checkSeedsOrientations(resultSeedVec))   { twoPB.readSeedOrderDiscordant++;   continue; }
        // check distance between to seeds
        if(checkSeedsDistance(resultSeedVec))       { twoPB.readAbnormalDistance++;      continue; }
        // check overlap length 
        if(checkOverlapLength(resultSeedVec))       { twoPB.readAbnormalOverlapLength++; continue; }
        // kmer threshold 20, gkmer threshold 3
        if( resultSeedVec.SeedInfoVec.size() < (isKmer ? 20 : 3) ) { twoPB.readSeedsThreshold++;  continue; }
        // check repeat by read index, old version. It can't detect every repeat seed.
        if(checkRepeatByRead(resultSeedVec,false))  { twoPB.readRepeatSeeds++;           continue; }
        
		checkRepeatByRead(resultSeedVec,true);
		
        std::cout << "kkmm\t"
                  << selfReadIter->second[0].firstContig   << "\t"
                  << selfReadIter->second[0].firstSide << "\t"
                  << (*readVecIter).first   << "\t"
                  << t->second[0].firstContig << "\t"
                  << t->second[0].firstSide << "\t"
                  << t->second[0].readIndex << "\t"
                  //<< tt->second[0].readIndex << "<\t"
                  //<< currentSeedVec.SeedInfoVec.size() << "\t"
                  << resultSeedVec.SeedInfoVec.size() << "\t"
                  //<< mostReadConnectContig(currentSeedVec) << "\t" 
                  << "\n";

        insertContigRelationHash((*readVecIter),resultSeedVec,inputContigPair);
    }
}

void ReassemblerPostProcess::layoutGraph(ContigRelationHashMap &inputContigPair)
{
    std::cout << "---------build graph by two read---------\n\n";
    std::cout<< "contig side number:\t" << inputContigPair.size() << "\n";
    
    SparseHashMap<int, bool, Int64Hasher> chimeraRecord;
    
    for(ContigRelationHashMap::iterator contigSideIter = inputContigPair.begin() ; contigSideIter!=inputContigPair.end() ; ++contigSideIter )
    {
        std::vector<OverlapRelation>::iterator iter ;
        std::string endContig = mostConnectContig(contigSideIter->second);

        // contigA have most read to connect contigB, so and contigB
        if( endContig=="" ) continue;
        else
        {
            std::string startContig;
            int tmpOverlapLength = -1;

            int ec_same = 0;
            int ec_rev  = 0;
            bool ec_vote = true;
            bool endContigNoOverlap = false;
            
            // find most connect direct
            for(std::vector<OverlapRelation>::iterator tmpIter = contigSideIter->second.begin(); tmpIter != contigSideIter->second.end(); ++tmpIter )
            {
                // check connect contig is correct
                if((*tmpIter).secondContig == endContig)
                {
                    if((*tmpIter).isSameStrand) ec_same++;
                    else ec_rev++;
                }
            }
            
            if(ec_same == ec_rev) continue;
            if(ec_same < ec_rev) ec_vote = false;
            
            // check each 
            for(std::vector<OverlapRelation>::iterator tmpIter = contigSideIter->second.begin(); tmpIter != contigSideIter->second.end(); ++tmpIter )
            {
                if( (*tmpIter).secondContig == endContig )
                {    
                    if((*tmpIter).isSameStrand != ec_vote) continue;
                    
                    if( ((*tmpIter).overlapLength > tmpOverlapLength) && ((*tmpIter).firstContig == startContig || tmpOverlapLength == -1) )
                    {
                        tmpOverlapLength = (*tmpIter).overlapLength;
                        
                        iter = tmpIter;
                        
                        std::string insertKey = (*tmpIter).secondContig+((*iter).secondSide ? "_Head":"_Tail");
                        
                        ContigRelationHashMap::iterator contigIterator = inputContigPair.find(insertKey);

                        if( contigIterator != inputContigPair.end() ) 
                            startContig = mostConnectContig(contigIterator->second);
                        else 
                            endContigNoOverlap = true;
                    }
                }
            }

            if( endContigNoOverlap )
            {
                // only one edge from contigA to contigB and contigB have no overlapping with others
            }
            else if( startContig != (*iter).firstContig ) continue;
        }
        
        std::string firstContig  = (*iter).firstContig;
        std::string secondContig = (*iter).secondContig;
        bool firstSide    = (*iter).firstSide;
        bool secondSide   = (*iter).secondSide;
        bool firstStrand  = (*iter).firstStrand;
        bool secondStrand = (*iter).secondStrand;

        //std::cout<<"beforeCheck" << "\t" << firstContig << "\t" << firstSide << "\t" << secondContig << "\t" << secondSide << "\n";
        
        // check two contigs side are not connect, store this pair in hash if not connect 
        if(checkConnect(firstContig,firstSide,secondContig,secondSide))continue;

        SparseHashMap<int, bool, Int64Hasher>::iterator chimeraIter1 = chimeraRecord.find((int)(*iter).firIdx);
        SparseHashMap<int, bool, Int64Hasher>::iterator chimeraIter2 = chimeraRecord.find((int)(*iter).secIdx);
        
        if( chimeraIter1 != chimeraRecord.end() )
        {
            if( chimeraIter1->second ) continue;
        }
        else
        {
            bool isChimera = checkChimera((*iter).firIdx);
            chimeraRecord.insert(std::make_pair((int)(*iter).firIdx,isChimera));
            if( isChimera )
            {    
                twoPB.chimera++;
                continue;
            }
        }
        
        if( chimeraIter2 != chimeraRecord.end() )
        {
            if( chimeraIter2->second ) continue;
        }
        else
        {
            bool isChimera = checkChimera((*iter).secIdx);
            chimeraRecord.insert(std::make_pair((int)(*iter).secIdx,isChimera));
            if( isChimera )
            {    
                twoPB.chimera++;
                continue;
            }
        }
        
		
        twoPB.sucsess++;
        
            std::cout<<"rree" <<"\t"
                     //<< contigSideIter->first << "\t"
                     //<< contigSideIter->second.size() << "\t"
                     //<< (connectAndOverlap-1 ? "->...<-" : "-><-   ") << "\t"
                     << (*iter).firIdx              << "\t"
                     //<< (*iter).readLength             << "\t"
                     //<< (*iter).overlapLength            << "\t"
                     << firstContig                      << "\t"
                     << (firstSide    ? "Head" : "Tail") << "\t"
                     << (firstStrand  ? "NonRC" : "RC")  << "\t"
                     //<< (*iter).firstContigLength        << "\t"
                     << (firstSide ? "anti" : "sense")   << "\t"
                     << (*iter).secIdx              << "\t"
                     << secondContig                     << "\t"
                     << (secondSide   ? "Head" : "Tail") << "\t"
                     << (secondStrand ? "NonRC" : "RC")  << "\t"
                     //<< (*iter).secondContigLength       << "\t"
                     << (secondSide ? "anti" : "sense")  << "\t"
                     << ((*iter).isSameStrand ? "same" : "reverse")      << "\t"
                     << "\n";
            
            SeqCoord firstSeqCoord((*iter).firContigStart,(*iter).firContigEnd,(*iter).firstContigLength); 
            SeqCoord secondSeqCoord((*iter).secContigStart,(*iter).secContigEnd,(*iter).secondContigLength);
            
            ContigVertex* pVL = ((ContigVertex*)m_pGraph->getVertex(firstContig));
            ContigVertex* pVR = ((ContigVertex*)m_pGraph->getVertex(secondContig));
            
            ContigEdge* edgeLR = new ContigEdge(pVR, firstSide ? ED_ANTISENSE : ED_SENSE, (*iter).isSameStrand ? EC_SAME : EC_REVERSE, firstSeqCoord, firstStrand, secondStrand, (*iter).readFragment, (*iter).overlapLength);
            ContigEdge* edgeRL = new ContigEdge(pVL, secondSide ? ED_ANTISENSE : ED_SENSE, (*iter).isSameStrand ? EC_SAME : EC_REVERSE, secondSeqCoord, secondStrand, firstStrand, (*iter).readFragment, (*iter).overlapLength);
            
            edgeLR->setTwin(edgeRL);
            edgeRL->setTwin(edgeLR);
     
            m_pGraph->addEdge(pVL,edgeLR);
            m_pGraph->addEdge(pVR,edgeRL);    
    }
}

        ///******** filter error read function **********///

void ReassemblerPostProcess::combineEveryCheckFunction(OverlapSeedVecInfo firVec, OverlapSeedVecInfo secVec, std::string firContig, std::string secContig, int readIndex, int section)
{
    std::string originRead  = backtrackRead( m_params.indices, readIndex );
    
    section == 1 ? onePB.totalCount++ : twoPB.totalCount;
    // each contig have three or more seeds on this read
    if( firVec.SeedInfoVec.size()<3 || ( section == 1 ? secVec.SeedInfoVec.size()<3 : false ) ) section == 1 ? onePB.seedsThreshold++ : twoPB.seedsThreshold;
    // filter read if it concurrently exist different strand seed and all seeds length are rather than 19
    else if( checkContigSeed(firVec,19) || ( section == 1 ? checkContigSeed(secVec,19) : false ) ) section == 1 ? onePB.seedDifferentStrand++ : twoPB.seedDifferentStrand;
    // check repeat by contig index.
    else if( checkRepeatByContig(firVec) || ( section == 1 ? checkRepeatByContig(secVec) : false ) ) section == 1 ? onePB.repeatSeeds++ : twoPB.repeatSeeds;
    // check seeds orientations concordant or discordant 
    else if( checkSeedsOrientations(firVec) || ( section == 1 ? checkSeedsOrientations(secVec) : false ) ) section == 1 ? onePB.seedOrderDiscordant++ : twoPB.seedOrderDiscordant++;
    // check distance between to seeds
    else if( checkSeedsDistance(firVec) || ( section == 1 ? checkSeedsDistance(secVec) : false ) ) section == 1 ? onePB.abnormalDistance++ : twoPB.abnormalDistance++;
    // check overlap length 
    else if( checkOverlapLength(firVec) || ( section == 1 ? checkOverlapLength(secVec) : false ) ) section == 1 ? onePB.abnormalOverlapLength++ : twoPB.abnormalOverlapLength++;
    // filter palindrome reads
    else if( checkPalindrome(originRead)) section == 1 ? onePB.palindrome++ : twoPB.palindrome++;
    // different section and each unique check function
    else if( section == 1 )
    {
        bool firstStrand         = firVec.SeedInfoVec[0].strand;
        bool secondStrand        = secVec.SeedInfoVec[0].strand;
        bool firstContigSide     = firVec.SeedInfoVec[0].contigSide;
        bool secondContigSide    = secVec.SeedInfoVec[0].contigSide;
        int firstContigLength    = firVec.SeedInfoVec[0].contigLength;
        int secondContigLength   = secVec.SeedInfoVec[0].contigLength;
        // check two contigs side and strand status
        if(checkContigRelation(firstContigSide,firstStrand,secondContigSide,secondStrand)) onePB.repeatStatus++; 
        else
        {
            /*
            std::cout<< readIndex << "\t"
                     << firContig << "\t"
                     << firstContigSide << "\t"
                     << firstContigLength << "\t"
                     << firVec.SeedInfoVec.size() << "\t"
                     << secContig << "\t"
                     << secondContigSide << "\t"
                     << secondContigLength << "\t"
                     << secVec.SeedInfoVec.size() << "\n";
            */
            onePB.passCount++;
            
            std::vector<OverlapSeedVecInfo> tmp;
            tmp.push_back(firVec);
            tmp.push_back(secVec);
            // build relation hash, use to check how many contig are connect the other contig
            insertContigRelationHash(firContig,firstContigSide,firstStrand,secContig,secondContigSide,secondStrand,firstContigLength,secondContigLength,tmp);
            insertContigRelationHash(secContig,secondContigSide,secondStrand,firContig,firstContigSide,firstStrand,secondContigLength,firstContigLength,tmp);
        
            //std::cout << connectContigHashMap.size() << "\n";
        }
                
    }
    else if( section == 2 )
    {
        ReadScanAndOverlapPosition singleOverlap;
        // check nonOverlap length
        if(checkNonOverlapPartialLength(firVec,singleOverlap)) twoPB.insufficientLength++;
        else
        {
            twoPB.passCount++;
            insertSingleRelationHash(firVec,singleOverlap);
        }
    }
}
        
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
    // filter read if it have different seed status
    int shortSeedCount = 0;
    int reverseSeed = 0;
    bool allSmall = true;
    // check every seeds are same strand
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
    
    bool dir = ((*last).readLocation > (*first).readLocation);
    
    while(last!=input.SeedInfoVec.end()-1)
    {    
        if(dir != ((*last).readLocation > (*first).readLocation) )
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
        int readDistanceBetweenTwoSeeds   = std::abs((*last).readLocation - (*first).readLocation);

        if(std::abs(readDistanceBetweenTwoSeeds-contigDistanceBetweenTwoSeeds)>200 &&
           std::abs(readDistanceBetweenTwoSeeds-contigDistanceBetweenTwoSeeds)>contigDistanceBetweenTwoSeeds*0.2)
                return true;

        first++;
        last++;
    }
    return false;
}
        
bool ReassemblerPostProcess::checkRepeatByRead(OverlapSeedVecInfo input, bool show)
{
    int slideWindow = 15;
    int totalKmerNumber = 0;
    int headCount = 0;
    int tailCount = 0;
    
    if(show)
    {
        for(SeedSequenceInfoVec::iterator iter = input.SeedInfoVec.begin(); iter!=input.SeedInfoVec.end() ; ++iter)
            for( int i = 0 ; i< (int)((*iter).seedString.length()) - slideWindow + 1 ; i++ )
                std::cout<< getFrequency(m_params.indices,(*iter).seedString.substr(i, slideWindow)) << " ";
        std::cout<<  "\n";
        /*
        for(SeedSequenceInfoVec::iterator iter = input.SeedInfoVec.begin(); iter!=input.SeedInfoVec.end() ; ++iter)
            for( int i = 0 ; i< (int)((*iter).seedString.length()) - slideWindow + 1 ; i++ )
                std::cout<< getFrequency(m_params.indices,(*iter).seedString.substr(i, slideWindow)) << " ";
        std::cout<<  "\n";
		*/
    }
    
    for(SeedSequenceInfoVec::iterator iter = input.SeedInfoVec.begin(); iter!=input.SeedInfoVec.end() ; ++iter)
    {
        // check if the first 10 kmers are within repeat
        if( totalKmerNumber >= 10 ) break;
        
        for( int i = 0 ; i< (int)((*iter).seedString.length()) - slideWindow + 1 ; i++ )
        {
            if( totalKmerNumber >= 10 ) break;
            
            totalKmerNumber++;
            
            std::string kmer = (*iter).seedString.substr(i, slideWindow);
            size_t kmerFreqs = getFrequency(m_params.indices,kmer);

            // repeat kmer
            if((int)kmerFreqs>=kmerSizeFreqByRead_15[seedCount_15/10*9])headCount++;
        }
    }
    
    if(show)std::cout<<  "\n";
    
    totalKmerNumber = 0;
    
    for(SeedSequenceInfoVec::iterator iter = --input.SeedInfoVec.end(); iter!=input.SeedInfoVec.begin() ; --iter)
    {
        // check if the last 10 kmers are within repeat
        if( totalKmerNumber >= 10 ) break;
        
        for( int i = 0 ; i< (int)((*iter).seedString.length()) - slideWindow + 1 ; i++ )
        {
            if( totalKmerNumber >= 10 ) break;
            
            totalKmerNumber++;
            
            std::string kmer = (*iter).seedString.substr(i, slideWindow);
            size_t kmerFreqs = getFrequency(m_params.indices,kmer);
 
            // repeat kmer
            if((int)kmerFreqs>=kmerSizeFreqByRead_15[seedCount_15/10*9])tailCount++;
        }
    }

    if(show)std::cout<<  "\n";
    
    if( headCount >= 3 || tailCount >= 3 )return true;

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
        size_t kmerFreqs = getFrequency(m_params.contigIndices,(*iter).seedString);
            
        //    contig1   ----------*-----------*--------*--------
        //                            --------*--------*-------------*------- contig2
        //    pacibo               -------------------------------------
        // freq by contig index   1           2        2             1 
        
        if(kmerFreqs<=2)uniqueOrOverlap++;
        
        //    contig1   ----------*--*---*--           ----*----*--*------- contig2
        //        pacibo        -------------------------------------
        // freq by contig index   1  1   1                 1    1  1 
        
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

int Cmp(const void *lhs, const void *rhs)
{
    return ((const int *)lhs)[1]-((const int *)rhs)[1];
}

bool ReassemblerPostProcess::checkPalindrome(std::string readSeq)
{
    SparseHashMap<std::string, KmerInfo, StringHasher> readKmerHashMap;
    size_t kmerSize = 11;
    
    for(size_t i = 0 ; i+kmerSize <= (size_t)readSeq.length() ; i++)
    {
        std::string kmer = readSeq.substr(i, kmerSize);
        size_t kmerFreqs = getFrequency(m_params.indices,kmer);
        // skip obvious repeat
        if( ((int)kmerFreqs > kmerSizeFreqByRead_11[seedCount_11/10*7]) ) continue;
        
        SparseHashMap<std::string, KmerInfo, StringHasher>::iterator findKmerIter = readKmerHashMap.find(kmer);
        // store current kmer and it's position
        if( findKmerIter == readKmerHashMap.end() ) 
        {
            KmerInfo tmpKmer;
            tmpKmer.position = (int)i;
            tmpKmer.secondPosition = 0;
            tmpKmer.isRepeat = false;
            readKmerHashMap.insert(std::make_pair(kmer,tmpKmer));
        }
        else
        {
            // if this kmer appear more than once
			findKmerIter->second.isRepeat = true;
        }
        
        findKmerIter = readKmerHashMap.find(reverseComplement(kmer));
        
        if( findKmerIter!=readKmerHashMap.end() )
        {
            // mark reverse kmer position
			if(!findKmerIter->second.isRepeat && findKmerIter->second.secondPosition == 0 )
                findKmerIter->second.secondPosition = (int)i;
            else if(!findKmerIter->second.isRepeat)
            {
                // if this kmer appear more than once
				findKmerIter->second.isRepeat = true;
            }
        }
    }

    int tmp = 0;
    int count = 0;

    int positionArray[readSeq.length()][2];
    
    for(SparseHashMap<std::string, KmerInfo, StringHasher>::iterator kmerIter = readKmerHashMap.begin() ; kmerIter != readKmerHashMap.end() ; ++kmerIter )
    {
        // ignore repeat noise
		if(!kmerIter->second.isRepeat && kmerIter->second.secondPosition != 0)
        {
            positionArray[tmp][0] = kmerIter->second.position;
            positionArray[tmp][1] = kmerIter->second.secondPosition;
            tmp++;
        }
    }
    
    qsort(positionArray,tmp,sizeof(positionArray[0]),Cmp);
    
	// check palindrome using continue reverse kmer
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

int Cmp2(const void *lhs, const void *rhs)
{
    //std::cout<< "(" <<((KmerInfo*)lhs)->queryPos << "," << ((KmerInfo*)rhs)->queryPos << ")";
	return ((KmerInfo*)lhs)->queryPos - ((KmerInfo*)rhs)->queryPos;
}

bool ReassemblerPostProcess::checkChimera(size_t queryIndex)
{
    size_t kmerSize = 15;
	// restore the string
    std::string querySeq  = backtrackRead( m_params.indices, queryIndex );
    // used to record the coordinates of the query and the corresponding sequence
    SparseHashMap<int, std::vector<KmerInfo>, Int64Hasher> readKmerHashMap;
    
	// collect all kmer and backtrack each interval
    for(size_t i = 0 ; i+kmerSize <= (size_t)querySeq.length() ; i++)
    {
        std::string kmer = querySeq.substr(i, kmerSize);
        size_t kmerFreqs = getFrequency(m_params.indices,kmer);
        
        BWTInterval interval  = BWTAlgorithms::findInterval(m_params.indices,kmer) ;
        BWTInterval rinterval = BWTAlgorithms::findInterval(m_params.indices,reverseComplement(kmer)) ;
        
        // in order to speed up the time, so filtered repeat kmer
        //if( kmerFreqs >= (size_t)kmerSizeFreqByRead_15[seedCount_15/10*9] ) continue;
        // in order to speed up the time, so filtered repeat kmer
        if( kmerFreqs >= 50 ) continue;
		
        for(size_t j = interval.lower ; j <= (size_t)interval.upper ; j++)
        {
            // filter self index
            if(interval.lower + 1 == interval.upper )continue;
            
            std::pair<size_t, size_t> currentSeed = BacktrackReadIdx(m_params.indices,j,true);
            // filter self index
            if( queryIndex == currentSeed.first )continue;
            
            SparseHashMap<int, std::vector<KmerInfo>, Int64Hasher>::iterator readKmerIter = readKmerHashMap.find(currentSeed.first);
            
            KmerInfo tmp;
            tmp.queryPos = i;
            tmp.sbjctPos = currentSeed.second;
            
            if(readKmerIter!=readKmerHashMap.end())
            {
                readKmerIter->second.push_back(tmp);
            }
            else
            {
                std::vector<KmerInfo> tmpVec;
                tmpVec.push_back(tmp);
                readKmerHashMap.insert(std::make_pair(currentSeed.first,tmpVec));
            }
            //std::cout<< currentSeed.first << "\t" << tmp.queryPos << "\t" << tmp.sbjctPos << "\n";
        }
         
        for(size_t j = rinterval.lower ; j <= (size_t)rinterval.upper ; j++)
        {
            if(interval.lower + 1 == interval.upper )continue;
            
            std::pair<size_t, size_t> currentSeed = BacktrackReadIdx(m_params.indices,j,true);
            
            if( queryIndex == currentSeed.first )continue;
            
            SparseHashMap<int, std::vector<KmerInfo>, Int64Hasher>::iterator readKmerIter = readKmerHashMap.find(currentSeed.first);
            
            KmerInfo tmp;
            tmp.queryPos = i;
            tmp.sbjctPos = currentSeed.second;
            
            if(readKmerIter!=readKmerHashMap.end())
            {
                readKmerIter->second.push_back(tmp);
            }
            else
            {
                std::vector<KmerInfo> tmpVec;
                tmpVec.push_back(tmp);
                readKmerHashMap.insert(std::make_pair(currentSeed.first,tmpVec));
            }
        }
    }
    
    //std::cout<< queryIndex << "checkChimera\n";
    //std::cout<< "query length   : " << querySeq.length() << "\n";
    // sbjct correspond to the query position, check chimeric read by line cover problem
	std::vector<KmerInfo> querysbjctPosPairVec;
    // if all sbjct fragment can't cover query, check each sbjct prominent length (unaligned length)
	std::vector<ReadScanAndOverlapPosition> sbjctStartEndPosPairVec ;
	
    int pairCount = -1;
    
    // check each overlap read by number of kmer and distance between two kmer
    for(SparseHashMap<int, std::vector<KmerInfo>, Int64Hasher>::iterator readKmerIter = readKmerHashMap.begin();readKmerIter != readKmerHashMap.end();readKmerIter++)
    {
        // filter most error read, because of repeat or sequencing error
        if( readKmerIter->second.size() < 10 )continue;
        
        std::string sbjct  = backtrackRead( m_params.indices, readKmerIter->first );
        //std::cout<< "query length   : " << querySeq.length() << "\n";
        //std::cout<< "sbjct length : " << sbjct.length() << "\n";
        //std::cout<< "kmer number    : " << readKmerIter->second.size() <<"\n";
        
        int queryPos  = -1;
        int sbjctPos = -1;
        bool error = false;
        
        for(std::vector<KmerInfo>::iterator vecIter = readKmerIter->second.begin(); vecIter != readKmerIter->second.end() ; vecIter++ )
        {
            if( queryPos == -1 || sbjctPos == -1 )
            {
                queryPos  = (*vecIter).queryPos;
                sbjctPos = (*vecIter).sbjctPos;
                continue;
            }
            
            int queryDis  = std::abs( queryPos  - (int)(*vecIter).queryPos  );
            int sbjctDis = std::abs( sbjctPos - (int)(*vecIter).sbjctPos );
            
            //std::cout<< queryDis << "\t" << sbjctDis << "\n";
            
            if( ( std::abs(sbjctDis-queryDis)>200 && std::abs(sbjctDis-queryDis)>queryDis*0.2 ) )
            {
                error = true;
                break;
            }
            
            queryPos  = (*vecIter).queryPos;
            sbjctPos = (*vecIter).sbjctPos;

        }

        std::vector<KmerInfo>::iterator beginIter = readKmerIter->second.begin();
        std::vector<KmerInfo>::iterator endIter   = --readKmerIter->second.end();

        float overlapQueryLength = std::abs( (int)(*beginIter).queryPos - (int)(*endIter).queryPos );
        float overlapSbjctLength = std::abs( (int)(*beginIter).sbjctPos - (int)(*endIter).sbjctPos );
        float lengthSimilarityRatio = overlapQueryLength / overlapSbjctLength ;
        
        if( error )continue;
        if( lengthSimilarityRatio > 1.25 || lengthSimilarityRatio < 0.8 ) continue;
        if( std::abs( (int)(*beginIter).queryPos - (int)(*endIter).queryPos ) < 400 ) continue;
        
        //std::cout<< "sbjct index  : " << readKmerIter->first << "\t";
        //std::cout<< (*beginIter).queryPos << " ~ " << (*endIter).queryPos << "\t";
        //std::cout<< (*beginIter).sbjctPos << " ~ " << (*endIter).sbjctPos << "\t( " << sbjct.length()  << " )" <<"\n";
        
        KmerInfo tmp1,tmp2;
        pairCount++;
		
        tmp1.pairNumber = pairCount;
        tmp2.pairNumber = pairCount;
        
        tmp1.queryPos = (*beginIter).queryPos;
        tmp2.queryPos = (*endIter).queryPos;
		
        querysbjctPosPairVec.push_back(tmp1);
        querysbjctPosPairVec.push_back(tmp2);
        
        ReadScanAndOverlapPosition tmp3;
        
        tmp3.scanReadStart     = (*beginIter).queryPos;
        tmp3.scanReadEnd       = (*endIter).queryPos;
        tmp3.readFragmentStart = (*beginIter).sbjctPos;
        tmp3.readFragmentEnd   = (*endIter).sbjctPos;
        tmp3.readLength        = sbjct.length();
        
        sbjctStartEndPosPairVec.push_back(tmp3);
		
    }

	for(std::vector<KmerInfo>::iterator iter_i = querysbjctPosPairVec.begin(); iter_i != querysbjctPosPairVec.end() ; iter_i++ )
	{
		for(std::vector<KmerInfo>::iterator iter_j = iter_i; iter_j != querysbjctPosPairVec.end() ; iter_j++ )
		{
			if( (*iter_i).queryPos > (*iter_j).queryPos )
			{
				int tmpPos  = (*iter_i).queryPos;
				int tmpPair = (*iter_i).pairNumber;
				
				(*iter_i).queryPos   = (*iter_j).queryPos;
				(*iter_i).pairNumber = (*iter_j).pairNumber;
				
				(*iter_j).queryPos   = tmpPos;
				(*iter_j).pairNumber = tmpPair;
			}
		}
	}

    //qsort(querysbjctPosPairVec.data(), querysbjctPosPairVec.size(),sizeof(KmerInfo),Cmp2);
    
    std::map<int,int> pairMap; 
    int startAnchor,searchAnchor;
    int fragment = 0;
    int fragmentStart;
    int fragmentEnd;
    
    startAnchor = 0;
    searchAnchor = 0;
    
    while( startAnchor < (int)querysbjctPosPairVec.size() )
    {
        fragment++;
        if( fragment > 1 ) break;
        
        //std::cout<< querysbjctPosPairVec[startAnchor].queryPos << "~";
        fragmentStart = querysbjctPosPairVec[startAnchor].queryPos;
        pairMap[querysbjctPosPairVec[startAnchor].pairNumber]=1;
        searchAnchor++;
        
        while( startAnchor != searchAnchor )
        {
            while( pairMap[querysbjctPosPairVec[startAnchor].pairNumber] != 2 )
            {
                pairMap[querysbjctPosPairVec[searchAnchor].pairNumber]++;
                searchAnchor++;
            }
            
            while( pairMap[querysbjctPosPairVec[startAnchor].pairNumber] != 1 )
            {
                startAnchor++;
                if( startAnchor == searchAnchor ) break;
            }
        }    
        //std::cout<< querysbjctPosPairVec[startAnchor-1].queryPos << "\n";
        fragmentEnd = querysbjctPosPairVec[startAnchor-1].queryPos;
    }
    
    if( fragment == 1 )
    {
        float overlappingRatio = (float)( fragmentEnd - fragmentStart ) / (float)querySeq.length();
        // prevent the query is too short
        int nonOverlapLength = querySeq.length() - fragmentEnd + fragmentStart - 1;
        if( !(overlappingRatio >= 0.9 || nonOverlapLength <= 1000) ) fragment = 2;
    }

    float chimericRead = 0;
    float normalRead = 0;
    
	if( fragment != 1 )
	{
		for( std::vector<ReadScanAndOverlapPosition>::iterator iter = sbjctStartEndPosPairVec.begin() ; iter != sbjctStartEndPosPairVec.end() ; iter++ )
		{
			int leftProminentLength  = -1;
			int rightProminentLength = -1;
			
			// left prominent length             |----|
			// query                             ----------------------------
			//                                        ||||||||
			// sbjct                         ----------------------
			// right prominent length                        |-----|
			
			// sbjct direction  --->
			if( (*iter).readFragmentStart < (*iter).readFragmentEnd )
			{
				leftProminentLength  = std::min( (*iter).scanReadStart , (*iter).readFragmentStart );
				rightProminentLength = std::min( (int)querySeq.length() - (*iter).scanReadEnd , (*iter).readLength - (*iter).readFragmentEnd );
			}
			// sbjct direction  <---
			else
			{
				leftProminentLength  = std::min( (*iter).scanReadStart , (*iter).readLength - (*iter).readFragmentStart );
				rightProminentLength = std::min( (int)querySeq.length() - (*iter).scanReadEnd , (*iter).readFragmentEnd );
			}
			
			assert( leftProminentLength != -1 && rightProminentLength != -1 );
			
			if( leftProminentLength > 1500 ) chimericRead++;
			else normalRead++;
			
			if( rightProminentLength > 1500 ) chimericRead++;
			else normalRead++;
		}
	}

    if( chimericRead / ( normalRead + normalRead ) >= 0.3 )
	{
		std::cout<< "chimera : " << queryIndex << "\n";
		return true;
	}
    return false;
}

bool ReassemblerPostProcess::checkOverlapLength(OverlapSeedVecInfo input)
{
    SeedSequenceInfoVec::iterator first = input.SeedInfoVec.begin();
    SeedSequenceInfoVec::iterator last  = input.SeedInfoVec.end();
    --last;
    
    if(std::abs((int)(*first).contigStartPos - (int)(*last).contigStartPos)<400) return true;
    if(std::abs((int)(*first).readLocation - (int)(*last).readLocation)<400) return true;
    
    return false;
}

bool ReassemblerPostProcess::checkNonOverlapPartialLength(OverlapSeedVecInfo input, ReadScanAndOverlapPosition &position)
{
    bool contigSide  = input.SeedInfoVec[0].contigSide;
    int readIndex  = input.SeedInfoVec[0].readIndex;
    int contigLength = input.SeedInfoVec[0].contigLength;
    std::string originRead = backtrackRead( m_params.indices, readIndex );
    
    ReadScanAndOverlapPosition reScanPosition = reScanContigOnReadRange( input.SeedInfoVec[0].contigPartialSeq, originRead );
    
    int thresholdLength = 400; 
    
    SeedSequenceInfoVec::iterator first = input.SeedInfoVec.begin();
    SeedSequenceInfoVec::iterator last  = input.SeedInfoVec.end();
    last--;

    // contig head side
    if(contigSide)
    {
        // read dir ---->
        if( (*first).readLocation < (*last).readLocation )
        {
            if( (int)(*first).readLocation - (int)(*first).contigStartPos < thresholdLength ) return true;
            position.scanReadStart = 0;
            position.scanReadEnd   = (*first).readLocation - (*first).contigStartPos;
            position.readFragmentStart = 0;
            position.readFragmentEnd   = (*first).readLocation;
            //position.cutContigToEndPointLength = (*first).contigStartPos;
            
            if(reScanPosition.kmerCount>10)
            {
                position.scanReadEnd     = reScanPosition.scanReadStart;
                position.readFragmentEnd = reScanPosition.scanReadStart;
            }
        }
        // read dir <----
        else
        {
            if( (int)originRead.length() - ( (*first).readLocation + (*first).contigStartPos ) < thresholdLength ) return true;
            position.scanReadStart = (*first).readLocation + (*first).contigStartPos;
            position.scanReadEnd   = originRead.length();
            position.readFragmentStart = (*first).readLocation;
            position.readFragmentEnd   = originRead.length();
            //position.cutContigToEndPointLength = (*first).contigStartPos;
            
            if(reScanPosition.kmerCount>10)
            {
                position.scanReadEnd     = reScanPosition.scanReadEnd;
                position.readFragmentEnd = reScanPosition.scanReadEnd;
            }
        }
    }
    // contig tail side
    else
    {
        int toEndLength = contigLength - (*last).contigStartPos;
        // read dir ---->
        if( (*first).readLocation < (*last).readLocation )
        {
            if( ( (int)originRead.length() - ( (*last).readLocation + toEndLength ) ) < thresholdLength ) return true;
            position.scanReadStart = (*last).readLocation + toEndLength;
            position.scanReadEnd   = originRead.length();
            position.readFragmentStart = (*last).readLocation;
            position.readFragmentEnd   = originRead.length();
            //position.cutContigToEndPointLength = toEndLength;
            
            if(reScanPosition.kmerCount>10)
            {
                position.scanReadEnd     = reScanPosition.scanReadEnd;
                position.readFragmentEnd = reScanPosition.scanReadEnd;
            }
        }
        // read dir <----
        else
        {
            if( (*last).readLocation - toEndLength < thresholdLength ) return true;
            position.scanReadStart = 0;
            position.scanReadEnd   = (*first).readLocation - toEndLength;
            position.readFragmentStart = 0;
            position.readFragmentEnd   = (*last).readLocation;
            //position.cutContigToEndPointLength = toEndLength;
            
            if(reScanPosition.kmerCount>10)
            {
                position.scanReadEnd     = reScanPosition.scanReadStart;
                position.readFragmentEnd = reScanPosition.scanReadStart;
            }
        }
    }
    
    //std::cout<<position.scanReadStart << "\t" << position.scanReadEnd <<"   range\n";
    return false;
}

        ///**************** tool function *****************///

void ReassemblerPostProcess::showDetailInformation()
{    
    std::cout << "\n";
    
    std::cout << "read's 11mer freq by read indices :\n"  
              << kmerSizeFreqByRead_11[seedCount_11/10*1] << "\t" << kmerSizeFreqByRead_11[seedCount_11/10*2] << "\t"
              << kmerSizeFreqByRead_11[seedCount_11/10*3] << "\t" << kmerSizeFreqByRead_11[seedCount_11/10*4] << "\t"
              << kmerSizeFreqByRead_11[seedCount_11/10*5] << "\t" << kmerSizeFreqByRead_11[seedCount_11/10*6] << "\t"
              << kmerSizeFreqByRead_11[seedCount_11/10*7] << "\t" << kmerSizeFreqByRead_11[seedCount_11/10*8] << "\t"
              << kmerSizeFreqByRead_11[seedCount_11/10*9] << "\t" << kmerSizeFreqByRead_11[seedCount_11-1]    << "\n";
    
    std::cout << "read's 15mer freq by read indices :\n"  
              << kmerSizeFreqByRead_15[seedCount_15/10*1] << "\t" << kmerSizeFreqByRead_15[seedCount_15/10*2] << "\t"
              << kmerSizeFreqByRead_15[seedCount_15/10*3] << "\t" << kmerSizeFreqByRead_15[seedCount_15/10*4] << "\t"
              << kmerSizeFreqByRead_15[seedCount_15/10*5] << "\t" << kmerSizeFreqByRead_15[seedCount_15/10*6] << "\t"
              << kmerSizeFreqByRead_15[seedCount_15/10*7] << "\t" << kmerSizeFreqByRead_15[seedCount_15/10*8] << "\t"
              << kmerSizeFreqByRead_15[seedCount_15/10*9] << "\t" << kmerSizeFreqByRead_15[seedCount_15-1]    << "\n";
                                        
    std::cout<<"\n";
    
    if( m_params.isThird )
    {    
        std::cout<< "read align one contig : "                       << twoPB.totalCount            << "\n"
                 << "pass read number      : "                       << twoPB.passCount             << "\n"
                 << "case 1. seeds num less than two : "             << twoPB.seedsThreshold        << "\n"
                 << "case 2. seeds different strand  : "             << twoPB.seedDifferentStrand   << "\n"
                 << "case 3. seeds order discordant  : "             << twoPB.seedOrderDiscordant   << "\n"
                 << "case 4. abnormal overlap length : "             << twoPB.abnormalOverlapLength << "\n"
                 << "case 5. abnormal distance between two seeds : " << twoPB.abnormalDistance      << "\n"
                 << "case 6. insufficient non overlap length     : " << twoPB.insufficientLength    << "\n"
                 << "case 7. palindrome reads : "                    << twoPB.palindrome            << "\n";
    }    
    if( m_params.isFirst || m_params.isSecond )
    {
        std::cout<< "read align two contig    : "                    << onePB.totalCount            << "\n"
                 << "align more than two contig : "                  << onePB.partialCount          << "\n"
                 << "pass read number      : "                       << onePB.passCount             << "\n"
                 << "case 1. two reads are same side : "             << onePB.repeatStatus          << "\n"
                 << "case 2. two flank are repeats   : "             << onePB.repeatSeeds           << "\n"
                 << "case 3. seeds num less than two : "             << onePB.seedsThreshold        << "\n"
                 << "case 4. seeds different strand  : "             << onePB.seedDifferentStrand   << "\n"
                 << "case 5. seeds order discordant  : "             << onePB.seedOrderDiscordant   << "\n"
                 << "case 6. abnormal overlap length : "             << onePB.abnormalOverlapLength << "\n"
                 << "case 7. abnormal distance between two seeds : " << onePB.abnormalDistance      << "\n" 
                 << "case 8. palindrome reads : "                    << onePB.palindrome            << "\n";
                 //<< "reassemble number : "                           << onePB.sucsess               << "\n";
    }
}
        
void ReassemblerPostProcess::insertContigRelationHash(std::string firstContig, bool firstSide, bool firstStrand, std::string secondContig, bool secondSide, bool secondStrand, int firstContigLength, int secondContigLength, std::vector<OverlapSeedVecInfo> contigInfoVec)
{
    OverlapRelation contigPair;
    std::vector<OverlapRelation> insertVec;
    std::string insertKey = firstContig+(firstSide ? "_Head":"_Tail");

    std::string originRead = backtrackRead( m_params.indices, (int)contigInfoVec[0].SeedInfoVec[0].readIndex );
    contigReadPosition f   = OnReadRange( contigInfoVec[0].SeedInfoVec , (int64_t)originRead.length() );
    contigReadPosition s   = OnReadRange( contigInfoVec[1].SeedInfoVec , (int64_t)originRead.length() );

    contigPair.firstContig   = firstContig;
    contigPair.secondContig  = secondContig;
    contigPair.firstSide     = firstSide;
    contigPair.secondSide    = secondSide;
    contigPair.firstStrand   = firstStrand;
    contigPair.secondStrand  = secondStrand;
    
    contigPair.readLength  = originRead.length();
    contigPair.readIndex   = contigInfoVec[0].SeedInfoVec[0].readIndex;
    contigPair.connect       = (f.maxReadPos < s.minReadPos || s.maxReadPos < f.minReadPos);
    
    contigPair.firstContigLength  = firstContigLength;
    contigPair.secondContigLength = secondContigLength;

    if( (f.firstReadStartPos < f.lastReadStartPos && s.firstReadStartPos < s.lastReadStartPos) || 
        (f.firstReadStartPos > f.lastReadStartPos && s.firstReadStartPos > s.lastReadStartPos))
        contigPair.isSameStrand = true;
    else
        contigPair.isSameStrand = false;

    //non overlap
    if(contigPair.connect)
    {
        if(f.minReadPos > s.maxReadPos)
            contigPair.readFragment = originRead.substr(s.maxReadPos,f.minReadPos - s.maxReadPos);
        else
            contigPair.readFragment = originRead.substr(f.maxReadPos,s.minReadPos - f.maxReadPos); 
        
        contigPair.overlapLength  = 0;
        contigPair.firContigStart = 0;
        contigPair.firContigEnd   = 0;
        contigPair.secContigStart = 0;
        contigPair.secContigEnd   = 0;
    }
    // overlap
    else
    {  
        contigPair.overlapLength = std::abs(f.contigEndPointOnReadPos-s.contigEndPointOnReadPos) + f.offsetOverlapLength + s.offsetOverlapLength;
        
        std::pair<int,int> firOverlap = overlapPosition(firstSide,contigPair.firstContigLength,contigPair.overlapLength);
        std::pair<int,int> secOverlap = overlapPosition(secondSide,contigPair.secondContigLength,contigPair.overlapLength);
        
        contigPair.firContigStart = firOverlap.first;
        contigPair.firContigEnd   = firOverlap.second;
        contigPair.secContigStart = secOverlap.first;
        contigPair.secContigEnd   = secOverlap.second;
    }

    // prevent different read build same relation edge  
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

void ReassemblerPostProcess::insertContigRelationHash(ReadKmerPair readVecIter, OverlapSeedVecInfo resultSeedVec, ContigRelationHashMap &inputContigPair)
{
    
    ReadConnectContigHashMap::iterator selfReadIter = singleContigHashMap.find(readVecIter.first);
    ReadConnectContigHashMap::iterator targetReadIter = singleContigHashMap.find(resultSeedVec.SeedInfoVec[0].readIndex);
    
    OverlapRelation connectResult;

    connectResult.firstContig        = selfReadIter->second[0].firstContig;
    connectResult.firstContigLength  = selfReadIter->second[0].firstContigLength;
    connectResult.firstSide          = selfReadIter->second[0].firstSide;
    connectResult.firstStrand        = selfReadIter->second[0].firstStrand;
        
    connectResult.secondContig       = targetReadIter->second[0].firstContig;
    connectResult.secondContigLength = targetReadIter->second[0].firstContigLength;
    connectResult.secondSide         = targetReadIter->second[0].firstSide;
    connectResult.secondStrand       = (resultSeedVec.SeedInfoVec[0].strand ?  targetReadIter->second[0].firstStrand : !targetReadIter->second[0].firstStrand );

    connectResult.firContigStart = 0;
    connectResult.firContigEnd   = 0;
    connectResult.secContigStart = 0;
    connectResult.secContigEnd   = 0;
    connectResult.overlapLength  = 0;
        
    connectResult.connect = true;
        
    connectResult.isSameStrand = checkSameStrand( readVecIter.first, resultSeedVec.SeedInfoVec[0].readIndex, resultSeedVec.SeedInfoVec);

    connectResult.firIdx = readVecIter.first;
    connectResult.secIdx = resultSeedVec.SeedInfoVec[0].readIndex;

    std::string firstReadFragment;
    std::string secondReadFragment = concordantNonOverlapFragment(resultSeedVec);
    std::string originRead  = backtrackRead( m_params.indices, readVecIter.first );
        
    if( selfReadIter->second[0].singleOverlap.scanReadStart == 0 ) 
    {
        SeedSequenceInfoVec::iterator first = resultSeedVec.SeedInfoVec.begin();
        SeedSequenceInfoVec::iterator last = resultSeedVec.SeedInfoVec.end();
        last--;
            
        int start = (*first).contigStartPos > (*last).contigStartPos ? (*first).contigStartPos : (*last).contigStartPos;
        int end   = selfReadIter->second[0].singleOverlap.scanReadEnd;
        firstReadFragment  = originRead.substr( start, end - start );
            
        //firstReadFragment  = originRead.substr( selfReadIter->second[0].singleOverlap.scanReadStart, selfReadIter->second[0].singleOverlap.scanReadEnd );
        connectResult.readFragment = secondReadFragment + firstReadFragment;
            
        //std::cout<< "start1\n";
        //std::cout<< (*first).contigStartPos << "\t" << (*last).contigStartPos <<"\n";
        //std::cout<< start << "\t" << end << "\t idx " << readVecIter.first << "\n";
    }
    else 
    {
        SeedSequenceInfoVec::iterator first = resultSeedVec.SeedInfoVec.begin();
        SeedSequenceInfoVec::iterator last = resultSeedVec.SeedInfoVec.end();
        last--;
            
        int start = selfReadIter->second[0].singleOverlap.scanReadStart;
        int end   = (*first).contigStartPos < (*last).contigStartPos ? (*first).contigStartPos : (*last).contigStartPos;
            
        firstReadFragment  = originRead.substr( start, end - start );
            
        //firstReadFragment  = originRead.substr( selfReadIter->second[0].singleOverlap.scanReadStart, originRead.length() - selfReadIter->second[0].singleOverlap.scanReadStart );
        connectResult.readFragment = firstReadFragment + secondReadFragment;
            
        //std::cout<< "end1\n";
        //std::cout<< (*first).contigStartPos << "\t" << (*last).contigStartPos <<"\n";
        //std::cout<< start << "\t" << end << "\t idx " << readVecIter.first << "\n";
    }
    
    std::string insertKey = selfReadIter->second[0].firstContig+(selfReadIter->second[0].firstSide ? "_Head":"_Tail");
    std::vector<OverlapRelation> insertVec;
    // prevent different read build same relation edge  
    ContigRelationHashMap::iterator contigIterator = inputContigPair.find(insertKey);
    //this contigSide already exist
    if( contigIterator != inputContigPair.end() ) 
    {
        // push this contigSide to hash
        contigIterator->second.push_back(connectResult);
    }
    //this contigSide not exist
    else
    {
        insertVec.push_back(connectResult);
        inputContigPair.insert(std::make_pair(insertKey,insertVec));
    }
}

void ReassemblerPostProcess::insertSingleRelationHash(OverlapSeedVecInfo input, ReadScanAndOverlapPosition scanRange)
{
    OverlapRelation contigPair;
    std::vector<OverlapRelation> insertVec;
    
    std::string originRead = backtrackRead( m_params.indices, (int)input.SeedInfoVec[0].readIndex );
    contigReadPosition f   = OnReadRange( input.SeedInfoVec, (int)originRead.length() );
    
    contigPair.firstContig   = input.contig;
    contigPair.firstSide     = input.SeedInfoVec[0].contigSide;
    contigPair.firstStrand   = input.SeedInfoVec[0].strand;
    contigPair.readLength  = originRead.length();
    contigPair.readIndex   = input.SeedInfoVec[0].readIndex;
    contigPair.firstContigLength  = input.SeedInfoVec[0].contigLength;
    //contigPair.scanStart = scanRange.first;
    //contigPair.scanEnd   = scanRange.second;
    contigPair.singleOverlap = scanRange;
    
    contigPair.readDirection = f.firstReadStartPos < f.lastReadStartPos;
    
    ReadConnectContigHashMap::iterator readIterator = singleContigHashMap.find(contigPair.readIndex);
    
    if( readIterator != singleContigHashMap.end() ) 
    {
        readIterator->second.push_back(contigPair);
    }
    else
    {
        insertVec.push_back(contigPair);
        singleContigHashMap.insert(std::make_pair(contigPair.readIndex,insertVec));
    }
    
    //std::cout << "(debug)contig : " << input.contig << "\tread : " << input.SeedInfoVec[0].readIndex << "\n";
}

contigReadPosition ReassemblerPostProcess::OnReadRange(SeedSequenceInfoVec inputSeedInfoVec, int64_t readLength)
{
    contigReadPosition result;
    
    int64_t firstContigStartPos = inputSeedInfoVec[0].contigStartPos;
    int64_t lastContigStartPos  = inputSeedInfoVec[0].contigStartPos;
    int64_t firstReadStartPos = inputSeedInfoVec[0].readLocation;
    int64_t lastReadStartPos  = inputSeedInfoVec[0].readLocation;
    int64_t contigEndPointOnReadPos;

    result.offsetOverlapLength = 0;
    
    std::string originRead = backtrackRead( m_params.indices, (int)inputSeedInfoVec[0].readIndex );
    ReadScanAndOverlapPosition reScanPosition = reScanContigOnReadRange( inputSeedInfoVec[0].contigPartialSeq, originRead );
    
    
    for(SeedSequenceInfoVec::iterator iter =inputSeedInfoVec.begin(); iter!=inputSeedInfoVec.end() ; ++iter)
    {
        if( (*iter).contigStartPos < firstContigStartPos )
        {
            firstContigStartPos = (*iter).contigStartPos;
            firstReadStartPos = (*iter).readLocation;
        }
        if( (*iter).contigStartPos + (*iter).seedLength > lastContigStartPos )
        {
            lastContigStartPos = (*iter).contigStartPos + (*iter).seedLength;
            lastReadStartPos = (*iter).readLocation + (*iter).seedLength;
        }
    }
    
    //head
    if(inputSeedInfoVec[0].contigSide)
    {
        //SENSE
        if( firstReadStartPos < lastReadStartPos )
        {
            contigEndPointOnReadPos = std::max( firstReadStartPos - firstContigStartPos, (int64_t)0 );
            
            if(reScanPosition.kmerCount>10)
                contigEndPointOnReadPos  = reScanPosition.scanReadStart;
            
            //if( firstContigStartPos > firstReadStartPos )
            //    result.offsetOverlapLength = std::abs(firstReadStartPos - firstContigStartPos);
        }
        //ANTISENSE
        else 
        {    
            contigEndPointOnReadPos = std::min( firstReadStartPos + firstContigStartPos, readLength );
            
            if(reScanPosition.kmerCount>=10)
                contigEndPointOnReadPos  = reScanPosition.scanReadEnd;
            
            //if( firstReadStartPos - firstContigStartPos > readLength)
            //    result.offsetOverlapLength = std::abs(firstReadStartPos + firstContigStartPos - readLength );
        }
    }
    //tail
    else
    {
        int64_t contigLength = inputSeedInfoVec[0].contigLength;
        //SENSE
        if( firstReadStartPos < lastReadStartPos ) 
        {
            contigEndPointOnReadPos = std::min( lastReadStartPos + (contigLength - lastContigStartPos), readLength );
            
            if(reScanPosition.kmerCount>=10)
                contigEndPointOnReadPos  = reScanPosition.scanReadEnd;
            
            //if( lastReadStartPos + (contigLength - lastContigStartPos) > readLength)
            //    result.offsetOverlapLength = std::abs( lastReadStartPos + (contigLength - lastContigStartPos) - readLength );
        }
        //ANTISENSE
        else     
        {
            contigEndPointOnReadPos = std::max( lastReadStartPos - (contigLength - lastContigStartPos) , (int64_t)0 );
            
            if(reScanPosition.kmerCount>10)
                contigEndPointOnReadPos  = reScanPosition.scanReadStart;
            
            //if( lastReadStartPos - (contigLength - lastContigStartPos) < 0)
            //    result.offsetOverlapLength = std::abs(lastReadStartPos - (contigLength - lastContigStartPos));
        }
    }

    result.minReadPos = std::min(lastReadStartPos,firstReadStartPos);
    result.minReadPos = std::min(result.minReadPos,contigEndPointOnReadPos);
    
    result.maxReadPos = std::max(lastReadStartPos,firstReadStartPos);
    result.maxReadPos = std::max(result.maxReadPos,contigEndPointOnReadPos);
    
    result.firstContigStartPos = firstContigStartPos;
    result.lastContigStartPos  = lastContigStartPos;
    result.firstReadStartPos = firstReadStartPos;
    result.lastReadStartPos  = lastReadStartPos;
    result.contigEndPointOnReadPos = contigEndPointOnReadPos;
    
    return result;
}

ReadScanAndOverlapPosition ReassemblerPostProcess::reScanContigOnReadRange(std::string partialContig, std::string originRead)
{
    ReadScanAndOverlapPosition result;
    
    SparseHashMap<std::string, KmerInfo, StringHasher> contigPartialKmerHashMap;
    std::vector<ReadScanAndOverlapPosition> rePositionVec;
    
    int minPosition = 0;
    
    size_t kmerSize = 11;
    
    result.kmerCount = 0;
    result.scanReadStart = 0;
    result.scanReadEnd   = 0;

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
    
    for(size_t i = 0 ; i+kmerSize <= originRead.length() ; i++)
    {
        std::string forwardkmer = originRead.substr(i, kmerSize);
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
    int contigOnReadPosition[partialContig.length()];
    
    for(SparseHashMap<std::string, KmerInfo, StringHasher>::iterator kmerIter = contigPartialKmerHashMap.begin() ; kmerIter!=contigPartialKmerHashMap.end() ; ++kmerIter)
    {
        if(!kmerIter->second.isRepeat && kmerIter->second.position != -1)
        {
            contigOnReadPosition[tmp] = kmerIter->second.position;
            tmp++;
        }
    }
    
    std::sort(contigOnReadPosition,contigOnReadPosition+tmp);
    
    for(int i =1;i<tmp;i++)
    {
        if( (std::abs(contigOnReadPosition[i] - contigOnReadPosition[i-1])) < 40 )
        {
            if(minPosition==0)
            {
                minPosition = contigOnReadPosition[i-1];
                ReadScanAndOverlapPosition tmpPosition;
                tmpPosition.scanReadStart = contigOnReadPosition[i-1];
                tmpPosition.scanReadEnd   = contigOnReadPosition[i];
                tmpPosition.kmerCount = 1;
                rePositionVec.push_back(tmpPosition);
            }
            else
            {
                std::vector<ReadScanAndOverlapPosition>::iterator tmpIter = rePositionVec.end();
                --tmpIter;
                (*tmpIter).scanReadEnd = contigOnReadPosition[i];
                (*tmpIter).kmerCount++;
            }
        }
        else minPosition = 0;    
    }
    
    for(std::vector<ReadScanAndOverlapPosition>::iterator tmpIter = rePositionVec.begin() ; tmpIter != rePositionVec.end() ; ++tmpIter )
    {
        if( (*tmpIter).kmerCount > result.kmerCount )
        {
            result.kmerCount = (*tmpIter).kmerCount;
            result.scanReadStart = (*tmpIter).scanReadStart;
            result.scanReadEnd = (*tmpIter).scanReadEnd;
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
        
        //std::cout << "mostConnect : " << (*iter).firstContig << "\t" <<(*iter).secondContig << "\n";
        
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
    //std::cout << "top two most seeds connect : " <<first << "\t" << second << "\t" << result << "\n";
    if(first==second) return "";
    return result;
}
  
int ReassemblerPostProcess::mostConnectRead(OverlapSeedVecInfo input, bool show)
{
    SparseHashMap<int, int, Int64Hasher> readAndCount;
    int mostconnectReadIndex = 0;
    int maxSeedCount = 0;
    
    for( SeedSequenceInfoVec::iterator tmpSeedIter = input.SeedInfoVec.begin() ; tmpSeedIter != input.SeedInfoVec.end() ; ++tmpSeedIter )
    {    
        SparseHashMap<int, int, Int64Hasher>::iterator findReadIndex = readAndCount.find((*tmpSeedIter).readIndex);
        
        if( findReadIndex != readAndCount.end() )
            findReadIndex->second++;
        else
            readAndCount.insert(std::make_pair((*tmpSeedIter).readIndex,1));
    }
    
    for( SparseHashMap<int, int, Int64Hasher>::iterator tmp = readAndCount.begin() ; tmp != readAndCount.end() ; tmp++ )
    {
        if(show)std::cout << tmp->first << "\t" << tmp->second << "\t" << input.SeedInfoVec.size() << "\n";
        
        if( tmp->second > maxSeedCount )
        {
            mostconnectReadIndex = tmp->first;
            maxSeedCount  = tmp->second;
        }
    }
    
    return mostconnectReadIndex;
}

std::string ReassemblerPostProcess::mostReadConnectContig(OverlapSeedVecInfo input)
{
    SparseHashMap<std::string, int, StringHasher> readConnectContig;
    std::string mostReadConnectContigString = "";
    int maxConnectCount = 0;
    
    
    for( SeedSequenceInfoVec::iterator tmpSeedIter = input.SeedInfoVec.begin() ; tmpSeedIter != input.SeedInfoVec.end() ; ++tmpSeedIter )
    {    
        ReadConnectContigHashMap::iterator t = singleContigHashMap.find((*tmpSeedIter).readIndex);
        SparseHashMap<std::string, int, StringHasher>::iterator findReadConnectContig = readConnectContig.find(t->second[0].firstContig);
        
        if( findReadConnectContig != readConnectContig.end() )
            findReadConnectContig->second++;
        else
            readConnectContig.insert(std::make_pair(t->second[0].firstContig,1));
        
        //std::cout << "(debug) " << t->second[0].firstContig << "\n";
    }

    for( SparseHashMap<std::string, int, StringHasher>::iterator tmp = readConnectContig.begin() ; tmp != readConnectContig.end() ; tmp++ )
    {
        if( tmp->second > maxConnectCount )
        {
            mostReadConnectContigString = tmp->first;
            maxConnectCount  = tmp->second;
        }
    }

    return mostReadConnectContigString;
}

bool ReassemblerPostProcess::mostSeedStrand(OverlapSeedVecInfo input)
{
    int forward  = 0;
    int backward = 0;
    
    for( SeedSequenceInfoVec::iterator tmpSeedIter = input.SeedInfoVec.begin() ; tmpSeedIter != input.SeedInfoVec.end() ; ++tmpSeedIter )
    {
        if((*tmpSeedIter).strand) forward++;
        else backward++;
    }
    
    if( forward > backward ) return true;
    return false;
}    

int ReassemblerPostProcess::connectVote(size_t readIdx, std::vector<OverlapRelation> inputRelationVec)
{
    // 0:same number of overlap and connect 
    // 1:overlap case 
    // 2:nonOverlap case
    
    size_t overlap = 0;
    size_t connect = 0;
    size_t readStatus = 0;
    
    for(std::vector<OverlapRelation>::iterator iter = inputRelationVec.begin(); iter != inputRelationVec.end(); ++iter )
    {
        if((*iter).readIndex == (int)readIdx)
        {
            if((*iter).connect) readStatus = 2;
            else readStatus = 1;
        }
        
        if((*iter).connect) connect++;
        else overlap++;
    }

    if( readStatus==1 && connect < overlap ) return 1;
    if( readStatus==2 && connect > overlap ) return 2;
    
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

void ReassemblerPostProcess::filterErrorStrand(OverlapSeedVecInfo input, OverlapSeedVecInfo &result, bool mostStrand)
{
    for( SeedSequenceInfoVec::iterator tmpSeedIter = input.SeedInfoVec.begin() ; tmpSeedIter != input.SeedInfoVec.end() ; ++tmpSeedIter )
    {
        if( (*tmpSeedIter).strand != mostStrand ) continue;
        result.SeedInfoVec.push_back((*tmpSeedIter));
    }
}

void ReassemblerPostProcess::filterErrorSeeds(OverlapSeedVecInfo input, OverlapSeedVecInfo &result, int mostRead)
{
    for( SeedSequenceInfoVec::iterator tmpSeedIter = input.SeedInfoVec.begin() ; tmpSeedIter != input.SeedInfoVec.end() ; ++tmpSeedIter )
    {
        if( (*tmpSeedIter).readIndex != mostRead ) continue;
        result.SeedInfoVec.push_back((*tmpSeedIter));
    }
}

void ReassemblerPostProcess::buildKmerHashAndReadKmerVec( std::string readConnectContig, bool connectContigSide, int connectContigLength, int readIndex, std::vector<ReassemblerSeedFeature> inputSeedVec, KmerHashMap &allKmerHash, ReadKmerVec &readVec)
{
    KmerInfoVec currentReadKmerVec;
    
    for( std::vector<ReassemblerSeedFeature>::iterator iter = inputSeedVec.begin(); iter != inputSeedVec.end(); ++iter )
    {    
        KmerInfo currentKmer;
        
        currentKmer.readConnectContigSide   = connectContigSide;
        currentKmer.readConnectContigID     = readConnectContig;
        currentKmer.readConnectContigLength = connectContigLength;
        currentKmer.readIndex = readIndex;
        currentKmer.kmerStr     = (*iter).seedStr;
        currentKmer.position    = (*iter).seedStartPos;
        
        currentReadKmerVec.push_back(currentKmer);
        
		// build all kmer information
		// so we can know that a kmer appears on those reads
		
        KmerHashMap::iterator findforwardKmerIter = allKmerHash.find((*iter).seedStr);
        
        if( findforwardKmerIter!=allKmerHash.end() )
        {
            currentKmer.kmerStrand = true;
            findforwardKmerIter->second.push_back(currentKmer);
        }
        else
        {
            std::vector<KmerInfo> tmpVec;
            currentKmer.kmerStrand = true;
            tmpVec.push_back(currentKmer);
            allKmerHash.insert(std::make_pair((*iter).seedStr,tmpVec));
        }
        
		currentKmer.kmerStr = reverseComplement((*iter).seedStr);
        KmerHashMap::iterator findreverseKmerIter = allKmerHash.find(currentKmer.kmerStr);
        
        if( findreverseKmerIter!=allKmerHash.end() )
        {
            currentKmer.kmerStrand = false;
            findreverseKmerIter->second.push_back(currentKmer);
        }
        else
        {
            std::vector<KmerInfo> tmpVec;
            currentKmer.kmerStrand = false;
            tmpVec.push_back(currentKmer);
            allKmerHash.insert(std::make_pair(reverseComplement((*iter).seedStr),tmpVec));
        }
    }
    readVec.push_back(std::make_pair(readIndex,currentReadKmerVec));
}

bool ReassemblerPostProcess::checkSameStrand(int startReadIndex, int targetReadIndex, SeedSequenceInfoVec brigdeSeedVec)
{
    bool strand;
    
    ReadConnectContigHashMap::iterator startReadIter  = singleContigHashMap.find(startReadIndex);
    ReadConnectContigHashMap::iterator targetReadIter = singleContigHashMap.find(targetReadIndex);
    
    std::vector<OverlapRelation>::iterator startFirst = startReadIter->second.begin();
    std::vector<OverlapRelation>::iterator targetFirst = targetReadIter->second.begin();


    if( (*startFirst).readDirection == (*targetFirst).readDirection ) strand = true;
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
    ReadConnectContigHashMap::iterator targetReadIter = singleContigHashMap.find(input.SeedInfoVec[0].readIndex);
    std::string secondReadRead = backtrackRead( m_params.indices, input.SeedInfoVec[0].readIndex );    
    std::string fragment;
    
    SeedSequenceInfoVec::iterator first = input.SeedInfoVec.begin();
    SeedSequenceInfoVec::iterator last = input.SeedInfoVec.end();
    last--;
    
    if( targetReadIter->second[0].singleOverlap.scanReadStart == 0 )
    {
        //int fragmentLength = targetReadIter->second[0].singleOverlap.scanReadEnd - overlapLength(input);
        //fragment = secondReadRead.substr( overlapLength(input), fragmentLength );
        int start = (*first).readLocation < (*last).readLocation ? (*first).readLocation : (*last).readLocation;
        int end   = targetReadIter->second[0].singleOverlap.scanReadEnd;
        fragment = secondReadRead.substr( start, end - start );
        
        //std::cout<< "start2\n";
        //std::cout<< (*first).readLocation << "\t" << (*last).readLocation <<"\n";
        //std::cout<< start << "\t" << end << "\t idx " << input.SeedInfoVec[0].readIndex << "\n";
    }
    else
    {
        //int fragmentLength = targetReadIter->second[0].singleOverlap.scanReadEnd - targetReadIter->second[0].singleOverlap.scanReadStart - overlapLength(input);
        //fragment = secondReadRead.substr( targetReadIter->second[0].singleOverlap.scanReadStart, fragmentLength );
        
        int start = targetReadIter->second[0].singleOverlap.scanReadStart;
        int end   = (*first).readLocation > (*last).readLocation ? (*first).readLocation : (*last).readLocation;
        fragment = secondReadRead.substr( start, end - start );
        
        //std::cout<< "end2\n";
        //std::cout<< (*first).readLocation << "\t" << (*last).readLocation <<"\n";
        //std::cout<< start << "\t" << end << "\t idx " << input.SeedInfoVec[0].readIndex << "\n";
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
        size_t kmerFreqsUsingRead   = getFrequency(m_params.indices,kmer);

        if(seedCount_15>=100000)break;
        kmerSizeFreqByRead_15[seedCount_15] = kmerFreqsUsingRead;
		
        seedCount_15++;
    }
    
    slideWindow = 11;
    
    for( int i = 0 ; i < (int)seed.length() - slideWindow + 1 ; i++ )
    {
        std::string kmer = seed.substr(i, slideWindow);
        size_t kmerFreqsUsingRead = getFrequency(m_params.indices,kmer);
		
        if(seedCount_11>=100000)break;
        kmerSizeFreqByRead_11[seedCount_11] = kmerFreqsUsingRead;
		
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

void SeedSequenceInfo::setContigSeq(std::string contigID, std::string contigPartialSeq, int contigLength, int contigStartPos, int contigEndPos, bool contigSide)
{
    this->contigID         = contigID;
    this->contigPartialSeq = contigPartialSeq;
    this->contigLength     = contigLength;
    this->contigStartPos   = contigStartPos;
    this->contigEndPos     = contigEndPos;
    this->contigSide       = contigSide;
    seedLength = contigEndPos - contigStartPos + 1;
}

void SeedSequenceInfo::setReadSeq(std::string seedString, int frequency, int readIndex, int readLocation, bool strand)
{
    this->seedString   = seedString;
    this->frequency    = frequency;
    this->readIndex    = readIndex;
    this->readLocation = readLocation;
    this->strand = strand;
} 

void SeedSequenceInfo::printSeedInfo()
{
    std::cout<< contigID       << "\t" 
             << contigLength   << "\t"
             << contigStartPos << "\t"
             << contigEndPos   << "\t"
             << seedLength     << "\t"
             << (contigSide ? "Head" : "Tail") << "\t"
             << readIndex    << "\t"
             << readLocation << "\t"
             << (strand ? "NonRC" : "RC") << "\t"
             << frequency      << "\t"
             //<< seedString     << "\t"
             << "\n";
}



        ///************************************************///
        ///******* Read Reassemble Basic Elements *********///
        ///************************************************///
        
ReadReassemblerBasicElements::ReadReassemblerBasicElements(const ReassemblerParameters params): m_params(params)
{
}
ReadReassemblerBasicElements::~ReadReassemblerBasicElements()
{
}

std::vector<ReassemblerSeedFeature> ReadReassemblerBasicElements::seedingByDynamicKmer(const std::string readSeq)
{
    return seedingByDynamicKmer( readSeq, 0, readSeq.length() );
}

std::vector<ReassemblerSeedFeature> ReadReassemblerBasicElements::seedingByDynamicKmer(const std::string readSeq, size_t scanStart, size_t scanEnd)
{
    std::vector<ReassemblerSeedFeature> seedVec;
    const size_t smallKmerSize = m_params.kmerLength;
    const size_t largeKmerSize = smallKmerSize/*+6*/;
    
    // set dynamic kmer as largest kmer initially, which will reduce size when no seeds can be found within a maximum interval
    size_t dynamicKmerSize = largeKmerSize;
    size_t kmerThreshold = m_params.seedKmerThreshold;
    
    // prevention of short reads
    if(readSeq.length() < dynamicKmerSize) return seedVec;

    // search for solid kmers and group consecutive solids kmers into one seed
    // reduce kmer size and kmerThreshold if no seeds can be found within maxSeedInterval
    for(size_t i = scanStart ; i+dynamicKmerSize <= scanEnd ; i++)
    {
        std::string kmer = readSeq.substr(i, dynamicKmerSize);
        size_t fwdKmerFreqs = BWTAlgorithms::countSequenceOccurrencesSingleStrand(kmer, m_params.indices);
        size_t rvcKmerFreqs = BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseComplement(kmer), m_params.indices);
        size_t kmerFreqs = fwdKmerFreqs + rvcKmerFreqs;
        
         //std::cout << i << ": " << kmer << "\t" << kmerFreqs <<":" << fwdKmerFreqs << ":" << rvcKmerFreqs << "\n";
        
        if(kmerFreqs >= kmerThreshold && fwdKmerFreqs>=1 && rvcKmerFreqs>=1)
        {    
            // skip low-complexity seeds, e.g., AAAAAAAATAAAA

            if(isLowComplexity(kmer, 0.7)) continue;
            
            size_t seedStartPos = i, 
            seedLen = 0;
            
            // Group consecutive solid kmers into one seed if possible
            size_t maxKmerFreq = kmerFreqs;
            for(i++ ; i+dynamicKmerSize <= scanEnd; i++)
            {
                kmer = readSeq.substr(i, dynamicKmerSize);

                fwdKmerFreqs = BWTAlgorithms::countSequenceOccurrencesSingleStrand(kmer, m_params.indices);
                rvcKmerFreqs = BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseComplement(kmer), m_params.indices);
                kmerFreqs = fwdKmerFreqs + rvcKmerFreqs;
        
                // std::cout << i << ": " << kmer << "\t" << kmerFreqs <<":" << fwdKmerFreqs << ":" << rvcKmerFreqs << "\n";

                if(isLowComplexity(kmer, 0.7)) break;

                maxKmerFreq = std::max(maxKmerFreq,kmerFreqs);

                if( kmerFreqs >= (size_t) kmerThreshold && fwdKmerFreqs>=1 && rvcKmerFreqs>=1) seedLen++;
                else break;
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
                seedVec.push_back(newSeed);
            }
            else
            {
                // push concatenated seeds into seed vector
                ReassemblerSeedFeature newSeed(seedStartPos, readSeq.substr(seedStartPos, seedEndPos-seedStartPos+1), false, dynamicKmerSize, kmerThreshold+10);
                seedVec.push_back(newSeed);
            }
            
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

std::vector<ReassemblerSeedFeature> ReadReassemblerBasicElements::collectKmer(const std::string readSeq, int scanStart, int scanEnd)
{
    std::vector<ReassemblerSeedFeature> seedVec;
    const size_t kmerSize = m_params.kmerLength;

    for(size_t i = scanStart ; i+kmerSize <= (size_t)scanEnd ; i++)
    {
        std::string kmer = readSeq.substr(i, kmerSize);
        size_t KmerFreqsUsingRead = getFrequency(m_params.indices,kmer);
        
        size_t fwdKmerFreqsUsingRead = BWTAlgorithms::countSequenceOccurrencesSingleStrand(kmer, m_params.indices);
        size_t rvcKmerFreqsUsingRead = BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseComplement(kmer), m_params.indices);

        if( fwdKmerFreqsUsingRead == 0 || rvcKmerFreqsUsingRead == 0 )continue;
        
        if(KmerFreqsUsingRead<2)continue;
        
        if(isLowComplexity(kmer, 0.7)) continue;
        
        ReassemblerSeedFeature newSeed(i, readSeq.substr(i, kmerSize), false, kmerSize, 0);
        seedVec.push_back(newSeed);    
    }
    return seedVec;
}

std::pair<size_t, size_t> ReadReassemblerBasicElements::refineRepeatSeed(const std::string readSeq, size_t& seedStartPos, size_t& seedEndPos)
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

bool ReadReassemblerBasicElements::isLowComplexity (std::string seq, float threshold)
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

std::pair<size_t, size_t> ReadReassemblerBasicElements::BacktrackReadIdx(const BWTIndexSet indices, size_t InputIdx, bool activeHash)
{
    size_t backtrackLength = 0;
    size_t idx = InputIdx;
    
    while(1)
    {
        char b = indices.pBWT->getChar(idx);
        size_t new_idx = indices.pBWT->getPC(b) + indices.pBWT->getOcc(b, idx - 1);
        
        if(activeHash)
        {        
            SparseHashMap<size_t, sizetPair, SizetHasher>::iterator findIter = recordBacktrackIndex.find(new_idx);
            if( findIter != recordBacktrackIndex.end() )
            {
                //std::cout <<"insert : " << InputIdx << " find : " << findIter->second.first << " $ is : " << findIter->first << " length : " << backtrackLength << "\n";

                SparseHashMap<size_t, sizetPair, SizetHasher>::iterator insertIter = recordBacktrackIndex.find(InputIdx);
                
                if( insertIter == recordBacktrackIndex.end() )
                {
                    sizetPair tmp = std::make_pair(findIter->first,backtrackLength+findIter->second.second);
                    recordBacktrackIndex.insert(std::make_pair(InputIdx,tmp));
                    return tmp;
                }
                return std::make_pair(findIter->second.first,findIter->second.second);
            }
        }
        
        if(b == '$')
        {
            idx = indices.pSSA->lookupLexoRank(new_idx);
            sizetPair tmp = std::make_pair(idx,backtrackLength);
            
            if(activeHash)
            {
                SparseHashMap<size_t, sizetPair, SizetHasher>::iterator insertIter = recordBacktrackIndex.find(InputIdx);
                
                if( insertIter == recordBacktrackIndex.end() )
                    recordBacktrackIndex.insert(std::make_pair(InputIdx,tmp));
                
                //std::cout <<"insert : " << InputIdx << " $ is : " << idx << " length : " << backtrackLength << "\n";
            }
            return tmp;
        }
        else idx = new_idx;
        
        backtrackLength++;
        if( idx == InputIdx )
        {
            return std::make_pair(0,backtrackLength);
        }        
    }
}

std::string ReadReassemblerBasicElements::backtrackRead(const BWTIndexSet indices, int64_t inputIdx)
{
    std::string backtrackRead;
            
    while(1)
    {
        char b = indices.pRBWT->getChar(inputIdx);
        int64_t new_idx = indices.pRBWT->getPC(b) + indices.pRBWT->getOcc(b, inputIdx - 1);
        
        if(b == '$')break;
        else inputIdx = new_idx;
        
        backtrackRead += b ;
    }
    return backtrackRead;
}

size_t ReadReassemblerBasicElements::getFrequency(const BWTIndexSet indices, std::string query)
{
    size_t fwdKmerFreqs = BWTAlgorithms::countSequenceOccurrencesSingleStrand(query, indices);
    size_t rvcKmerFreqs = BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseComplement(query), indices);
    return fwdKmerFreqs+rvcKmerFreqs;
}

void ReadReassemblerBasicElements::showReadFrequency(size_t queryIndex, int kmerSize)
{
    std::string readSeq  = backtrackRead( m_params.indices, queryIndex );
	
	for( int i = 0 ; i < (int)readSeq.length() - kmerSize + 1 ; i++ )
    {
        std::string kmer = readSeq.substr(i, kmerSize);
        std::cout<< getFrequency(m_params.indices,kmer) << " ";
    }
	std::cout<<"\n";
}

