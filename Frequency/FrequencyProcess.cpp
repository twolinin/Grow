///-----------------------------------------------
// 
// Written by Jyun-Hong Lin
// 
//-----------------------------------------------
//
// PacBioLongReadFrequency
//

#include "FrequencyProcess.h"
#include <set>
//#include <boost/function.hpp>
#include <boost/bind.hpp>

        ///************************************************///
        ///*********** Read Frequency Process ************///
        ///************************************************///

	
FrequencyProcess::FrequencyProcess(const FrequencyParameters params) 
: ReadFrequencyBasicElements(params)
{
}
FrequencyProcess::~FrequencyProcess()
{
}

FrequencyResult FrequencyProcess::PBFrequency(const SequenceWorkItem& workItem)
{    
    FrequencyResult result;
	std::string originSequence = workItem.read.seq.toString();
	
	
	static int countNumer = 1;
 
	//std::cout<< "sequence" << countNumer << "\n" ;
	/*
	if( checkChimera(originSequence) )
	{
		std::cout<< "\t is chimera\n";
	}
	else
	{
		std::cout << "\t no chimera\n";
	}
	*/
	countNumer++;
	

	repeatThreshold(15);
	//showReadFrequency(originSequence,15);
	//findRepeatRegion(originSequence,15,6);
	//showFrequency(originSequence);
	
    return result;
}

bool FrequencyProcess::checkChimera(size_t readIndex)
{
	std::string querySeq  = backtrackRead( m_params.indices, readIndex );
	return checkChimera(querySeq);
}

bool FrequencyProcess::checkChimera(std::string querySeq)
{
	return false;
	clock_t p1,p2;
	p1 = clock();
	
	size_t kmerSize = 15;
	bool showStepTime = false;
	size_t prominentLenghtThreshold = 3000;
	float chimeraRatio = 0.1;
	// restore the string
    //std::string querySeq  = backtrackRead( m_params.indices, queryIndex );
	size_t queryIndex = getSeqIdx( m_params.indices, querySeq );
	//std::cout<<queryIndex<<"\n";
    // used to record the coordinates of the query and the corresponding sequence
    SparseHashMap<size_t, std::vector<KmerInfo>, SizetHasher> readKmerHashMap;
    
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
        if( kmerFreqs > 50 ) continue;
		
        for(size_t j = interval.lower ; j <= (size_t)interval.upper ; j++)
        {
            // filter self index
            if(interval.lower + 1 == interval.upper )continue;
            
            std::pair<size_t, size_t> currentSeed = BacktrackReadIdx(m_params.indices,j,true);
            // filter self index
            if( queryIndex == currentSeed.first )continue;
            
            SparseHashMap<size_t, std::vector<KmerInfo>, SizetHasher>::iterator readKmerIter = readKmerHashMap.find(currentSeed.first);
            
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
            
            SparseHashMap<size_t, std::vector<KmerInfo>, SizetHasher>::iterator readKmerIter = readKmerHashMap.find(currentSeed.first);
            
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
    
	p2 = clock();
	if(showStepTime)std::cout<< (p2 - p1) / CLOCKS_PER_SEC <<" sec\t";
	
	
    //std::cout<< queryIndex << "checkChimera\n";
    //std::cout<< "query length   : " << querySeq.length() << "\n";
    // sbjct correspond to the query position, check chimeric read by line cover problem
	std::vector<KmerInfo> querysbjctPosPairVec;
    // if all sbjct fragment can't cover query, check each sbjct prominent length (unaligned length)
	std::vector<ReadScanAndOverlapPosition> sbjctStartEndPosPairVec ;
	
    int pairCount = -1;
    
    // check each overlap read by number of kmer and distance between two kmer
    for(SparseHashMap<size_t, std::vector<KmerInfo>, SizetHasher>::iterator readKmerIter = readKmerHashMap.begin();readKmerIter != readKmerHashMap.end();readKmerIter++)
    {
        // because of repeat and sequencing error, so need filter most error read
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
            
            int queryDis = std::abs( queryPos - (int)(*vecIter).queryPos );
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
        
        tmp3.scanQueryStart = (*beginIter).queryPos;
        tmp3.scanQueryEnd   = (*endIter).queryPos;
        tmp3.scanSbjctStart = (*beginIter).sbjctPos;
        tmp3.scanSbjctEnd   = (*endIter).sbjctPos;
        tmp3.sbjctLength    = sbjct.length();
        
        sbjctStartEndPosPairVec.push_back(tmp3);
		
    }
	
	// because qsort have some problem
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
	
    p2 = clock();
	if(showStepTime)std::cout<< (p2 - p1) / CLOCKS_PER_SEC <<" sec\t";
	
	
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
		assert( nonOverlapLength >= 0 );
        //if( !(overlappingRatio >= 0.9 || nonOverlapLength <= 1000) ) fragment = 2;
		if( overlappingRatio < 0.9 || (size_t)nonOverlapLength > prominentLenghtThreshold ) fragment = 2;
    }
	
	if( fragment != 1 )
	{
		float chimericRead = 0;
		float normalRead = 0;
		bool showDetail = false;
		
		for( std::vector<ReadScanAndOverlapPosition>::iterator iter = sbjctStartEndPosPairVec.begin() ; iter != sbjctStartEndPosPairVec.end() ; iter++ )
		{
			int leftProminentLength  = -1;
			int rightProminentLength = -1;
			
			// left prominent length             |----|
			// query                             ----------------------------
			//                                        ||||||||
			// sbjct                         ----------------------
			// right prominent length                        |-----|
			
			if(showDetail)std::cout<< "qry : " << (*iter).scanQueryStart << "\t" << (*iter).scanQueryEnd << "\t" << querySeq.length() <<"\n";
			if(showDetail)std::cout<< "sbj : " << (*iter).scanSbjctStart << "\t" << (*iter).scanSbjctEnd << "\t" << (*iter).sbjctLength <<"\n";
			
			// sbjct direction  --->
			if( (*iter).scanSbjctStart < (*iter).scanSbjctEnd )
			{
				leftProminentLength  = std::min( (*iter).scanQueryStart , (*iter).scanSbjctStart );
				rightProminentLength = std::min( (int)querySeq.length() - (*iter).scanQueryEnd , (*iter).sbjctLength - (*iter).scanSbjctEnd );
				
				if(showDetail)std::cout<<"-->\n";
				if(showDetail)std::cout<< "left  min ( " << (*iter).scanQueryStart << " , " << (*iter).scanSbjctStart << " ) = " << leftProminentLength << "\n";
				if(showDetail)std::cout<< "right min ( " << (int)querySeq.length() << " - " << (*iter).scanQueryEnd << " , " << (*iter).sbjctLength << " - " << (*iter).scanSbjctEnd <<" ) = " << rightProminentLength << "\n";
			}
			// sbjct direction  <---
			else
			{
				leftProminentLength  = std::min( (*iter).scanQueryStart , (*iter).sbjctLength - (*iter).scanSbjctStart );
				rightProminentLength = std::min( (int)querySeq.length() - (*iter).scanQueryEnd , (*iter).scanSbjctEnd );
			
				if(showDetail)std::cout<<"<--\n";
				if(showDetail)std::cout<< "left  min ( " << (*iter).scanQueryStart << " , " << (*iter).sbjctLength << " - " << (*iter).scanSbjctStart << " ) = " << leftProminentLength << "\n";
				if(showDetail)std::cout<< "right min ( " << (int)querySeq.length() << " - " << (*iter).scanQueryEnd << " , " << (*iter).scanSbjctEnd << " ) = " << rightProminentLength << "\n";
			}
			
			if(showDetail)std::cout<<"-----------------------------------------------\n";
			
			//std::cout<< "ProminentLength " << leftProminentLength  <<"\n";
			//std::cout<< "ProminentLength " << rightProminentLength <<"\n";
			
			assert( leftProminentLength >= 0 && rightProminentLength >= 0);
			
			if( (size_t)leftProminentLength > prominentLenghtThreshold ) chimericRead++;
			else normalRead++;
			
			if( (size_t)rightProminentLength > prominentLenghtThreshold ) chimericRead++;
			else normalRead++;
		}
		
		if( chimericRead + normalRead > 0)
		{
			if(showDetail)
				std::cout<< "chimericReadRatio " << chimericRead << "\t" << normalRead << "\t" << chimericRead / ( chimericRead + normalRead ) <<"\n";
		
			if( chimericRead / ( chimericRead + normalRead ) > chimeraRatio )
			{
				//std::cout<< "chimera : " << queryIndex << "\n";
				return true;
			}
		}
	}

    return false;
}

int FrequencyProcess::repeatThreshold(int kmerSize)
{
	KmerFreqHashMap kmerFreq;
	size_t totalSampleKmerCount = 0;
	size_t tmpCount = 0;
	
	// sample number of read
	for( size_t j = 0 ; j < 1000; j++ )
    {
        // get origin read
		std::string readSeq = backtrackRead( m_params.indices, j );
		// prevent sort read
		if( (int)readSeq.length() < kmerSize ) continue;
		// collect kemr freq and insert to hash
		for( int i = 0 ; i < (int)readSeq.length() - kmerSize + 1 ; i++, totalSampleKmerCount++ )
		{
			std::string kmer = readSeq.substr(i, kmerSize);
			size_t freq = getFrequency(m_params.indices,kmer);
			KmerFreqHashMap::iterator iter = kmerFreq.find(freq);
			
			if( iter != kmerFreq.end() )
			{			
				iter->second = iter->second + 1;
			}
			else
			{
				kmerFreq.insert(std::make_pair(freq,1));
			}
		}
		//std::cout<< std::flush << '\r' << j << "\n";
    }
	//std::cout << kmerFreq.size() << "\n";
	
	sizetVec freqVec;
	// in order to sort the count of frequencies ,so translate hash to vector
	for(KmerFreqHashMap::iterator iter = kmerFreq.begin(); iter != kmerFreq.end(); ++iter )
    {
		sizetPair tmp = std::make_pair(iter->first,iter->second);
		freqVec.push_back(tmp);
	}
	
	std::sort(freqVec.begin(), freqVec.end(), 
			  boost::bind(&std::pair<size_t, size_t>::first, _1) <
			  boost::bind(&std::pair<size_t, size_t>::first, _2));
	
	//for(sizetVec::iterator iter_i = freqVec.begin(); iter_i != freqVec.end() ; iter_i++ )
	//	std::cout<<( *iter_i).first << "\t" << (*iter_i).second << "\n";
	
	std::cout<< "total " << totalSampleKmerCount << "\n";
	std::cout<< "total*0.9 " << (int)(totalSampleKmerCount*0.9) << "\n";
	for(sizetVec::iterator iter_i = freqVec.begin(); iter_i != freqVec.end() ; iter_i++ )
	{
		tmpCount += (*iter_i).second;
		
		if( tmpCount > (int)(totalSampleKmerCount*0.9) )
		{
			std::cout<< "threshold " << (*iter_i).first << "\n";
			return (*iter_i).first;
			break;
		}
	}
	assert(false);
	return -1;
}

void FrequencyProcess::showFrequency(std::string originSequence)
{
	
	static int num = 0;
	
	std::cout << "sequence" << num << "\n";
	
    showReadFrequency(originSequence,13);
	showReadFrequency(originSequence,15);
	showReadFrequency(originSequence,17);
	showReadFrequency(originSequence,19);
	
	num++;
}


        ///************************************************///
        ///********* Read Frequency Post Process *********///
        ///************************************************///

FrequencyPostProcess::FrequencyPostProcess(const FrequencyParameters params) 
: ReadFrequencyBasicElements(params)
{

}
FrequencyPostProcess::~FrequencyPostProcess()
{
}
        
void FrequencyPostProcess::process(const SequenceWorkItem& item, /*const*/ FrequencyResult& result)
{    

}

        ///************************************************///
        ///******* Read Reassemble Basic Elements *********///
        ///************************************************///
        
ReadFrequencyBasicElements::ReadFrequencyBasicElements(const FrequencyParameters params): m_params(params)
{
}
ReadFrequencyBasicElements::~ReadFrequencyBasicElements()
{
}

std::pair<size_t, size_t> ReadFrequencyBasicElements::BacktrackReadIdx(const BWTIndexSet indices, size_t InputIdx, bool activeHash)
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
                size_t passLength = findIter->second.second;
				size_t targetIdx  = findIter->second.first;
				//std::cout <<"insert : " << InputIdx << " find : " << findIter->second.first << " $ is : " << findIter->first << " length : " << backtrackLength << "\n";
                SparseHashMap<size_t, sizetPair, SizetHasher>::iterator insertIter = recordBacktrackIndex.find(InputIdx);
                if( insertIter == recordBacktrackIndex.end() )
                {
                    sizetPair tmp = std::make_pair( targetIdx, passLength + backtrackLength );
                    recordBacktrackIndex.insert(std::make_pair(InputIdx,tmp));
                    return tmp;
                }
                return std::make_pair( targetIdx, passLength + backtrackLength );
            }
        }
        
        if(b == '$')
        {
            idx = indices.pSSA->lookupLexoRank(new_idx);
            sizetPair tmp = std::make_pair(idx,backtrackLength);
            //std::cout <<"insert : " << InputIdx << " $ is : " << idx << " length : " << backtrackLength << "\n";
            if(activeHash)
            {
                SparseHashMap<size_t, sizetPair, SizetHasher>::iterator insertIter = recordBacktrackIndex.find(InputIdx);
                
                if( insertIter == recordBacktrackIndex.end() )
                    recordBacktrackIndex.insert(std::make_pair(InputIdx,tmp));
            }
            return tmp;
        }
        else idx = new_idx;
        
        backtrackLength++;       
    }
}

size_t ReadFrequencyBasicElements::getSeqIdx(const BWTIndexSet indices, std::string query)
{
	size_t fwdKmerFreqs = BWTAlgorithms::countSequenceOccurrencesSingleStrand(query, indices);
    size_t rvcKmerFreqs = BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseComplement(query), indices);
	
	BWTInterval interval  = BWTAlgorithms::findInterval(indices,query) ;
    BWTInterval rinterval = BWTAlgorithms::findInterval(indices,query) ;
	
	if( fwdKmerFreqs + rvcKmerFreqs > 2 || fwdKmerFreqs + rvcKmerFreqs < 1 )
	{
		std::cout<< fwdKmerFreqs << "\t" << rvcKmerFreqs << "\n";
		
		std::string rvquery = reverseComplement(query);
		std::cout<< query.substr(0,30) << "\t" << rvquery.substr(query.length()-30,30) << "\n";
		std::cout<< getFrequency(indices,query) << "\t" << getFrequency(indices,rvquery) << "\n";
		std::cout<< interval.lower << "\t" << interval.upper << "\n";
		std::cout<< rinterval.lower << "\t" << rinterval.upper << "\n";
	
	}
	assert( fwdKmerFreqs + rvcKmerFreqs == 1 || fwdKmerFreqs + rvcKmerFreqs == 2 );
	return interval.lower;
}

std::string ReadFrequencyBasicElements::backtrackRead(const BWTIndexSet indices, int64_t inputIdx)
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

size_t ReadFrequencyBasicElements::getFrequency(const BWTIndexSet indices, std::string query)
{
    size_t fwdKmerFreqs = BWTAlgorithms::countSequenceOccurrencesSingleStrand(query, indices);
    size_t rvcKmerFreqs = BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseComplement(query), indices);
    return fwdKmerFreqs+rvcKmerFreqs;
}

void ReadFrequencyBasicElements::showReadFrequency(size_t queryIndex, int kmerSize)
{
    std::string readSeq  = backtrackRead( m_params.indices, queryIndex );
	
	for( int i = 0 ; i < (int)readSeq.length() - kmerSize + 1 ; i++ )
    {
        std::string kmer = readSeq.substr(i, kmerSize);
        std::cout<< getFrequency(m_params.indices,kmer) << " ";
    }
	std::cout<<"\n";
}

void ReadFrequencyBasicElements::showReadFrequency(std::string sequence, int kmerSize)
{
    for( int i = 0 ; i < (int)sequence.length() - kmerSize + 1 ; i++ )
    {
        std::string kmer = sequence.substr(i, kmerSize);
        std::cout<< getFrequency(m_params.indices,kmer) << " ";
    }
	std::cout<<"\n";
}
/*
void ReadFrequencyBasicElements::findRepeatRegion(std::string sequence, int kmerSize, int repeatThreshold)
{
    int repeatKmer = 0;
	
	for( int i = 0 ; i < (int)sequence.length() - 400 ; i++ )
    {
        std::string kmer = sequence.substr(i, kmerSize);
        //std::cout<< getFrequency(m_params.indices,kmer) << " ";
		if(i<400)
		{
			if( getFrequency(m_params.indices,kmer) > (size_t)repeatThreshold ) repeatKmer++;
			continue;
		}
		
		( getFrequency(m_params.indices,kmer) > (size_t)repeatThreshold ) ? repeatKmer++ : repeatKmer-- ;
		
		if( repeatKmer > 50 )
		{
			std::cout << "start "<< i << "~";
			while( repeatKmer > 50 )
			{
				i++;
				if( i >= (int)sequence.length() - 400 )break;
				kmer = sequence.substr(i, kmerSize);
				( getFrequency(m_params.indices,kmer) > (size_t)repeatThreshold ) ? repeatKmer++ : repeatKmer-- ;
			}
			std::cout << i <<"\n";
		}
		
    }
	//std::cout<<"\n";
}
*/

