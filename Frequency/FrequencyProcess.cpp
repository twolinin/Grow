///-----------------------------------------------
// 
// Written by Jyun-Hong Lin
// 
//-----------------------------------------------
//
// PacBioLongReadFrequency
//

#include "FrequencyProcess.h"


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
    
	static int num = 0;
	// get origin sequence
    std::string originSequence = workItem.read.seq.toString();
    
	std::cout << "sequence" << num++ << "\n";
	
    showReadFrequency(originSequence,13);
	showReadFrequency(originSequence,15);
	showReadFrequency(originSequence,17);
	showReadFrequency(originSequence,19);
	
	
    return result;
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
    /*for(SeedByReadHashMap::iterator currentPbIdx = result.seedHash.begin() ; currentPbIdx != result.seedHash.end() ; ++currentPbIdx)
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
	}*/
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

