

#include "ContigVertex.h"
#include "ContigEdge.h"


void ContigVertex::mergePacbio(Edge* sourceToTargetEdge)
{
	std::string PacBioString = ((ContigEdge*)sourceToTargetEdge)->m_pbFragment;
	
	// Merge the sequence
	// if label->getComp() == EC_REVERSE, this seq will reverseComplement.
	// label is target sequence
    DNAEncodedString targetSequence = sourceToTargetEdge->getLabel();
	
    size_t label_len = targetSequence.length();
	
	sourceToTargetEdge->updateSeqLen(m_seq.length() + label_len + PacBioString.length());
	
	Edge* targetToSourceEdge = sourceToTargetEdge->getTwin();
	/*
    std::cout<<  sourceToTargetEdge->getStart()->getID()  << "\t" 
	         << (sourceToTargetEdge->getDir() ? "head" : "tail") << "\t"
			 <<  targetToSourceEdge->getStart()->getID()  << "\t"
			 << (targetToSourceEdge->getDir() ? "head" : "tail") << "\t"
			 << (sourceToTargetEdge->getComp() ? "same" : "reverse")   << "\t"
			 << (((ContigEdge*)sourceToTargetEdge)->getV1PacbioStrand() ? "forwardPB" : "reversePB") <<"\t"
			 << (((ContigEdge*)sourceToTargetEdge)->getV2PacbioStrand() ? "forwardPB" : "reversePB") <<"\n";
    */
	if(sourceToTargetEdge->getComp() == EC_REVERSE)
	{
		//std::cout<<"reverse pacbio direct\n";
		
		((ContigEdge*)sourceToTargetEdge)->setV2Strand( !((ContigEdge*)sourceToTargetEdge)->getV2PacbioStrand() );
		/*
		std::cout<< sourceToTargetEdge->getStart()->getID()  << "\t" 
		         << (sourceToTargetEdge->getDir() ? "head" : "tail") << "\t"
			     <<  targetToSourceEdge->getStart()->getID()  << "\t"
			     << (targetToSourceEdge->getDir() ? "head" : "tail") << "\t"
			     << (sourceToTargetEdge->getComp() ? "same" : "reverse")   << "\t"
			     << (((ContigEdge*)sourceToTargetEdge)->getV1PacbioStrand() ? "forwardPB" : "reversePB") <<"\t"
			     << (((ContigEdge*)sourceToTargetEdge)->getV2PacbioStrand() ? "forwardPB" : "reversePB") <<"\n";
		*/
	}
	else
	{
		//std::cout<<"keep pacbio direct\n";
	}
	
    if(sourceToTargetEdge->getDir() == ED_SENSE)
    {   //tail , v1 contig append v2 contig
		//std::cout<< "source + ";
		if(!((ContigEdge*)sourceToTargetEdge)->getV1PacbioStrand())
		{   // pacbio antisense mapping v1 contig, reversePB
			//std::cout<<"reversePB + ";
			PacBioString = reverseComplement(PacBioString);
		}
		//else std::cout<<"forwardPB + ";
			
		//std::cout<< "target\n";
		
		m_seq.append(PacBioString);
		m_seq.append(targetSequence);
    }
    else
    {   //head , v2 contig append v1 contig
		//std::cout<< "target + ";
		if(!((ContigEdge*)sourceToTargetEdge)->getV2PacbioStrand())
		{   // pacbio antisense mapping v2 contig, reversePB
		    //std::cout<<"reversePB + ";
			PacBioString = reverseComplement(PacBioString);
		}
		//else std::cout<<"forwardPB + ";
			
		//std::cout<< "source\n";
		
		targetSequence.append(PacBioString);
		targetSequence.append(m_seq);
        std::swap(m_seq, targetSequence);
    }

	//std::cout<<"\n";
	
    // Update the coverage value of the vertex
    m_coverage += sourceToTargetEdge->getEnd()->getCoverage();

#ifdef VALIDATE
    VALIDATION_WARNING("Vertex::merge")
    validate();
#endif

}
