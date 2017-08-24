

#include "ContigVertex.h"
#include "ContigEdge.h"


void ContigVertex::mergePacbio(Edge* pEdge)
{
	std::string PacBioString = ((ContigEdge*)pEdge)->m_pbFragment;
	
	// Merge the sequence
	// if label->getComp() == EC_REVERSE, this seq will reverseComplement.
	// label is v2 sequence
    DNAEncodedString label = pEdge->getLabel();
	
    size_t label_len = label.length();
	
	pEdge->updateSeqLen(m_seq.length() + label_len + PacBioString.length());
	
	/*
    std::cout<< pEdge->getStart()->getID()  << "\t" 
	         << (pEdge->getDir() ? "ED_ANTISENSE" : "ED_SENSE") << "\t"
			 << (pEdge->getComp() ? "EC_SAME" : "EC_REVERSE")   << "\t"
			 << (((ContigEdge*)pEdge)->getV1PacbioStrand() ? "NonRC" : "RC") <<"\t"
			 << (((ContigEdge*)pEdge)->getV2PacbioStrand() ? "NonRC" : "RC") <<"\n";
    */
	if(pEdge->getComp() == EC_REVERSE)
	{
		((ContigEdge*)pEdge)->setV2Strand( !((ContigEdge*)pEdge)->getV2PacbioStrand() );
		/*
		std::cout<< pEdge->getStart()->getID()  << "\t" 
	         << (pEdge->getDir() ? "ANTISENSE" : "SENSE") << "\t\t"
			 << (((ContigEdge*)pEdge)->getV1PacbioStrand() ? "NonRC" : "RC") <<"\t"
			 << (((ContigEdge*)pEdge)->getV2PacbioStrand() ? "NonRC" : "RC") <<"\n";
		*/
	}
	
	
    if(pEdge->getDir() == ED_SENSE)
    {   //tail , v1 contig append v2 contig
		
		if(!((ContigEdge*)pEdge)->getV1PacbioStrand())
		{   // pacbio antisense mapping v1 contig, RC
			//std::cout<<"rvPB1\n";
			PacBioString = reverseComplement(PacBioString);
		}
		
		m_seq.append(PacBioString);
		m_seq.append(label);
    }
    else
    {   //head , v2 contig append v1 contig

		if(!((ContigEdge*)pEdge)->getV2PacbioStrand())
		{   // pacbio antisense mapping v2 contig, RC
		    //std::cout<<"rvPB2\n";
			PacBioString = reverseComplement(PacBioString);
		}
		
		label.append(PacBioString);
		label.append(m_seq);
        std::swap(m_seq, label);
    }

    // Update the coverage value of the vertex
    m_coverage += pEdge->getEnd()->getCoverage();

#ifdef VALIDATE
    VALIDATION_WARNING("Vertex::merge")
    validate();
#endif

}
