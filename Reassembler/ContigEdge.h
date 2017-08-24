

#ifndef CONTIGEDGE_H
#define CONTIGEDGE_H

#include "Edge.h"
#include "ContigVertex.h"

class ContigEdge : public Edge
{
	public: 
	ContigEdge(	Vertex* end, EdgeDir dir, EdgeComp comp, SeqCoord m, bool v1_pacbioStrand ,bool v2_pacbioStrand, std::string pacbioFragment, int overlapLength, GraphColor c= GC_WHITE )
    {
		m_firstStrand  = v1_pacbioStrand;
		m_secondStrand  = v2_pacbioStrand;
		m_overlap    = overlapLength;
		m_pbFragment = pacbioFragment;
		
		m_pEnd       = end;
		m_pTwin      = NULL;
		m_matchCoord = m;
		m_color      = c ;
		
		m_edgeData.setDir(dir);
        m_edgeData.setComp(comp);
    }

	
	bool getV1PacbioStrand(){return m_firstStrand;}
	bool getV2PacbioStrand(){return m_secondStrand;}
	
	void reverseV1Strand(){ m_firstStrand = !m_firstStrand;}
	void reverseV2Strand(){ m_secondStrand = !m_secondStrand;}
	
	void setV1Strand(bool inputStrand){ m_firstStrand = inputStrand;}
	void setV2Strand(bool inputStrand){ m_secondStrand = inputStrand;}
	
	std::string m_pbFragment;
	bool m_firstStrand;
	bool m_secondStrand;

	int m_overlap;
	
};

#endif