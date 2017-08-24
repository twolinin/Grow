

#ifndef CONTIGVERTEX_H
#define CONTIGVERTEX_H

#include "Vertex.h"
#include "ContigEdge.h"


class ContigVertex : public Vertex
{
	public:

		void mergePacbio(Edge* pEdge);

};

#endif