
#ifndef CONTIGBIGRAPH_H
#define CONTIGBIGRAPH_H

#include "Bigraph.h"
#include "ContigEdge.h"
#include "ContigVertex.h"

class ContigGraph;
typedef bool(*ContigVertexVisitFunction)(ContigGraph*, ContigVertex*);

typedef SparseHashMap<std::string, ContigVertex*, StringHasher> ContigVertexPtrMap;

typedef ContigVertexPtrMap::iterator ContigVertexPtrMapIter;
typedef ContigVertexPtrMap::const_iterator ContigVertexPtrMapConstIter;

class ContigGraph : public Bigraph
{
	public:
		
		ContigGraph():Bigraph()
		{
		};

        void Simplify();
		void pacbioSimplify();
		
		size_t Simplify(Vertex* pV, EdgeDir dir);
		size_t pacbioSimplify(Vertex* pV, EdgeDir dir);
	
		// Merge vertices that are joined by the specified edge
        void mergePacbio(Vertex* pV1, Edge* pEdge);
		void mergeOrigin(Vertex* pV1, Edge* pEdge);

		
};

namespace CGUtil
{

	// Load a string graph from a fasta file.
	// Returns a graph where each sequence in the fasta is a vertex but there are no edges in the graph.
	ContigGraph* loadFASTA(const std::string& filename);

};

#endif