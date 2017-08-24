
#include "SeqReader.h"
#include "ContigGraph.h"
#include "ContigEdge.h"

//*
// Simplify the graph by compacting singular edges
void ContigGraph::Simplify()
{
	assert(!hasContainment());
	size_t mergeCount = 0 ;
	
	VertexPtrMapIter iter = m_vertices.begin();

	while(iter != m_vertices.end())
	{
		mergeCount += Simplify(iter->second, ED_SENSE);
		mergeCount += Simplify(iter->second, ED_ANTISENSE);
		++iter;
	}

	if(mergeCount>0)
		std::cout << "<Simplify> Merge Vertices : "  <<  mergeCount << std::endl;

}

size_t ContigGraph::Simplify(Vertex* pV, EdgeDir dir)
{
	size_t mergeCount = 0 ;
	
	// Get the edges for this direction
	EdgePtrVec edges = pV->getEdges(dir);
	
	while(edges.size() == 1)
	{
		assert(!edges.front()->isSelf());
		Edge* pSingle = edges.front();
		
		Edge* pTwin = pSingle->getTwin();
		Vertex* pW = pSingle->getEnd();
		
		// Check that the twin edge dir of pW is simple as well
		if(pW->countEdges(pTwin->getDir()) == 1)
		{
			//merge pW seq into pV and move its transitive edges to pV
			
			//mergePacbio(pV, pSingle);
			merge(pV, pSingle);
			
			mergeCount++ ;
			
			//remove self edges produced by V->W->V=V<->V
			edges = pV->getEdges(dir);
			size_t selfEdgeCount = 0 ;
			for(EdgePtrVecIter edge_iter = edges.begin(); edge_iter != edges.end(); ++edge_iter)
			{
				Edge* pEdge = *edge_iter;
				//self edge found
				if (pEdge->isSelf())
				{
					pV->deleteEdge(pEdge->getTwin());
					pV->deleteEdge(pEdge);
					selfEdgeCount++;
				}
			}
			//retrieve edges again if updated by selfedges
			if(selfEdgeCount >0)
				edges = pV->getEdges(dir);
		}//end of if pW single edge
		else
			break;
	}//end of if pV single edge

	
	return mergeCount ;
}

void ContigGraph::pacbioSimplify()
{
	assert(!hasContainment());
	size_t mergeCount = 0 ;
	
	VertexPtrMapIter iter = m_vertices.begin();

	while(iter != m_vertices.end())
	{
		mergeCount += pacbioSimplify(iter->second, ED_SENSE);
		mergeCount += pacbioSimplify(iter->second, ED_ANTISENSE);
		++iter;
	}
	
	if(mergeCount>0)
		std::cout << "<Simplify> Merge Vertices : "  <<  mergeCount << std::endl;

}

// merge unipaths from pV to pW in EdgeDir
size_t ContigGraph::pacbioSimplify(Vertex* pV, EdgeDir dir)
{
	size_t mergeCount = 0 ;
	
	// Get the edges for this direction
	EdgePtrVec edges = pV->getEdges(dir);

	while(edges.size() == 1)
	{
		assert(!edges.front()->isSelf());
		Edge* pSingle = edges.front();
		
		Edge* pTwin = pSingle->getTwin();
		Vertex* pW = pSingle->getEnd();
		
		// Check that the twin edge dir of pW is simple as well
		if(pW->countEdges(pTwin->getDir()) == 1)
		{
			//merge pW seq into pV and move its transitive edges to pV
			
			mergePacbio(pV, pSingle);

			mergeCount++ ;
			
			//remove self edges produced by V->W->V=V<->V
			edges = pV->getEdges(dir);
			size_t selfEdgeCount = 0 ;
			for(EdgePtrVecIter edge_iter = edges.begin(); edge_iter != edges.end(); ++edge_iter)
			{
				Edge* pEdge = *edge_iter;
				//self edge found
				if (pEdge->isSelf())
				{
					pV->deleteEdge(pEdge->getTwin());
					pV->deleteEdge(pEdge);
					selfEdgeCount++;
				}
			}
			//retrieve edges again if updated by selfedges
			if(selfEdgeCount >0)
				edges = pV->getEdges(dir);
		}//end of if pW single edge
		else
			break;
	}//end of if pV single edge

	
	return mergeCount ;
} 

void ContigGraph::mergePacbio(Vertex* pV1, Edge* pEdge)
{
	Vertex* pV2 = pEdge->getEnd();
	//std::cout << "Merging " << pV1->getID() << " with " << pV2->getID() << "\n";
	
	// Merge the data
	((ContigVertex*)pV1)->mergePacbio(pEdge);

	// Get the twin edge (the edge in v2 that points to v1)
	Edge* pTwin = pEdge->getTwin();

	// Ensure v2 has the twin edge
	assert(pV2->hasEdge(pTwin));

	//
	assert(pV1->hasEdge(pEdge));
	size_t transLength = pV2->getOriginLength(!pTwin->getDir());
	pV1->setOriginLength(transLength,pEdge->getDir());

	// Get the edge set opposite of the twin edge (which will be the new edges in this direction for V1)
	EdgePtrVec transEdges = pV2->getEdges(!pTwin->getDir());
	
	// Move the edges from pV2 to pV1
	for(EdgePtrVecIter iter = transEdges.begin(); iter != transEdges.end(); ++iter)
	{
		Edge* pTransEdge = *iter;
		
		//std::cout<< pTransEdge->getStart()->getID() << "\t" << pTransEdge->getEnd()->getID() << "\n";

		// pacbio direction 
		if(pEdge->getDir() == ED_SENSE && pEdge->getComp() == EC_REVERSE)
		{   // reverse
			((ContigEdge*)pTransEdge)->setV1Strand( !((ContigEdge*)pTransEdge)->getV1PacbioStrand() );
		}
		else if(pEdge->getDir() == ED_ANTISENSE && pEdge->getComp() == EC_REVERSE)
		{
			//std::cout<<"debug case.\n";
		}

		
		// Remove the edge from V2, this does not destroy the edge
		pV2->removeEdge(pTransEdge);

		// Join pEdge to the start of transEdge
		// This updates the starting point of pTransEdge to be V1
		// This calls Edge::extend on the twin edge
		pTransEdge->join(pEdge);
		
		assert(pTransEdge->getDir() == pEdge->getDir());
		
		
		pV1->addEdge(pTransEdge); // add to V1

		// Notify the edges they have been updated
		pTransEdge->update();
		pTransEdge->getTwin()->update();
	}

	// Remove the edge from pV1 to pV2
	pV1->removeEdge(pEdge);
	delete pEdge;
	pEdge = 0;

	// Remove the edge from pV2 to pV1
	pV2->removeEdge(pTwin);
	delete pTwin;
	pEdge = 0;

	// Remove V2
	// It is guarenteed to not be connected
	removeIslandVertex(pV2);
	//validate();
}

//*/
// Load a graph (with no edges) from a fasta file
ContigGraph* CGUtil::loadFASTA(const std::string& filename)
{
	ContigGraph* pGraph = new ContigGraph;
	SeqReader reader(filename);
	SeqRecord record;

	while(reader.get(record))
	{
		Vertex* pVertex = new(pGraph->getVertexAllocator()) Vertex(record.id, record.seq.toString());
		pGraph->addVertex(pVertex);
	}
	return pGraph;
}
