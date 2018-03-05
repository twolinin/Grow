
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
	// check each sequence has links to other sequences
	while(iter != m_vertices.end())
	{
		// check if the tail of this sequence has other sequences
		mergeCount += Simplify(iter->second, ED_SENSE);
		// check if the head of this sequence has other sequences
		mergeCount += Simplify(iter->second, ED_ANTISENSE);
		++iter;
	}

	if(mergeCount>0)
		std::cout << "<Simplify> Merge Vertices : "  <<  mergeCount << std::endl;

}

size_t ContigGraph::Simplify(Vertex* sourceVertex, EdgeDir dir)
{
	size_t mergeCount = 0 ;
	
	// Get the edges for this direction
	EdgePtrVec edges = sourceVertex->getEdges(dir);
	
	while(edges.size() == 1)
	{
		assert(!edges.front()->isSelf());
		Edge* pSingle = edges.front();
		
		Edge* pTwin = pSingle->getTwin();
		Vertex* targetVertex = pSingle->getEnd();
		
		// Check that the twin edge dir of targetVertex is simple as well
		if(targetVertex->countEdges(pTwin->getDir()) == 1)
		{
			//merge targetVertex seq into sourceVertex and move its transitive edges to sourceVertex
			
			//mergePacbio(sourceVertex, pSingle);
			merge(sourceVertex, pSingle);
			
			mergeCount++ ;
			
			//remove self edges produced by V->W->V=V<->V
			edges = sourceVertex->getEdges(dir);
			size_t selfEdgeCount = 0 ;
			for(EdgePtrVecIter edge_iter = edges.begin(); edge_iter != edges.end(); ++edge_iter)
			{
				Edge* pEdge = *edge_iter;
				//self edge found
				if (pEdge->isSelf())
				{
					sourceVertex->deleteEdge(pEdge->getTwin());
					sourceVertex->deleteEdge(pEdge);
					selfEdgeCount++;
				}
			}
			//retrieve edges again if updated by selfedges
			if(selfEdgeCount >0)
				edges = sourceVertex->getEdges(dir);
		}//end of if targetVertex single edge
		else
			break;
	}//end of if sourceVertex single edge

	
	return mergeCount ;
}

void ContigGraph::pacbioSimplify()
{
	assert(!hasContainment());
	size_t mergeCount = 0 ;
	
	VertexPtrMapIter iter = m_vertices.begin();
	// check each sequence has links to other sequences
	while(iter != m_vertices.end())
	{
		// check if the tail of this sequence has other sequences
		mergeCount += pacbioSimplify(iter->second, ED_SENSE);
		// check if the head of this sequence has other sequences
		mergeCount += pacbioSimplify(iter->second, ED_ANTISENSE);
		++iter;
	}
	
	if(mergeCount>0)
		std::cout << "<Simplify> Merge Vertices : "  <<  mergeCount << std::endl;

}

// merge unipaths from sourceVertex to targetVertex in EdgeDir
size_t ContigGraph::pacbioSimplify(Vertex* sourceVertex, EdgeDir dir)
{
	size_t mergeCount = 0 ;
	
	// Get the edges for this direction
	EdgePtrVec sourceEdges = sourceVertex->getEdges(dir);

	while(sourceEdges.size() == 1)
	{
		assert(!sourceEdges.front()->isSelf());
		Edge* sourceToTargetEdge = sourceEdges.front();
		
		Edge* pTwin = sourceToTargetEdge->getTwin();
		Vertex* targetVertex = sourceToTargetEdge->getEnd();
		
		// Check that the twin edge dir of targetVertex is simple as well
		if(targetVertex->countEdges(pTwin->getDir()) == 1)
		{
			//merge targetVertex seq into sourceVertex and move its transitive edges to sourceVertex
			
			mergePacbio(sourceVertex, sourceToTargetEdge);

			mergeCount++ ;
			
			//remove self edges produced by source -> target ->source = source <-> source
			sourceEdges = sourceVertex->getEdges(dir);
			size_t selfEdgeCount = 0 ;
			for(EdgePtrVecIter edge_iter = sourceEdges.begin(); edge_iter != sourceEdges.end(); ++edge_iter)
			{
				Edge* pEdge = *edge_iter;
				//self edge found
				if (pEdge->isSelf())
				{
					sourceVertex->deleteEdge(pEdge->getTwin());
					sourceVertex->deleteEdge(pEdge);
					selfEdgeCount++;
				}
			}
			//retrieve edges again if updated by selfedges
			if(selfEdgeCount >0)
				sourceEdges = sourceVertex->getEdges(dir);
		}//end of if targetVertex single edge
		else
			break;
	}//end of if sourceVertex single edge

	
	return mergeCount ;
} 

void ContigGraph::mergePacbio(Vertex* sourceVertex, Edge* sourceToTargetEdge)
{
	Vertex* targetVertex = sourceToTargetEdge->getEnd();
	//std::cout << "Merging " << sourceVertex->getID() << " with " << targetVertex->getID() << "\n";
	
	// Merge the data
	((ContigVertex*)sourceVertex)->mergePacbio(sourceToTargetEdge);

	// Get the twin edge (the edge in targetVertex that points to sourceVertex)
	Edge* pTwin = sourceToTargetEdge->getTwin();

	// Ensure targetVertex has the twin edge
	assert(targetVertex->hasEdge(pTwin));

	//
	assert(sourceVertex->hasEdge(sourceToTargetEdge));
	size_t transLength = targetVertex->getOriginLength(!pTwin->getDir());
	sourceVertex->setOriginLength(transLength,sourceToTargetEdge->getDir());

	// Get the edge set opposite of the twin edge (which will be the new edges in this direction for sourceVertex)
	EdgePtrVec transEdges = targetVertex->getEdges(!pTwin->getDir());
	
	// Move the edges from targetVertex to sourceVertex
	for(EdgePtrVecIter iter = transEdges.begin(); iter != transEdges.end(); ++iter)
	{
		Edge* pTransEdge = *iter;
		
		//std::cout<< "trans edge\n";
		//std::cout<< pTransEdge->getStart()->getID() << "\t" << pTransEdge->getEnd()->getID() << "\n";

		// pacbio direction
		// Dir  : sense is connect using tail of sequence, anti is using head
		// Comp : same, two sequences are same strand
		if(sourceToTargetEdge->getDir() == ED_SENSE && sourceToTargetEdge->getComp() == EC_REVERSE)
		{   // reverse
			((ContigEdge*)pTransEdge)->setV1Strand( !((ContigEdge*)pTransEdge)->getV1PacbioStrand() );
		}
		else if(sourceToTargetEdge->getDir() == ED_ANTISENSE && sourceToTargetEdge->getComp() == EC_REVERSE)
		{
			// it's need to fix some bug.
			//std::cout<<"debug case.\n";
		}

		
		// Remove the edge from targetVertex, this does not destroy the edge
		targetVertex->removeEdge(pTransEdge);

		// Join sourceToTargetEdge to the start of transEdge
		// This updates the starting point of pTransEdge to be sourceVertex
		// This calls Edge::extend on the twin edge
		pTransEdge->join(sourceToTargetEdge);
		
		/*if(pTransEdge->getDir() != sourceToTargetEdge->getDir())
		{
			std::cout << "Merging " << sourceVertex->getID() << " with " << targetVertex->getID() << "\n";
			std::cout<< pTransEdge->getStart()->getID() << "\t" << pTransEdge->getEnd()->getID() << "\n";
		}*/
		
		// unclear assert meaning 
		//assert(pTransEdge->getDir() == sourceToTargetEdge->getDir());
		
		
		sourceVertex->addEdge(pTransEdge); // add to sourceVertex

		// Notify the edges they have been updated
		pTransEdge->update();
		pTransEdge->getTwin()->update();
	}

	// Remove the edge from sourceVertex to targetVertex
	sourceVertex->removeEdge(sourceToTargetEdge);
	delete sourceToTargetEdge;
	sourceToTargetEdge = 0;

	// Remove the edge from targetVertex to sourceVertex
	targetVertex->removeEdge(pTwin);
	delete pTwin;
	sourceToTargetEdge = 0;

	// Remove V2
	// It is guarenteed to not be connected
	removeIslandVertex(targetVertex);
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
		Vertex* sourceVertexertex = new(pGraph->getVertexAllocator()) Vertex(record.id, record.seq.toString());
		pGraph->addVertex(sourceVertexertex);
	}
	return pGraph;
}
