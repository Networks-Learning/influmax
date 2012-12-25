#include "stdafx.h"
#include "cascinf.h"

int MetaId = 1e9;

void TDominatorGraph::ComputeDomTree(const TDomMethod& method) {
	PNGraph GroundTruthNoSources = new TNGraph(*GroundTruth);

	// remove sources from the graph before computing the postdominator tree
	for (int i=0; i<Sources.Len(); i++) {
		if (GroundTruthNoSources->IsNode(Sources[i])) {GroundTruthNoSources->DelNode(Sources[i]);};
	}

	// we convert the graph to the format needed by third-party code that compute dominator tree
	// we rather use consecutive id's for the third-party code starting at id=1
	TIntH HIds;
	int p = 1;
	for (TNGraph::TNodeI NI = GroundTruthNoSources->BegNI();
				NI < GroundTruthNoSources->EndNI();
				NI++) {
		HIds.AddDat(NI.GetId()) = p++;
	}

	DominatorGraph DG;

	int *arclist = new int [2*GroundTruthNoSources->GetEdges()];
	p = 0;
	for (TNGraph::TEdgeI EI = GroundTruthNoSources->BegEI();
			EI < GroundTruthNoSources->EndEI();
			EI++) {
		arclist[p++] = HIds.GetKeyId(EI.GetDstNId())+1;
		arclist[p++] = HIds.GetKeyId(EI.GetSrcNId())+1;
	}

	DG.buildGraph (GroundTruthNoSources->GetNodes(),
				   GroundTruthNoSources->GetEdges(),
				   HIds.GetKeyId(Sink[0])+1,
				   arclist,
				   true);

	delete [] arclist;

	int *idom = new int [GroundTruthNoSources->GetNodes()+1];

	switch(method) {
	case SNCA:
		DG.snca(HIds.GetKeyId(Sink[0])+1, idom);
		break;
	case SLT:
		DG.slt(HIds.GetKeyId(Sink[0])+1, idom);
		break;
	case LT:
		DG.lt(HIds.GetKeyId(Sink[0])+1, idom);
		break;
	default:
		break;
	}

	DominatorTree = TNGraph::New();

	// add nodes to the dominator tree
	for (int i=1; i<GroundTruthNoSources->GetNodes()+1; i++) {
		if (idom[i]!=0) { DominatorTree->AddNode(HIds.GetKey(i-1)); }
	}

	// add edges to the dominator tree
	for (int i=1; i<GroundTruthNoSources->GetNodes()+1; i++) {
		if (idom[i]!=0 && i!=idom[i]) { DominatorTree->AddEdge(HIds.GetKey(i-1), HIds.GetKey(idom[i]-1)); }
	}

	// be sure the GroundTruthNoSources memory is released
	GroundTruthNoSources.Clr();

	// remove dominator tree vector
	delete [] idom;

	return;
}

void TDominatorGraph::GetDomGamma(TIntV& Gamma, TIntV& M) {
	// compute Gamma (boundary of the source)
	for (int i=0; i<Sources.Len(); i++) {
		if (!GroundTruth->IsNode(Sources[i])) {continue;}
		TNGraph::TNodeI NI = GroundTruth->GetNI(Sources[i]);
		for (int j=0; j<NI.GetOutDeg(); j++) {
			if (!Sources.IsIn(NI.GetOutNId(j))) { Gamma.AddUnique(NI.GetOutNId(j)); }
		}
	}

	// compute M (boundary of the source with respect to the dominator relation)
	for (int i=0; i<Gamma.Len(); i++) {
		if (!DominatorTree->IsNode(Gamma[i])) continue;
		TIntV Nu; GetNu(Gamma[i], Nu);
		int isc=0; for (int k=0;k<Nu.Len();k++) { if (Gamma.IsIn(Nu[k])) { isc++; } }
		if (isc==1) { M.Add(Gamma[i]); }
	}
}

void TDominatorGraph::GetK(TIntV& K) {
	for (TNGraph::TNodeI NI = GroundTruth->BegNI(); NI<GroundTruth->EndNI(); NI++) {
		if (!DominatorTree->IsNode(NI.GetId()) && !Sources.IsIn(NI.GetId())) {
			K.Add(NI.GetId()); }
	}
}

void TDominatorGraph::SaveDomGraph(const TStr& OutFNm) {
	TFOut FOut(OutFNm);

	// write nodes to file
	for (TNGraph::TNodeI NI = DominatorTree->BegNI(); NI < DominatorTree->EndNI(); NI++) {
		FOut.PutStr(TStr::Fmt("%d\n", NI.GetId())); // nodes
	}

	FOut.PutStr("\n");

	// write edges to file
	for (TNGraph::TEdgeI EI = DominatorTree->BegEI(); EI < DominatorTree->EndEI(); EI++) {
		FOut.PutStr(TStr::Fmt("%d,%d\n", EI.GetSrcNId(), EI.GetDstNId())); // edges
	}
}

double TSources::GetAllNodeProb(const int& Source, THash<TInt, TNodeInfo>& NodeInfoH) {
    double P = 0.0;

    if (Source==-1)
    	return P;

    P += (UpdateInf(Source, NodeInfoH) - CurGain); // marginal gain

	return P;
}

void TSources::SetReachableGroundTruthFromSources(const int& Sink, TIntV &ReachableNodesV) {
	// only consider reachable nodes (forward) from the sources
	TIntH ReachableNodes;
	TIntH CurReachNodes;

	bool verbose = false;

	// remove ingoing edges to sources
	TIntPrV EdgesToRemove;
	for (int i=0; i<NIdV.Len(); i++) {
		TNGraph::TNodeI NI = GroundTruthReachable->GetNI(NIdV[i]);
		for (int j=0; j<NI.GetInDeg(); j++) { 
			EdgesToRemove.Add(TIntPr(NI.GetInNId(j), NI.GetId())); 
			if (verbose) { printf("Removing %d->%d\n", NI.GetInNId(j), NI.GetId()); }	
		}
	}

	// remove ingoing edges to potential source
	if (PNId != -1) {
		TNGraph::TNodeI NI = GroundTruthReachable->GetNI(PNId);
		for (int j=0; j<NI.GetInDeg(); j++) {
			EdgesToRemove.Add(TIntPr(NI.GetInNId(j), NI.GetId()));
			if (verbose) { printf("Removing %d->%d\n", NI.GetInNId(j), NI.GetId()); }
		}
	}

	// remove outgoing edges from sink
	TNGraph::TNodeI NI = GroundTruthReachable->GetNI(Sink);
	for (int j=0; j<NI.GetOutDeg(); j++) {
		EdgesToRemove.Add(TIntPr(Sink, NI.GetOutNId(j)));
		if (verbose) { printf("Removing %d->%d\n", Sink, NI.GetOutNId(j)); }
	}

	for (int j=0; j<EdgesToRemove.Len(); j++) { GroundTruthReachable->DelEdge(EdgesToRemove[j].Val1, EdgesToRemove[j].Val2); }

	for (int i=0; i<NIdV.Len(); i++) {
		ReachableNodes.AddDat(NIdV[i]) = 0;
		if (MxReachable > -1) { TSnap::GetShortPath(GroundTruthReachable, NIdV[i], CurReachNodes, true, MxReachable); }
		else { TSnap::GetShortPath(GroundTruthReachable, NIdV[i], CurReachNodes, true); }

		for (int j=0; j<CurReachNodes.Len(); j++) {
			if (!ReachableNodes.IsKey(CurReachNodes.GetKey(j))) { ReachableNodes.AddDat(CurReachNodes.GetKey(j))  = TInt::Mx; }
			ReachableNodes.GetDat(CurReachNodes.GetKey(j)) = TMath::Mn(ReachableNodes.GetDat(CurReachNodes.GetKey(j)), CurReachNodes[j]);
			if (verbose) { printf("Node %d is reachable from %d\n", CurReachNodes.GetKey(j).Val, NIdV[i].Val); }
		}
	}

	if (PNId != -1) {
		ReachableNodes.AddDat(PNId) = 0;
		if (MxReachable > -1) { TSnap::GetShortPath(GroundTruthReachable, PNId, CurReachNodes, true, MxReachable); }
		else { TSnap::GetShortPath(GroundTruthReachable, PNId, CurReachNodes, true); }
		for (int j=0; j<CurReachNodes.Len(); j++) {
			if (!ReachableNodes.IsKey(CurReachNodes.GetKey(j))) { ReachableNodes.AddDat(CurReachNodes.GetKey(j))  = TInt::Mx; }
			ReachableNodes.GetDat(CurReachNodes.GetKey(j)) = TMath::Mn(ReachableNodes.GetDat(CurReachNodes.GetKey(j)), CurReachNodes[j]);
			if (verbose) { printf("Node %d is reachable from %d\n", CurReachNodes.GetKey(j).Val, PNId.Val); }
		}
	}

	ReachableNodes.GetKeyV(ReachableNodesV);
	if (verbose) { printf("Reachable nodes:%d\n", ReachableNodesV.Len()); }
	GroundTruthReachable = TSnap::GetSubGraph(GroundTruthReachable, ReachableNodesV, false);
}

void TSources::SetReachableGroundTruthFromSink(const int& Sink, TIntV &ReachableNodesV) {
	TIntH ReachableNodes;

	bool verbose = false;

	// only consider nodes that are reachable backwards from Sink
	if (MxReachable > -1) { TSnap::GetShortPathReverse(GroundTruthReachable, Sink, ReachableNodes, true, MxReachable); }
	else { TSnap::GetShortPathReverse(GroundTruthReachable, Sink, ReachableNodes, true); }

	ReachableNodes.GetKeyV(ReachableNodesV);
	if (verbose) { printf("Reachable nodes:%d\n", ReachableNodesV.Len()); }

	GroundTruthReachableCurrent = TSnap::GetSubGraph(GroundTruthReachable, ReachableNodesV, false);// remove outgoing edges from the sink
}

// init average number of infected nodes before T when the source set is empty
void TSources::InitInf() {
    Clr();
}

// update the source set given a new node Source
double TSources::UpdateInf(const int& Source, THash<TInt, TNodeInfo>& NodeInfoH) {
	if (NIdV.IsIn(Source)) { return CurGain; }

	PNId = Source;

	return GetInf(NodeInfoH);
}

// return average number of infected nodes before T (i.e., probability of infection for every node)
double TSources::GetInf(THash<TInt, TNodeInfo>& NodeInfoH) {
    double P = 0;
    double PSink = 0;

    bool verbose = true;

    if (verbose && PNId.Val==-1) {
    	printf("Evaluating source set of size %d...\n", NIdV.Len());
    	NodeProbsH.AddDat(NodeProbsH.Len()) = TIntFltH();
    }

	if (verbose && PNId.Val!=-1) { printf("Evaluating source %d...\n", PNId.Val); }

	// if source set contains one element and doesn't have outgoing links, return 1!
	if ( (Len() == 0) && PNId.Val!=-1 && (GroundTruth->GetNI(PNId).GetOutDeg() == 0) ) {
		P = 1.0;
		printf("....::Sigma_A=%f, Delta_%d=%f::\n", P, PNId.Val, P - CurGain);
		return P;
	}

    // from closer to further in a BFS sense
	TIntV NodesId; GroundTruth->GetNIdV(NodesId);
    for (int i=0; i<NodesId.Len(); i++) {
    	TInt Sink = NodesId[i];
    	if (PNId != -1) {
    		P += GetInf(Sink, NodeInfoH);
    	} else {
    		PSink = GetInf(Sink);
    		P += PSink;
    		NodeProbsH.GetDat(NodeProbsH.Len()-1).AddDat(Sink) = PSink;
    	}
    }

    if (verbose && PNId.Val==-1) { printf("....::Sigma_A=%f::\n", P); }
    if (verbose && PNId.Val!=-1) { printf("....::Sigma_A=%f, Delta_%d=%f::\n", P, PNId.Val, P - CurGain); }

    return P;
}

// compute P(t_n < T)
double TSources::GetInf(const int& Sink, THash<TInt, TNodeInfo>& NodeInfoH) {
	double Prob = 0.0;

	bool verbose = false;

	// if the sink node is a source, trivially infected with prob 1.0
	if (IsNode(Sink) || Sink==PNId) {
		if (!NodeProbsH.IsKey(PNId)) { NodeProbsH.AddDat(PNId) = TIntFltH(); }
		if (!NodeProbsH.GetDat(PNId).IsKey(Sink)) { NodeProbsH.GetDat(PNId).AddDat(Sink) = 0.0; }
		NodeProbsH.GetDat(PNId).GetDat(Sink) = 1.0;

		if (verbose) { printf("%d.", Sink); fflush(stdout); }

		return 1.0;
	}

	TIntV ReachableNodesV;

	GroundTruthReachable = new TNGraph(*GroundTruth);
	SetReachableGroundTruthFromSources(Sink, ReachableNodesV);
	if (!GroundTruthReachable->IsNode(Sink)) { return 0.0; }

	printf("%d.", Sink); fflush(stdout);

	// consider only nodes that are reachable from the sink (backwards)
	GroundTruthReachableCurrent.Clr(); SetReachableGroundTruthFromSink(Sink, ReachableNodesV);

	if (GroundTruthReachableCurrent->IsNode(Sink) && Bfs) { return 1.0; }

	// get DAG
	TPair<PNGraph, TIntPrFltH> DAG(PNGraph::New(), TIntPrFltH());
	GetDAG(Sink, DAG);

	// if the sink node is not connected to sources
	if (DAG.Val2.Len() == 0)
		return 0.0;

	// return probability P(t > T) = [1 0 0 ...]'*ET*[1 1 ... 1]
	Prob = ExpM(DAG, ExpMethod);

	if (verbose) {printf("P(t_%d <= %f) = %f\n", Sink, Window.Val, (1-Prob)); }

	// if we are going to store node info, add infection probability to current infection probs
	if (!NodeProbsH.IsKey(PNId)) { NodeProbsH.AddDat(PNId) = TIntFltH(); }
	if (!NodeProbsH.GetDat(PNId).IsKey(Sink)) { NodeProbsH.GetDat(PNId).AddDat(Sink) = 0.0; }
	NodeProbsH.GetDat(PNId).GetDat(Sink) = (1-Prob);

	// be sure we release the memory of the DAG
	DAG.Val1.Clr();
	DAG.Val2.Clr();

	// return P(t <= T)
	return (1-Prob);
}

double TSources::GetInf(const int& Sink) {
	double Prob = 0.0;

	bool verbose = true;

	// if the sink node is a source, trivially infected with prob 1.0
	if (IsNode(Sink) || Sink == PNId) {
		if (verbose) { printf("%d.", Sink); fflush(stdout); }
		return 1.0;
	}

	TIntV ReachableNodesV;

	GroundTruthReachable = new TNGraph(*GroundTruth);
	SetReachableGroundTruthFromSources(Sink, ReachableNodesV);
	if (!GroundTruthReachable->IsNode(Sink)) { return 0.0; }

	// consider only nodes that are reachable from the sink (backwards)
	GroundTruthReachableCurrent.Clr(); SetReachableGroundTruthFromSink(Sink, ReachableNodesV);

	if (GroundTruthReachableCurrent->IsNode(Sink) && Bfs) { return 1.0; }

	// get DAG
	TPair<PNGraph, TIntPrFltH> DAG(PNGraph::New(), TIntPrFltH());
	GetDAG(Sink, DAG);

	// if the sink node is not connected to sources
	if (DAG.Val2.Len() == 0)
		return 0.0;

	// return probability P(t > T) = [1 0 0 ...]'*ET*[1 1 ... 1]
	Prob = ExpM(DAG, ExpMethod);

	// if (verbose) { printf("P(t_%d <= %f) = %f\n", Sink, Window.Val, (1-Prob)); }

	// be sure we release the memory of the DAG
	DAG.Val1.Clr();
	DAG.Val2.Clr();

	if (verbose) { printf("%d.", Sink); fflush(stdout); }

	// return P(t <= T)
	return (1-Prob);
}

double TSources::GetInf(const int& Sink, TFltPrV& Cdf, const double& MinT, const double& MaxT, const double& Step) {
	bool verbose = true;
	double Prob = 0.0;
	TFltV Times;

	for (double T=MinT; T<=MaxT; T+=Step) { Times.Add(T); }

	// if the sink node is a source, trivially infected with prob 1.0
	if (IsNode(Sink)) {
		for (int i=0; i<Times.Len(); i++) { Cdf.Add(TFltPr(Times[i], 1.0)); }
		return 1.0;
	}

	TIntV ReachableNodesV;

	GroundTruthReachable = new TNGraph(*GroundTruth);
	SetReachableGroundTruthFromSources(Sink, ReachableNodesV);
	if (!GroundTruthReachable->IsNode(Sink)) { return 0.0; }

	// consider only nodes that are reachable from the sink (backwards)
	GroundTruthReachableCurrent.Clr(); SetReachableGroundTruthFromSink(Sink, ReachableNodesV);

	if (GroundTruthReachableCurrent->IsNode(Sink) && Bfs) { return 1.0; }

	// get DAG
	TPair<PNGraph, TIntPrFltH> DAG(PNGraph::New(), TIntPrFltH());
	GetDAG(Sink, DAG);

	// if the sink node is not connected to sources
	if (DAG.Val2.Len() == 0) {
		for (int i=0; i<Times.Len(); i++) { Cdf.Add(TFltPr(Times[i], 0.0)); }
		return 0.0;
	}

	// compute probability P(t <= T) = [1 0 0 ...]'*ET*[1 1 ... 1] for T \in (MinT, MaxT) with Step
	for (int i=0; i<Times.Len(); i++) {
		SetWindow(Times[i]);
		Prob = 1-ExpM(DAG, ExpMethod);
		Cdf.Add(TFltPr(Times[i], Prob));
	}

	return Cdf[Cdf.Len()-1].Val2;
}

void TSources::GetDAG(const int& Sink, TPair<PNGraph, TIntPrFltH> &DAG) {
	TIntV S, T;

	bool verbose = false;

	// add all sources & potential source if any that can reach the Sink
	for (int i=0; i<Len(); i++) { if (GroundTruthReachableCurrent->IsNode(GetNode(i))) { S.Add(GetNode(i)); } }
	if ( (PNId != -1) && (!IsNode(PNId)) ) { if (GroundTruthReachableCurrent->IsNode(PNId)) { S.Add(PNId); } }

	// add sink node
	T.Add(Sink);

	// get metanodes (sets)
	TVec<TIntFltH> Sets, Connectivity;
	if (verbose) {printf("Computing minimal cuts...\n");}
	GetSets(S, T, Sets, 0);

	// we could discard sink nodes not connected to sources in a more efficient way
	if (Sets.Len() == 0) { return; }

	// add full set
	Sets.Add(TIntFltH());
	for (int i=0; i<Sets.Len()-1; i++) {
		for (int j=0;j<Sets[i].Len();j++) {
			if (!Sets[Sets.Len()-1].IsKey(Sets[i].GetKey(j))) { Sets[Sets.Len()-1].AddDat(Sets[i].GetKey(j)) = 0.0; }
		}
	}

	Sets[Sets.Len()-1].AddDat(Sink) = 0.0;

	if (verbose) {
		for (int i=0; i<Sets.Len(); i++) {
			printf("Set %d: ", i);
			for (int j=0; j<Sets[i].Len(); j++) { printf("%d ", Sets[i].GetKey(j).Val); }
			printf("\n");
		}
	}

	// compute connectivity of Sets
	if (verbose) {printf("Computing connectivity between %d minimal cuts...\n", Sets.Len());}
	for (int i=0; i<Sets.Len()-1; i++) { Connectivity.Add(TIntFltH());}
	ComputeEdges(Sink, Sets, Connectivity);

	// add metanodes to the dag
	for (int i=0; i<Sets.Len(); i++) { DAG.Val1->AddNode(i); }

	// add edges to build the DAG of metanodes
	for (int i=0; i<Connectivity.Len(); i++) {
		for (TIntFltH::TIter NI = Connectivity[i].BegI(); NI<Connectivity[i].EndI(); NI++) {
			DAG.Val1->AddEdge(i, NI.GetKey());
			DAG.Val2.AddDat(TIntPr(i, NI.GetKey())) = NI.GetDat();
		}
	}

	if (verbose) {
		for (TIntPrFltH::TIter EI = DAG.Val2.BegI(); EI<DAG.Val2.EndI(); EI++) {
			printf("%d->%d:%f\n", EI.GetKey().Val1.Val, EI.GetKey().Val2.Val, EI.GetDat().Val);
		}

		printf("Total number of edges:%d\n", DAG.Val1->GetEdges());
	}

	return;
}

void TSources::ComputeEdges(const int& Sink, TVec<TIntFltH> &Sets, TVec<TIntFltH> &Connectivity) {
	// we want only nodes in the set of final disable nodes
	TIntV SupersetNotInSubset;
	for (int s=1; s<Sets.Len(); s++) {
		for (int i=0; i<s; i++) {
			// if Sets[s] < Sets[i], there is no edge
			if (Sets[i].Len() >= Sets[s].Len()) { continue; }

			// if Sets[i] \notin Sets[s], there is no edge
			// going backwards will be much faster, if there is a single element different, it will be at the end
			bool NotValid=false; for (int j=Sets[i].Len()-1;j>=0;j--) { if (!Sets[s].IsKey(Sets[i].GetKey(j))) {NotValid=true; break;} };
			if (NotValid) { continue; }

			// if Sets[s] = S(Sets[i] U {v}) for some v \in Sets[s] \ Sets[i], there is edge
			SupersetNotInSubset.Clr();
			int nt = 0;
			for (int j=Sets[s].Len()-1; j>=0; j--) {
				// if ( (j+nt+1) == Sets[i].Len() ) { break; } // don't continue, all the remaining elements belong to the subset
				if (Sets[i].IsKey(Sets[s].GetKey(j))) { nt++; continue; }
				SupersetNotInSubset.Add(Sets[s].GetKey(j));
			}

			for (int j=0;j<SupersetNotInSubset.Len(); j++) {
				int nn = 0;
				// check that {v} has an incoming link from Sets[i]
				bool NotReachable = true;
				for (int n=0; n<Sets[i].Len(); n++) {
					if (GroundTruthReachableCurrent->IsEdge(Sets[i].GetKey(n), SupersetNotInSubset[j]) &&
							!SupersetNotInSubset.IsIn(Sets[i].GetKey(n)))
					{ NotReachable = false; break; }
				}

				if (NotReachable) { continue; }

				// check if adding {v} dominates only the elements in SupersetNotInSubset
				TNGraph::TNodeI NI = GroundTruthReachableCurrent->GetNI(SupersetNotInSubset[j]);

				if (SupersetNotInSubset[j]!=Sink) {
					GroundTruthMetanode = new TNGraph(*GroundTruthReachableCurrent);
					nn = CheckBlocks(Sink, Sets[i], SupersetNotInSubset[j], SupersetNotInSubset);
					GroundTruthMetanode.Clr();
				}

				// build meta-node made up of {v} U Sets[i] and see if it dominates all elements in Sets[s] \ Sets[i]
				if (nn==SupersetNotInSubset.Len() || SupersetNotInSubset[j]==Sink) {
					for (int no=0; no<Sets[i].Len(); no++) {
						if (Alphas.IsKey(TIntPr(Sets[i].GetKey(no), SupersetNotInSubset[j])) &&
								!SupersetNotInSubset.IsIn(Sets[i].GetKey(no))) {
							if (!Connectivity[i].IsKey(s)) { Connectivity[i].AddDat(s) = 0.0; }
							Connectivity[i].GetDat(s) = Connectivity[i].GetDat(s) + Alphas.GetDat(TIntPr(Sets[i].GetKey(no), SupersetNotInSubset[j]));
						}
					}
					break;
				}
			}
		}
	}
}

int TSources::CheckBlocks(const int& Sink, TIntFltH &Set, const int &v, TIntV& SupersetNotInSubset) {
	int nn = 1;

	if (!GroundTruthMetanode->IsNode(MetaId)) {
		GroundTruthMetanode->AddNode(MetaId);

		for (int i=0; i<Set.Len(); i++) {
			if (!GroundTruthMetanode->IsNode(Set.GetKey(i))) { continue; }

			TNGraph::TNodeI NI = GroundTruthMetanode->GetNI(Set.GetKey(i));

			for (int j=0; j<NI.GetOutDeg(); j++) {
				if (Set.IsKey(NI.GetOutNId(j))) { continue; }
				GroundTruthMetanode->AddEdge(MetaId, NI.GetOutNId(j));
			}

			for (int j=0; j<NI.GetInDeg(); j++) {
				if (Set.IsKey(NI.GetInNId(j))) { continue; }
				GroundTruthMetanode->AddEdge(NI.GetInNId(j), MetaId);
			}
		}
	}

	PNGraph GroundTruthMetanodeCurrent = new TNGraph(*GroundTruthMetanode);

	TNGraph::TNodeI NIv = GroundTruthMetanodeCurrent->GetNI(v);
	for (int j=0; j<NIv.GetOutDeg(); j++) {
		if (Set.IsKey(NIv.GetOutNId(j))) { continue; }
		GroundTruthMetanodeCurrent->AddEdge(MetaId, NIv.GetOutNId(j));
	}

	for (int j=0; j<NIv.GetInDeg(); j++) {
		if (Set.IsKey(NIv.GetInNId(j))) { continue; }
		GroundTruthMetanodeCurrent->AddEdge(NIv.GetInNId(j), MetaId);
	}

	for (int i=0; i<Set.Len(); i++)	{ if (GroundTruthMetanodeCurrent->IsNode(Set.GetKey(i))) { GroundTruthMetanodeCurrent->DelNode(Set.GetKey(i)); } }
	GroundTruthMetanodeCurrent->DelNode(v);

	TDominatorGraph DG;
	DG.AddGroundTruth(GroundTruthMetanodeCurrent);
	TIntV Sk; Sk.Add(Sink); DG.SetSink(Sk);
	DG.ComputeDomTree(SNCA);

	PNGraph BTree = TSnap::GetBfsTree(DG.DominatorTree, MetaId, false, true);

	TIntV DomNI; BTree->GetNIdV(DomNI);

	// for (int i=0; i<DomNI.Len(); i++) { printf("%d ", DomNI[i].Val); }
	// printf("\n");

	if (BTree->GetNodes()!=SupersetNotInSubset.Len()) {
		BTree.Clr();
		GroundTruthMetanodeCurrent.Clr();
		DG.Clr();
		return nn;
	}

	BTree.Clr(); // release BTree memory

	for (int i=0; i<DomNI.Len(); i++) { if (SupersetNotInSubset.IsIn(DomNI[i])) { nn++; } else if (DomNI[i]!=MetaId) { break; } }

	// be sure we release the memory
	GroundTruthMetanodeCurrent.Clr();
	DG.Clr();

	return nn;
}

// compute metanodes (minimal cuts) using Provan&Shier (Paradigm for listing (s,t)-cuts in graphs)
void TSources::GetSets(TIntV& S, TIntV& T, TVec<TIntFltH> &Sets, int NumInfected) {
	TPair<TInt, TIntV> NodePivot(TInt(-1), TIntV());

	// if we limit the number of hops for an infection to take place (it is approx. but major speed-up),
	// check length of S and stop if necessary
	if (MxLengthPath!=-1 && NumInfected>MxLengthPath) { return; }

	GetPivot(S, T, NodePivot, Sets.Len());

	TIntV Sp(S);

	if (NodePivot.Val2.Len() == 0) {
		// add set to the collection of metanodes (nodes per set)
		TIntFltH Set;
		for (int i=0; i<Sp.Len(); i++) { Set.AddDat(Sp[i]) = 0.0; }
		if (Set.Len() > 0) { Sets.Add(Set); }

		return;
	}

	TIntV Tp(T);
	Tp.Add(NodePivot.Val1);

	for (int i=0; i<NodePivot.Val2.Len(); i++) {
		if (!Sp.IsIn(NodePivot.Val2[i])) {Sp.Add(NodePivot.Val2[i]);}
	}

	GetSets(S, Tp, Sets, NumInfected);
	GetSets(Sp, T, Sets, ++NumInfected);
}

void TSources::GetPivot(TIntV& S, TIntV& T, TPair<TInt, TIntV>& NodePivot, const int& NumSets) {
	TIntV Gamma, M;
	TDominatorGraph DGraph(S, T);

	// add ground truth
	DGraph.AddGroundTruth(GroundTruthReachableCurrent);

	// compute dominator tree using SNCA
	DGraph.ComputeDomTree(SNCA);

	// compute minimal elements of gamma with respect to dominator tree
	DGraph.GetDomGamma(Gamma, M);

	// try to find at least one pivot
	for (int i=0; i<M.Len(); i++) {
		TIntV Nu, Iv, K;

		DGraph.GetNu(M[i], Nu);
		DGraph.GetK(K);

		PNGraph IGraph = TSnap::GetSubGraph(DGraph.GroundTruth, Nu);

		bool NoValid = false;
		for (int j=0; j<Gamma.Len(); j++) {
			if (!IGraph->IsNode(Gamma[j]))
				continue;

			TIntV Ivv;
			PNGraph BTree = TSnap::GetBfsTree(IGraph, Gamma[j], false, true);
			BTree->GetNIdV(Ivv);
			BTree.Clr(); // release BTree memory

			for (int k=0; k<Ivv.Len(); k++) {
				Iv.AddUnique(Ivv[k]);
				if (T.IsIn(Ivv[k]) && !K.IsIn(Ivv[k])) { NoValid=true; break; }
			}
			if (NoValid) {break;}
		}

		if (NoValid) {continue;}

		Nu.AddV(K);
		IGraph = TSnap::GetSubGraph(DGraph.GroundTruth, Nu);

		for (int j=0; j<Gamma.Len(); j++) {
			if (!IGraph->IsNode(Gamma[j]))
				continue;

			TIntV Ivv;
			PNGraph BTree = TSnap::GetBfsTree(IGraph, Gamma[j], false, true);
			BTree->GetNIdV(Ivv);
			BTree.Clr(); // release BTree memory

			for (int k=0; k<Ivv.Len(); k++) { Iv.AddUnique(Ivv[k]); }
		}

		NodePivot.Val1 = M[i];
		NodePivot.Val2.AddV(Iv);

		DGraph.Clr();

		return;
	}

	DGraph.Clr();
}

double TSources::ExpM(TPair<PNGraph, TIntPrFltH>& DAG, const TExpMethod& ExpMethod) {
	TVec<TIntFltH> A, ET;
	double Prob = 0.0, MinA = TFlt::Mx;
	TFltV ET1stRow;

	bool verbose = false;

	// we do not need the sink node
	for (int i=0; i<DAG.Val1->GetNodes()-1; i++) { A.Add(TIntFltH()); A[i].AddDat(i) = 0.0; ET.Add(TIntFltH()); }

	for (TIntPrFltH::TIter EI = DAG.Val2.BegI(); EI < DAG.Val2.EndI(); EI++) {
		if (EI.GetKey().Val1==DAG.Val1->GetNodes()-1) continue; // in principle, it cannot happen

		if (EI.GetKey().Val2!=DAG.Val1->GetNodes()-1) { A[EI.GetKey().Val1].AddDat(EI.GetKey().Val2) = EI.GetDat(); }

		A[EI.GetKey().Val1].GetDat(EI.GetKey().Val1) = A[EI.GetKey().Val1].GetDat(EI.GetKey().Val1) - EI.GetDat();
	}

	int nzmax = 0;
	for (int i=0; i<DAG.Val1->GetNodes()-1; i++) { if (A[i].GetDat(i) <= MinA) { MinA = A[i].GetDat(i); } nzmax += (A[i].Len()); }
	if (verbose) { printf("q:%f\n", -MinA); }

	// Parlett method to compute exponential matrix of upper triangular (DAG) matrix
	// we need to compute the full upper triangular matrix -- we need only the first row, but all others
	// are used to compute the 1st (backtracking)!
	if (ExpMethod==PARLETT) {
			if (verbose) {printf("EXPM:PARLETT\n");}
			for (int i=0; i<A.Len(); i++) {
				// avoid having two terms in the diagonal equal by summing up small quantity
				ET[i].AddDat(i) = exp(A[i].GetDat(i).Val*Window.Val);
			}

			for (int i=A.Len()-2; i>=0; i--) {
				for (int j=i+1; j<A.Len(); j++) {
					double aux = 0.0;
					for (int k=0; k<(j-i); k++) {
						if (ET[i].IsKey(i+k) && A[i+k].IsKey(j)) {
							aux = aux + (ET[i].GetDat(i+k).Val * A[i+k].GetDat(j).Val); }
						if (ET[j-k].IsKey(j) && A[i].IsKey(j-k)) {
							aux = aux - (ET[j-k].GetDat(j).Val * A[i].GetDat(j-k).Val); }
					}

					if (fabs(aux) > 1e-9 && fabs(A[i].GetDat(i).Val-A[j].GetDat(j).Val)>1e-9)
						ET[i].AddDat(j) = aux/(A[i].GetDat(i).Val-A[j].GetDat(j).Val);
					else {
						ET[i].AddDat(j) = 0.0;
					}
				}
			}

			ET[0].GetDatV(ET1stRow);
			for (int i=0; i<ET1stRow.Len(); i++) { Prob = Prob + ET1stRow[i]; }

	} else if (ExpMethod==KRYLOV || (ExpMethod==MIXED && A.Len()>=500)) {
		if (nzmax==1) {
			Prob = exp(Window.Val*A[0][0]);
			if (verbose) { printf("1x1 matrix\n"); }
			return Prob; }

		if (verbose) {printf("EXPM:KRYLOV\n");}
		// variables for Krylov using expokit
		int nmax = A.Len();
		int n = nmax;
		int m = TInt::GetMn(30, n-1);
		int nz = nzmax;
		int mmax = m; // maximum size for the krylov basis

		int lwsp = (int)(nmax*(mmax+2)+5*pow(mmax+2,2)+7);
		int liwsp = TInt::GetMx(nmax, m+2);

		double *H = new double[nzmax]; // matrix
		int *ia = new int[nzmax]; // for indexing
		int *ja = new int[nzmax]; // for indexing
		double *v = new double[nmax];
		double *w = new double[nmax];
		double *wsp = new double[lwsp];
		int *iwsp = new int[liwsp];

		int iflag = 0;
		int itrace = 0;
		double tol = 0;
		double anorm = 0;
		double t = Window.Val;

		int aa = 0;
		for (int i=0; i<A.Len(); i++) {
			for (int j=0; j<A[i].Len(); j++, aa++) {
				ia[aa] = i+1; // row indices
				ja[aa] = A[i].GetKey(j).Val+1; // col indices
				H[aa] = A[i][j].Val;
			}
		}

		for (int i=0; i<n; i++) { wsp[i] = 0.0; v[i] = 1.0; }
		for (int i=0; i<nz; i++) { wsp[ia[i]-1] = wsp[ia[i]-1] + fabs( H[i] ); };
		anorm = wsp[0];
		for (int i=1; i<n; i++) { if (anorm < wsp[i]) { anorm = wsp[i]; } }

		if (verbose) {
			printf("%d < %d\n", lwsp, (int)(n*(m+2)+5*pow(m+2,2)+7) );
			printf("%d < %d\n", liwsp, m+2);
			printf("%d >= %d\n", m, n);
		}

		// call fortran routine (expokit library)
		dgexpv_ ( &n, &nz, &m, &t, v, w, &tol, &anorm, wsp, &lwsp, iwsp, &liwsp, &itrace, H, ia, ja, &iflag);

		if (iflag != 0) { printf("WARNING: There is a problem computing a matrix exponential: %d!\n", iflag); }

		Prob = w[0];

		delete [] H;
		delete [] ia; // for indexing
		delete [] ja; // for indexing
		delete [] iwsp;
		delete [] v;
		delete [] w;
		delete [] wsp;
	} else if (ExpMethod==HESSENBERG || (ExpMethod==MIXED && A.Len()<500)) {
		if (verbose) {printf("EXPM:HESSENBERG\n");}
        int m = A.Len();
        double t = Window.Val;
        double *col = new double[A.Len()];
        double *wsp = new double[2*m*(m+2)];

        // allocation
        double *H = new double[m*m];

		// fortran organize matrix as (j,i) instead of (i,j)
		for (int i=0; i<A.Len(); i++) {
			col[i] = 1.0;
			for (int j=0; j<A.Len(); j++) {
				if (A[i].IsKey(j)) { H[j*A.Len()+i] = A[i].GetDat(j).Val; }
				else { H[j*A.Len()+i] = 0.0; }
			}
		}

		// call fortran routine (expokit library)
		dnchbv_( &m, &t, H, &m, col, wsp );
		Prob = col[0];

		delete [] col;
	    delete [] wsp;
	    delete [] H;
	}

	return Prob;
}

// load ground truth graph from txt file
void TMaxInfBs::LoadGroundTruthTxt(TSIn& SIn) {
	GroundTruth = TNGraph::New(); TStr Line;

	// add nodes
	SIn.GetNextLn(Line);
	while (!SIn.Eof() && Line != "") {
		TStrV NIdV; Line.SplitOnAllCh(',', NIdV);
		GroundTruth->AddNode(NIdV[0].GetInt());
		if (NIdV.Len()>2) {
			Positions.AddDat(NIdV[0].GetInt()) = TFltPr(NIdV[2].GetFlt(), NIdV[3].GetFlt());
		}
		SIn.GetNextLn(Line); }

	// add edges
	while (!SIn.Eof()) {
		SIn.GetNextLn(Line);
		TStrV NIdV; Line.SplitOnAllCh(',', NIdV);
		// printf("%f %f %f\n", NIdV[0].GetFlt(), NIdV[1].GetFlt(), NIdV[2].GetFlt());
		GroundTruth->AddEdge((int)NIdV[0].GetFlt(), (int)NIdV[1].GetFlt());
		if (NIdV.Len()>2) {
			Alphas.AddDat(TIntPr((int)NIdV[0].GetFlt(), (int)NIdV[1].GetFlt())) = NIdV[2].GetFlt();
		}
	}

	printf("groundtruth nodes:%d edges:%d\n", GroundTruth->GetNodes(), GroundTruth->GetEdges());
}

void TMaxInfBs::GenerateGroundTruth(const int& TNetwork, const int& NNodes, const int& NEdges, const TStr& NetworkParams) {
	  TKronMtx SeedMtx;
	  TStr MtxNm;

	  switch (TNetwork) {
	  // 2-dimension kronecker network
	  case 0:
		  printf("Kronecker graph for Ground Truth\n");
		  SeedMtx = TKronMtx::GetMtx(NetworkParams.CStr()); // 0.5,0.5,0.5,0.5

		  printf("\n*** Seed matrix:\n");
		  SeedMtx.Dump();

		  GroundTruth = TKronMtx::GenFastKronecker(SeedMtx, (int)TMath::Log2(NNodes+1), NEdges, true, 0);

		  break;

	  // forest fire network
	  case 1:
		  printf("Forest Fire graph for Ground Truth\n");
		  TStrV NetworkParamsV; NetworkParams.SplitOnAllCh(';', NetworkParamsV);

		  TFfGGen FF(true, // BurnExpFireP
					 NetworkParamsV[0].GetInt(), // StartNNodes (1)
					 NetworkParamsV[1].GetFlt(), // ForwBurnProb (0.2)
					 NetworkParamsV[2].GetFlt(), // BackBurnProb (0.17)
					 NetworkParamsV[3].GetInt(), // DecayProb (1)
					 NetworkParamsV[4].GetInt(), // Take2AmbasPrb (0)
					 NetworkParamsV[5].GetInt()); // OrphanPrb (0)

		  FF.GenGraph(NNodes, false);
		  GroundTruth = FF.GetGraph();

		  break;
	  }
}

void TMaxInfBs::GenerateAlphas(const int& DAlphas, const TStr& RAlphas) {
	int NumDistributions = 1, CurrentDistribution = 0;
	TStrV RAlphasV; RAlphas.SplitOnAllCh(';', RAlphasV);

	if (GroundTruth->GetNodes() == 0) {
			printf("Ground truth must be generated before running GenerateAlphas!\n");
			return;
	}

	// if several distributions of a kind, assign to one uniformly at random
	if (DAlphas==0 || DAlphas==1) { NumDistributions = RAlphasV.Len()/2; }
	else if (DAlphas==2) { NumDistributions = RAlphasV.Len(); }
	else { FailR("Unknown alpha distribution!"); }

	for (TNGraph::TEdgeI EI = GroundTruth->BegEI(); EI < GroundTruth->EndEI(); EI++) {
		if (NumDistributions > 1) { CurrentDistribution = TFlt::Rnd.GetUniDevInt(0, NumDistributions-1); }

		if (DAlphas==0) { // uniform (RAlphasV[0] is mean, RAlphasV[1] is range)
			IAssert(RAlphasV[CurrentDistribution*2].GetFlt()>RAlphasV[CurrentDistribution*2+1].GetFlt()/2);

			Alphas.AddDat(TIntPr(EI.GetSrcNId(), EI.GetDstNId())) =
						RAlphasV[CurrentDistribution*2].GetFlt() +
						(TFlt::Rnd.GetUniDev()-0.5) * RAlphasV[CurrentDistribution*2+1].GetFlt();

		} else if (DAlphas==1) { // gaussian (RAlphasV[0] is mean, RAlphasV[1] is std dev, we limit values to be over 0)
			Alphas.AddDat(TIntPr(EI.GetSrcNId(), EI.GetDstNId())) =
						  TFlt::Rnd.GetNrmDev(RAlphasV[CurrentDistribution*2].GetFlt(),
								  	  	  	  RAlphasV[CurrentDistribution*2+1].GetFlt(),
								  	  	  	  0,
								  	  	  	  TFlt::Mx);
		} else { // rayleigh (RAlphasV[0] is sigma)
			Alphas.AddDat(TIntPr(EI.GetSrcNId(), EI.GetDstNId())) =
						  TFlt::Rnd.GetRayleigh(RAlphasV[CurrentDistribution].GetFlt());
		}
	}
}


void TMaxInfBs::Init() {
    // reset vector
	Sources.CurGain = 0.0;
	Sources.SetGroundTruth(GroundTruth);
	Sources.SetAlphas(Alphas);
	Sources.GroundTruthReachable = TNGraph::New();
	NodeGainV.Clr();

    // add nodes
    for (TNGraph::TNodeI NI = GroundTruth->BegNI(); NI < GroundTruth->EndNI(); NI++)
    	NodeGainV.Add(TFltIntPr(TFlt::Mx, NI.GetId()));
}

int TMaxInfBs::GetBestNode(double& LastGain, bool& msort, int &attempts) {
	TInt BestN;
	TVec<TInt> KeysV;
	TFltIntPrV NodeGainCopyToSortV;
	TIntV NodeZero;
	double BestGain = TFlt::Mn;
	int BestGainIndex = -1;

	bool verbose = false;

    if (msort) {
    	for (int i=0; i<TMath::Mn(attempts-1, NodeGainV.Len()); i++)
    	    NodeGainCopyToSortV.Add(NodeGainV[i]);

    	if (verbose) {
    		printf("Sorting sublist of size %d of marginal gains!\n", NodeGainCopyToSortV.Len()); }

    	// sort this list
    	NodeGainCopyToSortV.Sort(false);

    	if (verbose) {
    		printf("Sublist sorted!\n"); }

    	// clever way of resorting without need to copy (google interview question! :-))
    	for (int i=0, ii=0, j=0; ii < NodeGainCopyToSortV.Len(); j++) {
    		if ( (i+NodeGainCopyToSortV.Len() < NodeGainV.Len()) && (NodeGainCopyToSortV[ii].Val1 < NodeGainV[i+NodeGainCopyToSortV.Len()].Val1) ) {
    			NodeGainV[j] = NodeGainV[i+NodeGainCopyToSortV.Len()];
    			i++;
    		} else {
    			NodeGainV[j] = NodeGainCopyToSortV[ii];
    			ii++;
    		}
    	}
    }

    attempts = 0;

	for (int e = 0; e < NodeGainV.Len(); e++) {
	  const TInt& Node = NodeGainV[e].Val2;
	  if (Sources.IsNode(Node)) { continue; } // if node was already included in the sources set

	  const double NProb = Sources.GetAllNodeProb(Node, NodeInfoH);

	  NodeGainV[e].Val1 = NProb; // update marginal gain

	  if (BestGain < NProb) {
		BestGain = NProb;
		BestGainIndex = e;
		BestN = Node;
	  }

	  // if we only update one weight, we don't need to sort the list
	  attempts++;

	  // keep track of zero nodes after sorting once the full list
	  if (!Sources.IsNode(Node) && Sources.Len() > 1) {
		  if (NProb == 0)
			  NodeZero.Add(e);
	  }

	  // lazy evaluation
	  if (e+1 == NodeGainV.Len() || BestGain >= NodeGainV[e+1].Val1) {
		  Sources.CurGain += BestGain;
		  Sources.Add(BestN);

		  printf("sources:%d -> sigma:%f (+source:%d -> delta:%f, %d updates)\n", Sources.NIdV.Len(), Sources.CurGain.Val, BestN.Val, BestGain, attempts);

		  if (BestGain == 0)
			  return (-1);

		  // if we are going to store node info, add infection probability
		  NodeGainV.Del(BestGainIndex);

		  // we know the edges in 0 will be in sorted order, so we start from the biggest
		  for (int i=NodeZero.Len()-1; i>=0; i--) {
			  if (NodeZero[i] > BestGainIndex)
				  NodeGainV.Del(NodeZero[i]-1);
			  else
				  NodeGainV.Del(NodeZero[i]);
		  }

		  if (NodeZero.Len() > 2)
			  attempts -= (NodeZero.Len()-1);

		  msort = (attempts > 1);

		  LastGain = BestGain;

		  // average # of infected nodes for MAXINF solution
		  NumInfectedV.Add(NumInfectedV.Len()>0? NumInfectedV[NumInfectedV.Len()-1]+LastGain : LastGain);

		  // compute bound if needed (not efficiently implemented)
		  if (BoundOn) { NumInfectedBoundV.Add(GetBound(NumInfectedV[NumInfectedV.Len()-1].Val)); }

		  return BestN.Val;
	  }
	}

	printf("Nodes exhausted!\n");
	return (-1);
}

void TMaxInfBs::GetBestNodeCompleteSearch(const int &NumSources, int &IndexOptimalSourceSet, TIntV& OptimalSourceSet) {
	TIntIntVH CompleteCandidateSets;
	double BestGain = TFlt::Mn;
	int BestGainIndex = -1;

	OptimalSourceSet.Clr();
	Sources.NodeProbsH.Clr();

	// Generate all candidate sets for the # of sources (it will be massive if the network is large or the # of sources is large)
	GenerateCompleteCandidateSets(NumSources, CompleteCandidateSets);

	printf("%d candidate sets of length %d\n", CompleteCandidateSets.Len(), NumSources);

	for (int e = 0; e < CompleteCandidateSets.Len(); e++) {
	  const TIntV& CandidateSources = CompleteCandidateSets[e];
	  // add sources
	  Sources.Clr();
	  Sources.PNId = -1;
	  for (int i=0; i<CandidateSources.Len(); i++) { Sources.Add(CandidateSources[i]); }

	  const double NProb = Sources.GetInf(NodeInfoH);

	  if (BestGain < NProb) {
		BestGain = NProb;
		BestGainIndex = e;
	  }
	}

	for (int i=0; i<CompleteCandidateSets[BestGainIndex].Len(); i++) {
		 OptimalSourceSet.Add(CompleteCandidateSets[BestGainIndex][i]);
	}

	IndexOptimalSourceSet = BestGainIndex;

	// average # of infected nodes for MAXINF solution
	NumInfectedV.Add(BestGain);
}

void TMaxInfBs::GenerateCompleteCandidateSets(const int &NumSources, TIntIntVH& CompleteCandidateSets) {
	int n = GroundTruth->GetNodes();
	int m = NumSources;
	int *a = new int[n];
	int *p = new int[m];
	bool done = false;
	bool move_found = false;
	int i = 0; int j = 0; int k = 0;

	for (TNGraph::TNodeI NI = GroundTruth->BegNI(); NI < GroundTruth->EndNI(); NI++) { a[i++] = NI.GetId(); }

	// This is an adhoc algo that keeps m pointers to the next valid combination
	for (i=0; i<m; i++) { p[i]=i; } // the p[.] contain indices to the a vector whose elements constitute next combination

	while (!done) {
	    CompleteCandidateSets.AddDat(CompleteCandidateSets.Len()) = TIntV();

	    for (i=0; i<m; i++) { CompleteCandidateSets[CompleteCandidateSets.Len()-1].Add(a[p[i]]); }

	    // update combination
	    // method: start with p[m-1]. try to increment it. if it is already at the end, then try moving p[m-2] ahead.
	    // if this is possible, then reset p[m-1] to 1 more than (the new) p[m-2].
	    // if p[m-2] can not also be moved, then try p[m-3]. move that ahead. then reset p[m-2] and p[m-1].
	    // repeat all the way down to p[0]. if p[0] can not also be moved, then we have generated all combinations.
	    j = m-1;
	    i = 1;
	    move_found=false;
	    while ((j>=0) && !move_found)
	    {
	        if (p[j]<(n-i))
	        {
	            move_found=true;
	            p[j]++; // point p[j] to next index
	            for (k=j+1;k<m;k++)
	            {
	                p[k]=p[j]+(k-j);
	            }
	        }
	        else
	        {
	            j--;
	            i++;
	        }
	    }
	    if (!move_found) done=true;
	}
}

// computes bound using lazy evaluation concept to speed up
double TMaxInfBs::GetBound(const double& CurGain) {
	double Bound = 0;
	int k = 0;
	int attempts = 0;
	TFltV Bounds;

	// bound could be computed faster (using lazy evaluation, as in the optimization procedure)
	for (int e=0; e < NodeGainV.Len(); e++) {
		k = 0;
		const TInt& NN = NodeGainV[e].Val2;

		if (!Sources.IsNode(NN)) {
			const double NProb = Sources.GetAllNodeProb(NN, NodeInfoH);

			attempts++;

			for (int i=0; i<Bounds.Len(); i++) { if (e+1==NodeGainV.Len() || Bounds[i] >= NodeGainV[e+1].Val1) k++; }

			if (k == Sources.Len()) break;

			if (Bounds.Len() == 0) { Bounds.Add(NProb); continue; }

			Bounds.Sort(false);

			if (NProb > Bounds[Bounds.Len()-1]) { Bounds[Bounds.Len()-1] = NProb; }
			else if (Bounds.Len() < Sources.Len()) { Bounds.Add(NProb); }
		}
	}

	Bound += CurGain;
	for (int i=0; i<Sources.Len() && i<Bounds.Len(); i++) Bound += Bounds[i];

	printf("sources:%d -> sigma:%f / bound:%f (%d updates for bound)\n", Sources.Len(), Sources.CurGain.Val, Bound, attempts);

	return Bound;
}

// not scalable, only works with small networks
void TMaxInfBs::CompleteSearchOpt(const int &MxNodes) {
	double LastGain = TFlt::Mx;
	int attempts = 0;
	bool msort = false;
	int IndexOptimalSourceSet = -1;

	TIntV OptimalSourceSet;

	printf("nodes:%d number of sources:%d\nRunning complete search for influence maximization...\n", GroundTruth->GetNodes(), MxNodes);

	for (int k = 0; k < MxNodes && k < GroundTruth->GetNodes(); k++) {
		GetBestNodeCompleteSearch(k+1, IndexOptimalSourceSet, OptimalSourceSet);

		TIntFltH &Probs = Sources.NodeProbsH.GetDat(IndexOptimalSourceSet); // get probs of the optimal set of size k+1

		// save node info in hash
		for (int i=0; i<Probs.Len(); i++) { // TIntFltH
			if (!NodeInfoH.IsKey(Probs.GetKey(i))) { NodeInfoH.AddDat(Probs.GetKey(i)) = TNodeInfo(); }
			NodeInfoH.GetDat(Probs.GetKey(i)).Source = (OptimalSourceSet.IsIn(Probs.GetKey(i)))? k+1 : 0;
			NodeInfoH.GetDat(Probs.GetKey(i)).Prob.AddDat(k+1) = Probs[i];
		}
	}
}

void TMaxInfBs::GreedyOpt(const int& MxNodes) {
    double LastGain = TFlt::Mx;
    int attempts = 0;
    bool msort = false;

	printf("nodes:%d number of sources:%d\nRunning MAXINF...\n", GroundTruth->GetNodes(), MxNodes);

    Sources.CurGain = Sources.GetAllNodeProb(-1, NodeInfoH);

    for (int k = 0; k < MxNodes && NodeGainV.Len() > 0; k++) {
      double LastGain = Sources.CurGain;

      int BestN = -1;
      TExeTm tt;
      BestN = GetBestNode(LastGain, msort, attempts);
      if (TimingOn) { TimingV.Add(tt.GetSecs()); }

      if (BestN == -1) // if we cannot add more nodes, we stop
    	  break;

      // save node info in hash
      TIntFltH &Probs = Sources.NodeProbsH.GetDat(BestN);
      for (int i=0; i<Probs.Len(); i++) { // TIntFltH
    	  if (!NodeInfoH.IsKey(Probs.GetKey(i))) { NodeInfoH.AddDat(Probs.GetKey(i)) = TNodeInfo(); }
    	  NodeInfoH.GetDat(Probs.GetKey(i)).Source = (Probs.GetKey(i)==BestN)? k+1 : NodeInfoH.GetDat(Probs.GetKey(i)).Source.Val;
    	  NodeInfoH.GetDat(Probs.GetKey(i)).Prob.AddDat(k+1) = Probs[i];
      }
    }
}

void TMaxInfBs::OutDegreeOpt(const int &MxNodes) {
	TIntH NIdH;

	printf("nodes:%d number of sources:%d\nRunning outdegree for influence maximization...\n", GroundTruth->GetNodes(), MxNodes);

	// choose nodes with highest outdegree
	for (TNGraph::TNodeI NI = GroundTruth->BegNI(); NI < GroundTruth->EndNI(); NI++) { NIdH.AddDat(NI.GetId()) = NI.GetOutDeg(); }
	NIdH.SortByDat(false);

	// add nodes sequentially to store information about them
	for (int j=0; j<MxNodes; j++) {
		Sources.PNId = NIdH.GetKey(j);
		NumInfectedV.Add(Sources.GetInf(NodeInfoH));
		Sources.Add(Sources.PNId);

		// save node info in hash if needed, average over iters on the fly
		TIntFltH &Probs = Sources.NodeProbsH.GetDat(Sources.PNId);
		for (int i=0; i<Probs.Len(); i++) { // TIntFltH
			if (!NodeInfoH.IsKey(Probs.GetKey(i))) { NodeInfoH.AddDat(Probs.GetKey(i)) = TNodeInfo(); }
			NodeInfoH.GetDat(Probs.GetKey(i)).Source = (Probs.GetKey(i)==Sources.PNId)? j+1 : NodeInfoH.GetDat(Probs.GetKey(i)).Source.Val;
			if (!NodeInfoH.GetDat(Probs.GetKey(i)).Prob.IsKey(j+1)) { NodeInfoH.GetDat(Probs.GetKey(i)).Prob.AddDat(j+1) = 0.0; }
			NodeInfoH.GetDat(Probs.GetKey(i)).Prob.GetDat(j+1) = Probs[i];
		}
	}
}

void TMaxInfBs::RandomOpt(const int &MxNodes, const int &ItersRandom) {
	printf("nodes:%d number of sources:%d\nFinding %d random source sets for influence maximization...\n", GroundTruth->GetNodes(), MxNodes, ItersRandom);

	for (int reps=0; reps<ItersRandom; reps++) {
		Sources.Clr();
		TIntV NIdV; TSnap::GetRndSubGraph(GroundTruth, MxNodes)->GetNIdV(NIdV);

		// add nodes sequentially to store information about them
		for (int j=0; j<MxNodes; j++) {
			Sources.PNId = NIdV[j];
			if (reps==0) { NumInfectedV.Add(Sources.GetInf(NodeInfoH)); }
			else { NumInfectedV[j] = NumInfectedV[j] + Sources.GetInf(NodeInfoH); }
			Sources.Add(Sources.PNId);

			// save node info in hash if needed, average over iters on the fly
			if (Sources.NodeProbsH.IsKey(Sources.PNId)) {
			TIntFltH &Probs = Sources.NodeProbsH.GetDat(Sources.PNId);
			for (int i=0; i<Probs.Len(); i++) { // TIntFltH
				if (!NodeInfoH.IsKey(Probs.GetKey(i))) { NodeInfoH.AddDat(Probs.GetKey(i)) = TNodeInfo(); }
				NodeInfoH.GetDat(Probs.GetKey(i)).Source = (Probs.GetKey(i)==Sources.PNId)? j+1 : NodeInfoH.GetDat(Probs.GetKey(i)).Source.Val;
				if (!NodeInfoH.GetDat(Probs.GetKey(i)).Prob.IsKey(j+1)) { NodeInfoH.GetDat(Probs.GetKey(i)).Prob.AddDat(j+1) = 0.0; }
				NodeInfoH.GetDat(Probs.GetKey(i)).Prob.GetDat(j+1) = Probs[i];
			}
			}
		}

		SaveInfluence(TStr::Fmt("%s-rep-%d", OutFNm.CStr(), reps));
		Sources.NodeProbsH.Clr();
		NodeInfoH.Clr();
	}

	for (int j=0; j<MxNodes; j++) { NumInfectedV[j] = NumInfectedV[j]/ItersRandom; }
}

void TMaxInfBs::SourceListOpt(const TStr& ListSources) {
	TIntV NIdV;

	TStrV NIdStrV; ListSources.SplitOnAllCh(';', NIdStrV);
	for (int i=0; i<NIdStrV.Len(); i++) {
		NIdV.Add(NIdStrV[i].GetInt());
	}

	printf("nodes:%d number of sources:%d\nComputing influence for a given source list...\n", GroundTruth->GetNodes(), NIdV.Len());

	// add nodes sequentially to store information about them
	for (int j=0; j<NIdV.Len(); j++) {
		Sources.PNId = NIdV[j];
		NumInfectedV.Add(Sources.GetInf(NodeInfoH));
		Sources.Add(Sources.PNId);

		// save node info in hash if needed, average over iters on the fly
		TIntFltH &Probs = Sources.NodeProbsH.GetDat(Sources.PNId);
		for (int i=0; i<Probs.Len(); i++) { // TIntFltH
			if (!NodeInfoH.IsKey(Probs.GetKey(i))) { NodeInfoH.AddDat(Probs.GetKey(i)) = TNodeInfo(); }
			NodeInfoH.GetDat(Probs.GetKey(i)).Source = (Probs.GetKey(i)==Sources.PNId)? j+1 : NodeInfoH.GetDat(Probs.GetKey(i)).Source.Val;
			if (!NodeInfoH.GetDat(Probs.GetKey(i)).Prob.IsKey(j+1)) { NodeInfoH.GetDat(Probs.GetKey(i)).Prob.AddDat(j+1) = 0.0; }
			NodeInfoH.GetDat(Probs.GetKey(i)).Prob.GetDat(j+1) = Probs[i];
		}
	}
}

void TMaxInfBs::SourceFileOpt(const TStr& FileSources) {
	TIntV NIdV;

	TFIn FIn(FileSources);
	TStr Line; FIn.GetNextLn(Line); // first line contain number of sources

	printf("nodes:%d number of sources:%d\nComputing influence for a given source list...\n", GroundTruth->GetNodes(), Line.GetInt());

	while (!FIn.Eof()) {
		FIn.GetNextLn(Line);
		TStrV NIdStrV; Line.SplitOnAllCh('\t', NIdStrV);
		NIdV.Add(NIdStrV[0].GetInt());
	}

	// add nodes sequentially to store information about them
	for (int j=0; j<NIdV.Len(); j++) {
		Sources.PNId = NIdV[j];
		NumInfectedV.Add(Sources.GetInf(NodeInfoH));
		Sources.Add(Sources.PNId);

		// save node info in hash if needed, average over iters on the fly
		TIntFltH &Probs = Sources.NodeProbsH.GetDat(Sources.PNId);
		for (int i=0; i<Probs.Len(); i++) { // TIntFltH
			if (!NodeInfoH.IsKey(Probs.GetKey(i))) { NodeInfoH.AddDat(Probs.GetKey(i)) = TNodeInfo(); }
			NodeInfoH.GetDat(Probs.GetKey(i)).Source = (Probs.GetKey(i)==Sources.PNId)? j+1 : NodeInfoH.GetDat(Probs.GetKey(i)).Source.Val;
			if (!NodeInfoH.GetDat(Probs.GetKey(i)).Prob.IsKey(j+1)) { NodeInfoH.GetDat(Probs.GetKey(i)).Prob.AddDat(j+1) = 0.0; }
			NodeInfoH.GetDat(Probs.GetKey(i)).Prob.GetDat(j+1) = Probs[i];
		}
	}
}

void TMaxInfBs::SaveGroundTruth(const TStr& OutFNm) {
	TFOut FOut(OutFNm);

	// write nodes to file
	for (TNGraph::TNodeI NI = GroundTruth->BegNI(); NI < GroundTruth->EndNI(); NI++) {
		FOut.PutStr(TStr::Fmt("%d,%d\n", NI.GetId(), NI.GetId())); // nodes
	}

	FOut.PutStr("\n");

	// write edges to file (not allowing self loops in the network)
	for (TNGraph::TEdgeI EI = GroundTruth->BegEI(); EI < GroundTruth->EndEI(); EI++) {
		// not allowing self loops in the Kronecker network
		if (EI.GetSrcNId() != EI.GetDstNId()) {
			if (Alphas.IsKey(TIntPr(EI.GetSrcNId(), EI.GetDstNId())))
				FOut.PutStr(TStr::Fmt("%d,%d,%f\n", EI.GetSrcNId(), EI.GetDstNId(), Alphas.GetDat(TIntPr(EI.GetSrcNId(), EI.GetDstNId())).Val));
			else
				FOut.PutStr(TStr::Fmt("%d,%d,1\n", EI.GetSrcNId(), EI.GetDstNId()));
		}
	}
}

void TMaxInfBs::SaveGroundTruthForKdd10(const TStr& OutFNm) {
	TFOut FOut(OutFNm);

	THash<TIntPr, TFltPr> UEdges;

	// write edges to file (not allowing self loops in the network)
	for (TNGraph::TEdgeI EI = GroundTruth->BegEI(); EI < GroundTruth->EndEI(); EI++) {
		if (!UEdges.IsKey(TIntPr(EI.GetSrcNId(), EI.GetDstNId())) && !UEdges.IsKey(TIntPr(EI.GetDstNId(), EI.GetSrcNId()))) {
			UEdges.AddDat(TIntPr(EI.GetSrcNId(), EI.GetDstNId())) = TFltPr(1.0, 0.0);
		} else if (UEdges.IsKey(TIntPr(EI.GetSrcNId(), EI.GetDstNId()))) {
			UEdges.GetDat(TIntPr(EI.GetSrcNId(), EI.GetDstNId())).Val2 = 1.0;
		} else {
			UEdges.GetDat(TIntPr(EI.GetDstNId(), EI.GetSrcNId())).Val2 = 1.0;
		}
	}

	FOut.PutStr(TStr::Fmt("%d %d\n", GroundTruth->GetNodes(), UEdges.Len()));

	// there is a bug in the kdd10 code and deal badly with id's that are 0
	for (THash<TIntPr, TFltPr>::TIter EI = UEdges.BegI(); EI < UEdges.EndI(); EI++) {
		FOut.PutStr(TStr::Fmt("%d %d %f %f\n",
							  EI.GetKey().Val1.Val,
							  EI.GetKey().Val2.Val,
							  EI.GetDat().Val1.Val,
							  EI.GetDat().Val2.Val));

		FOut.PutStr(TStr::Fmt("%d %d %f %f\n",
							  EI.GetKey().Val2.Val,
							  EI.GetKey().Val1.Val,
							  EI.GetDat().Val2.Val,
							  EI.GetDat().Val1.Val));
	}
}

void TMaxInfBs::SavePajek(const TStr& OutFNm) {
    FILE *F = fopen(OutFNm.CStr(), "wt");
    fprintf(F, "*Vertices %d\n", GroundTruth->GetNodes());
	bool ZeroId = GroundTruth->IsNode(0);
    for (TNGraph::TNodeI NI = GroundTruth->BegNI(); NI < GroundTruth->EndNI(); NI++) {
      fprintf(F, "%d \"%d\" ic Blue\n", NI.GetId()+(int)ZeroId, NI.GetId()+(int)ZeroId);
    }
    fprintf(F, "*Arcs\n");
    for (TNGraph::TEdgeI EI = GroundTruth->BegEI(); EI < GroundTruth->EndEI(); EI++) {
      fprintf(F, "%d %d %f\n", EI.GetSrcNId()+(int)ZeroId, EI.GetDstNId()+(int)ZeroId, Alphas.GetDat(TIntPr(EI.GetSrcNId(), EI.GetDstNId())).Val);
    }
    fclose(F);
}

void TMaxInfBs::SaveBound(const TStr& OutFNm) {
	TFOut FOut(OutFNm);

	if (NumInfectedBoundV.Len() == NumInfectedV.Len()) { FOut.PutStr("% Sigma/Sigma bound\n"); }
	else { FOut.PutStr("% Sigma\n"); }

	for (int i=0; i<NumInfectedV.Len(); i++) {
		if (NumInfectedBoundV.Len() == NumInfectedV.Len()) {
			FOut.PutStr(TStr::Fmt("%f/%f\n", NumInfectedV[i].Val, NumInfectedBoundV[i].Val));
		} else {
			FOut.PutStr(TStr::Fmt("%f\n", NumInfectedV[i].Val));
		}
	}
}

void TMaxInfBs::SaveTiming(const TStr& OutFNm) {
	TFOut FOut(OutFNm);

	for (int i=0; i<TimingV.Len(); i++) {
		FOut.PutStr(TStr::Fmt("%f\n", TimingV[i].Val));
	}
}

void TMaxInfBs::SaveInfluence(const TStr& OutFNm) {
	TFOut FOut(OutFNm);

	FOut.PutStr("% Total average number of infected nodes\n");

	double Sigma = 0.0;
	for (THash<TInt, TNodeInfo>::TIter NI = NodeInfoH.BegI(); NI < NodeInfoH.EndI(); NI++) {
		Sigma += NI.GetDat().GetProbability();
	}

	FOut.PutStr(TStr::Fmt("%f\n", Sigma));

	FOut.PutStr("% ID,Infection Probability,Source\n");
	// write short info about each node
	for (TNGraph::TNodeI NI = GroundTruth->BegNI(); NI < GroundTruth->EndNI(); NI++) {
		if (!NodeInfoH.IsKey(NI.GetId())) { FOut.PutStr(TStr::Fmt("%d,0.0,0\n", NI.GetId())); continue; }

		FOut.PutStr(TStr::Fmt("%d,%f,%d\n",
							  NI.GetId(),
							  NodeInfoH.GetDat(NI.GetId()).GetProbability(),
							  NodeInfoH.GetDat(NI.GetId()).IsSource()? 1 : 0));
	}

	FOut.PutStr("\n% ID,Iter/Prob/Source;Iter/Prob/Source\n");
	// write long info about each node
	for (TNGraph::TNodeI NI = GroundTruth->BegNI(); NI < GroundTruth->EndNI(); NI++) {
		if (!NodeInfoH.IsKey(NI.GetId())) { FOut.PutStr(TStr::Fmt("%d,-1/0.0/0\n", NI.GetId())); continue; }

		FOut.PutStr(TStr::Fmt("%d", NI.GetId()));
		for (TIntFltH::TIter SI = NodeInfoH.GetDat(NI.GetId()).Prob.BegI(); SI < NodeInfoH.GetDat(NI.GetId()).Prob.EndI(); SI++) {
			FOut.PutStr(TStr::Fmt(",%d/%f/%d",
							  SI.GetKey().Val,
							  SI.GetDat().Val,
							  NodeInfoH.GetDat(NI.GetId()).Source.Val));
		}

		FOut.PutStr("\n");
	}
}

void TMaxInfBs::SaveInfluencePajek(const TStr& OutFNm) {
	FILE *F = fopen(OutFNm.CStr(), "wt");
	fprintf(F, "*Vertices %d\n", GroundTruth->GetNodes());
	bool ZeroId = GroundTruth->IsNode(0);
	for (TNGraph::TNodeI NI = GroundTruth->BegNI(); NI < GroundTruth->EndNI(); NI++) {
		TStr Position;
		if (Positions.Len() > 0) { Position = TStr::Fmt("%f %f 0.5", Positions.GetDat(NI.GetId()).Val1.Val, Positions.GetDat(NI.GetId()).Val2.Val); }
		else { Position = ""; }
		if (!NodeInfoH.IsKey(NI.GetId())) {
			fprintf(F, "%d \"%d\" %s ic White bc Black\n", NI.GetId()+(int)ZeroId, NI.GetId()+(int)ZeroId, Position.CStr());
		} else if (NodeInfoH.GetDat(NI.GetId()).IsSource()) {
			fprintf(F, "%d \"%d\" %s ic Black bc Red\n", NI.GetId()+(int)ZeroId, NI.GetId()+(int)ZeroId, Position.CStr());
		} else {
			int Prob = ((int)(100*NodeInfoH.GetDat(NI.GetId()).GetProbability()/5)) * 5;
			if (Prob > 95) { Prob = 95; }
			if (Prob < 5) { Prob = 5; }
			fprintf(F, "%d \"%d\" %s ic Gray%s%d bc Black\n", NI.GetId()+(int)ZeroId, NI.GetId()+(int)ZeroId, Position.CStr(), (Prob==5)? "0" : "",Prob);
		}
	}
	fprintf(F, "*Arcs\n");
	for (TNGraph::TEdgeI EI = GroundTruth->BegEI(); EI < GroundTruth->EndEI(); EI++) {
		fprintf(F, "%d %d %f\n", EI.GetSrcNId()+(int)ZeroId, EI.GetDstNId()+(int)ZeroId, Alphas.GetDat(TIntPr(EI.GetSrcNId(), EI.GetDstNId())).Val);
	}
	fclose(F);
}

void TMaxInfBs::SaveInfluenceSeqPajek(const int& Iters, const TStr& OutFNm) {
	for (int i=1; i<=Iters; i++) {
		FILE *F = fopen(TStr::Fmt("%s-%d-sources.net", OutFNm.CStr(), i).CStr(), "wt");
		fprintf(F, "*Vertices %d\n", GroundTruth->GetNodes());
		bool ZeroId = GroundTruth->IsNode(0);
		for (TNGraph::TNodeI NI = GroundTruth->BegNI(); NI < GroundTruth->EndNI(); NI++) {
			TStr Position;
			if (Positions.Len() > 0) { Position = TStr::Fmt("%f %f 0.5", Positions.GetDat(NI.GetId()).Val1.Val, Positions.GetDat(NI.GetId()).Val2.Val); }
			else { Position = ""; }
			if (!NodeInfoH.IsKey(NI.GetId()) || !NodeInfoH.GetDat(NI.GetId()).Prob.IsKey(i)) {
				fprintf(F, "%d \"%d\" %s ic White bc Black\n", NI.GetId()+(int)ZeroId, NI.GetId()+(int)ZeroId, Position.CStr());
			} else if (NodeInfoH.GetDat(NI.GetId()).Source.Val <= i && NodeInfoH.GetDat(NI.GetId()).Source.Val > 0) {
				fprintf(F, "%d \"%d\" %s ic Black bc Red\n", NI.GetId()+(int)ZeroId, NI.GetId()+(int)ZeroId, Position.CStr());
			} else {
				int Prob = ((int)(100*NodeInfoH.GetDat(NI.GetId()).Prob.GetDat(i)/5)) * 5;
				if (Prob > 95) { Prob = 95; }
				if (Prob < 5) { Prob = 5; }
				fprintf(F, "%d \"%d\" %s ic Gray%s%d bc Black\n", NI.GetId()+(int)ZeroId, NI.GetId()+(int)ZeroId, Position.CStr(), (Prob==5)? "0" : "",Prob);
			}
		}

		// edges are always the same for every set of sources
		fprintf(F, "*Arcs\n");
		for (TNGraph::TEdgeI EI = GroundTruth->BegEI(); EI < GroundTruth->EndEI(); EI++) {
			fprintf(F, "%d %d %f\n", EI.GetSrcNId()+(int)ZeroId, EI.GetDstNId()+(int)ZeroId, Alphas.GetDat(TIntPr(EI.GetSrcNId(), EI.GetDstNId())).Val);
		}
		fclose(F);
	}
}
