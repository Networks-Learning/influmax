#ifndef snap_cascinf_h
#define snap_cascinf_h

#include "Snap.h"
#include "dgraph.h"

typedef enum {
	SNCA, // snca algorithm
	SLT, // sophisticated Lengauer-Tarjan algorithm (with balanced link-eval)
	LT // simplified Lengauer-Tarjan algorithm (with standard link-eval)
} TDomMethod;

typedef enum {
	PARLETT, // Parlett's method (u_ii != u_jj for any ii, jj)
	KRYLOV, // Krylov's method for sparse matrices (using expokit fortran code)
	HESSENBERG, // Chebyshev's method for Hessenberg matrices (using expokit fortran code)
	MIXED // KRYLOV for big matrices, HESSENBERG for small ones
} TExpMethod;

extern"C" {
	void dnchbv_( int *m, double *t, double *H, int *ldh, double *y, double *wsp );
	void dgexpv_( int *n, int *nz, int *m, double *t, double *v, double *w, double *tol, double *anorm,
				  double *wsp, int *lwsp, int *iwsp, int *liwsp, int *itrace, double *H, int *ia, int *ja,
				  int *iflag );
}

class TDominatorGraph {
public:
	TIntV Sources, Sink;
	PNGraph GroundTruth;
	PNGraph DominatorTree;

public:
	TDominatorGraph() : Sources(), Sink() { }
	TDominatorGraph(TIntV& Ss, TIntV& Sk) : Sources(Ss), Sink(Sk) { }
	TDominatorGraph(TSIn& SIn) : Sources(SIn), Sink(SIn), DominatorTree(SIn) { }
	void Save(TSOut& SOut) const { Sources.Save(SOut); Sink.Save(SOut); DominatorTree.Save(SOut); }
	void Clr() { Sources.Clr(); Sink.Clr(); DominatorTree.Clr(); }

	void AddGroundTruth(PNGraph& gt) { GroundTruth = gt; }
	void SetSink(TIntV& Sk) { Sink = Sk; }
	void SetSources(TIntV& Ss) { Sources = Ss; }

	void ComputeDomTree(const TDomMethod& method);

	PNGraph GetDomTree() { return DominatorTree; }

	void GetNu(const int &root, TIntV &Nu) {
		PNGraph BTree = TSnap::GetBfsTree(DominatorTree, root, false, true);
		BTree->GetNIdV(Nu);
		BTree.Clr(); // release BTree memory
	}
	void GetDomGamma(TIntV& Gamma, TIntV &M);
	void GetK(TIntV& K);

	void SaveDomGraph(const TStr& OutFNm);
};

// Node info (name and probability of infection)
class TNodeInfo {
public:
  TStr Name;
  TIntFltH Prob;
  TInt Source;
public:
  TNodeInfo() { Source = 0; }
  TNodeInfo(const TStr& NodeNm) : Name(NodeNm), Source(0) { }
  TNodeInfo(TSIn& SIn) : Name(SIn), Prob(SIn), Source(SIn) { }
  void Save(TSOut& SOut) const { Name.Save(SOut); Prob.Save(SOut); Source.Save(SOut); }
  bool IsSource() { return (Source>0); }
  void SetSource(const int& Iter) { Source = Iter; }
  double GetProbability() { return Prob[Prob.Len()-1].Val; }
  void AddProb(const int& Iter, const int& Probability) { Prob.AddDat(Iter) = Probability; }
};

// Sources
class TSources {
public:
  TIntV NIdV; // optimal sources
  TInt PNId; // potential source
  TFlt CurGain; // current optimal objective function

  PNGraph GroundTruth; // groundtruth
  TIntPrFltH Alphas; // tx rates groundtruth

  PNGraph GroundTruthReachable, GroundTruthReachableCurrent, GroundTruthMetanode; // reachable groundtruth from sources/sink

  TInt MxLengthPath, MxReachable; // max length path / max min length path
  TFlt Window; // time window
  TExpMethod ExpMethod; // exponential method

  TBool Bfs;

  THash<TInt, TIntFltH> NodeProbsH; // current iteration Node Probabilities
public:
  TSources() : NIdV(), PNId(-1), MxLengthPath(-1), MxReachable(-1), CurGain(0.0), Window(10.0), ExpMethod(HESSENBERG), Bfs(false) { }
  TSources(const float& Win) : NIdV(), PNId(-1), MxLengthPath(-1), MxReachable(-1), CurGain(0.0), Window(Win), ExpMethod(HESSENBERG), Bfs(false) { }
  TSources(TSIn& SIn) : NIdV(SIn), PNId(SIn), MxLengthPath(SIn), MxReachable(SIn), CurGain(SIn), Window(SIn), ExpMethod(HESSENBERG), Bfs(false) { }
  void Save(TSOut& SOut) const  { NIdV.Save(SOut); CurGain.Save(SOut); Window.Save(SOut); }
  void Clr() { NIdV.Clr(); PNId = -1; CurGain = 0.0; }
  int Len() const { return NIdV.Len(); }
  int GetNode(const int& Elem) const { return NIdV[Elem]; }
  void Add(const int& NId) { NIdV.Add(NId); }
  void AddPotential(const int& NId) { PNId = NId; }
  bool IsNode(const int& NId) const { return NIdV.IsIn(NId); }

  void SetWindow(const float& Win) { Window = Win; }
  void SetGroundTruth(PNGraph& gt) { GroundTruth = gt; }
  void SetMxLengthPath(const int& Mxlp) { MxLengthPath = Mxlp; }
  void SetMxReachable(const int& Mxr) { MxReachable = Mxr; }
  void SetReachableGroundTruthFromSources(const int& Sink, TIntV &ReachableNodesV);
  void SetReachableGroundTruthFromSink(const int& Sink, TIntV &ReachableNodesV);
  void SetBfs(const bool& BfsOn) { Bfs = BfsOn; }

  void SetAlphas(TIntPrFltH& Alphs) { Alphas = Alphs; }
  void SetExpMethod(TExpMethod em) { ExpMethod = em; }

  double GetAllNodeProb(const int& Source, THash<TInt, TNodeInfo>& NodeInfoH);
  void InitInf();
  double UpdateInf(const int& Source, THash<TInt, TNodeInfo>& NodeInfoH);
  double GetInf(THash<TInt, TNodeInfo>& NodeInfoH);
  double GetInf(const int& Sink, THash<TInt, TNodeInfo>& NodeInfoH);
  double GetInf(const int& Sink);

  double GetInf(const int& Sink, TFltPrV& Cdf, const double& MinT, const double& MaxT, const double& Step);

  void GetDAG(const int& Sink, TPair<PNGraph, TIntPrFltH> &DAG);
  void ComputeEdges(const int& Sink, TVec<TIntFltH> &Sets, TVec<TIntFltH> &Connectivity);
  int CheckBlocks(const int& Sink, TIntFltH &Set, const int &v, TIntV& SupersetNotInSubset);
  void GetSets(TIntV& S, TIntV& T, TVec<TIntFltH> &Sets, int NumInfected);
  void GetPivot(TIntV& S, TIntV& T, TPair<TInt, TIntV>& NodePivot, const int& NumSets);

  double ExpM(TPair<PNGraph, TIntPrFltH>& DAG, const TExpMethod& ExpMethod);
};

// MAXINF algorithm class
class TMaxInfBs {
public:
  PNGraph GroundTruth;
  TIntPrFltH Alphas;
  THash<TInt, TFltPr> Positions;

  TFltIntPrV NodeGainV; // node gains
  THash<TInt, TNodeInfo> NodeInfoH; // nodes info
  TSources Sources; // sources object
  TFltV NumInfectedV, NumInfectedBoundV; // ave number of infected and bound

  TStr OutFNm;

  bool BoundOn, EqualRates, RandomSources, Bfs, TimingOn;

  TFltV TimingV;

public:
  TMaxInfBs( ) { BoundOn = false; TimingOn = false; }
  TMaxInfBs(TSIn& SIn) : GroundTruth(SIn), Alphas(SIn), NodeGainV(SIn), Sources(SIn) { }
  void Save(TSOut& SOut) const { GroundTruth.Save(SOut); Alphas.Save(SOut); NodeGainV.Save(SOut); Sources.Save(SOut); }

  void LoadGroundTruthTxt(TSIn& SIn);

  void AddGroundTruth(PNGraph& gt) { GroundTruth = gt; }
  void GenerateGroundTruth(const int& TNetwork, const int& NNodes, const int& NEdges, const TStr& NetworkParams);
  void GenerateAlphas(const int& DAlphas, const TStr& RAlphas);
  void SetWindow(const float& Win) { Sources.SetWindow(Win); }
  void SetBound(const float& Bound) { BoundOn = Bound; }
  void SetMxLengthPath(const int& Mxlp) { Sources.SetMxLengthPath(Mxlp); }
  void SetMxReachable(const int& Mxr) { Sources.SetMxReachable(Mxr); }
  void SetExpMethod(TExpMethod em) { Sources.SetExpMethod(em); }
  void SetEqualRates(const bool& EqRates) { EqualRates = EqRates; }
  void SetRandomSources(const bool& RandomSrcs) { RandomSources = RandomSrcs; }
  void SetBfs(const bool& BfsOn) { Bfs = BfsOn; Sources.SetBfs(BfsOn); }
  void SetTiming(const int& ComputeTiming) { TimingOn = (ComputeTiming==1); }

  void Init();
  int GetBestNode(double& LastGain, bool& msort, int &attempts);
  int GetBestNodeDfs(double& LastGain, bool& msort, int &attempts);
  void GetBestNodeCompleteSearch(const int& NumSources, int &IndexOptimalSourceSet, TIntV& OptimalSourceSet);
  double GetBound(const double& CurInf);

  void GenerateCompleteCandidateSets(const int &NumSources, TIntIntVH& CompleteCandidateSets);

  void CompleteSearchOpt(const int &MxNodes); // not scalable, only works with small networks
  void GreedyOpt(const int& MxNodes);
  void GreedyEqualAlphasOpt(const int &MxNodes, const double &AlphaEquals);
  void GreedyBfsOpt(const int& MxNodes);
  void OutDegreeOpt(const int &MxNodes);
  void OutRateOpt(const int &MxNodes);
  void RandomOpt(const int &MxNodes, const int &ItersRandom);
  void SourceListOpt(const TStr& ListSources);
  void SourceFileOpt(const TStr& FileSources);

  void SaveGroundTruth(const TStr& OutFNm);
  void SaveGroundTruthForKdd10(const TStr& OutFNm);
  void SavePajek(const TStr& OutFNm);
  void SaveBound(const TStr& OutFNm);
  void SaveTiming(const TStr& OutFNm);
  void SaveInfluence(const TStr& OutFNm);
  void SaveInfluencePajek(const TStr& OutFNm);
  void SaveInfluenceSeqPajek(const int& Iters, const TStr& OutFNm);
};

#endif
