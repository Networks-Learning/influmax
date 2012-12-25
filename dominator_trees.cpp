#include "stdafx.h"

int main(int argc, char* argv[]) {
  Env = TEnv(argc, argv, TNotify::StdNotify);
  Env.PrepArgs(TStr::Fmt("\nDomTree. build: %s, %s. Time: %s", __TIME__, __DATE__, TExeTm::GetCurTm()));
  TExeTm ExeTm;
  Try
  const TStr GroundTruthFNm = Env.GetIfArgPrefixStr("-n:", "example-network.txt", "Input ground-truth network (one file)");
  const TStr OutFNm  = Env.GetIfArgPrefixStr("-o:", "network", "Output file name(s) prefix");

  const TStr Sources  = Env.GetIfArgPrefixStr("-s:", "-1;-1", "Sources");
  const int Sink = Env.GetIfArgPrefixInt("-k:", 1, "Sink");

  TMaxInfBs NIB;

  // load ground truth network
  TFIn FInG(GroundTruthFNm);
  NIB.LoadGroundTruthTxt(FInG);

  TDominatorGraph DomGraph;

  TIntV SinksI, SourcesI;
  SinksI.Add(Sink);

  TStrV SourcesV; Sources.SplitOnAllCh(';', SourcesV);
  for (int i=0; i<SourcesV.Len();i++) { if (SourcesV[i].GetInt()!=-1) { SourcesI.Add(SourcesV[i].GetInt()); } }

  DomGraph.SetSink(SinksI);
  DomGraph.SetSources(SourcesI);
  DomGraph.AddGroundTruth(NIB.GroundTruth);

  // compute dominator tree
  DomGraph.ComputeDomTree(SNCA);
  PNGraph DomTree = DomGraph.GetDomTree();

  // compute Gamma and M
  TIntV Gamma, M;
  DomGraph.GetDomGamma(Gamma, M);

  // print out Gamma
  printf("Gamma: ");
  for (int i=0; i<Gamma.Len(); i++) { printf("%d ", Gamma[i].Val); }
  printf("\n");

  // print out M
  printf("M: ");
  for (int i=0; i<M.Len(); i++) { printf("%d ", M[i].Val); }
  printf("\n");

  // compute K
  TIntV K;
  DomGraph.GetK(K);

  // print out K
  printf("K: ");
  for (int i=0; i<K.Len(); i++) { printf("%d ", K[i].Val); }
  printf("\n");

  // compute nu for each element in M
  for (int i=0; i<M.Len(); i++) {
	  TIntV Nu; DomGraph.GetNu(M[i].Val, Nu);
	  printf("Nu(%d): ", M[i].Val);
	  for (int j=0; j<Nu.Len(); j++) { printf("%d ", Nu[j].Val); }
	  printf("\n");
  }

  // save dominator tree
  DomGraph.SaveDomGraph(TStr::Fmt("%s-dominator-tree.txt", OutFNm.CStr()));

  Catch
  printf("\nrun time: %s (%s)\n", ExeTm.GetTmStr(), TSecTm::GetCurTm().GetTmStr().CStr());
  return 0;
}
