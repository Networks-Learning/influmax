#include "stdafx.h"

int main(int argc, char* argv[]) {
  Env = TEnv(argc, argv, TNotify::StdNotify);
  Env.PrepArgs(TStr::Fmt("\nINFLUMAX. build: %s, %s. Time: %s", __TIME__, __DATE__, TExeTm::GetCurTm()));
  TExeTm ExeTm;
  Try
  const TStr GroundTruthFNm = Env.GetIfArgPrefixStr("-n:", "example-network.txt", "Input ground-truth network (one file)");
  const TStr OutFNm  = Env.GetIfArgPrefixStr("-o:", "network", "Output file name(s) prefix");
  const int NSources  = Env.GetIfArgPrefixInt("-s:", 1, "Number of sources");
  const double T = Env.GetIfArgPrefixFlt("-t:", 10.0, "Time window");
  const int ExpMethod = Env.GetIfArgPrefixInt("-e:", 3, "Matrix exponential. 0:Parlett, 1:Krylov, 2:Hessenberg, 3:Krylov(big)+Hessenberg(small)\n");
  const int MxReachable = Env.GetIfArgPrefixInt("-mr:", -1, "Maximum length for minimum reachable path from source(s)\n");
  const int MxLPath = Env.GetIfArgPrefixInt("-mp:", -1, "Maximum length for diffusion path from source(s)\n");
  const int Baselines = Env.GetIfArgPrefixInt("-b:", 0, "Baselines\n\
    0:INFLUMAX, 1:random sources, 2:outdegree, 3:complete search, 4:given source list, 5:given source file\n");
  const int ItersRandom = Env.GetIfArgPrefixInt("-ir:", 100, "Repetitions random sources selection if -b:1\n");
  const TStr ListSources = Env.GetIfArgPrefixStr("-ls:", "", "List of sources (id1;id2;id3) if -b:4 or source filename if -b:5\n");
      
  TMaxInfBs MIB;

  // load ground truth network
  TFIn FInG(GroundTruthFNm);

  bool ComputeBound = false; bool ComputeCompleteSearch = false;

  MIB.LoadGroundTruthTxt(FInG);
  MIB.SetWindow(T);
  MIB.SetBound(ComputeBound);
  MIB.SetMxLengthPath(MxLPath);
  MIB.SetMxReachable(MxReachable);
  MIB.SetExpMethod(ExpMethod==0? PARLETT : ((ExpMethod==1)? KRYLOV : ((ExpMethod==2)? HESSENBERG : MIXED)));
  MIB.SetBfs(false);

  if (MxLPath < MxReachable) { printf("WARNING: MxLPath:%d < MxReachable:%d!\n", MxLPath, MxReachable); }

  MIB.OutFNm = OutFNm;

  MIB.Init();

  if (Baselines==0) { // maxinf
  	  MIB.GreedyOpt(NSources);
  } else if (Baselines==1) { // equal alphas
	  MIB.RandomOpt(NSources, ItersRandom);
  } else if (Baselines==2) { // outdegree
	  MIB.OutDegreeOpt(NSources);
  } else if (Baselines==3) { // complete search (if the network is big, this is not scalable!)
	  MIB.CompleteSearchOpt(NSources);
  } else if (Baselines==4) {
	  MIB.SourceListOpt(ListSources);
  } else if (Baselines==5) {
	  MIB.SourceFileOpt(ListSources);
  } else {
	  FailR("Bad -s: parameter.");
  }

  // plot objective function
  MIB.SaveBound(TStr::Fmt("influence-average-%s.txt", OutFNm.CStr()));

  // save node info
  MIB.SaveInfluence(TStr::Fmt("influence-info-%s.txt", OutFNm.CStr()));

  Catch
  printf("\nrun time: %s (%s)\n", ExeTm.GetTmStr(), TSecTm::GetCurTm().GetTmStr().CStr());
  return 0;
}
