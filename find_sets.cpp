#include "stdafx.h"

int main(int argc, char* argv[]) {
  Env = TEnv(argc, argv, TNotify::StdNotify);
  Env.PrepArgs(TStr::Fmt("\nFindSets. build: %s, %s. Time: %s", __TIME__, __DATE__, TExeTm::GetCurTm()));
  TExeTm ExeTm;
  Try
  const TStr GroundTruthFNm = Env.GetIfArgPrefixStr("-n:", "example-network.txt", "Input ground-truth network (one file)");

  const TStr Sources  = Env.GetIfArgPrefixStr("-s:", "-1;-1", "Sources");
  const int Sink = Env.GetIfArgPrefixInt("-k:", 1, "Sink");

  const double T = Env.GetIfArgPrefixFlt("-t:", 10.0, "Time window");

  const int MxReachable = Env.GetIfArgPrefixInt("-mr:", -1, "Maximum length for minimum reachable path from source(s)\n");
  const int MxLPath = Env.GetIfArgPrefixInt("-mp:", -1, "Maximum length for diffusion path from source(s)\n");

  const int ExpMethod = Env.GetIfArgPrefixInt("-e:", 3, "Matrix exponential. 0:Parlett, 1:Krylov, 2:Hessenberg, 3:Krylov(big)+Hessenberg(small)\n");

  const int Output = Env.GetIfArgPrefixInt("-a:", 0, "0:infection probability, 1: state size, 2:all\n");
  TMaxInfBs MIB;
  TIntV SourcesI;

  // load ground truth network in NIB and OSources
  TFIn FInG(GroundTruthFNm);
  MIB.LoadGroundTruthTxt(FInG);

  MIB.SetWindow(T);
  MIB.SetMxLengthPath(MxLPath);
  MIB.SetMxReachable(MxReachable);
  MIB.Init();
  MIB.SetExpMethod(ExpMethod==0? PARLETT : ((ExpMethod==1)? KRYLOV : ((ExpMethod==2)? HESSENBERG : MIXED)));

  // add sources and sink
  TStrV SourcesV; Sources.SplitOnAllCh(';', SourcesV);
  for (int i=0; i<SourcesV.Len();i++) {
	  if (SourcesV[i].GetInt()!=-1) { SourcesI.Add(SourcesV[i].GetInt()); MIB.Sources.Add(SourcesI[SourcesI.Len()-1]); }
  }

  // compute P(t_sink <= T)
  double Prob = MIB.Sources.GetInf(Sink, MIB.NodeInfoH, true);

  if (Output==0 || Output==2) {
	  printf("\nP(t_%d <= %f) = %f\n", Sink, T, Prob);
  }

  if (Output==1 || Output==2) {
	  printf("\nState size -> Nodes:%d, Edges:%d", MIB.Sources.StateSize[0].Val1.Val, MIB.Sources.StateSize[0].Val2.Val);
  }

  Catch
  printf("\nrun time: %s (%s)\n", ExeTm.GetTmStr(), TSecTm::GetCurTm().GetTmStr().CStr());
  return 0;
}
