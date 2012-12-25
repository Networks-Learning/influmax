#include "stdafx.h"
#include "cascinf.h"

int main(int argc, char* argv[]) {
  Env = TEnv(argc, argv, TNotify::StdNotify);
  Env.PrepArgs(TStr::Fmt("\nGenerate different synthetic networks. build: %s, %s. Time: %s", __TIME__, __DATE__, TExeTm::GetCurTm()));
  TExeTm ExeTm;
  Try

  const int TNetwork = Env.GetIfArgPrefixInt("-t:", 0,
		  "Type of network to generate\n0:kronecker, 1:forest fire, 2:given (default:0)\n");

  const TStr NetworkParams = Env.GetIfArgPrefixStr("-g:", TStr("0.9 0.5; 0.5 0.9"),
		  "Parameters for the network (default:0.9 0.5; 0.5 0.9)\n");

  const int NNodes = Env.GetIfArgPrefixInt("-n:", 128, "Number of nodes (default:128)\n");
  const int NEdges = Env.GetIfArgPrefixInt("-e:", 128, "Number of edges (default:128)\n");

  const int DAlphas = Env.GetIfArgPrefixInt("-ad:", 0, "Distributions for alpha\n0:uniform, 1:gaussian, 2:rayleigh");

  // up to 2 distributions
  const TStr RAlphas = Env.GetIfArgPrefixStr("-ar:", TStr("1;0.5"),
		  "Sigma (-ad:2) or mean and either range (-ad:0) or std dev (-ad:1) for alpha (default:1;0.5)\n");

  const TStr FileName = Env.GetIfArgPrefixStr("-f:", TStr("example"), "Name of the network (default:example)\n");

  const TStr GroundTruthFileName = Env.GetIfArgPrefixStr("-n:", TStr("input"), "Name of the input network (default:input)\n");

  TMaxInfBs MIB;

  // Generate GroundTruth
  if (TNetwork < 2) { MIB.GenerateGroundTruth(TNetwork, NNodes, NEdges, NetworkParams); }
  else {
	  TFIn GFIn(GroundTruthFileName);
	  MIB.LoadGroundTruthTxt(GFIn);
  }

  // Generate Alphas
  MIB.GenerateAlphas(DAlphas, RAlphas);

  // Save GroundTruth in txt and pajek format
  MIB.SaveGroundTruth(TStr::Fmt("%s-network.txt", FileName.CStr()));

  Catch
  printf("\nrun time: %s (%s)\n", ExeTm.GetTmStr(), TSecTm::GetCurTm().GetTmStr().CStr());
  return 0;
}
