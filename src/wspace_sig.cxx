#include "wspace_sig.hpp"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <initializer_list>
#include <vector>
#include <string>
#include <stdlib.h>
#include <ctime>
#include <sys/stat.h>

#include <unistd.h>
#include <getopt.h>

#include "TString.h"
#include "TSystem.h"
#include "TDirectory.h"

#include "bin.hpp"
#include "process.hpp"
#include "utilities.hpp"
#include "systematic.hpp"
#include "cut.hpp"
#include "workspace_generator.hpp"


using namespace std;

namespace{
  double lumi = 20;
  double sig_strength = 0.;
  BlindLevel blind_level = BlindLevel::blinded;
  bool no_kappa = false;
  bool do_syst = true;
  bool use_r4 = true;
  unsigned n_toys = 0;
  string outfolder = "out/";
  string nb_bins("TML");
  string sigfile = "/net/cms27/cms27r0/babymaker/2016_04_29/mc/merged_higloose/*TChiHH_mChi-400*.root";
}

int main(int argc, char *argv[]){
  time_t begtime, endtime;
  time(&begtime);
  cout << fixed << setprecision(2);
  GetOptions(argc, argv);
  if(sigfile==""){
    cout<<endl<<"You need to specify the input file with -f. Exiting"<<endl<<endl;
    return 1;
  }
  string hostname = execute("echo $HOSTNAME");
  string basefolder("/net/cms2/cms2r0/babymaker/");
  if(Contains(hostname, "lxplus")) basefolder = "/afs/cern.ch/user/m/manuelf/work/";
  string foldermc(basefolder+"babies/2016_04_29/mc/merged_higloose/"); 
  string folderdata(basefolder+"babies/2016_04_29/data/merged_abcd/"); // Pointing to a random folder

  
  //Define processes. Try to minimize splitting
  string stitch_cuts("stitch&&pass");
  Process ttbar{"ttbar", {
      {foldermc+"/*_TTJets*Lept*.root/tree"},
	{foldermc+"/*_TTJets*HT*.root/tree"}
    },stitch_cuts};

  Process other{"other", {
      {foldermc+"/*DYJetsToLL*.root/tree"},
	{foldermc+"/*_WWTo*.root/tree"},
	  {foldermc+"/*_TTTT*.root/tree"},
	    {foldermc+"/*_WZ*.root/tree"},
	      {foldermc+"/*QCD_HT*.root/tree"},
		{foldermc+"/*_WJetsToLNu*.root/tree"},
		  {foldermc+"/*_TTWJets*.root/tree"},
		    {foldermc+"/*_TTZTo*.root/tree"},
		      {foldermc+"/*_ST_*.root/tree"},
			{foldermc+"/*_TTGJets*.root/tree"},
			  {foldermc+"/*ttHJetTobb*.root/tree"}
    },stitch_cuts};
  Process signal{"signal", {
      //{sigfile+"/tree"}
      {foldermc+"/*TChiHH_mChi-400*.root/tree"}
    },"pass", false, true};

  string data_cuts("(trig[4]||trig[8])&&pass");

  Process data{"data", {
      {folderdata+"/*.root/tree"}
    }, data_cuts, true};

  //Make list of all backgrounds. Backgrounds assumed to be orthogonal
  set<Process> backgrounds{ttbar, other};

  //Baseline selection applied to all bins and processes
  Cut baseline{"hig_drmax<2.2&&ntks==0"}; //njets>=4&&njets<=5&&!low_dphi&&nvleps==0 already in the skim

  string cut2b="nbt==2&&nbm==2", cut3b="nbt>=2&&nbm==3&&nbl==3", cut4b="nbt>=2&&nbm>=3&&nbl>=4";
  if(nb_bins=="TTL"){
    cut2b = "nbt==2";
    cut3b = "nbt==3&&nbl==3";
    cut4b = "nbt>=3&&nbl>=4";
  }
  if(nb_bins=="MMM"){
    cut2b = "nbm==2";
    cut3b = "nbm==3";
    cut4b = "nbm>=4";
  }
  if(nb_bins=="TMM"){
    cut2b = "nbt==2&&nbm==2";
    cut3b = "nbm==3";
    cut4b = "nbm>=4";
  }
  string cutme0="&&met>250&&met<=300", cutmet1="&&met>300";
  string cutsig="hig_am>100&&hig_am<140&&hig_dm<40", cutsbd="!("+cutsig+")";

  Bin sbd_2b_met0{"sbd_2b_met0", cut2b+"&&"+cutsbd+cutme0, blind_level>=BlindLevel::blinded};
  Bin sig_2b_met0{"sig_2b_met0", cut2b+"&&"+cutsig+cutme0, blind_level>=BlindLevel::blinded};
  Bin sbd_3b_met0{"sbd_3b_met0", cut3b+"&&"+cutsbd+cutme0, blind_level>=BlindLevel::blinded};
  Bin sig_3b_met0{"sig_3b_met0", cut3b+"&&"+cutsig+cutme0, blind_level>=BlindLevel::blinded};
  Bin sbd_4b_met0{"sbd_4b_met0", cut4b+"&&"+cutsbd+cutme0, blind_level>=BlindLevel::blinded};
  Bin sig_4b_met0{"sig_4b_met0", cut4b+"&&"+cutsig+cutme0, blind_level>=BlindLevel::blinded};

  Bin sbd_2b_met1{"sbd_2b_met1", cut2b+"&&"+cutsbd+cutmet1, blind_level>=BlindLevel::blinded};
  Bin sig_2b_met1{"sig_2b_met1", cut2b+"&&"+cutsig+cutmet1, blind_level>=BlindLevel::blinded};
  Bin sbd_3b_met1{"sbd_3b_met1", cut3b+"&&"+cutsbd+cutmet1, blind_level>=BlindLevel::blinded};
  Bin sig_3b_met1{"sig_3b_met1", cut3b+"&&"+cutsig+cutmet1, blind_level>=BlindLevel::blinded};
  Bin sbd_4b_met1{"sbd_4b_met1", cut4b+"&&"+cutsbd+cutmet1, blind_level>=BlindLevel::blinded};
  Bin sig_4b_met1{"sig_4b_met1", cut4b+"&&"+cutsig+cutmet1, blind_level>=BlindLevel::blinded};

  //// Defining the 2x3 ABCD in bins of met
  set<Block> blocks_abcd;

  blocks_abcd = {
    {"met0", {{sbd_2b_met0, sbd_3b_met0, sbd_4b_met0},
	      {sig_2b_met0, sig_3b_met0, sig_4b_met0}}},
    {"met1", {{sbd_2b_met1, sbd_3b_met1, sbd_4b_met1},
	      {sig_2b_met1, sig_3b_met1, sig_4b_met1}}}
  };

  //// Creating workspaces
  Cut *pbaseline(&baseline);
  set<Block> *pblocks(&blocks_abcd);

  string sysfolder = "txt/systematics/";
  string sysfile(sysfolder+"/sys_SMS-TChiHH_4b_mChi-400.txt");
  
  // If systematic file does not exist, complain
  struct stat buffer;   
  if(stat (sysfile.c_str(), &buffer) != 0) {
    cout<<endl<<"WARNING: "<<sysfile<<" does not exist. Using ";
    sysfile = "txt/systematics/sys_SMS-TChiHH_4b_mChi-400.txt";
    cout<<sysfile<<" instead"<<endl<<endl;
  }

  gSystem->mkdir(outfolder.c_str(), kTRUE);
  string sig_s = "_sig"+to_string(sig_strength);
  ReplaceAll(sig_s, ".000000","");
  string outname(outfolder+"/wspace_SMS-TChiHH_4b_mChi-400_nb"+nb_bins+sig_s+".root");
  if(!use_r4) ReplaceAll(outname, "wspace_","wspace_nor4_");

  float rmax = 20.;
  WorkspaceGenerator wgNom(*pbaseline, *pblocks, backgrounds, signal, data, sysfile, use_r4, sig_strength, 1.);
  wgNom.SetRMax(rmax);
  wgNom.SetKappaCorrected(!no_kappa);
  wgNom.SetDoSystematics(do_syst);
  wgNom.SetLuminosity(lumi);
  wgNom.SetDoSystematics(do_syst);
  wgNom.AddToys(n_toys);
  wgNom.WriteToFile(outname);

  time(&endtime); 
  cout<<"Finding workspaces took "<<fixed<<setprecision(0)<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;  
}



void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"sigfile", required_argument, 0, 'f'},
      {"lumi", required_argument, 0, 'l'},
      {"nb_bins", no_argument, 0, 'b'},
      {"unblind", required_argument, 0, 'u'},
      {"no_syst", no_argument, 0, 0},
      {"nokappa", no_argument, 0, 'k'},
      {"no_r4", no_argument, 0, '4'},
      {"toys", required_argument, 0, 0},
      {"sig_strength", required_argument, 0, 'g'},
      {"outfolder", required_argument, 0, 'o'},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "f:l:b:u:k4g:o:", long_options, &option_index);
    if( opt == -1) break;

    string optname;
    switch(opt){
    case 'l':
      lumi = atof(optarg);
      break;
    case 'b':
      nb_bins = optarg;
      break;
    case 'g':
      sig_strength = atof(optarg);
      break;
    case 'u':
      if(string(optarg)=="all"){
        blind_level = BlindLevel::unblinded;
      }else if(string(optarg)=="sideband"){
        blind_level = BlindLevel::r4_blinded;
      }else if(string(optarg)=="1b"){
        blind_level = BlindLevel::unblind_1b;
      }else{
        blind_level = BlindLevel::blinded;
      }
      break;
    case 'o':
      outfolder = optarg;
      break;
    case 'f':
      sigfile = optarg;
      break;
    case 'k':
      no_kappa = true;
      break;
    case '4':
      use_r4 = false;
      break;
    case 0:
      optname = long_options[option_index].name;
      if(optname == "no_syst"){
        do_syst = false;
      }else if(optname == "toys"){
        n_toys = atoi(optarg);
      }else{
        printf("Bad option! Found option name %s\n", optname.c_str());
      }
      break;
    default:
      printf("Bad option! getopt_long returned character code 0%o\n", opt);
      break;
    }
  }
}
