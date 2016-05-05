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
  string sigfile = "/net/cms27/cms27r0/babymaker/2016_04_29/mc/merged_met100nb2nj4nl0/mergedbaby__SMS-TChiHH_mChi-400_madgraphMLM-pythia8__met100nb2nj4nl0_nfiles_1.root";
  string outfolder = "out/";
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
  string foldermc(basefolder+"babies/2016_04_29/mc/merged_met100nb2nj4nl0/"); 
  string folderdata(basefolder+"babies/2016_04_29/data/merged_abcd/"); // Pointing to a random folder

  
  //Define processes. Try to minimize splitting
  string stitch_cuts("stitch&&pass");
  Process ttbar{"ttbar", {
      {foldermc+"/*_TTJets*.root/tree"}
    },stitch_cuts};

  Process other{"other", {
      {foldermc+"/*_WJetsToLNu*.root/tree"},
        {foldermc+"/*_TTWJets*.root/tree"},
          {foldermc+"/*_TTZTo*.root/tree"},
            {foldermc+"/*_ST_*.root/tree"},
              {foldermc+"/*DYJetsToLL*.root/tree"},
                {foldermc+"/*QCD_HT*.root/tree"},
                  {foldermc+"/*_WWTo*.root/tree"},
                    {foldermc+"/*_TTGJets*.root/tree"},
		      {foldermc+"/*_TTTT*.root/tree"},
			{foldermc+"/*_WZ*.root/tree"},
			  {foldermc+"/*ttHJetTobb*.root/tree"}
    },stitch_cuts};
  Process signal{"signal", {
      {sigfile+"/tree"}
    },stitch_cuts, false, true};

  string data_cuts("(trig[4]||trig[8])&&pass");

  Process data{"data", {
      {folderdata+"/*.root/tree"}
    }, data_cuts, true};

  //Make list of all backgrounds. Backgrounds assumed to be orthogonal
  set<Process> backgrounds{ttbar, other};

  //Baseline selection applied to all bins and processes
  Cut baseline{"njets>=4&&njets<=5&&!low_dphi"}; //MET and nvleps==0 already in the skim

  Bin sb_2b_met0{ "sb_2b_met0",  "hig_bin==21&&met>150&&met<=300", blind_level>=BlindLevel::blinded};
  Bin sig_2b_met0{"sig_2b_met0", "hig_bin==22&&met>150&&met<=300", blind_level>=BlindLevel::blinded};
  Bin sb_3b_met0{ "sb_3b_met0",  "hig_bin==31&&met>150&&met<=300", blind_level>=BlindLevel::blinded};
  Bin sig_3b_met0{"sig_3b_met0", "hig_bin==32&&met>150&&met<=300", blind_level>=BlindLevel::blinded};
  Bin sb_4b_met0{ "sb_4b_met0",  "hig_bin==41&&met>150&&met<=300", blind_level>=BlindLevel::blinded};
  Bin sig_4b_met0{"sig_4b_met0", "hig_bin==42&&met>150&&met<=300", blind_level>=BlindLevel::blinded};

  Bin sb_2b_met1{ "sb_2b_met1",  "hig_bin==21&&met>300", blind_level>=BlindLevel::blinded};
  Bin sig_2b_met1{"sig_2b_met1", "hig_bin==22&&met>300", blind_level>=BlindLevel::blinded};
  Bin sb_3b_met1{ "sb_3b_met1",  "hig_bin==31&&met>300", blind_level>=BlindLevel::blinded};
  Bin sig_3b_met1{"sig_3b_met1", "hig_bin==32&&met>300", blind_level>=BlindLevel::blinded};
  Bin sb_4b_met1{ "sb_4b_met1",  "hig_bin==41&&met>300", blind_level>=BlindLevel::blinded};
  Bin sig_4b_met1{"sig_4b_met1", "hig_bin==42&&met>300", blind_level>=BlindLevel::blinded};

  //// Defining the 2x3 ABCD in bins of met
  set<Block> blocks_abcd;

  blocks_abcd = {
    {"met0", {{sb_2b_met0, sb_3b_met0, sb_4b_met0},
	      {sig_2b_met0, sig_3b_met0, sig_4b_met0}}},
    {"met1", {{sb_2b_met1, sb_3b_met1, sb_4b_met1},
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
  string outname(outfolder+"/wspace_SMS-TChiHH_4b_mChi-400_xsecNom.root");
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
    opt = getopt_long(argc, argv, "l:u:k4g:f:o:", long_options, &option_index);
    if( opt == -1) break;

    string optname;
    switch(opt){
    case 'l':
      lumi = atof(optarg);
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
