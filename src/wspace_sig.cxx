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
#include "cross_sections.hpp"


using namespace std;

namespace{
  double lumi = 35.9;
  double sig_strength = 0.;
  BlindLevel blind_level = BlindLevel::unblinded;
  bool no_kappa = true;
  bool do_syst = true;
  bool use_r4 = true;
  bool deep = true;
  unsigned n_toys = 0;
  string outfolder = "";
  string nb_bins("TTML");
  string sigfile = "/net/cms27/cms27r0/babymaker/babies/2017_02_26/TChiHH/merged_higmc_higtight/mergedbaby__SMS-TChiHH_mGluino-400_mLSP-1_Tune_skim_higtight_higmc_nfiles_1.root";
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
  string basefolder("/net/cms29/cms29r0/babymaker/");
  if(Contains(hostname, "lxplus")) basefolder = "/afs/cern.ch/user/m/manuelf/work/";

  string foldermc(basefolder+"babies/2017_01_27/mc/merged_higmc_higtight/"); 
  string folderdata(basefolder+"babies/2017_02_14/data/merged_higdata_higloose/"); // Pointing to a random folder for now

  
  //Define processes. Try to minimize splitting
  string stitch_cuts("stitch_met&&pass");
  Process ttbar{"ttbar", {
      {foldermc+"/*_TTJets*Lep*.root/tree"}
    },stitch_cuts};

  Process other{"other", {
      {foldermc+"/*_TTW*.root/tree"},
	{foldermc+"/*_TTZ*.root/tree"},
	  {foldermc+"/*_TTGJets*.root/tree"},
	    {foldermc+"/*ttHJetTobb*.root/tree"},
	      {foldermc+"/*_TTTT*.root/tree"},
		{foldermc+"/*_ZJet*.root/tree"},
		  {foldermc+"/*_WJetsToLNu*.root/tree"},
		    {foldermc+"/*DYJetsToLL*.root/tree"},
		      {foldermc+"/*_ST_*.root/tree"},
			{foldermc+"/*QCD_HT*0_Tune*.root/tree"},
			  {foldermc+"/*QCD_HT*Inf_Tune*.root/tree"},
			    {foldermc+"/*_WH_HToBB*.root/tree"},
			      {foldermc+"/*_ZH_HToBB*.root/tree"},
				{foldermc+"/*_WWTo*.root/tree"},
				  {foldermc+"/*_WZ*.root/tree"},
				    {foldermc+"/*_ZZ*.root/tree"}
    },stitch_cuts};
    
  Process signal{"signal", {
      {sigfile+"/tree"}
    },"pass_goodv&&pass_ecaldeadcell&&pass_hbhe&&pass_hbheiso&&pass_fsmet", false, true};

  string data_cuts("(trig[13]||trig[33]||trig[14]||trig[15]||trig[30]||trig[31])&&pass");


  Process data{"data", {
      {folderdata+"/*.root/tree"}
    }, data_cuts, true};

  //Make list of all backgrounds. Backgrounds assumed to be orthogonal
  set<Process> backgrounds{ttbar, other};

  //Baseline selection applied to all bins and processes
  string basestr = "hig_drmax<2.2&&ntks==0&&njets>=4&&njets<=5&&!low_dphi&&nvleps==0&&pass_ra2_badmu&&met/met_calo<5";
  if (deep) basestr = "higd_drmax<2.2&&ntks==0&&njets>=4&&njets<=5&&!low_dphi&&nvleps==0&&pass_ra2_badmu&&met/met_calo<5";
  Cut baseline{basestr}; 

  string cut2b="nbt==2&&nbm==2";
  string cut3b="nbt>=2&&nbm==3&&nbl==3";
  string cut4b="nbt>=2&&nbm>=3&&nbl>=4";
  if (deep) {
    cut2b="nbdt==2&&nbdm==2";
    cut3b="nbdt>=2&&nbdm==3&&nbdl==3";
    cut4b="nbdt>=2&&nbdm>=3&&nbdl>=4";
  }

  string cutmet0="&&met>150&&met<=200", cutmet1="&&met>200&&met<=300", cutmet2="&&met>300&&met<=450";
  string cutmet3="&&met>450";

  //cutmet2="&&met>300";
  string cuthig="hig_am>100&&hig_am<140&&hig_dm<40";
  string cutsbd="!(hig_am>100&&hig_am<140)&&hig_dm<40&&hig_am<200";
  if (deep) {
    cuthig="higd_am>100&&higd_am<140&&higd_dm<40";
    cutsbd="!(higd_am>100&&higd_am<140)&&higd_dm<40&&higd_am<200";
  }

  Bin sbd_2b_met0{"sbd_2b_met0", cut2b+"&&"+cutsbd+cutmet0, blind_level>=BlindLevel::blinded};
  Bin hig_2b_met0{"hig_2b_met0", cut2b+"&&"+cuthig+cutmet0, blind_level>=BlindLevel::blinded};
  Bin sbd_3b_met0{"sbd_3b_met0", cut3b+"&&"+cutsbd+cutmet0, blind_level>=BlindLevel::blinded};
  Bin hig_3b_met0{"hig_3b_met0", cut3b+"&&"+cuthig+cutmet0, blind_level>=BlindLevel::blinded};
  Bin sbd_4b_met0{"sbd_4b_met0", cut4b+"&&"+cutsbd+cutmet0, blind_level>=BlindLevel::blinded};
  Bin hig_4b_met0{"hig_4b_met0", cut4b+"&&"+cuthig+cutmet0, blind_level>=BlindLevel::blinded};

  Bin sbd_2b_met1{"sbd_2b_met1", cut2b+"&&"+cutsbd+cutmet1, blind_level>=BlindLevel::blinded};
  Bin hig_2b_met1{"hig_2b_met1", cut2b+"&&"+cuthig+cutmet1, blind_level>=BlindLevel::blinded};
  Bin sbd_3b_met1{"sbd_3b_met1", cut3b+"&&"+cutsbd+cutmet1, blind_level>=BlindLevel::blinded};
  Bin hig_3b_met1{"hig_3b_met1", cut3b+"&&"+cuthig+cutmet1, blind_level>=BlindLevel::blinded};
  Bin sbd_4b_met1{"sbd_4b_met1", cut4b+"&&"+cutsbd+cutmet1, blind_level>=BlindLevel::blinded};
  Bin hig_4b_met1{"hig_4b_met1", cut4b+"&&"+cuthig+cutmet1, blind_level>=BlindLevel::blinded};

  Bin sbd_2b_met2{"sbd_2b_met2", cut2b+"&&"+cutsbd+cutmet2, blind_level>=BlindLevel::blinded};
  Bin hig_2b_met2{"hig_2b_met2", cut2b+"&&"+cuthig+cutmet2, blind_level>=BlindLevel::blinded};
  Bin sbd_3b_met2{"sbd_3b_met2", cut3b+"&&"+cutsbd+cutmet2, blind_level>=BlindLevel::blinded};
  Bin hig_3b_met2{"hig_3b_met2", cut3b+"&&"+cuthig+cutmet2, blind_level>=BlindLevel::blinded};
  Bin sbd_4b_met2{"sbd_4b_met2", cut4b+"&&"+cutsbd+cutmet2, blind_level>=BlindLevel::blinded};
  Bin hig_4b_met2{"hig_4b_met2", cut4b+"&&"+cuthig+cutmet2, blind_level>=BlindLevel::blinded};

  Bin sbd_2b_met3{"sbd_2b_met3", cut2b+"&&"+cutsbd+cutmet3, blind_level>=BlindLevel::blinded};
  Bin hig_2b_met3{"hig_2b_met3", cut2b+"&&"+cuthig+cutmet3, blind_level>=BlindLevel::blinded};
  Bin sbd_3b_met3{"sbd_3b_met3", cut3b+"&&"+cutsbd+cutmet3, blind_level>=BlindLevel::blinded};
  Bin hig_3b_met3{"hig_3b_met3", cut3b+"&&"+cuthig+cutmet3, blind_level>=BlindLevel::blinded};
  Bin sbd_4b_met3{"sbd_4b_met3", cut4b+"&&"+cutsbd+cutmet3, blind_level>=BlindLevel::blinded};
  Bin hig_4b_met3{"hig_4b_met3", cut4b+"&&"+cuthig+cutmet3, blind_level>=BlindLevel::blinded};

  //// Defining the 2x3 ABCD in bins of met
  set<Block> blocks_abcd;

  blocks_abcd = {
    {"met0", {{sbd_2b_met0, sbd_3b_met0, sbd_4b_met0},
	      {hig_2b_met0, hig_3b_met0, hig_4b_met0}}},
    {"met1", {{sbd_2b_met1, sbd_3b_met1, sbd_4b_met1},
	      {hig_2b_met1, hig_3b_met1, hig_4b_met1}}},
    {"met2", {{sbd_2b_met2, sbd_3b_met2, sbd_4b_met2},
	      {hig_2b_met2, hig_3b_met2, hig_4b_met2}}},
    {"met3", {{sbd_2b_met3, sbd_3b_met3, sbd_4b_met3},
    	      {hig_2b_met3, hig_3b_met3, hig_4b_met3}}}
  };

  //// Parsing the gluino and LSP masses
  int mglu, mlsp;
  parseMasses(sigfile, mglu, mlsp);
  string glu_lsp("mGluino-"+to_string(mglu)+"_mLSP-"+to_string(mlsp));

  cout<<"Creating workpaces"<<endl;
  //// Creating workspaces
  Cut *pbaseline(&baseline);
  set<Block> *pblocks(&blocks_abcd);

  string sysfolder = "/net/cms27/cms27r0/babymaker/sys/2017_02_26/TChiHH/";
  string sysfile(sysfolder+"/sys_SMS-TChiHH_mGluino-"+to_string(mglu)+"_mLSP-1_35p9ifb.txt");
  
  // If systematic file does not exist, complain
  struct stat buffer;   
  if(stat (sysfile.c_str(), &buffer) != 0) {
    cout<<endl<<"WARNING: "<<sysfile<<" does not exist. Using ";
    sysfile = "txt/systematics/sys_TChiHH.txt";
    cout<<sysfile<<" instead"<<endl<<endl;
  }

  // Cross sections
  float xsec, xsec_unc;
  xsec::higgsinoCrossSection(mglu, xsec, xsec_unc);

  gSystem->mkdir(outfolder.c_str(), kTRUE);

  int digits=0;
  if(lumi-floor(lumi)>0) digits = 1;
  TString lumi_s = "_lumi"+RoundNumber(lumi,digits);  lumi_s.ReplaceAll(".","p");
  digits=0;
  if(sig_strength-floor(sig_strength)>0) digits = 1;
  TString sig_s = "_sig"+RoundNumber(sig_strength,digits); sig_s.ReplaceAll(".","p");
  TString outname(outfolder+"wspace_TChiHH_"+glu_lsp+"_xsecNom_nb"+nb_bins+lumi_s+sig_s+".root");
  if(!use_r4) outname.ReplaceAll("wspace_","wspace_nor4_");

  float rmax = 20.;
  WorkspaceGenerator wgNom(*pbaseline, *pblocks, backgrounds, signal, data, sysfile, use_r4, sig_strength, 1.);
  wgNom.SetRMax(rmax);
  wgNom.SetKappaCorrected(!no_kappa);
  wgNom.SetDoSystematics(do_syst);
  wgNom.SetLuminosity(lumi);
  wgNom.SetDoSystematics(do_syst);
  wgNom.AddToys(n_toys);
  wgNom.WriteToFile(outname.Data());

  outname.ReplaceAll("Nom", "Up");
  WorkspaceGenerator wgUp(*pbaseline, *pblocks, backgrounds, signal, data, sysfile, use_r4, sig_strength, 1+xsec_unc);
  wgUp.SetRMax(rmax);
  wgUp.SetKappaCorrected(!no_kappa);
  wgUp.SetDoSystematics(do_syst);
  wgUp.SetLuminosity(lumi);
  wgUp.SetDoSystematics(do_syst);
  wgUp.AddToys(n_toys);
  wgUp.WriteToFile(outname.Data());

  outname.ReplaceAll("Up", "Down");
  WorkspaceGenerator wgDown(*pbaseline, *pblocks, backgrounds, signal, data, sysfile, use_r4, sig_strength, 1-xsec_unc);
  wgDown.SetRMax(rmax);
  wgDown.SetKappaCorrected(!no_kappa);
  wgDown.SetDoSystematics(do_syst);
  wgDown.SetLuminosity(lumi);
  wgDown.SetDoSystematics(do_syst);
  wgDown.AddToys(n_toys);
  wgDown.WriteToFile(outname.Data());

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
      {"csv", no_argument, 0, 0},
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
      }else if(optname == "csv"){
        deep = false;
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
