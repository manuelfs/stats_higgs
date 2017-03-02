#include "yield_manager.hpp"

#include <iostream>
#include <sstream>
#include <array>

#include "bin.hpp"
#include "process.hpp"
#include "cut.hpp"
#include "utilities.hpp"

using namespace std;

map<YieldKey, GammaParams> YieldManager::yields_ = map<YieldKey, GammaParams>();
const double YieldManager::store_lumi_ = 4.;

YieldManager::YieldManager(double lumi):
  local_lumi_(lumi),
  verbose_(false){
}

GammaParams YieldManager::GetYield(const YieldKey &key) const{
  if(!HaveYield(key)) ComputeYield(key);

  double factor = local_lumi_/store_lumi_;
  if(GetProcess(key).IsData()) factor = 1.;

  return factor*yields_.at(key);
}

GammaParams YieldManager::GetYield(const Bin &bin,
                                   const Process &process,
                                   const Cut &cut) const{
  return GetYield(YieldKey(bin, process, cut));
}

const double & YieldManager::Luminosity() const{
  return local_lumi_;
}

double & YieldManager::Luminosity(){
  return local_lumi_;
}

bool YieldManager::HaveYield(const YieldKey &key) const{
  return yields_.find(key) != yields_.end();
}

void YieldManager::ComputeYield(const YieldKey &key) const{
  const Bin &bin = GetBin(key);
  const Process &process = GetProcess(key);
  const Cut &cut = GetCut(key);

  GammaParams gps;

  if(HaveYield(key)){
    if(verbose_){
      cout << "Using known yield for " << key << endl;
    }
    gps = GetYield(key);
  }else if(process.GetEntries() == 0){
    if(verbose_){
      cout << "No entries found for " << key << endl;
    }
    gps.SetNEffectiveAndWeight(0., 0.);
  }else{
    if(verbose_){
      cout << "Computing yield for " << key << endl;
    }
    ostringstream oss;
    oss << local_lumi_ << flush;
    Cut lumi_weight;
    
    //string totweight = oss.str()+"*weight*((type>=7000&&type<8000)*((met> 100 && met<= 125)*0.092+(met> 125 && met<= 150)*0.217+(met> 150 && met<= 175)*0.383+(met> 175 && met<= 200)*0.535+(met> 200 && met<= 225)*0.642+(met> 225 && met<= 250)*0.711+(met> 250 && met<= 275)*0.758+(met> 275 && met<= 300)*0.779+(met> 300 && met<=9999)*0.815) + (type<7000||type>=8000)*((met> 100 && met<= 125)*0.121+(met> 125 && met<= 150)*0.331+(met> 150 && met<= 175)*0.602+(met> 175 && met<= 200)*0.798+(met> 200 && met<= 225)*0.898+(met> 225 && met<= 250)*0.943+(met> 250 && met<= 275)*0.966+(met> 275 && met<= 300)*0.975+(met> 300 && met<=9999)*0.985))";

    string totweight = oss.str()+"*weight/w_btag*w_bhig_deep*eff_trig";



    //totweight = oss.str()+"*weight";

    if(local_lumi_ > 3){ lumi_weight = process.IsData() ? Cut() : 
	(Contains(process.Name(), "sig")?Cut(totweight):Cut(totweight));}

    else{ lumi_weight = process.IsData() ? Cut() : 
	(Contains(process.Name(), "sig")?Cut(totweight):Cut(totweight));}


    array<Cut, 5> cuts;
    cuts.at(0) = lumi_weight*(cut && bin.Cut() && process.Cut());
    cuts.at(1) = lumi_weight*(cut && process.Cut());
    cuts.at(2) = lumi_weight*(process.Cut());
    cuts.at(3) = lumi_weight;
    cuts.at(4) = Cut();

    for(size_t icut = 0; icut < cuts.size() && gps.Weight()<=0.; ++icut){
      if(icut > 0 && !process.CountZeros()){
        gps.SetNEffectiveAndWeight(0., 0.);
        break;
      }
      Cut &this_cut = cuts.at(icut);
      if(verbose_){
        cout << "Trying cut " << this_cut << endl;
      }
      GammaParams temp_gps = process.GetYield(this_cut);
      //// Averaging signal yields cutting on met and met_tru, as prescripted by SUSY group
      //// https://twiki.cern.ch/twiki/bin/viewauth/CMS/SUSRecommendationsICHEP16#Special_treatment_of_MET_uncerta
      if(Contains(process.Name(), "sig")){
        string mettru_s = this_cut.GetCut(); 
        ReplaceAll(mettru_s, "met_calo", "XXXYYYZZZ_calo");
        ReplaceAll(mettru_s, "pass_fsmet", "pass_fsZZBBBZZ");
        ReplaceAll(mettru_s, "met", "met_tru");
        ReplaceAll(mettru_s, "XXXYYYZZZ_calo", "met_calo");
        ReplaceAll(mettru_s, "pass_fsZZBBBZZ", "pass_fsmet");
        Cut mettru_cut(mettru_s); 
        GammaParams mettru_gps = process.GetYield(mettru_cut);
        if(verbose_) cout<<"Yields: met "<<temp_gps.Yield()<<", met_tru "<<mettru_gps.Yield();
        temp_gps.SetYieldAndUncertainty(0.5*(temp_gps.Yield()+mettru_gps.Yield()),
          max(temp_gps.Uncertainty(), mettru_gps.Uncertainty()));
        // temp_gps.SetYieldAndUncertainty(mettru_gps.Yield(), mettru_gps.Uncertainty());
        if(verbose_) cout<<", average "<<temp_gps.Yield()<<" for bin "<<bin.Name()<<endl;
      } // If it is signal
      if(icut == 0) gps = temp_gps;
      else gps.SetNEffectiveAndWeight(0., temp_gps.Weight());
    }
  }

  if(verbose_){
    cout << "Found yield=" << gps << '\n' << endl;
  }
  double factor = store_lumi_/local_lumi_;
  if(process.IsData()) factor = 1.;
  yields_[key] = factor*gps;
}
