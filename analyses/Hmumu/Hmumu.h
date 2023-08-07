#ifndef FCCANALYZER_HMUMU_H
#define FCCANALYZER_HMUMU_H

#include "defines.h"

namespace FCCAnalyses {

Vec_rp FilterObjects(Vec_rp total, Vec_rp remove)
    { Vec_rp result;
     int ok = 0; 
     for ( int i = 0; i <= total.size(); i++)
     { ok = 0 ;
      auto &p = total[i];
      TLorentzVector p1;
    p1.SetXYZM(total[i].momentum.x, total[i].momentum.y, total[i].momentum.z, total[i].mass);
      for (int j = 0; j<= remove.size(); j++)
         { TLorentzVector p2;
    p2.SetXYZM(remove[j].momentum.x, remove[j].momentum.y, remove[j].momentum.z, remove[j].mass);
             if ( p1.DeltaR(p2) <= 0.001 )
             { ok = 1;
             break; }
         }
      if (ok ==0)
      {
          result.push_back(p);
             }
         }
     return result;
    }

 
float inv_mass(Vec_rp in){
    Vec_tlv tlv = makeLorentzVectors(in);
    TLorentzVector sum;
    sum.SetXYZM(0, 0, 0, 0);
    for (auto & particle: tlv){
        sum = sum + particle;
    }
    return sum.M();
}

    
}

#endif