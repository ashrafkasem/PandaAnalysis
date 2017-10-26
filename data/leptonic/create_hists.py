#!/usr/bin/env python

from sys import argv
argv = []

import ROOT as root 
from array import array 

'''
double weightEWKWZCorr(float mass){
  double value = 1.0;
  if     (mass<100) value = 1.006;
  else if(mass<160) value = 1.006;
  else if(mass<180) value = 0.998;
  else if(mass<200) value = 0.993;
  else if(mass<220) value = 0.993;
  else if(mass<250) value = 0.993;
  else if(mass<300) value = 0.992;
  else if(mass<350) value = 0.990;
  else if(mass<430) value = 0.988;
  else if(mass<500) value = 0.984;
  else if(mass<600) value = 0.980;
  else              value = 0.975;
  return value;
}
''' 
bins = array('f', [0, 100, 160, 180, 200, 220, 250, 300, 350, 430, 500, 600, 999])
vals = array('f', [1.006, 1.006, 0.998, 0.993, 0.993, 0.993, 0.992, 0.990, 0.988, 0.984, 0.980, 0.975])
h1 = root.TH1D('hEWKWZCorr', '', len(bins)-1, bins)
for ib, val in enumerate(vals):
    h1.SetBinContent(ib+1, val)


'''
float kfactor_qqZZ_qcd_M(float GENmassZZ, int finalState = 2)
{

    // finalState=1 : 4e/4mu/4tau
    // finalState=2 : 2e2mu/2mutau/2e2tau

    float k=0.0;

    if (finalState==1) {
        k+=1.23613311013*(abs(GENmassZZ)>0.0&&abs(GENmassZZ)<=25.0);
        k+=1.17550314639*(abs(GENmassZZ)>25.0&&abs(GENmassZZ)<=50.0);
        k+=1.17044565911*(abs(GENmassZZ)>50.0&&abs(GENmassZZ)<=75.0);
        k+=1.03141209689*(abs(GENmassZZ)>75.0&&abs(GENmassZZ)<=100.0);
        k+=1.05285574912*(abs(GENmassZZ)>100.0&&abs(GENmassZZ)<=125.0);
        k+=1.11287217794*(abs(GENmassZZ)>125.0&&abs(GENmassZZ)<=150.0);
        k+=1.13361441158*(abs(GENmassZZ)>150.0&&abs(GENmassZZ)<=175.0);
        k+=1.10355603327*(abs(GENmassZZ)>175.0&&abs(GENmassZZ)<=200.0);
        k+=1.10053981637*(abs(GENmassZZ)>200.0&&abs(GENmassZZ)<=225.0);
        k+=1.10972676811*(abs(GENmassZZ)>225.0&&abs(GENmassZZ)<=250.0);
        k+=1.12069120525*(abs(GENmassZZ)>250.0&&abs(GENmassZZ)<=275.0);
        k+=1.11589101635*(abs(GENmassZZ)>275.0&&abs(GENmassZZ)<=300.0);
        k+=1.13906170314*(abs(GENmassZZ)>300.0&&abs(GENmassZZ)<=325.0);
        k+=1.14854594271*(abs(GENmassZZ)>325.0&&abs(GENmassZZ)<=350.0);
        k+=1.14616229031*(abs(GENmassZZ)>350.0&&abs(GENmassZZ)<=375.0);
        k+=1.14573157789*(abs(GENmassZZ)>375.0&&abs(GENmassZZ)<=400.0);
        k+=1.13829430515*(abs(GENmassZZ)>400.0&&abs(GENmassZZ)<=425.0);
        k+=1.15521193686*(abs(GENmassZZ)>425.0&&abs(GENmassZZ)<=450.0);
        k+=1.13679822698*(abs(GENmassZZ)>450.0&&abs(GENmassZZ)<=475.0);
        k+=1.13223956942*(abs(GENmassZZ)>475.0);
    }

    if (finalState==2) {
        k+=1.25094466582*(abs(GENmassZZ)>0.0&&abs(GENmassZZ)<=25.0);
        k+=1.22459455362*(abs(GENmassZZ)>25.0&&abs(GENmassZZ)<=50.0);
        k+=1.19287368979*(abs(GENmassZZ)>50.0&&abs(GENmassZZ)<=75.0);
        k+=1.04597506451*(abs(GENmassZZ)>75.0&&abs(GENmassZZ)<=100.0);
        k+=1.08323413771*(abs(GENmassZZ)>100.0&&abs(GENmassZZ)<=125.0);
        k+=1.09994968030*(abs(GENmassZZ)>125.0&&abs(GENmassZZ)<=150.0);
        k+=1.16698455800*(abs(GENmassZZ)>150.0&&abs(GENmassZZ)<=175.0);
        k+=1.10399053155*(abs(GENmassZZ)>175.0&&abs(GENmassZZ)<=200.0);
        k+=1.10592664340*(abs(GENmassZZ)>200.0&&abs(GENmassZZ)<=225.0);
        k+=1.10690381480*(abs(GENmassZZ)>225.0&&abs(GENmassZZ)<=250.0);
        k+=1.11194928918*(abs(GENmassZZ)>250.0&&abs(GENmassZZ)<=275.0);
        k+=1.13522586553*(abs(GENmassZZ)>275.0&&abs(GENmassZZ)<=300.0);
        k+=1.11895090244*(abs(GENmassZZ)>300.0&&abs(GENmassZZ)<=325.0);
        k+=1.13898508615*(abs(GENmassZZ)>325.0&&abs(GENmassZZ)<=350.0);
        k+=1.15463977506*(abs(GENmassZZ)>350.0&&abs(GENmassZZ)<=375.0);
        k+=1.17341664594*(abs(GENmassZZ)>375.0&&abs(GENmassZZ)<=400.0);
        k+=1.20093349763*(abs(GENmassZZ)>400.0&&abs(GENmassZZ)<=425.0);
        k+=1.18915554919*(abs(GENmassZZ)>425.0&&abs(GENmassZZ)<=450.0);
        k+=1.18546007375*(abs(GENmassZZ)>450.0&&abs(GENmassZZ)<=475.0);
        k+=1.12864505708*(abs(GENmassZZ)>475.0);
    }

    if (k==0.0) return 1.1;
    else return k; // if something goes wrong return inclusive k-factor

}
'''

xbins = array('f', [0.5, 1.5, 2.5])
ybins = array('f', [x*25 for x in xrange(20)])
data = [[1.23613311013, 1.17550314639, 1.17044565911, 1.03141209689, 1.05285574912, 1.11287217794, 
         1.13361441158, 1.10355603327, 1.10053981637, 1.10972676811, 1.12069120525, 1.11589101635, 
         1.13906170314, 1.14854594271, 1.14616229031, 1.14573157789, 1.13829430515, 1.15521193686, 
         1.13679822698, 1.13223956942],
        [1.25094466582, 1.22459455362, 1.19287368979, 1.04597506451, 1.08323413771, 1.09994968030, 
         1.16698455800, 1.10399053155, 1.10592664340, 1.10690381480, 1.11194928918, 1.13522586553, 
         1.11895090244, 1.13898508615, 1.15463977506, 1.17341664594, 1.20093349763, 1.18915554919, 
         1.18546007375, 1.12864505708]]
h2 = root.TH2D('hqqZZKfactor', '', len(xbins)-1, xbins, len(ybins)-1, ybins)
for ix, d in enumerate(data):
    for iy, v in enumerate(d):
        h2.SetBinContent(ix+1, iy+1, v)

fout = root.TFile('data/leptonic/data.root', 'RECREATE')
fout.WriteTObject(h1)
fout.WriteTObject(h2)
fout.Close()