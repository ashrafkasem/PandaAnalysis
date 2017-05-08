from PandaCore.Tools.Misc import *
from re import sub


metTrigger='(trigger&1)!=0'
eleTrigger='(trigger&2)!=0'
phoTrigger='(trigger&4)!=0'


cuts = {}
weights = {}
triggers = {
    'met':'(trigger&1)!=0',
    'ele':'(trigger&2)!=0',
    'pho':'(trigger&4)!=0',
}

#baseline = 'metFilter==1 && nFatjet==1 && fj1Pt>200 && nTau==0 && fj1ECFN_2_3_20/pow(fj1ECFN_1_2_20,2.00)<0.113'
#baseline = 'metFilter==1 && nFatjet==1 && fj1Pt>200 && nTau==0 && N2DDT<0.0 && Sum$(jetPt>30 && jetIso)<2'
baseline = 'metFilter==1 && nFatjet==1 && fj1Pt>200 && nTau==0 && Sum$(jetPt>30 && jetIso)<2 && N2DDT<0'
#baseline = 'metFilter==1 && nFatjet==1 && fj1Pt>200 && nTau==0 && Sum$(jetPt>30 && jetIso)<2 && fj1MSD_corr>100 && fj1MSD_corr<150'

#cuts for specific regions
cuts['signal'] = tAND(baseline,'nLooseLep==0 && nLooseElectron==0 && nLoosePhoton==0 && puppimet>200 && isojetNBtags==0 && fj1MSD_corr>100 && fj1MSD_corr<150 && dphipuppimet>0.4 && fj1DoubleCSV>0.75')
cuts['tm'] = tAND(baseline,'(nLooseElectron+nLoosePhoton+nTau)==0 && nLooseMuon==1 && looseLep1IsTight==1 && puppiUWmag>200 && isojetNBtags==1 && fj1MSD_corr>100 && fj1MSD_corr<150 && dphipuppiUW>0.4 && fj1DoubleCSV>0.75')
cuts['te'] = tAND(baseline,'(nLooseMuon+nLoosePhoton+nTau)==0 && nLooseElectron==1 && looseLep1IsTight==1 && puppiUWmag>200 && isojetNBtags==1 && fj1MSD_corr>100 && fj1MSD_corr<150 && dphipuppiUW>0.4 && puppimet>50 && fj1DoubleCSV>0.75')
cuts['wmn'] = tAND(baseline,'(nLooseElectron+nLoosePhoton+nTau)==0 && nLooseMuon==1 && looseLep1IsTight==1 && puppiUWmag>200 && isojetNBtags==0 && fj1MSD_corr>50 && fj1MSD_corr<150 && dphipuppiUW>0.4 && fj1DoubleCSV>0.75')
cuts['wen'] = tAND(baseline,'(nLooseMuon+nLoosePhoton+nTau)==0 && nLooseElectron==1 && looseLep1IsTight==1 && puppiUWmag>200 && isojetNBtags==0 && fj1MSD_corr>50 && fj1MSD_corr<150 && dphipuppiUW>0.4 && puppimet>50 && fj1DoubleCSV>0.75')
cuts['wmn_fail'] = tAND(baseline,'(nLooseElectron+nLoosePhoton+nTau)==0 && nLooseMuon==1 && looseLep1IsTight==1 && puppiUWmag>200 && isojetNBtags==0 && fj1MSD_corr>50 && fj1MSD_corr<150 && dphipuppiUW>0.4 && fj1DoubleCSV<0.75')
cuts['wen_fail'] = tAND(baseline,'(nLooseMuon+nLoosePhoton+nTau)==0 && nLooseElectron==1 && looseLep1IsTight==1 && puppiUWmag>200 && isojetNBtags==0 && fj1MSD_corr>50 && fj1MSD_corr<150 && dphipuppiUW>0.4 && puppimet>50 && fj1DoubleCSV<0.75')
cuts['zmm'] = tAND(baseline,'(nLooseElectron+nLoosePhoton+nTau)==0 && nLooseMuon==2 && looseLep1IsTight==1 && puppiUZmag>200 && isojetNBtags==0 && fj1MSD_corr>50 && fj1MSD_corr<250 && dphipuppiUZ>0.4 && diLepMass>60 && diLepMass<120 && fj1DoubleCSV>0.75')
cuts['zee'] = tAND(baseline,'(nLooseMuon+nLoosePhoton+nTau)==0 && nLooseElectron==2 && looseLep1IsTight==1 && puppiUZmag>200 && isojetNBtags==0 && fj1MSD_corr>50 && fj1MSD_corr<250 && dphipuppiUZ>0.4 && diLepMass>60 && diLepMass<120 && fj1DoubleCSV>0.75')
cuts['zmm_fail'] = tAND(baseline,'(nLooseElectron+nLoosePhoton+nTau)==0 && nLooseMuon==2 && looseLep1IsTight==1 && puppiUZmag>200 && isojetNBtags==0 && fj1MSD_corr>50 && fj1MSD_corr<250 && dphipuppiUZ>0.4 && diLepMass>60 && diLepMass<120 && fj1DoubleCSV<0.75')
cuts['zee_fail'] = tAND(baseline,'(nLooseMuon+nLoosePhoton+nTau)==0 && nLooseElectron==2 && looseLep1IsTight==1 && puppiUZmag>200 && isojetNBtags==0 && fj1MSD_corr>50 && fj1MSD_corr<250 && dphipuppiUZ>0.4 && diLepMass>60 && diLepMass<120 && fj1DoubleCSV<0.75')
cuts['pho'] = tAND(baseline,'(nLooseMuon+nLooseElectron+nTau)==0 && nLoosePhoton==1 && loosePho1IsTight==1 && puppiUAmag>200 && isojetNBtags==0 && fj1MSD_corr>50 && fj1MSD_corr<250 && dphipuppiUA>0.4 && fj1DoubleCSV>0.75')

#weights for specific regions
weights['base'] = 'normalizedWeight*sf_pu*sf_ewkV*sf_qcdV'
weights['signal'] = tTIMES(weights['base'],'%f*sf_metTrig*sf_tt*sf_btag0')
weights['tm'] = tTIMES(weights['base'],'%f*sf_metTrig*sf_lepID*sf_lepIso*sf_lepTrack*sf_tt*sf_btag1')
weights['te'] = tTIMES(weights['base'],'%f*sf_eleTrig*sf_lepID*sf_lepTrack*sf_tt*sf_btag1')
weights['wmn'] = tTIMES(weights['base'],'%f*sf_metTrig*sf_lepID*sf_lepIso*sf_lepTrack*sf_tt*sf_btag0')
weights['wmn_fail'] = tTIMES(weights['base'],'%f*sf_metTrig*sf_lepID*sf_lepIso*sf_lepTrack*sf_tt*sf_btag0')
weights['wen'] = tTIMES(weights['base'],'%f*sf_eleTrig*sf_lepID*sf_lepIso*sf_lepTrack*sf_tt*sf_btag0')
weights['wen_fail'] = tTIMES(weights['base'],'%f*sf_eleTrig*sf_lepID*sf_lepIso*sf_lepTrack*sf_tt*sf_btag0')
weights['zmm'] = tTIMES(weights['base'],'%f*sf_metTrig*sf_lepID*sf_lepIso*sf_lepTrack*sf_tt*sf_btag0')
weights['zee'] = tTIMES(weights['base'],'%f*sf_eleTrig*sf_lepID*sf_lepIso*sf_lepTrack*sf_tt*sf_btag0')
weights['zmm_fail'] = tTIMES(weights['base'],'%f*sf_metTrig*sf_lepID*sf_lepIso*sf_lepTrack*sf_tt*sf_btag0')
weights['zee_fail'] = tTIMES(weights['base'],'%f*sf_eleTrig*sf_lepID*sf_lepIso*sf_lepTrack*sf_tt*sf_btag0')
weights['pho'] = tTIMES(weights['base'],'%f*sf_phoTrig*sf_btag0*1.0')


for r in ['signal','wmn','wen','tm','te','zmm','zee','zmm_fail','zee_fail','wmn_fail','wen_fail','pho']:
#    print r
    for shift in ['BUp','BDown','MUp','MDown']:
        for cent in ['sf_btag']:
            if 'btag0' in weights[r]:
                weights[r+'_'+cent+shift] = sub(cent+'0',cent+'0'+shift,weights[r])
#                print weights[r+'_'+cent+shift]
            if 'btag1' in weights[r]:
                weights[r+'_'+cent+shift] = sub(cent+'1',cent+'1'+shift,weights[r])
#                print weights[r+'_'+cent+shift]


for r in ['signal','wmn','wen','tm','te','zmm','zee','pho']:
    for shift in ['ScaleUp','ScaleDown','PDFUp','PDFDown']:
        if shift == 'ScaleUp':
            weights[r+'_'+shift] = weights[r] + "*(scale[3]+1)"
        if shift == 'ScaleDown':
            weights[r+'_'+shift] = weights[r] + "*(scale[5]+1)"
        if shift == 'PDFUp':
            weights[r+'_'+shift] = weights[r] + "*pdfUp"
        if shift == 'PDFDown':
            weights[r+'_'+shift] = weights[r] + "*pdfDown"
'''
   
for r in ['signal','wmn','wen','tm','te','zmm','zee','pho','zmm_fail','zee_fail']:
    for proc in ['','ttbar','vh','wjets','zll','zvv']:
        for shift in ['ScaleUp','ScaleDown','PDFUp','PDFDown']:
            if shift == 'ScaleUp':
                weights[r+'_'+proc+shift] = weights[r] + "*(scale[3]+1)"
            if shift == 'ScaleDown':
                weights[r+'_'+proc+shift] = weights[r] + "*(scale[5]+1)"
            if shift == 'PDFUp':
                weights[r+'_'+proc+shift] = weights[r] + "*pdfUp"
            if shift == 'PDFDown':
                weights[r+'_'+proc+shift] = weights[r] + "*pdfDown"

'''
     
            
weights['zmm_fail_PDFUp'] = tTIMES(weights['base'],'%f*sf_metTrig*sf_lepID*sf_lepIso*sf_lepTrack*sf_tt*sf_btag0*pdfUp')
weights['zmm_fail_PDFDown'] = tTIMES(weights['base'],'%f*sf_metTrig*sf_lepID*sf_lepIso*sf_lepTrack*sf_tt*sf_btag0*pdfDown')
weights['zmm_fail_ScaleUp'] = tTIMES(weights['base'],'%f*sf_metTrig*sf_lepID*sf_lepIso*sf_lepTrack*sf_tt*sf_btag0*(scale[3]+1)')
weights['zmm_fail_ScaleDown'] = tTIMES(weights['base'],'%f*sf_metTrig*sf_lepID*sf_lepIso*sf_lepTrack*sf_tt*sf_btag0*(scale[5]+1)')
weights['zee_fail_PDFUp'] = tTIMES(weights['base'],'%f*sf_eleTrig*sf_lepID*sf_lepIso*sf_lepTrack*sf_tt*sf_btag0*pdfUp')
weights['zee_fail_PDFDown'] = tTIMES(weights['base'],'%f*sf_eleTrig*sf_lepID*sf_lepIso*sf_lepTrack*sf_tt*sf_btag0*pdfDown')
weights['zee_fail_ScaleUp'] = tTIMES(weights['base'],'%f*sf_eleTrig*sf_lepID*sf_lepIso*sf_lepTrack*sf_tt*sf_btag0*(scale[3]+1)')
weights['zee_fail_ScaleDown'] = tTIMES(weights['base'],'%f*sf_eleTrig*sf_lepID*sf_lepIso*sf_lepTrack*sf_tt*sf_btag0*(scale[5]+1)')
weights['wmn_fail_PDFUp'] = tTIMES(weights['base'],'%f*sf_metTrig*sf_lepID*sf_lepIso*sf_lepTrack*sf_tt*sf_btag0*pdfUp')
weights['wmn_fail_PDFDown'] = tTIMES(weights['base'],'%f*sf_metTrig*sf_lepID*sf_lepIso*sf_lepTrack*sf_tt*sf_btag0*pdfDown')
weights['wmn_fail_ScaleUp'] = tTIMES(weights['base'],'%f*sf_metTrig*sf_lepID*sf_lepIso*sf_lepTrack*sf_tt*sf_btag0*(scale[3]+1)')
weights['wmn_fail_ScaleDown'] = tTIMES(weights['base'],'%f*sf_metTrig*sf_lepID*sf_lepIso*sf_lepTrack*sf_tt*sf_btag0*(scale[5]+1)')
weights['wen_fail_PDFUp'] = tTIMES(weights['base'],'%f*sf_eleTrig*sf_lepID*sf_lepIso*sf_lepTrack*sf_tt*sf_btag0*pdfUp')
weights['wen_fail_PDFDown'] = tTIMES(weights['base'],'%f*sf_eleTrig*sf_lepID*sf_lepIso*sf_lepTrack*sf_tt*sf_btag0*pdfDown')
weights['wen_fail_ScaleUp'] = tTIMES(weights['base'],'%f*sf_eleTrig*sf_lepID*sf_lepIso*sf_lepTrack*sf_tt*sf_btag0*(scale[3]+1)')
weights['wen_fail_ScaleDown'] = tTIMES(weights['base'],'%f*sf_eleTrig*sf_lepID*sf_lepIso*sf_lepTrack*sf_tt*sf_btag0*(scale[5]+1)')


