#include "../interface/GeneralTree.h"

#define NJET 20
#define NSUBJET 2

GeneralTree::GeneralTree() {
//STARTCUSTOMCONST
  for (unsigned iS=0; iS!=6; ++iS) {
    scale[iS] = 1;
  }
  
  for (auto ibeta : ibetas) {
    for (auto N : Ns) {
      for (auto order : orders) {
        ECFParams p;
        p.ibeta = ibeta;
        p.N = N;
        p.order = order;
        ecfParams.push_back(p);
        fj1ECFNs[p] = -1;
      }
    }
  }

  for (unsigned iShift=0; iShift!=bNShift; ++iShift) {
    for (unsigned iJet=0; iJet!=bNJet; ++iJet) {
      for (unsigned iTags=0; iTags!=bNTags; ++iTags) {
        BTagParams p;
        p.jet = (BTagJet)iJet;
        p.tag = (BTagTags)iTags;
        p.shift = (BTagShift)iShift;
        btagParams.push_back(p);
        sf_btags[p] = 1;
      }
    }
  }

  for (unsigned int iSJ=0; iSJ!=NSUBJET; ++iSJ) {
    fj1sjPt[iSJ] = -1;
    fj1sjEta[iSJ] = -1;
    fj1sjPhi[iSJ] = -1;
    fj1sjM[iSJ] = -1;
    fj1sjCSV[iSJ] = -1;
    fj1sjQGL[iSJ] = -1;
  }
  for (unsigned int iJ=0; iJ!=NJET; ++iJ) {
    jetPt[iJ] = -1;
    jetEta[iJ] = -1;
    jetPhi[iJ] = -1;
    jetE[iJ] = -1;
    jetCSV[iJ] = -1;
    jetIso[iJ] = -1;
    jetQGL[iJ] = -1;
  }

//ENDCUSTOMCONST
}

GeneralTree::~GeneralTree() {
//STARTCUSTOMDEST
//ENDCUSTOMDEST
}

void GeneralTree::Reset() {
//STARTCUSTOMRESET
  for (unsigned iS=0; iS!=6; ++iS) {
    scale[iS] = 1;
  }

  for (auto p : ecfParams) { 
    fj1ECFNs[p] = -1;
  }

  for (auto p : btagParams) { 
    sf_btags[p] = 1;
  }

  for (unsigned int iSJ=0; iSJ!=NSUBJET; ++iSJ) {
    fj1sjPt[iSJ] = -99;
    fj1sjEta[iSJ] = -99;
    fj1sjPhi[iSJ] = -99;
    fj1sjM[iSJ] = -99;
    fj1sjCSV[iSJ] = -99;
    fj1sjQGL[iSJ] = -99;
  }
  for (unsigned int iJ=0; iJ!=NJET; ++iJ) {
    jetPt[iJ] = -99;
    jetEta[iJ] = -99;
    jetPhi[iJ] = -99;
    jetE[iJ] = -99;
    jetCSV[iJ] = -99;
    jetIso[iJ] = -99;
    jetQGL[iJ] = -99;
  }

//ENDCUSTOMRESET
    fj1SubMaxCSV = -1;
    looseLep1IsHLTSafe = 0;
    looseLep2IsHLTSafe = 0;
    runNumber = 0;
    lumiNumber = 0;
    eventNumber = 0;
    npv = 0;
    pu = 0;
    mcWeight = -1;
    trigger = 0;
    metFilter = 0;
    egmFilter = 0;
    filter_maxRecoil = -1;
    filter_whichRecoil = -1;
    sf_ewkV = 1;
    sf_qcdV = 1;
    sf_ewkV2j = 1;
    sf_qcdV2j = 1;
    sf_qcdTT = 1;
    sf_lepID = 1;
    sf_lepIso = 1;
    sf_lepTrack = 1;
    sf_pho = 1;
    sf_eleTrig = 1;
    sf_phoTrig = 1;
    sf_metTrig = 1;
    sf_pu = 1;
    sf_npv = 1;
    sf_tt = 1;
    sf_tt_ext = 1;
    sf_tt_bound = 1;
    sf_tt8TeV = 1;
    sf_tt8TeV_ext = 1;
    sf_tt8TeV_bound = 1;
    sf_phoPurity = 1;
    pfmet = -1;
    pfmetphi = -1;
    pfmetnomu = -1;
    puppimet = -1;
    puppimetphi = -1;
    calomet = -1;
    calometphi = -1;
    pfcalobalance = -1;
    sumET = -1;
    trkmet = -1;
    puppiUWmag = -1;
    puppiUWphi = -1;
    puppiUZmag = -1;
    puppiUZphi = -1;
    puppiUAmag = -1;
    puppiUAphi = -1;
    puppiUperp = -1;
    puppiUpara = -1;
    puppiUmag = -1;
    puppiUphi = -1;
    pfUWmag = -1;
    pfUWphi = -1;
    pfUZmag = -1;
    pfUZphi = -1;
    pfUAmag = -1;
    pfUAphi = -1;
    pfUperp = -1;
    pfUpara = -1;
    pfUmag = -1;
    pfUphi = -1;
    dphipfmet = -1;
    dphipuppimet = -1;
    dphipuppiUW = -1;
    dphipuppiUZ = -1;
    dphipuppiUA = -1;
    dphipfUW = -1;
    dphipfUZ = -1;
    dphipfUA = -1;
    dphipuppiU = -1;
    dphipfU = -1;
    trueGenBosonPt = -1;
    genBosonPt = -1;
    genBosonEta = -1;
    genBosonMass = -1;
    genBosonPhi = -1;
    genWPlusPt = -1;
    genWMinusPt = -1;
    genWPlusEta = -1;
    genWMinusEta = -1;
    genTopPt = -1;
    genTopIsHad = 0;
    genTopEta = -1;
    genAntiTopPt = -1;
    genAntiTopIsHad = 0;
    genAntiTopEta = -1;
    genTTPt = -1;
    genTTEta = -1;
    nJet = 0;
    nIsoJet = 0;
    jet1Flav = 0;
    jet1Phi = -1;
    jet1Pt = -1;
    jet1GenPt = -1;
    jet1Eta = -1;
    jet1CSV = -1;
    jet1IsTight = 0;
    jet2Flav = 0;
    jet2Phi = -1;
    jet2Pt = -1;
    jet2GenPt = -1;
    jet2Eta = -1;
    jet2CSV = -1;
    isojet1Pt = -1;
    isojet1CSV = -1;
    isojet1Flav = 0;
    isojet2Pt = -1;
    isojet2CSV = -1;
    isojet2Flav = 0;
    jet12Mass = -1;
    jet12DEta = -1;
    jetNBtags = 0;
    isojetNBtags = 0;
    nFatjet = 0;
    fj1Tau32 = -1;
    fj1Tau21 = -1;
    fj1Tau32SD = -1;
    fj1Tau21SD = -1;
    fj1MSD = -1;
    fj1MSDScaleUp = -1;
    fj1MSDScaleDown = -1;
    fj1MSDSmeared = -1;
    fj1MSDSmearedUp = -1;
    fj1MSDSmearedDown = -1;
    fj1MSD_corr = -1;
    fj1Pt = -1;
    fj1PtScaleUp = -1;
    fj1PtScaleDown = -1;
    fj1PtSmeared = -1;
    fj1PtSmearedUp = -1;
    fj1PtSmearedDown = -1;
    fj1Phi = -1;
    fj1Eta = -1;
    fj1M = -1;
    fj1MaxCSV = -1;
    fj1MinCSV = -1;
    fj1DoubleCSV = -1;
    fj1GenPt = -1;
    fj1GenSize = -1;
    fj1IsMatched = 0;
    fj1GenWPt = -1;
    fj1GenWSize = -1;
    fj1IsWMatched = 0;
    fj1HighestPtGen = 0;
    fj1HighestPtGenPt = -1;
    fj1IsTight = 0;
    fj1IsLoose = 0;
    fj1RawPt = -1;
    fj1IsHF = 0;
    fj1HTTMass = -1;
    fj1HTTFRec = -1;
    fj1IsClean = 0;
    fj1NConst = 0;
    fj1NSDConst = 0;
    fj1EFrac100 = -1;
    fj1SDEFrac100 = -1;
    isHF = 0;
    nLoosePhoton = 0;
    nTightPhoton = 0;
    loosePho1IsTight = 0;
    loosePho1Pt = -1;
    loosePho1Eta = -1;
    loosePho1Phi = -1;
    nLooseLep = 0;
    nLooseElectron = 0;
    nLooseMuon = 0;
    nTightLep = 0;
    nTightElectron = 0;
    nTightMuon = 0;
    looseLep1PdgId = 0;
    looseLep2PdgId = 0;
    looseLep1IsTight = 0;
    looseLep2IsTight = 0;
    looseLep1Pt = -1;
    looseLep1Eta = -1;
    looseLep1Phi = -1;
    looseLep2Pt = -1;
    looseLep2Eta = -1;
    looseLep2Phi = -1;
    diLepMass = -1;
    nTau = 0;
    mT = -1;
    hbbpt = -1;
    hbbeta = -1;
    hbbphi = -1;
    hbbm = -1;
    scaleUp = -1;
    scaleDown = -1;
    pdfUp = -1;
    pdfDown = -1;
}

void GeneralTree::WriteTree(TTree *t) {
  treePtr = t;
//STARTCUSTOMWRITE
  if (monohiggs) {
    treePtr->Branch("jetPt",jetPt,"jetPt[nJet]/F");
    treePtr->Branch("jetEta",jetEta,"jetEta[nJet]/F");
    treePtr->Branch("jetPhi",jetPhi,"jetPhi[nJet]/F");
    treePtr->Branch("jetE",jetE,"jetE[nJet]/F");
    treePtr->Branch("jetCSV",jetCSV,"jetCSV[nJet]/F");
    treePtr->Branch("jetIso",jetIso,"jetIso[nJet]/F");
    treePtr->Branch("jetQGL",jetQGL,"jetQGL[nJet]/F");
    treePtr->Branch("fj1sjPt",fj1sjPt,"fj1sjPt[2]/F");
    treePtr->Branch("fj1sjPhi",fj1sjPhi,"fj1sjPhi[2]/F");
    treePtr->Branch("fj1sjEta",fj1sjEta,"fj1sjEta[2]/F");
    treePtr->Branch("fj1sjM",fj1sjM,"fj1sjM[2]/F");
    treePtr->Branch("fj1sjCSV",fj1sjCSV,"fj1sjCSV[2]/F");
    treePtr->Branch("fj1sjQGL",fj1sjQGL,"fj1sjQGL[2]/F");
    treePtr->Branch("fj1MSD_corr",&fj1MSD_corr,"fj1MSD_corr/F");
    treePtr->Branch("fj1DoubleCSV",&fj1DoubleCSV,"fj1DoubleCSV/F");
    treePtr->Branch("hbbpt",&hbbpt,"hbbpt/F");
    treePtr->Branch("hbbeta",&hbbeta,"hbbeta/F");
    treePtr->Branch("hbbphi",&hbbphi,"hbbphi/F");
    treePtr->Branch("hbbm",&hbbm,"hbbm/F");
    treePtr->Branch("hbbjtidx",hbbjtidx,"hbbjtidx[2]/I");
  }
  treePtr->Branch("scale",scale,"scale[6]/F");

  for (auto p : ecfParams) { 
    TString ecfn(makeECFString(p));
    treePtr->Branch("fj1"+ecfn,&(fj1ECFNs[p]),"fj1"+ecfn+"/F");
  }

  for (auto p : btagParams) {
    TString btagn(makeBTagSFString(p));
    treePtr->Branch(btagn,&(sf_btags[p]),btagn+"/F");
  }
//ENDCUSTOMWRITE
    treePtr->Branch("fj1SubMaxCSV",&fj1SubMaxCSV,"fj1SubMaxCSV/F");
    treePtr->Branch("looseLep1IsHLTSafe",&looseLep1IsHLTSafe,"looseLep1IsHLTSafe/I");
    treePtr->Branch("looseLep2IsHLTSafe",&looseLep2IsHLTSafe,"looseLep2IsHLTSafe/I");
    treePtr->Branch("runNumber",&runNumber,"runNumber/I");
    treePtr->Branch("lumiNumber",&lumiNumber,"lumiNumber/I");
    treePtr->Branch("eventNumber",&eventNumber,"eventNumber/l");
    treePtr->Branch("npv",&npv,"npv/I");
    treePtr->Branch("pu",&pu,"pu/I");
    treePtr->Branch("mcWeight",&mcWeight,"mcWeight/F");
    treePtr->Branch("trigger",&trigger,"trigger/I");
    treePtr->Branch("metFilter",&metFilter,"metFilter/I");
    treePtr->Branch("egmFilter",&egmFilter,"egmFilter/I");
    treePtr->Branch("filter_maxRecoil",&filter_maxRecoil,"filter_maxRecoil/F");
    treePtr->Branch("filter_whichRecoil",&filter_whichRecoil,"filter_whichRecoil/F");
    treePtr->Branch("sf_ewkV",&sf_ewkV,"sf_ewkV/F");
    treePtr->Branch("sf_qcdV",&sf_qcdV,"sf_qcdV/F");
    treePtr->Branch("sf_ewkV2j",&sf_ewkV2j,"sf_ewkV2j/F");
    treePtr->Branch("sf_qcdV2j",&sf_qcdV2j,"sf_qcdV2j/F");
    treePtr->Branch("sf_qcdTT",&sf_qcdTT,"sf_qcdTT/F");
    treePtr->Branch("sf_lepID",&sf_lepID,"sf_lepID/F");
    treePtr->Branch("sf_lepIso",&sf_lepIso,"sf_lepIso/F");
    treePtr->Branch("sf_lepTrack",&sf_lepTrack,"sf_lepTrack/F");
    treePtr->Branch("sf_pho",&sf_pho,"sf_pho/F");
    treePtr->Branch("sf_eleTrig",&sf_eleTrig,"sf_eleTrig/F");
    treePtr->Branch("sf_phoTrig",&sf_phoTrig,"sf_phoTrig/F");
    treePtr->Branch("sf_metTrig",&sf_metTrig,"sf_metTrig/F");
    treePtr->Branch("sf_pu",&sf_pu,"sf_pu/F");
    treePtr->Branch("sf_npv",&sf_npv,"sf_npv/F");
    treePtr->Branch("sf_tt",&sf_tt,"sf_tt/F");
    treePtr->Branch("sf_tt_ext",&sf_tt_ext,"sf_tt_ext/F");
    treePtr->Branch("sf_tt_bound",&sf_tt_bound,"sf_tt_bound/F");
    treePtr->Branch("sf_tt8TeV",&sf_tt8TeV,"sf_tt8TeV/F");
    treePtr->Branch("sf_tt8TeV_ext",&sf_tt8TeV_ext,"sf_tt8TeV_ext/F");
    treePtr->Branch("sf_tt8TeV_bound",&sf_tt8TeV_bound,"sf_tt8TeV_bound/F");
    treePtr->Branch("sf_phoPurity",&sf_phoPurity,"sf_phoPurity/F");
    treePtr->Branch("pfmet",&pfmet,"pfmet/F");
    treePtr->Branch("pfmetphi",&pfmetphi,"pfmetphi/F");
    treePtr->Branch("pfmetnomu",&pfmetnomu,"pfmetnomu/F");
    treePtr->Branch("puppimet",&puppimet,"puppimet/F");
    treePtr->Branch("puppimetphi",&puppimetphi,"puppimetphi/F");
    treePtr->Branch("calomet",&calomet,"calomet/F");
    treePtr->Branch("calometphi",&calometphi,"calometphi/F");
    treePtr->Branch("pfcalobalance",&pfcalobalance,"pfcalobalance/F");
    treePtr->Branch("sumET",&sumET,"sumET/F");
    treePtr->Branch("trkmet",&trkmet,"trkmet/F");
    treePtr->Branch("puppiUWmag",&puppiUWmag,"puppiUWmag/F");
    treePtr->Branch("puppiUWphi",&puppiUWphi,"puppiUWphi/F");
    treePtr->Branch("puppiUZmag",&puppiUZmag,"puppiUZmag/F");
    treePtr->Branch("puppiUZphi",&puppiUZphi,"puppiUZphi/F");
    treePtr->Branch("puppiUAmag",&puppiUAmag,"puppiUAmag/F");
    treePtr->Branch("puppiUAphi",&puppiUAphi,"puppiUAphi/F");
    treePtr->Branch("puppiUperp",&puppiUperp,"puppiUperp/F");
    treePtr->Branch("puppiUpara",&puppiUpara,"puppiUpara/F");
    treePtr->Branch("puppiUmag",&puppiUmag,"puppiUmag/F");
    treePtr->Branch("puppiUphi",&puppiUphi,"puppiUphi/F");
    treePtr->Branch("pfUWmag",&pfUWmag,"pfUWmag/F");
    treePtr->Branch("pfUWphi",&pfUWphi,"pfUWphi/F");
    treePtr->Branch("pfUZmag",&pfUZmag,"pfUZmag/F");
    treePtr->Branch("pfUZphi",&pfUZphi,"pfUZphi/F");
    treePtr->Branch("pfUAmag",&pfUAmag,"pfUAmag/F");
    treePtr->Branch("pfUAphi",&pfUAphi,"pfUAphi/F");
    treePtr->Branch("pfUperp",&pfUperp,"pfUperp/F");
    treePtr->Branch("pfUpara",&pfUpara,"pfUpara/F");
    treePtr->Branch("pfUmag",&pfUmag,"pfUmag/F");
    treePtr->Branch("pfUphi",&pfUphi,"pfUphi/F");
    treePtr->Branch("dphipfmet",&dphipfmet,"dphipfmet/F");
    treePtr->Branch("dphipuppimet",&dphipuppimet,"dphipuppimet/F");
    treePtr->Branch("dphipuppiUW",&dphipuppiUW,"dphipuppiUW/F");
    treePtr->Branch("dphipuppiUZ",&dphipuppiUZ,"dphipuppiUZ/F");
    treePtr->Branch("dphipuppiUA",&dphipuppiUA,"dphipuppiUA/F");
    treePtr->Branch("dphipfUW",&dphipfUW,"dphipfUW/F");
    treePtr->Branch("dphipfUZ",&dphipfUZ,"dphipfUZ/F");
    treePtr->Branch("dphipfUA",&dphipfUA,"dphipfUA/F");
    treePtr->Branch("dphipuppiU",&dphipuppiU,"dphipuppiU/F");
    treePtr->Branch("dphipfU",&dphipfU,"dphipfU/F");
    treePtr->Branch("trueGenBosonPt",&trueGenBosonPt,"trueGenBosonPt/F");
    treePtr->Branch("genBosonPt",&genBosonPt,"genBosonPt/F");
    treePtr->Branch("genBosonEta",&genBosonEta,"genBosonEta/F");
    treePtr->Branch("genBosonMass",&genBosonMass,"genBosonMass/F");
    treePtr->Branch("genBosonPhi",&genBosonPhi,"genBosonPhi/F");
    treePtr->Branch("genWPlusPt",&genWPlusPt,"genWPlusPt/F");
    treePtr->Branch("genWMinusPt",&genWMinusPt,"genWMinusPt/F");
    treePtr->Branch("genWPlusEta",&genWPlusEta,"genWPlusEta/F");
    treePtr->Branch("genWMinusEta",&genWMinusEta,"genWMinusEta/F");
    treePtr->Branch("genTopPt",&genTopPt,"genTopPt/F");
    treePtr->Branch("genTopIsHad",&genTopIsHad,"genTopIsHad/I");
    treePtr->Branch("genTopEta",&genTopEta,"genTopEta/F");
    treePtr->Branch("genAntiTopPt",&genAntiTopPt,"genAntiTopPt/F");
    treePtr->Branch("genAntiTopIsHad",&genAntiTopIsHad,"genAntiTopIsHad/I");
    treePtr->Branch("genAntiTopEta",&genAntiTopEta,"genAntiTopEta/F");
    treePtr->Branch("genTTPt",&genTTPt,"genTTPt/F");
    treePtr->Branch("genTTEta",&genTTEta,"genTTEta/F");
    treePtr->Branch("nJet",&nJet,"nJet/I");
    treePtr->Branch("nIsoJet",&nIsoJet,"nIsoJet/I");
    treePtr->Branch("jet1Flav",&jet1Flav,"jet1Flav/I");
    treePtr->Branch("jet1Phi",&jet1Phi,"jet1Phi/F");
    treePtr->Branch("jet1Pt",&jet1Pt,"jet1Pt/F");
    treePtr->Branch("jet1GenPt",&jet1GenPt,"jet1GenPt/F");
    treePtr->Branch("jet1Eta",&jet1Eta,"jet1Eta/F");
    treePtr->Branch("jet1CSV",&jet1CSV,"jet1CSV/F");
    treePtr->Branch("jet1IsTight",&jet1IsTight,"jet1IsTight/I");
    treePtr->Branch("jet2Flav",&jet2Flav,"jet2Flav/I");
    treePtr->Branch("jet2Phi",&jet2Phi,"jet2Phi/F");
    treePtr->Branch("jet2Pt",&jet2Pt,"jet2Pt/F");
    treePtr->Branch("jet2GenPt",&jet2GenPt,"jet2GenPt/F");
    treePtr->Branch("jet2Eta",&jet2Eta,"jet2Eta/F");
    treePtr->Branch("jet2CSV",&jet2CSV,"jet2CSV/F");
    treePtr->Branch("isojet1Pt",&isojet1Pt,"isojet1Pt/F");
    treePtr->Branch("isojet1CSV",&isojet1CSV,"isojet1CSV/F");
    treePtr->Branch("isojet1Flav",&isojet1Flav,"isojet1Flav/I");
    treePtr->Branch("isojet2Pt",&isojet2Pt,"isojet2Pt/F");
    treePtr->Branch("isojet2CSV",&isojet2CSV,"isojet2CSV/F");
    treePtr->Branch("isojet2Flav",&isojet2Flav,"isojet2Flav/I");
    treePtr->Branch("jet12Mass",&jet12Mass,"jet12Mass/F");
    treePtr->Branch("jet12DEta",&jet12DEta,"jet12DEta/F");
    treePtr->Branch("jetNBtags",&jetNBtags,"jetNBtags/I");
    treePtr->Branch("isojetNBtags",&isojetNBtags,"isojetNBtags/I");
    treePtr->Branch("nFatjet",&nFatjet,"nFatjet/I");
    treePtr->Branch("fj1Tau32",&fj1Tau32,"fj1Tau32/F");
    treePtr->Branch("fj1Tau21",&fj1Tau21,"fj1Tau21/F");
    treePtr->Branch("fj1Tau32SD",&fj1Tau32SD,"fj1Tau32SD/F");
    treePtr->Branch("fj1Tau21SD",&fj1Tau21SD,"fj1Tau21SD/F");
    treePtr->Branch("fj1MSD",&fj1MSD,"fj1MSD/F");
    treePtr->Branch("fj1MSDScaleUp",&fj1MSDScaleUp,"fj1MSDScaleUp/F");
    treePtr->Branch("fj1MSDScaleDown",&fj1MSDScaleDown,"fj1MSDScaleDown/F");
    treePtr->Branch("fj1MSDSmeared",&fj1MSDSmeared,"fj1MSDSmeared/F");
    treePtr->Branch("fj1MSDSmearedUp",&fj1MSDSmearedUp,"fj1MSDSmearedUp/F");
    treePtr->Branch("fj1MSDSmearedDown",&fj1MSDSmearedDown,"fj1MSDSmearedDown/F");
    treePtr->Branch("fj1MSD_corr",&fj1MSD_corr,"fj1MSD_corr/F");
    treePtr->Branch("fj1Pt",&fj1Pt,"fj1Pt/F");
    treePtr->Branch("fj1PtScaleUp",&fj1PtScaleUp,"fj1PtScaleUp/F");
    treePtr->Branch("fj1PtScaleDown",&fj1PtScaleDown,"fj1PtScaleDown/F");
    treePtr->Branch("fj1PtSmeared",&fj1PtSmeared,"fj1PtSmeared/F");
    treePtr->Branch("fj1PtSmearedUp",&fj1PtSmearedUp,"fj1PtSmearedUp/F");
    treePtr->Branch("fj1PtSmearedDown",&fj1PtSmearedDown,"fj1PtSmearedDown/F");
    treePtr->Branch("fj1Phi",&fj1Phi,"fj1Phi/F");
    treePtr->Branch("fj1Eta",&fj1Eta,"fj1Eta/F");
    treePtr->Branch("fj1M",&fj1M,"fj1M/F");
    treePtr->Branch("fj1MaxCSV",&fj1MaxCSV,"fj1MaxCSV/F");
    treePtr->Branch("fj1MinCSV",&fj1MinCSV,"fj1MinCSV/F");
    treePtr->Branch("fj1DoubleCSV",&fj1DoubleCSV,"fj1DoubleCSV/F");
    treePtr->Branch("fj1GenPt",&fj1GenPt,"fj1GenPt/F");
    treePtr->Branch("fj1GenSize",&fj1GenSize,"fj1GenSize/F");
    treePtr->Branch("fj1IsMatched",&fj1IsMatched,"fj1IsMatched/I");
    treePtr->Branch("fj1GenWPt",&fj1GenWPt,"fj1GenWPt/F");
    treePtr->Branch("fj1GenWSize",&fj1GenWSize,"fj1GenWSize/F");
    treePtr->Branch("fj1IsWMatched",&fj1IsWMatched,"fj1IsWMatched/I");
    treePtr->Branch("fj1HighestPtGen",&fj1HighestPtGen,"fj1HighestPtGen/I");
    treePtr->Branch("fj1HighestPtGenPt",&fj1HighestPtGenPt,"fj1HighestPtGenPt/F");
    treePtr->Branch("fj1IsTight",&fj1IsTight,"fj1IsTight/I");
    treePtr->Branch("fj1IsLoose",&fj1IsLoose,"fj1IsLoose/I");
    treePtr->Branch("fj1RawPt",&fj1RawPt,"fj1RawPt/F");
    treePtr->Branch("fj1IsHF",&fj1IsHF,"fj1IsHF/I");
    treePtr->Branch("fj1HTTMass",&fj1HTTMass,"fj1HTTMass/F");
    treePtr->Branch("fj1HTTFRec",&fj1HTTFRec,"fj1HTTFRec/F");
    treePtr->Branch("fj1IsClean",&fj1IsClean,"fj1IsClean/I");
    treePtr->Branch("fj1NConst",&fj1NConst,"fj1NConst/I");
    treePtr->Branch("fj1NSDConst",&fj1NSDConst,"fj1NSDConst/I");
    treePtr->Branch("fj1EFrac100",&fj1EFrac100,"fj1EFrac100/F");
    treePtr->Branch("fj1SDEFrac100",&fj1SDEFrac100,"fj1SDEFrac100/F");
    treePtr->Branch("isHF",&isHF,"isHF/I");
    treePtr->Branch("nLoosePhoton",&nLoosePhoton,"nLoosePhoton/I");
    treePtr->Branch("nTightPhoton",&nTightPhoton,"nTightPhoton/I");
    treePtr->Branch("loosePho1IsTight",&loosePho1IsTight,"loosePho1IsTight/I");
    treePtr->Branch("loosePho1Pt",&loosePho1Pt,"loosePho1Pt/F");
    treePtr->Branch("loosePho1Eta",&loosePho1Eta,"loosePho1Eta/F");
    treePtr->Branch("loosePho1Phi",&loosePho1Phi,"loosePho1Phi/F");
    treePtr->Branch("nLooseLep",&nLooseLep,"nLooseLep/I");
    treePtr->Branch("nLooseElectron",&nLooseElectron,"nLooseElectron/I");
    treePtr->Branch("nLooseMuon",&nLooseMuon,"nLooseMuon/I");
    treePtr->Branch("nTightLep",&nTightLep,"nTightLep/I");
    treePtr->Branch("nTightElectron",&nTightElectron,"nTightElectron/I");
    treePtr->Branch("nTightMuon",&nTightMuon,"nTightMuon/I");
    treePtr->Branch("looseLep1PdgId",&looseLep1PdgId,"looseLep1PdgId/I");
    treePtr->Branch("looseLep2PdgId",&looseLep2PdgId,"looseLep2PdgId/I");
    treePtr->Branch("looseLep1IsTight",&looseLep1IsTight,"looseLep1IsTight/I");
    treePtr->Branch("looseLep2IsTight",&looseLep2IsTight,"looseLep2IsTight/I");
    treePtr->Branch("looseLep1Pt",&looseLep1Pt,"looseLep1Pt/F");
    treePtr->Branch("looseLep1Eta",&looseLep1Eta,"looseLep1Eta/F");
    treePtr->Branch("looseLep1Phi",&looseLep1Phi,"looseLep1Phi/F");
    treePtr->Branch("looseLep2Pt",&looseLep2Pt,"looseLep2Pt/F");
    treePtr->Branch("looseLep2Eta",&looseLep2Eta,"looseLep2Eta/F");
    treePtr->Branch("looseLep2Phi",&looseLep2Phi,"looseLep2Phi/F");
    treePtr->Branch("diLepMass",&diLepMass,"diLepMass/F");
    treePtr->Branch("nTau",&nTau,"nTau/I");
    treePtr->Branch("mT",&mT,"mT/F");
    treePtr->Branch("scaleUp",&scaleUp,"scaleUp/F");
    treePtr->Branch("scaleDown",&scaleDown,"scaleDown/F");
    treePtr->Branch("pdfUp",&pdfUp,"pdfUp/F");
    treePtr->Branch("pdfDown",&pdfDown,"pdfDown/F");
}

