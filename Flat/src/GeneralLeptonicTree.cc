#include "../interface/GeneralLeptonicTree.h"

GeneralLeptonicTree::GeneralLeptonicTree() {

  for (unsigned iS=0; iS!=6; ++iS) {
    scale[iS] = 1;
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

}

GeneralLeptonicTree::~GeneralLeptonicTree() {

}

void GeneralLeptonicTree::Reset() {

  for (unsigned iS=0; iS!=6; ++iS) {
    scale[iS] = 1;
  }

  for (auto p : btagParams) { 
    sf_btags[p] = 1;
  }

  for (auto iter=signal_weights.begin(); iter!=signal_weights.end(); ++iter) {
    signal_weights[iter->first] = 1; // does pair::second return a reference?
  }

    runNumber = -1;
    lumiNumber = -1;
    eventNumber = -1;
    npv = -1;
    pu = -1;
    mcWeight = 1;
    trigger = 0;
    metFilter = 0;
    egmFilter = 0;

    nLooseLep = 0;
    looseGenLep1PdgId = 0;
    looseGenLep2PdgId = 0;
    looseGenLep3PdgId = 0;
    looseGenLep4PdgId = 0;
    looseLep1PdgId = -1;
    looseLep2PdgId = -1;
    looseLep3PdgId = -1;
    looseLep4PdgId = -1;
    looseLep1SelBit = 0;
    looseLep2SelBit = 0;
    looseLep3SelBit = 0;
    looseLep4SelBit = 0;
    looseLep1Pt = -1;
    looseLep2Pt = -1;
    looseLep3Pt = -1;
    looseLep4Pt = -1;
    looseLep1Eta = -1;
    looseLep2Eta = -1;
    looseLep3Eta = -1;
    looseLep4Eta = -1;
    looseLep1Phi = -1;
    looseLep2Phi = -1;
    looseLep3Phi = -1;
    looseLep4Phi = -1;

    nJet = 0;
    jetNLBtags = 0;
    jetNMBtags = 0;
    jetNTBtags = 0;

    jet1Pt    = -1;
    jet2Pt    = -1;
    jet3Pt    = -1;
    jet4Pt    = -1;
    jet1Eta   = -1;
    jet2Eta   = -1;
    jet3Eta   = -1;
    jet4Eta   = -1;
    jet1Phi   = -1;
    jet2Phi   = -1;
    jet3Phi   = -1;
    jet4Phi   = -1;
    jet1BTag  = -1;
    jet2BTag  = -1;
    jet3BTag  = -1;
    jet4BTag  = -1;
    jet1GenPt = -1;
    jet2GenPt = -1;
    jet3GenPt = -1;
    jet4GenPt = -1;
    jet1Flav    = -1;
    jet2Flav    = -1;
    jet3Flav    = -1;
    jet4Flav    = -1;
    jet1SelBit  = 0;
    jet2SelBit  = 0;
    jet3SelBit  = 0;
    jet4SelBit  = 0;

    jet1PtUp    = -1;
    jet2PtUp    = -1;
    jet3PtUp    = -1;
    jet4PtUp    = -1;
    jet1PtDown  = -1;
    jet2PtDown  = -1;
    jet3PtDown  = -1;
    jet4PtDown  = -1;
    jet1EtaUp   = -1;
    jet2EtaUp   = -1;
    jet3EtaUp   = -1;
    jet4EtaUp   = -1;
    jet1EtaDown = -1;
    jet2EtaDown = -1;
    jet3EtaDown = -1;
    jet4EtaDown = -1;

    pfmet = -1;
    pfmetphi = -1;
    pfmetRaw = -1;
    pfmetUp = -1;
    pfmetDown = -1;
    pfmetnomu = -1;
    puppimet = -1;
    puppimetphi = -1;
    calomet = -1;
    calometphi = -1;
    trkmet = -1;
    trkmetphi = -1;
    dphipfmet = -1;
    dphipuppimet = -1;

    genLep1Pt = -1;
    genLep2Pt = -1;
    genLep1Eta = -1;
    genLep2Eta = -1;
    genLep1Phi = -1;
    genLep2Phi = -1;
    genLep1PdgId = 0;
    genLep2PdgId = 0;

    nTau = 0;
    pdfUp = 1;
    pdfDown = 1;

    nLoosePhoton = 0;
    loosePho1Pt  = -1;
    loosePho1Eta = -1;
    loosePho1Phi = -1;
   
    sf_pu = 1.0;
    sf_puUp = 1.0;
    sf_puDown = 1.0;
    sf_zz = 1.0;
    sf_zzUnc = 1.0;
    sf_wz = 1.0;
    sf_zh = 1.0;
    sf_zhUp = 1.0;
    sf_zhDown = 1.0;
    sf_tt = 1.0;

    sf_trk1 = 1.0;
    sf_trk2 = 1.0;
    sf_trk3 = 1.0;
    sf_trk4 = 1.0;
    sf_loose1 = 1.0;
    sf_loose2 = 1.0;
    sf_loose3 = 1.0;
    sf_loose4 = 1.0;
    sf_medium1 = 1.0;
    sf_medium2 = 1.0;
    sf_medium3 = 1.0;
    sf_medium4 = 1.0;
    sf_tight1 = 1.0;
    sf_tight2 = 1.0;
    sf_tight3 = 1.0;
    sf_tight4 = 1.0;

    sf_unc1 = 1.0;
    sf_unc2 = 1.0;
    sf_unc3 = 1.0;
    sf_unc4 = 1.0;

}

void GeneralLeptonicTree::WriteTree(TTree *t) {
  treePtr = t;

  for (auto iter=signal_weights.begin(); iter!=signal_weights.end(); ++iter) {
    Book("rw_"+iter->first,&(signal_weights[iter->first]),"rw_"+iter->first+"/F");
  }

  Book("scale",scale,"scale[6]/F");

  for (auto p : btagParams) {
    TString btagn(makeBTagSFString(p));
    Book(btagn,&(sf_btags[p]),btagn+"/F");
  }

  Book("runNumber",&runNumber,"runNumber/I");
  Book("lumiNumber",&lumiNumber,"lumiNumber/I");
  Book("eventNumber",&eventNumber,"eventNumber/l");
  Book("npv",&npv,"npv/I");
  Book("pu",&pu,"pu/I");
  Book("mcWeight",&mcWeight,"mcWeight/F");
  Book("trigger",&trigger,"trigger/I");
  Book("metFilter",&metFilter,"metFilter/I");
  Book("egmFilter",&egmFilter,"egmFilter/I");

  Book("nLooseLep",&nLooseLep,"nLooseLep/I");
  Book("looseGenLep1PdgId",&looseGenLep1PdgId,"looseGenLep1PdgId/I");
  Book("looseGenLep2PdgId",&looseGenLep2PdgId,"looseGenLep2PdgId/I");
  Book("looseGenLep3PdgId",&looseGenLep3PdgId,"looseGenLep3PdgId/I");
  Book("looseGenLep4PdgId",&looseGenLep4PdgId,"looseGenLep4PdgId/I");
  Book("looseLep1PdgId",&looseLep1PdgId,"looseLep1PdgId/I");
  Book("looseLep2PdgId",&looseLep2PdgId,"looseLep2PdgId/I");
  Book("looseLep3PdgId",&looseLep3PdgId,"looseLep3PdgId/I");
  Book("looseLep4PdgId",&looseLep4PdgId,"looseLep4PdgId/I");
  Book("looseLep1SelBit",&looseLep1SelBit,"looseLep1SelBit/I");
  Book("looseLep2SelBit",&looseLep2SelBit,"looseLep2SelBit/I");
  Book("looseLep3SelBit",&looseLep3SelBit,"looseLep3SelBit/I");
  Book("looseLep4SelBit",&looseLep4SelBit,"looseLep4SelBit/I");
  Book("looseLep1Pt",&looseLep1Pt,"looseLep1Pt/F");
  Book("looseLep2Pt",&looseLep2Pt,"looseLep2Pt/F");
  Book("looseLep3Pt",&looseLep3Pt,"looseLep3Pt/F");
  Book("looseLep4Pt",&looseLep4Pt,"looseLep4Pt/F");
  Book("looseLep1Eta",&looseLep1Eta,"looseLep1Eta/F");
  Book("looseLep2Eta",&looseLep2Eta,"looseLep2Eta/F");
  Book("looseLep3Eta",&looseLep3Eta,"looseLep3Eta/F");
  Book("looseLep4Eta",&looseLep4Eta,"looseLep4Eta/F");
  Book("looseLep1Phi",&looseLep1Phi,"looseLep1Phi/F");
  Book("looseLep2Phi",&looseLep2Phi,"looseLep2Phi/F");
  Book("looseLep3Phi",&looseLep3Phi,"looseLep3Phi/F");
  Book("looseLep4Phi",&looseLep4Phi,"looseLep4Phi/F");

  Book("nJet",&nJet,"nJet/I");
  Book("jetNLBtags",&jetNLBtags,"jetNLBtags/I");
  Book("jetNMBtags",&jetNMBtags,"jetNMBtags/I");
  Book("jetNTBtags",&jetNTBtags,"jetNTBtags/I");

  Book("jet1Pt",&jet1Pt,"jet1Pt/F");
  Book("jet2Pt",&jet2Pt,"jet2Pt/F");
  Book("jet3Pt",&jet3Pt,"jet3Pt/F");
  Book("jet4Pt",&jet4Pt,"jet4Pt/F");
  Book("jet1Eta",&jet1Eta,"jet1Eta/F");
  Book("jet2Eta",&jet2Eta,"jet2Eta/F");
  Book("jet3Eta",&jet3Eta,"jet3Eta/F");
  Book("jet4Eta",&jet4Eta,"jet4Eta/F");
  Book("jet1Phi",&jet1Phi,"jet1Phi/F");
  Book("jet2Phi",&jet2Phi,"jet2Phi/F");
  Book("jet3Phi",&jet3Phi,"jet3Phi/F");
  Book("jet4Phi",&jet4Phi,"jet4Phi/F");
  Book("jet1BTag",&jet1BTag,"jet1BTag/F");
  Book("jet2BTag",&jet2BTag,"jet2BTag/F");
  Book("jet3BTag",&jet3BTag,"jet3BTag/F");
  Book("jet4BTag",&jet4BTag,"jet4BTag/F");
  Book("jet1GenPt",&jet1GenPt,"jet1GenPt/F");
  Book("jet2GenPt",&jet2GenPt,"jet2GenPt/F");
  Book("jet3GenPt",&jet3GenPt,"jet3GenPt/F");
  Book("jet4GenPt",&jet4GenPt,"jet4GenPt/F");
  Book("jet1Flav",&jet1Flav,"jet1Flav/I");
  Book("jet2Flav",&jet2Flav,"jet2Flav/I");
  Book("jet3Flav",&jet3Flav,"jet3Flav/I");
  Book("jet4Flav",&jet4Flav,"jet4Flav/I");
  Book("jet1SelBit",&jet1SelBit,"jet1SelBit/I");
  Book("jet2SelBit",&jet2SelBit,"jet2SelBit/I");
  Book("jet3SelBit",&jet3SelBit,"jet3SelBit/I");
  Book("jet4SelBit",&jet4SelBit,"jet4SelBit/I");

  Book("jet1PtUp",&jet1PtUp,"jet1PtUp/F");
  Book("jet2PtUp",&jet2PtUp,"jet2PtUp/F");
  Book("jet3PtUp",&jet3PtUp,"jet3PtUp/F");
  Book("jet4PtUp",&jet4PtUp,"jet4PtUp/F");
  Book("jet1PtDown",&jet1PtDown,"jet1PtDown/F");
  Book("jet2PtDown",&jet2PtDown,"jet2PtDown/F");
  Book("jet3PtDown",&jet3PtDown,"jet3PtDown/F");
  Book("jet4PtDown",&jet4PtDown,"jet4PtDown/F");
  Book("jet1EtaUp",&jet1EtaUp,"jet1EtaUp/F");
  Book("jet2EtaUp",&jet2EtaUp,"jet2EtaUp/F");
  Book("jet3EtaUp",&jet3EtaUp,"jet3EtaUp/F");
  Book("jet4EtaUp",&jet4EtaUp,"jet4EtaUp/F");
  Book("jet1EtaDown",&jet1EtaDown,"jet1EtaDown/F");
  Book("jet2EtaDown",&jet2EtaDown,"jet2EtaDown/F");
  Book("jet3EtaDown",&jet3EtaDown,"jet3EtaDown/F");
  Book("jet4EtaDown",&jet4EtaDown,"jet4EtaDown/F");

  Book("pfmet",&pfmet,"pfmet/F");
  Book("pfmetphi",&pfmetphi,"pfmetphi/F");
  Book("pfmetRaw",&pfmetRaw,"pfmetRaw/F");
  Book("pfmetUp",&pfmetUp,"pfmetUp/F");
  Book("pfmetDown",&pfmetDown,"pfmetDown/F");
  Book("pfmetnomu",&pfmetnomu,"pfmetnomu/F");
  Book("puppimet",&puppimet,"puppimet/F");
  Book("puppimetphi",&puppimetphi,"puppimetphi/F");
  Book("calomet",&calomet,"calomet/F");
  Book("calometphi",&calometphi,"calometphi/F");
  Book("trkmet",&trkmet,"trkmet/F");
  Book("trkmetphi",&trkmetphi,"trkmetphi/F");
  Book("dphipfmet",&dphipfmet,"dphipfmet/F");
  Book("dphipuppimet",&dphipuppimet,"dphipuppimet/F");

  Book("genLep1Pt",&genLep1Pt,"genLep1Pt/F");
  Book("genLep2Pt",&genLep2Pt,"genLep2Pt/F");
  Book("genLep1Eta",&genLep1Eta,"genLep1Eta/F");
  Book("genLep2Eta",&genLep2Eta,"genLep2Eta/F");
  Book("genLep1Phi",&genLep1Phi,"genLep1Phi/F");
  Book("genLep2Phi",&genLep2Phi,"genLep2Phi/F");
  Book("genLep1PdgId",&genLep1PdgId,"genLep1PdgId/I");
  Book("genLep2PdgId",&genLep2PdgId,"genLep2PdgId/I");

  Book("nTau",&nTau,"nTau/I");
  Book("pdfUp",&pdfUp,"pdfUp/F");
  Book("pdfDown",&pdfDown,"pdfDown/F");

  Book("nLoosePhoton",&nLoosePhoton,"nLoosePhoton/I");
  Book("loosePho1Pt",&loosePho1Pt,"loosePho1Pt/F");
  Book("loosePho1Eta",&loosePho1Eta,"loosePho1Eta/F");
  Book("loosePho1Phi",&loosePho1Phi,"loosePho1Phi/F");

  Book("sf_pu",&sf_pu,"sf_pu/F");
  Book("sf_puUp",&sf_puUp,"sf_puUp/F");
  Book("sf_puDown",&sf_puDown,"sf_puDown/F");
  Book("sf_zz",&sf_zz,"sf_zz/F");
  Book("sf_zzUnc",&sf_zzUnc,"sf_zzUnc/F");
  Book("sf_wz",&sf_wz,"sf_wz/F");
  Book("sf_zh",&sf_zh,"sf_zh/F");
  Book("sf_zhUp",&sf_zhUp,"sf_zhUp/F");
  Book("sf_zhDown",&sf_zhDown,"sf_zhDown/F");
  Book("sf_tt",&sf_tt,"sf_tt/F");

  Book("sf_trk1",&sf_trk1,"sf_trk1/F");
  Book("sf_trk2",&sf_trk2,"sf_trk2/F");
  Book("sf_trk3",&sf_trk3,"sf_trk3/F");
  Book("sf_trk4",&sf_trk4,"sf_trk4/F");
  Book("sf_loose1",&sf_loose1,"sf_loose1/F");
  Book("sf_loose2",&sf_loose2,"sf_loose2/F");
  Book("sf_loose3",&sf_loose3,"sf_loose3/F");
  Book("sf_loose4",&sf_loose4,"sf_loose4/F");
  Book("sf_medium1",&sf_medium1,"sf_medium1/F");
  Book("sf_medium2",&sf_medium2,"sf_medium2/F");
  Book("sf_medium3",&sf_medium3,"sf_medium3/F");
  Book("sf_medium4",&sf_medium4,"sf_medium4/F");
  Book("sf_tight1",&sf_tight1,"sf_tight1/F");
  Book("sf_tight2",&sf_tight2,"sf_tight2/F");
  Book("sf_tight3",&sf_tight3,"sf_tight3/F");
  Book("sf_tight4",&sf_tight4,"sf_tight4/F");
  Book("sf_unc1",&sf_unc1,"sf_unc1/F");
  Book("sf_unc2",&sf_unc2,"sf_unc2/F");
  Book("sf_unc3",&sf_unc3,"sf_unc3/F");
  Book("sf_unc4",&sf_unc4,"sf_unc4/F");

}
