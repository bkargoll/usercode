import ROOT
from DataFormats.FWLite import Events, Handle
import math
import array
from tdrstyle import setTDRStyle


def deltaPhi(phi1, phi2):
     PHI = abs(phi1-phi2)
     if (PHI<=3.14159265):
         return PHI
     else:
         return 2*3.14159265-PHI

def deltaR(eta1, phi1, eta2, phi2) :
    deta = eta1-eta2
    dphi = deltaPhi(phi1,phi2)
    return math.sqrt(deta*deta + dphi*dphi)




#xbin = [10,15,18,20,22,25,30,35,40,45,50,60,75,100,160]
xbin = [5,10,15,17,18,19,20,25,30,35,40,45,50,55,60,70,80,90,100,120,140,160,200]
xbins = array.array('d',xbin)


def histograms(fileName, labelOffline, L2PtCut = 25., MTCut = 40.):

    events = Events(fileName)
    # vertex
    VertexHandle = Handle('std::vector<reco::Vertex>')
    vertexLabel = "offlinePrimaryVertices"
    # vertexLabel = "hltPixelVertices"

    # taus
    TauHandle  = Handle('std::vector<reco::PFTau>')
    JetHandle = Handle('std::vector<reco::CaloJet>')
    METHandle = Handle('std::vector<reco::PFMET>')
    PFJetHandle = Handle('std::vector<reco::PFJet>')
    MuonHandle = Handle('std::vector<reco::Muon>')
    
    # Histograms
    mTH = ROOT.TH1F("mTH","",100,0,100)
    numeratorJet = ROOT.TH1F("numeratorJet","",len(xbin)-1,xbins)
    denominatorJet = ROOT.TH1F("denominatorJet","",len(xbin)-1,xbins)
    numeratorJetEta = ROOT.TH1F("numeratorJetEta","",60,-3.0,3.0)
    denominatorJetEta = ROOT.TH1F("denominatorJetEta","",60,-3.0,3.0)
    

    # xbinEta = [0.,0.5,1.,1.5,2.5]
    # xbinsEta = array.array('d',xbinEta)
    # numeratorEta = ROOT.TH1F("numeratorEta","",len(xbinEta)-1,xbinsEta)
    # denominatorEta = ROOT.TH1F("denominatorEta","",len(xbinEta)-1,xbinsEta)
    # numeratorL1Vertex = ROOT.TH1F("numeratorL1Vertex","",10,0,40)
    # numeratorL2Vertex = ROOT.TH1F("numeratorL2Vertex","",10,0,40)
    # numeratorL2p5Vertex = ROOT.TH1F("numeratorL2p5Vertex","",10,0,40)
    # numeratorL3Vertex = ROOT.TH1F("numeratorL3Vertex","",10,0,40)
    # denominatorVertex = ROOT.TH1F("denominatorVertex","",10,0,40)
    # ROOT.SetOwnership(resolutionPt,True)
    
    NEvt = 0

    for event in events:
        NEvt += 1
        if NEvt % 1000 == 0: print "Process Event ", NEvt
        #if NEvt > 100: break
        # getting the handle from the event
        event.getByLabel(labelOffline, TauHandle)
        offlineTaus = TauHandle.product()

        event.getByLabel("ak5PFJetsL1L2L3Residual", PFJetHandle)
        offlineJets = PFJetHandle.product()
        
        event.getByLabel("hltCaloJetL1FastJetCorrected",JetHandle)
        l2Jets = JetHandle.product()

        event.getByLabel("pfMet",METHandle)
        met = METHandle.product().at(0)
        
        event.getByLabel("offlineSelectedMuons",MuonHandle)
        muons= MuonHandle.product()
        
        event.getByLabel(vertexLabel,VertexHandle)
        offlinevertex = VertexHandle.product().at(0).position()
        vertexSize = VertexHandle.product().size()

        # matching variables
        matchingConeL1 = 0.5
        matchingConeHLT = 0.3
        foundL1 = False
        foundL2 = False
        foundL2p5 = False
        foundL3 = False
        # making MT variable
        MT = 0.
        
        myMuons = []
        for muon in muons:
             if muon.pt()<10: continue
             if abs(muon.eta()) > 2.1: continue
             #apply vertex cuts which have been missing before
             if abs(muon.innerTrack().dxy(offlinevertex)) > 0.045: continue
             if abs(muon.innerTrack().dz(offlinevertex)) > 0.1: continue
             myMuons.append((muon.pt(), muon))
        myMuons.sort()
        # get at least 1 muon
        
        if len(myMuons) < 1: continue 
        muon = myMuons.pop()[1]
        MT = (muon.pt()+met.pt())*(muon.pt()+met.pt())-(muon.px()+met.px())*(muon.px()+met.px()) - (muon.py()+met.py())*(muon.py()+met.py())
        MT = math.sqrt(MT)
        # Select W+Jets
        mTH.Fill(MT)        
        if MT < MTCut: continue
        
        # taking the offline PFJet
        myJets = []
        for jet in offlineJets:
             if jet.pt() < 10: continue
             if abs( jet.eta() ) > 3.0: continue
             matched = False
             for tau in offlineTaus:
                  dr = deltaR(tau.eta(),tau.phi(),jet.eta(),jet.phi())
                  if dr <0.5:
                       matched = True
                       break
             if not(matched):
                  myJets.append((jet.pt(),jet))

        myJets.sort()
        if len(myJets) < 1: continue    
        myJet = myJets.pop()[1]
        if deltaR(myJet.eta(),myJet.phi(), muon.eta(),muon.phi()) < 0.3:
             if len(myJets) > 0:
                  myJet = myJets.pop()[1]
                  denominatorJet.Fill(myJet.pt())
                  denominatorJetEta.Fill(myJet.eta())
             else:
                  continue
        else:
             denominatorJet.Fill(myJet.pt())
        for jet in l2Jets:
             if jet.pt() < L2PtCut: continue
             if abs( jet.eta() )> 3.0: continue
             if deltaR(jet.eta(),jet.phi(), myJet.eta(), myJet.phi()) < 0.5:
                  numeratorJet.Fill(myJet.pt())
                  numeratorJetEta.Fill(myJet.eta())
                  break
    

       
    numeratorJet.Sumw2()
    denominatorJet.Sumw2()
    numeratorJetEta.Sumw2()
    denominatorJetEta.Sumw2()

    return mTH, numeratorJet, denominatorJet, numeratorJetEta, denominatorJetEta


# running the macro

# just to avoid opening windows
ROOT.gROOT.SetBatch()

setTDRStyle()

'''
## Run2012A
dataset = "Run2012A"
folder = "dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/user/bkargoll/TriggerEfficiencies/SingleMuDataset/MuTauSkim/HLTrerunDoubleTauJet_v3/"
files = ["hltDoubleTauJet_offlineSelectedMuTau_3_1_n7X.root",
         "hltDoubleTauJet_offlineSelectedMuTau_5_1_Dzk.root",
         "hltDoubleTauJet_offlineSelectedMuTau_4_1_rwb.root",
         "hltDoubleTauJet_offlineSelectedMuTau_7_1_xgB.root",
         "hltDoubleTauJet_offlineSelectedMuTau_2_1_Aey.root",
         "hltDoubleTauJet_offlineSelectedMuTau_1_1_ETO.root",
         "hltDoubleTauJet_offlineSelectedMuTau_6_1_I6O.root",
         "hltDoubleTauJet_offlineSelectedMuTau_8_1_p0U.root"]

## Run2012B FirstPart
dataset = "FirstPartRun2012B"
folder = "dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/user/bkargoll/TriggerEfficiencies/SingleMuDataset/MuTauSkim/Run2012B/FirstPart/HLTrerunDoubleTauJet_v3/"
files = ["hltDoubleTauJet_offlineSelectedMuTau_3_1_zKI.root",
         "hltDoubleTauJet_offlineSelectedMuTau_6_1_xJO.root", 
         "hltDoubleTauJet_offlineSelectedMuTau_5_1_HEZ.root", 
         "hltDoubleTauJet_offlineSelectedMuTau_7_1_sAD.root", 
         "hltDoubleTauJet_offlineSelectedMuTau_2_1_LgX.root", 
         "hltDoubleTauJet_offlineSelectedMuTau_4_1_rur.root", 
         "hltDoubleTauJet_offlineSelectedMuTau_1_1_d6N.root", 
         "hltDoubleTauJet_offlineSelectedMuTau_8_1_UDb.root"] 
'''
## Run2012B SecondPart
dataset = "SecondPartRun2012B"
folder = "dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/user/bkargoll/TriggerEfficiencies/SingleMuDataset/MuTauSkim/Run2012B/SecondPart/HLTrerunDoubleTauJet_v3/"
files = ["hltDoubleTauJet_offlineSelectedMuTau_3_1_bxj.root",
         "hltDoubleTauJet_offlineSelectedMuTau_1_1_wS5.root", 
         "hltDoubleTauJet_offlineSelectedMuTau_5_1_PN8.root", 
         "hltDoubleTauJet_offlineSelectedMuTau_4_1_SJ5.root",
         "hltDoubleTauJet_offlineSelectedMuTau_2_1_f1K.root"] 


input = []
for file in files:
    input.append(folder+file)

offlineTauLabel = "offlineSelectedTausXXMisodbMelecrej"
hltPt = 30
MT = 60

myhistos = histograms(input,offlineTauLabel,hltPt, MT)
filename = files[0].split("/")[-1]
filename = filename.split(".")[0]
filename = filename.split("_")[0]+filename.split("_")[1]
newfile = ROOT.TFile("/user/kargoll/results/turnOnJetQCDCorrected_"+dataset+"_"+filename+"_Pt"+str(hltPt)+"_MT"+str(MT)+".root","RECREATE")
for i in range(5):
     myhistos[i].Write()
newfile.Close()     


