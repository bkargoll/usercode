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
xbin = [5,10,15,17,18,19,20,25,30,35,40,45,50,55,60,120,200]
xbins = array.array('d',xbin)


def histograms(fileName, labelOffline, labelL1, labelL2,labelL2p5, labelL3, L2PtCut = 25., L3PtCut = 25., MTCut = 40., nProngsCut = 5, dZCut = 0.2, offlineProngs = 0, maxInvM = 10000.):

    events = Events(fileName)
    # vertex
    VertexHandle = Handle('std::vector<reco::Vertex>')
    vertexLabel = "offlinePrimaryVertices"
    # vertexLabel = "hltPixelVertices"

    # taus
    TauHandle  = Handle('std::vector<reco::PFTau>')
    JetHandle = Handle('std::vector<reco::CaloJet>')
    L1Handle = Handle('std::vector<l1extra::L1JetParticle>')
    METHandle = Handle('std::vector<reco::PFMET>')
    PFJetHandle = Handle('std::vector<reco::PFJet>')
    MuonHandle = Handle('std::vector<reco::Muon>')
    
    # Histograms
    deltaZTauTau = ROOT.TH1F("deltaZTauTau","",200,0,10.)
    mTH = ROOT.TH1F("mTH","",100,0,100)
    invM = ROOT.TH1F("invM","",200,0,200)
    invMSel = ROOT.TH1F("invMSel","",200,0,200)
    numeratorL1 = ROOT.TH1F("numeratorL1","",len(xbin)-1,xbins)
    numeratorL2 = ROOT.TH1F("numeratorL2","",len(xbin)-1,xbins)
    numeratorL1L2 = ROOT.TH1F("numeratorL1L2","",len(xbin)-1,xbins)
    numeratorL2p5 = ROOT.TH1F("numeratorL2p5","",len(xbin)-1,xbins)
    numeratorL1L2L2p5 = ROOT.TH1F("numeratorL1L2L2p5","",len(xbin)-1,xbins)
    numeratorL3 = ROOT.TH1F("numeratorL3","",len(xbin)-1,xbins)
    numeratorL1_HLT = ROOT.TH1F("numeratorL1_HLT","",len(xbin)-1,xbins)
    numeratorJet = ROOT.TH1F("numeratorJet","",len(xbin)-1,xbins)
    
    denominator = ROOT.TH1F("denominator","",len(xbin)-1,xbins)
    denominatorJet = ROOT.TH1F("denominatorJet","",len(xbin)-1,xbins)
    
    nTaus = ROOT.TH1I("NTaus","",20,-0.5,19.5)
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
    NEvtMuon = 0
    NEvtMuonMT = 0
    NEvtMuonMTTau = 0
    NMatchedL1 = 0 
    NMatchedL2 = 0 
    NMatchedL2p5 = 0 
    NMatchedL3 = 0      

    for event in events:
        NEvt += 1
        if NEvt % 1000 == 0: print "Process Event ", NEvt
        if NEvt > 2000: break
         
        # getting the handle from the event
        event.getByLabel(labelOffline, TauHandle)
        offlineTaus = TauHandle.product()
        
        event.getByLabel("ak5PFJets", PFJetHandle)
        offlineJets = PFJetHandle.product()
        
        event.getByLabel(labelL1,"Tau",L1Handle)
        l1Taus = L1Handle.product()             

        event.getByLabel(labelL1,"Central",L1Handle)
        l1Jets = L1Handle.product()             

        event.getByLabel(labelL2,JetHandle)
        l2Taus = JetHandle.product()

        #event.getByLabel("hltCaloJetL1FastJetCorrected",JetHandle)
        #l2Jets = JetHandle.product()

        event.getByLabel(labelL2p5,JetHandle)
        l2p5Taus = JetHandle.product()             

        event.getByLabel(labelL3,TauHandle)
        l3Taus = TauHandle.product()

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
             #apply vertex cuts which have been missing in muon selection before
             if abs(muon.innerTrack().dxy(offlinevertex)) > 0.045: continue
             if abs(muon.innerTrack().dz(offlinevertex)) > 0.1: continue
             myMuons.append((muon.pt(), muon))
        myMuons.sort()
        # get at least 1 muon
        if len(myMuons) < 1: continue 
        # additional muon veto
        if len(myMuons) > 1: continue
        NEvtMuon += 1
        muon = myMuons.pop()[1]
        MT = (muon.pt()+met.pt())*(muon.pt()+met.pt())-(muon.px()+met.px())*(muon.px()+met.px()) - (muon.py()+met.py())*(muon.py()+met.py())
        MT = math.sqrt(MT)
        # Remove W+Jets
        mTH.Fill(MT)        
        if MT > MTCut: continue
        NEvtMuonMT += 1

   
        #select best tau
        cntTaus = 0
        tau = None
        bestZmass = 999999.0
        for iTau in offlineTaus:
            cntTaus += 1
            if iTau.pt() < 5: continue
            if abs(iTau.eta()) > 2.1: continue
            # muon and tau with opposite charge
            if iTau.charge() != -1 * muon.charge(): continue
            # require seperation between muon and tau
            dR_mutau = deltaR(iTau.eta(),iTau.phi(),muon.eta(),muon.phi())
            if dR_mutau < 0.5 : continue
            # make sure tau has a leading track, which is needed later for the dZ cut
            if not iTau.leadPFChargedHadrCand().trackRef().isNonnull(): continue
            # offline taus: only 1prongs, only 3prongs or all?
            if not offlineProngs == 0:
                if not iTau.signalPFChargedHadrCands().size() == offlineProngs: continue
            NEvtMuonMTTau += 1
            Zmass = 2*(iTau.mass()*muon.mass() + muon.energy()*iTau.energy() - muon.px()*iTau.px() - muon.py()*iTau.py() - muon.pz()*iTau.pz())
            Zmass = math.sqrt(Zmass)
            invM.Fill(Zmass)
            if tau == None:
                bestZmass = Zmass
                tau = iTau
                continue
            if iTau.pt() > tau.pt():
                bestZmass = Zmass
                tau = iTau
            
        #skip Z->mumu events, where 1 mu fakes a tau
        if bestZmass > maxInvM: continue
        
        #fill denominator
        nTaus.Fill(cntTaus)
        if tau == None: continue
        invMSel.Fill(bestZmass)
        denominator.Fill(tau.pt())

        # loop on numerators
        myL1Taus = []
        myL1Jets = []
        myL2Taus = []
        myL2p5Taus = []
        myL3Taus = []
        for l1Tau in l1Taus:
             if l1Tau.pt() < 44: continue
             if abs(l1Tau.eta()) > 2.1: continue
             dr = deltaR(tau.eta(),tau.phi(),l1Tau.eta(),l1Tau.phi())
             myL1Taus.append((dr, l1Tau))            
        for l1Jet in l1Jets:
             if l1Jet.pt() < 64: continue
             if abs(l1Jet.eta()) > 2.1: continue
             dr = deltaR(tau.eta(),tau.phi(),l1Jet.eta(),l1Jet.phi())
             myL1Jets.append((dr,l1Jet))
        for l2Tau in l2Taus:
             if l2Tau.pt() < L2PtCut: continue
             if abs(l2Tau.eta()) > 2.1: continue
             dr = deltaR(tau.eta(),tau.phi(),l2Tau.eta(),l2Tau.phi())
             myL2Taus.append((dr, l2Tau))            
        for l2p5Tau in l2p5Taus:
             if l2p5Tau.pt() < L2PtCut: continue
             if abs(l2p5Tau.eta()) > 2.1: continue
             dr = deltaR(tau.eta(),tau.phi(),l2p5Tau.eta(),l2p5Tau.phi())
             myL2p5Taus.append((dr, l2p5Tau))            
        for l3Tau in l3Taus:
             if l3Tau.pt() < L3PtCut: continue
             if abs(l3Tau.eta()) > 2.1: continue
             dr = deltaR(tau.eta(),tau.phi(),l3Tau.eta(),l3Tau.phi())
             myL3Taus.append((dr, l3Tau))            

        # Filling L1 condition
        if len(myL1Taus) > 0:
             myL1Taus.sort()
             if myL1Taus[0][0] < matchingConeL1:
                  numeratorL1.Fill(tau.pt())
                  foundL1 = True
                  NMatchedL1 += 1
        if len(myL1Jets) > 0:
             myL1Jets.sort()
             if myL1Jets[0][0] < matchingConeL1:
                  numeratorL1.Fill(tau.pt())
                  foundL1 = True
                  NMatchedL1 += 1
        # Filling L2 condition
        if len(myL2Taus) >0:
             myL2Taus.sort()
             if myL2Taus[0][0] < matchingConeHLT:
                  numeratorL2.Fill(tau.pt())
                  foundL2 = True
                  NMatchedL2 += 1
        if len(myL2p5Taus) >0:
             myL2p5Taus.sort()
             if myL2p5Taus[0][0] < matchingConeHLT:
                  numeratorL2p5.Fill(tau.pt())
                  foundL2p5 = True
                  NMatchedL2p5 += 1
        if len(myL3Taus) >0:
             myL3Taus.sort()
             if myL3Taus[0][0] < matchingConeHLT:
                  if myL3Taus[0][1].leadPFChargedHadrCand().isNonnull():
                       if myL3Taus[0][1].signalPFChargedHadrCands().size() < nProngsCut:
                           if myL3Taus[0][1].leadPFChargedHadrCand().trackRef().isNonnull():
                               dZ = abs(tau.leadPFChargedHadrCand().trackRef().dz(offlinevertex) - myL3Taus[0][1].leadPFChargedHadrCand().trackRef().dz(offlinevertex))
                               deltaZTauTau.Fill(dZ)
                               if dZ < dZCut:
                                   foundL3 = True
                                   numeratorL3.Fill(tau.pt())
                                   NMatchedL3 += 1

        
        if foundL1 and foundL2 :
            numeratorL1L2.Fill(tau.pt())
            if foundL2p5:
                numeratorL1L2L2p5.Fill(tau.pt())
                if foundL3:
                    numeratorL1_HLT.Fill(tau.pt())
           

    numeratorL1.Sumw2()
    numeratorL2.Sumw2()
    numeratorL2p5.Sumw2()
    numeratorL3.Sumw2()
    numeratorL1L2.Sumw2()
    numeratorL1L2L2p5.Sumw2()
    numeratorL1_HLT.Sumw2()
    denominator.Sumw2()
    numeratorJet.Sumw2()
    denominatorJet.Sumw2()
    
    print "Events: ", NEvt, ", with Muon: ", NEvtMuon, ", with good MT: ", NEvtMuonMT, ", with tau: ", NEvtMuonMTTau
    print "Matched Taus:   L1: ", NMatchedL1, ", L2: ", NMatchedL2, ", L2.5: ", NMatchedL2p5, ", L3: ", NMatchedL3

    return denominator, numeratorL1, numeratorL2, numeratorL2p5, numeratorL3, numeratorL1L2, numeratorL1L2L2p5, numeratorL1_HLT, deltaZTauTau, mTH, numeratorJet, denominatorJet, invM, invMSel, nTaus


# running the macro

# just to avoid opening windows
ROOT.gROOT.SetBatch()

setTDRStyle()
frame = ROOT.TH1F("frame","",len(xbin)-1,xbins)
frame.SetMinimum(0.01)
frame.SetMaximum(1.01)

### define input
dataset = "Whole2012"
input = []

if dataset == "Whole2012" or dataset == "Run2012A":
    print "Process Run2012A..."
    folder = "dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/user/bkargoll/TriggerEfficiencies/SingleMuDataset/MuTauSkim/HLTrerunDoubleTauJet_v3/"
    files = ["hltDoubleTauJet_offlineSelectedMuTau_3_1_n7X.root",
             "hltDoubleTauJet_offlineSelectedMuTau_5_1_Dzk.root",
             "hltDoubleTauJet_offlineSelectedMuTau_4_1_rwb.root",
             "hltDoubleTauJet_offlineSelectedMuTau_7_1_xgB.root",
             "hltDoubleTauJet_offlineSelectedMuTau_2_1_Aey.root",
             "hltDoubleTauJet_offlineSelectedMuTau_1_1_ETO.root",
             "hltDoubleTauJet_offlineSelectedMuTau_6_1_I6O.root",
             "hltDoubleTauJet_offlineSelectedMuTau_8_1_p0U.root"]
    for file in files:
        input.append(folder+file)

if dataset == "Whole2012" or dataset == "Run2012B":
    print "Process Run2012B Part 1..."
    folder = "dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/user/bkargoll/TriggerEfficiencies/SingleMuDataset/MuTauSkim/Run2012B/FirstPart/HLTrerunDoubleTauJet_v3/"
    files = ["hltDoubleTauJet_offlineSelectedMuTau_3_1_zKI.root",
             "hltDoubleTauJet_offlineSelectedMuTau_6_1_xJO.root", 
             "hltDoubleTauJet_offlineSelectedMuTau_5_1_HEZ.root", 
             "hltDoubleTauJet_offlineSelectedMuTau_7_1_sAD.root", 
             "hltDoubleTauJet_offlineSelectedMuTau_2_1_LgX.root", 
             "hltDoubleTauJet_offlineSelectedMuTau_4_1_rur.root", 
             "hltDoubleTauJet_offlineSelectedMuTau_1_1_d6N.root", 
             "hltDoubleTauJet_offlineSelectedMuTau_8_1_UDb.root"] 
    for file in files:
        input.append(folder+file)     
    print "Process Run2012B Part 2..."
    folder = "dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/user/bkargoll/TriggerEfficiencies/SingleMuDataset/MuTauSkim/Run2012B/SecondPart/HLTrerunDoubleTauJet_v3/"
    files = ["hltDoubleTauJet_offlineSelectedMuTau_3_1_bxj.root",
             "hltDoubleTauJet_offlineSelectedMuTau_1_1_wS5.root", 
             "hltDoubleTauJet_offlineSelectedMuTau_5_1_PN8.root", 
             "hltDoubleTauJet_offlineSelectedMuTau_4_1_SJ5.root",
             "hltDoubleTauJet_offlineSelectedMuTau_2_1_f1K.root"]
    for file in files:
        input.append(folder+file)
    print "Process Run2012B Part 3..."
    folder = "dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/user/bkargoll/TriggerEfficiencies/SingleMuDataset/MuTauSkim/Run2012B/ThirdPart/HLTrerunDoubleTauJet_v3/"
    files = ["hltDoubleTauJet_offlineSelectedMuTau_1_1_dfk.root"]
    for file in files:
        input.append(folder+file)


offlineTauLabel = "offlineSelectedTausXXMisodbMelecrej" #"offlineSelectedTausH2TauJetMVAEle" #offlineSelectedTausMuTauBaseline
l1Label = "hltL1extraParticles"
l2Label = "hltL2TauJets"
l2p5Label = "hltL2TauJetsIso"
l3Label = "hltSelectedPFTausTrackPt5MediumIsolation"
l2Pt = 30
l3Pt = 30
MT = 40
nProngs = 5
dZ = 0.2
offlProngs = 0 # 0 = all, 1 =1prongs, 3=3prongs
invMCut = 9999999



myhistos = histograms(input,offlineTauLabel,l1Label,l2Label,l2p5Label,l3Label,l2Pt,l3Pt, MT, nProngs, dZ, offlProngs, invMCut)
# c1= ROOT.TCanvas()
# c1.Divide(1,5)
filename = files[0].split("/")[-1]
filename = filename.split(".")[0]
filename = filename.split("_")[0]+filename.split("_")[1]
if offlProngs == 0: offlProngs = "all"
else: offlProngs = str(offlProngs)
newfile = ROOT.TFile("/user/kargoll/results/TurnOns/turnOn_"+dataset+"_"+filename+"_"+offlineTauLabel+"_l2Pt"+str(l2Pt)+"_l3Pt"+str(l3Pt)+"_TrackPt5MediumIso_dZ"+str(dZ)+"_nProngs"+str(nProngs)+"_noWJetsMT"+str(MT)+"_"+offlProngs+"OfflineProngs_maxInvM"+str(invMCut)+".root","RECREATE")
# myhistos[i+1].Divide(myhistos[i+1],myhistos[0],1,1,"B")
# frame.DrawCopy()
# myhistos[i+1].DrawCopy("pE2same")
# myhistos[i+1].Write()
for i in range(15):
     myhistos[i].Write()
newfile.Close()     

# c1.SaveAs("turnOn_"+filename+".png")

