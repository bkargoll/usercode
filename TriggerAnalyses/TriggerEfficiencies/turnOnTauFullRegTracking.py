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


def histograms(fileName, labelOffline, labelL1, labelL2,labelL2p5, labelL3Full, labelL3Reg, L2PtCut = 25., L3PtCut = 25., MTCut = 40., dZCut = 0.2, offlineProngs = 0):

    events = Events(fileName)
    # vertex
    VertexHandle = Handle('std::vector<reco::Vertex>')
    vertexLabel = "offlinePrimaryVertices"

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
    numeratorL3Full = ROOT.TH1F("numeratorL3Full","",len(xbin)-1,xbins)
    numeratorL3Reg = ROOT.TH1F("numeratorL3Reg","",len(xbin)-1,xbins)
    numeratorL1_HLTFull = ROOT.TH1F("numeratorL1_HLTFull","",len(xbin)-1,xbins)
    numeratorL1_HLTReg = ROOT.TH1F("numeratorL1_HLTReg","",len(xbin)-1,xbins)
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
    NMatchedL3Full = 0   
    NMatchedL3Reg = 0   

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

        event.getByLabel(labelL3Full,TauHandle)
        l3TausFull = TauHandle.product()
        
        event.getByLabel(labelL3Reg,TauHandle)
        l3TausReg = TauHandle.product()       

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
        foundL3Full = False
        foundL3Reg = False
        
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
        mTau = 1.77682
        #mZ = 91.188
        bestZmass = 999999.0
        for iTau in offlineTaus:
            cntTaus += 1
            if iTau.pt() < 5: continue
            if abs(iTau.eta()) > 2.1: continue
            # muon and tau with opposite charge
            if iTau.charge() != -1 * muon.charge(): continue
            # offline taus: only 1prongs, only 3prongs or all?
            if not offlineProngs == 0:
                if not iTau.signalPFChargedHadrCands().size() == offlineProngs: continue
            NEvtMuonMTTau += 1
            Zmass = 2*(mTau*mTau + muon.energy()*iTau.energy() - muon.px()*iTau.px() - muon.py()*iTau.py() - muon.pz()*iTau.pz())
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
        #if bestZmass > maxInvM: continue
        
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
        myL3TausFull = []
        myL3TausReg = []
        for l1Tau in l1Taus:
             if l1Tau.pt() < 44: continue
             dr = deltaR(tau.eta(),tau.phi(),l1Tau.eta(),l1Tau.phi())
             myL1Taus.append((dr, l1Tau))            
        for l1Jet in l1Jets:
             if l1Jet.pt() < 64: continue
             dr = deltaR(tau.eta(),tau.phi(),l1Jet.eta(),l1Jet.phi())
             myL1Jets.append((dr,l1Jet))
        for l2Tau in l2Taus:
             if l2Tau.pt() < L2PtCut: continue
             dr = deltaR(tau.eta(),tau.phi(),l2Tau.eta(),l2Tau.phi())
             myL2Taus.append((dr, l2Tau))            
        for l2p5Tau in l2p5Taus:
             if l2p5Tau.pt() < L2PtCut: continue
             dr = deltaR(tau.eta(),tau.phi(),l2p5Tau.eta(),l2p5Tau.phi())
             myL2p5Taus.append((dr, l2p5Tau))            
        for l3TauFull in l3TausFull:
             if l3TauFull.pt() < L3PtCut: continue
             dr = deltaR(tau.eta(),tau.phi(),l3TauFull.eta(),l3TauFull.phi())
             myL3TausFull.append((dr, l3TauFull))
        for l3TauReg in l3TausReg:
             if l3TauReg.pt() < L3PtCut: continue
             dr = deltaR(tau.eta(),tau.phi(),l3TauReg.eta(),l3TauReg.phi())
             myL3TausReg.append((dr, l3TauReg))                

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
        # Filling L3 condition                  
        if len(myL3TausFull) >0:
             myL3TausFull.sort()
             if myL3TausFull[0][0] < matchingConeHLT:
                 if myL3TausFull[0][1].leadPFChargedHadrCand().isNonnull():
                       #if myL3TausFull[0][1].signalPFChargedHadrCands().size() < nProngsCut:
                           if tau.leadPFChargedHadrCand().trackRef().isNonnull():
                               if myL3TausFull[0][1].leadPFChargedHadrCand().trackRef().isNonnull():
                                   dZ = abs(tau.leadPFChargedHadrCand().trackRef().dz(offlinevertex) - myL3TausFull[0][1].leadPFChargedHadrCand().trackRef().dz(offlinevertex))
                                   deltaZTauTau.Fill(dZ)
                                   if dZ < dZCut:
                                       foundL3Full = True
                                       numeratorL3Full.Fill(tau.pt())
                                       NMatchedL3Full += 1
                   

        if len(myL3TausReg) >0:
             myL3TausReg.sort()
             if myL3TausReg[0][0] < matchingConeHLT:
                 if myL3TausReg[0][1].leadPFChargedHadrCand().isNonnull():
                       #if myL3TausReg[0][1].signalPFChargedHadrCands().size() < nProngsCut:
                           if tau.leadPFChargedHadrCand().trackRef().isNonnull():
                               if myL3TausReg[0][1].leadPFChargedHadrCand().trackRef().isNonnull():
                                   dZ = abs(tau.leadPFChargedHadrCand().trackRef().dz(offlinevertex) - myL3TausReg[0][1].leadPFChargedHadrCand().trackRef().dz(offlinevertex))
                                   deltaZTauTau.Fill(dZ)
                                   if dZ < dZCut:
                                       foundL3Reg = True
                                       numeratorL3Reg.Fill(tau.pt())
                                       NMatchedL3Reg += 1

        # fill plots for accumulated efficiencies
        if foundL1 and foundL2 :
            numeratorL1L2.Fill(tau.pt())
            if foundL2p5:
                numeratorL1L2L2p5.Fill(tau.pt())
                if foundL3Full:
                    numeratorL1_HLTFull.Fill(tau.pt())
                if foundL3Reg:
                    numeratorL1_HLTReg.Fill(tau.pt())
           

    numeratorL1.Sumw2()
    numeratorL2.Sumw2()
    numeratorL2p5.Sumw2()
    numeratorL3Full.Sumw2()
    numeratorL3Reg.Sumw2()
    numeratorL1L2.Sumw2()
    numeratorL1L2L2p5.Sumw2()
    numeratorL1_HLTFull.Sumw2()
    numeratorL1_HLTReg.Sumw2()
    denominator.Sumw2()
    numeratorJet.Sumw2()
    denominatorJet.Sumw2()
    
    print "Events: ", NEvt, ", with Muon: ", NEvtMuon, ", with good MT: ", NEvtMuonMT, ", with tau: ", NEvtMuonMTTau
    print "Matched Taus:   L1: ", NMatchedL1, ", L2: ", NMatchedL2, ", L2.5: ", NMatchedL2p5, ", L3Full: ", NMatchedL3Full, ", L3Reg: ", NMatchedL3Reg

    return denominator, numeratorL1, numeratorL2, numeratorL2p5, numeratorL3Full, numeratorL3Reg, numeratorL1L2, numeratorL1L2L2p5, numeratorL1_HLTFull, numeratorL1_HLTReg, deltaZTauTau, mTH, numeratorJet, denominatorJet, invM, invMSel, nTaus


# running the macro

# just to avoid opening windows
ROOT.gROOT.SetBatch()

setTDRStyle()
frame = ROOT.TH1F("frame","",len(xbin)-1,xbins)
frame.SetMinimum(0.01)
frame.SetMaximum(1.01)

### input files
## Run2012A
dataset = "Whole2012"
folder = "dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/user/bkargoll/TriggerEfficiencies/SingleMuDataset/MuTauSkim/HLTRerunDoubleTauForMichal/Run2012A/"
files = ["hltDiTauForMichal_offlineSelectedMuTau_6_1_KHx.root",
         "hltDiTauForMichal_offlineSelectedMuTau_4_1_Q0W.root",
         "hltDiTauForMichal_offlineSelectedMuTau_1_1_grf.root",
         "hltDiTauForMichal_offlineSelectedMuTau_5_1_B2E.root",
         "hltDiTauForMichal_offlineSelectedMuTau_8_1_dOY.root",
         "hltDiTauForMichal_offlineSelectedMuTau_9_1_auc.root",
         "hltDiTauForMichal_offlineSelectedMuTau_2_1_6NG.root",
         "hltDiTauForMichal_offlineSelectedMuTau_3_1_1in.root",
         "hltDiTauForMichal_offlineSelectedMuTau_7_1_Mvt.root", 
         "hltDiTauForMichal_offlineSelectedMuTau_10_1_Gyg.root"]

input = []
for file in files:
    input.append(folder+file)
    
## Run2012B FirstPart
folder = "dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/user/bkargoll/TriggerEfficiencies/SingleMuDataset/MuTauSkim/HLTRerunDoubleTauForMichal/FirstPart2012B/"
files = ["hltDiTauForMichal_offlineSelectedMuTau_3_1_p5y.root",
         "hltDiTauForMichal_offlineSelectedMuTau_4_1_niX.root", 
         "hltDiTauForMichal_offlineSelectedMuTau_1_1_lFD.root", 
         "hltDiTauForMichal_offlineSelectedMuTau_7_1_b2h.root", 
         "hltDiTauForMichal_offlineSelectedMuTau_2_1_v4D.root", 
         "hltDiTauForMichal_offlineSelectedMuTau_6_1_pQW.root", 
         "hltDiTauForMichal_offlineSelectedMuTau_5_1_INC.root", 
         "hltDiTauForMichal_offlineSelectedMuTau_8_1_Vq4.root"] 

for file in files:
    input.append(folder+file)

## Run2012B SecondPart
folder = "dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/user/bkargoll/TriggerEfficiencies/SingleMuDataset/MuTauSkim/HLTRerunDoubleTauForMichal/SecondPart2012B/"
files = ["hltDiTauForMichal_offlineSelectedMuTau_4_1_Tgn.root",
         "hltDiTauForMichal_offlineSelectedMuTau_5_1_jyX.root", 
         "hltDiTauForMichal_offlineSelectedMuTau_3_1_3tw.root", 
         "hltDiTauForMichal_offlineSelectedMuTau_1_1_urs.root",
         "hltDiTauForMichal_offlineSelectedMuTau_8_1_sAu.root",
         "hltDiTauForMichal_offlineSelectedMuTau_2_1_1hH.root",
         "hltDiTauForMichal_offlineSelectedMuTau_6_1_TCC.root",
         "hltDiTauForMichal_offlineSelectedMuTau_7_1_fSX.root",
         "hltDiTauForMichal_offlineSelectedMuTau_9_1_vtE.root"] 
for file in files:
    input.append(folder+file)         

## Run2012B ThirdPart
folder = "dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/user/bkargoll/TriggerEfficiencies/SingleMuDataset/MuTauSkim/HLTRerunDoubleTauForMichal/ThirdPart2012B/"
files = ["hltDiTauForMichal_offlineSelectedMuTau_1_1_Aja.root",
         "hltDiTauForMichal_offlineSelectedMuTau_6_1_Geo.root",
         "hltDiTauForMichal_offlineSelectedMuTau_7_1_SxM.root",
         "hltDiTauForMichal_offlineSelectedMuTau_5_1_wIJ.root",
         "hltDiTauForMichal_offlineSelectedMuTau_9_1_MKO.root",
         "hltDiTauForMichal_offlineSelectedMuTau_2_1_aJp.root",
         "hltDiTauForMichal_offlineSelectedMuTau_3_1_TWX.root",
         "hltDiTauForMichal_offlineSelectedMuTau_8_1_pk3.root",         
         "hltDiTauForMichal_offlineSelectedMuTau_4_1_EIH.root",        
         "hltDiTauForMichal_offlineSelectedMuTau_10_1_82M.root"]
for file in files:
    input.append(folder+file)

### define parameters
offlineTauLabel = "offlineSelectedTausXXMisodbMelecrej" 
l1Label = "hltL1extraParticles"
l2Label = "hltL2TauJets"
l2p5Label = "hltL2TauJetsIso"
l3LabelFull = "hltSelectedPFTausTrackPt1MediumIsolationProng4"
l3LabelReg = "hltSelectedPFTausTrackPt1MediumIsolationProng4Reg"
l2Pt = 35
l3Pt = 35
MT = 40
nProngs = "Prong4"
dZ = 0.2
offlProngs = 0 # 0 = all, 1 =1prongs, 3=3prongs

myhistos = histograms(input,offlineTauLabel,l1Label,l2Label,l2p5Label,l3LabelFull, l3LabelReg, l2Pt, l3Pt, MT, dZ, offlProngs)

filename = files[0].split("/")[-1]
filename = filename.split(".")[0]
filename = filename.split("_")[0]+filename.split("_")[1]
if offlProngs == 0: offlProngs = "all"
else: offlProngs = str(offlProngs)
newfile = ROOT.TFile("/user/kargoll/results/TurnOns/DoubleTauForMichal/TESTturnOn_"+dataset+"_"+filename+"_"+offlineTauLabel+"_l2Pt"+str(l2Pt)+"_l3Pt"+str(l3Pt)+"_TrackPt1MediumIso_dZ"+str(dZ)+"_"+nProngs+"_noWJetsMT"+str(MT)+"_"+offlProngs+"OfflineProngs.root","RECREATE")

for i in range(17):
     myhistos[i].Write()
newfile.Close()     


