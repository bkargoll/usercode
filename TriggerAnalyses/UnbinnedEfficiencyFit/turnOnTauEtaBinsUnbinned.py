import ROOT
from DataFormats.FWLite import Events, Handle
import math
import array
import imp
from tdrstyle import setTDRStyle
from time import time


def deltaPhi(phi1, phi2):
     PHI = abs(phi1 - phi2)
     if (PHI <= 3.14159265):
         return PHI
     else:
         return 2 * 3.14159265 - PHI
     
def deltaPhiWithSign(phi1, phi2):
    PI = 3.14159265
    if phi1 < 0 :  phi1 = phi1 + 2 * PI
    if phi2 < 0 :  phi2 = phi2 + 2 * PI
    deltaphi = phi1 - phi2
    if deltaphi > 2 * PI : deltaphi = deltaphi - 2 * PI
    if deltaphi < -2 * PI : deltaphi = deltaphi + 2 * PI
    if deltaphi > PI : deltaphi = -(2 * PI - deltaphi)
    if deltaphi < -PI: deltaphi = (2 * PI + deltaphi)
    return deltaphi;


def deltaR(eta1, phi1, eta2, phi2) :
    deta = eta1 - eta2
    dphi = deltaPhi(phi1, phi2)
    return math.sqrt(deta * deta + dphi * dphi)


eventCount = 0


def histograms(fileName, labelOffline, labelL1, labelL2, labelL2p5, labelL3, L1TauPtCut=0, L1JetPtCut=0, L2PtCut=25., L3PtCut=25., MTCut=40., invertMT=False,
               nProngsCut=5, dZCut=0.2, offlineProngs=0, maxInvM=10000., offlinePtCut=40., etaBinList=[]):

    events = Events(fileName)
    # vertex
    VertexHandle = Handle('std::vector<reco::Vertex>')
    vertexLabel = "offlinePrimaryVertices"
    # vertexLabel = "hltPixelVertices"

    # taus
    TauHandle = Handle('std::vector<reco::PFTau>')
    JetHandle = Handle('std::vector<reco::CaloJet>')
    L1Handle = Handle('std::vector<l1extra::L1JetParticle>')
    METHandle = Handle('std::vector<reco::PFMET>')
    PFJetHandle = Handle('std::vector<reco::PFJet>')
    MuonHandle = Handle('std::vector<reco::Muon>')
    
    # parse eta bin information
    if etaBinList == []: etaBinList = [[0.0, 2.1, "Total"]]
    
    # book histograms
    deltaZTauTau = ROOT.TH1F("deltaZTauTau", "", 200, 0, 10.)
    mTH = ROOT.TH1F("mTH", "", 100, 0, 100)
    invM = ROOT.TH1F("invM", "", 200, 0, 200)
    invMSel = ROOT.TH1F("invMSel", "", 200, 0, 200)
    nTaus = ROOT.TH1I("NTaus", "", 20, -0.5, 19.5)
    deltaRL1L3 = ROOT.TH1F("deltaRL1L3", "", 200, 0., 5.)
    tauPt = ROOT.TH1F("tauPt", "", 500, 0., 500.)
    tauEta = ROOT.TH1F("tauEta", "", 500, -2.5, 2.5)
    # efficiency graphs
    effL1 = []
    effL2 = []
    effL1L2 = []
    effL2p5 = []
    effL1L2L2p5 = []
    effL2L2p5 = []
    effL3 = []
    effL1_HLT = []
    # efficiencies over eta
    effEtaL1 = []
    effEtaL2 = []
    effEtaL1L2 = []
    effEtaL2p5 = []
    effEtaL1L2L2p5 = []
    effEtaL2L2p5 = []
    effEtaL3 = []
    effEtaL1_HLT = []
    
    # resolution histos
    resEtaL1 = []
    resPhiL1 = []
    resRL1 = []
    resPtL1 = []
    resRelptL1 = []
    resEtaL2 = []
    resPhiL2 = []
    resRL2 = []
    resPtL2 = []
    resRelptL2 = []
    resEtaL2p5 = []
    resPhiL2p5 = []
    resRL2p5 = []
    resPtL2p5 = []
    resRelptL2p5 = []
    resEtaL3 = []
    resPhiL3 = []
    resRL3 = []
    resPtL3 = []
    resRelptL3 = []
    
    NEvt = 0
    NEvtMuon = 0
    NEvtMuonMT = 0
    NEvtMuonMTTau = 0
    NTauFilled = 0
    
    for bin in etaBinList:
        arrayLength = 1000
        
        tmpL1 = ROOT.TGraph(arrayLength)
        tmpL1.SetName("effL1_" + bin[2])
        effL1.append(tmpL1)
        
        tmpL2 = ROOT.TGraph(arrayLength)
        tmpL2.SetName("effL2_" + bin[2])
        effL2.append(tmpL2)
        
        tmpL1L2 = ROOT.TGraph(arrayLength)
        tmpL1L2.SetName("effL1L2_" + bin[2])
        effL1L2.append(tmpL1L2)
        
        tmpL2p5 = ROOT.TGraph(arrayLength)
        tmpL2p5.SetName("effL2p5_" + bin[2])
        effL2p5.append(tmpL2p5)
        
        tmpL1L2L2p5 = ROOT.TGraph(arrayLength)
        tmpL1L2L2p5.SetName("effL1L2L2p5_" + bin[2])
        effL1L2L2p5.append(tmpL1L2L2p5)
        
        tmpL2L2p5 = ROOT.TGraph(arrayLength)
        tmpL2L2p5.SetName("effL2L2p5_" + bin[2])
        effL2L2p5.append(tmpL2L2p5)
        
        tmpL3 = ROOT.TGraph(arrayLength)
        tmpL3.SetName("effL3_" + bin[2])
        effL3.append(tmpL3)
        
        tmpL1_HLT = ROOT.TGraph(arrayLength)
        tmpL1_HLT.SetName("effL1L2L2p5L3_" + bin[2])
        effL1_HLT.append(tmpL1_HLT)
        
        # over eta
        tmpEtaL1 = ROOT.TGraph(arrayLength)
        tmpEtaL1.SetName("effEtaL1_" + bin[2])
        effEtaL1.append(tmpEtaL1)
        
        tmpEtaL2 = ROOT.TGraph(arrayLength)
        tmpEtaL2.SetName("effEtaL2_" + bin[2])
        effEtaL2.append(tmpEtaL2)
        
        tmpEtaL1L2 = ROOT.TGraph(arrayLength)
        tmpEtaL1L2.SetName("effEtaL1L2_" + bin[2])
        effEtaL1L2.append(tmpEtaL1L2)
        
        tmpEtaL2p5 = ROOT.TGraph(arrayLength)
        tmpEtaL2p5.SetName("effEtaL2p5_" + bin[2])
        effEtaL2p5.append(tmpEtaL2p5)
        
        tmpEtaL1L2L2p5 = ROOT.TGraph(arrayLength)
        tmpEtaL1L2L2p5.SetName("effEtaL1L2L2p5_" + bin[2])
        effEtaL1L2L2p5.append(tmpEtaL1L2L2p5)
        
        tmpEtaL2L2p5 = ROOT.TGraph(arrayLength)
        tmpEtaL2L2p5.SetName("effEtaL2L2p5_" + bin[2])
        effEtaL2L2p5.append(tmpEtaL2L2p5)
        
        tmpEtaL3 = ROOT.TGraph(arrayLength)
        tmpEtaL3.SetName("effEtaL3_" + bin[2])
        effEtaL3.append(tmpEtaL3)
        
        tmpEtaL1_HLT = ROOT.TGraph(arrayLength)
        tmpEtaL1_HLT.SetName("effEtaL1L2L2p5L3_" + bin[2])
        effEtaL1_HLT.append(tmpEtaL1_HLT)
        

        
        resEtaL1.append(ROOT.TH1F("resEtaL1_" + bin[2], "", 400, -1.0, 1.0))   
        resPhiL1.append(ROOT.TH1F("resPhiL1_" + bin[2], "", 400, -1., 1.))   
        resRL1.append(ROOT.TH1F("resRL1_" + bin[2], "", 400, 0., 1.))      
        resPtL1.append(ROOT.TH1F("resPtL1_" + bin[2], "", 1000, -100., 100.))     
        resRelptL1.append(ROOT.TH1F("resRelptL1_" + bin[2], "", 1000, -5., 5.))  
        resEtaL2.append(ROOT.TH1F("resEtaL2_" + bin[2], "", 400, -1.0, 1.0)) 
        resPhiL2 .append(ROOT.TH1F("resPhiL2_" + bin[2], "", 400, -1., 1.))   
        resRL2.append(ROOT.TH1F("resRL2_" + bin[2], "", 400, 0., 1.))       
        resPtL2.append(ROOT.TH1F("resPtL2_" + bin[2], "", 1000, -100., 100.))      
        resRelptL2.append(ROOT.TH1F("resRelptL2_" + bin[2], "", 1000, -5., 5.))   
        resEtaL2p5.append(ROOT.TH1F("resEtaL2p5_" + bin[2], "", 400, -1.0, 1.0))  
        resPhiL2p5.append(ROOT.TH1F("resPhiL2p5_" + bin[2], "", 400, -1., 1.))  
        resRL2p5.append(ROOT.TH1F("resRL2p5_" + bin[2], "", 400, 0., 1.))     
        resPtL2p5.append(ROOT.TH1F("resPtL2p5_" + bin[2], "", 1000, -100., 100.))    
        resRelptL2p5.append(ROOT.TH1F("resRelptL2p5_" + bin[2], "", 1000, -5., 5.)) 
        resEtaL3.append(ROOT.TH1F("resEtaL3_" + bin[2], "", 400, -1.0, 1.0))    
        resPhiL3.append(ROOT.TH1F("resPhiL3_" + bin[2], "", 400, -1., 1.))    
        resRL3.append(ROOT.TH1F("resRL3_" + bin[2], "", 400, 0., 1.))       
        resPtL3.append(ROOT.TH1F("resPtL3_" + bin[2], "", 1000, -100., 100.))      
        resRelptL3.append(ROOT.TH1F("resRelptL3_" + bin[2], "", 1000, -5., 5.))   
        
   

    print "Start event loop..."
    for event in events:
        if NEvt == 0: print "    ...started"
        NEvt += 1
        if NEvt % 5000 == 0: print "Process Event ", NEvt
        #if NEvt > 2000: break
         
        # getting the handle from the event
        event.getByLabel(labelOffline, TauHandle)
        offlineTaus = TauHandle.product()
        
        event.getByLabel("ak5PFJets", PFJetHandle)
        offlineJets = PFJetHandle.product()
        
        event.getByLabel(labelL1, "Tau", L1Handle)
        l1Taus = L1Handle.product()             

        event.getByLabel(labelL1, "Central", L1Handle)
        l1Jets = L1Handle.product()             

        event.getByLabel(labelL2, JetHandle)
        l2Taus = JetHandle.product()

        # event.getByLabel("hltCaloJetL1FastJetCorrected",JetHandle)
        # l2Jets = JetHandle.product()

        event.getByLabel(labelL2p5, JetHandle)
        l2p5Taus = JetHandle.product()             

        event.getByLabel(labelL3, TauHandle)
        l3Taus = TauHandle.product()

        event.getByLabel("pfMet", METHandle)
        met = METHandle.product().at(0)
        
        event.getByLabel("offlineSelectedMuons", MuonHandle)
        muons = MuonHandle.product()
        
        event.getByLabel(vertexLabel, VertexHandle)
        offlinevertex = VertexHandle.product().at(0).position()
        vertexSize = VertexHandle.product().size()

        # matching variables
        matchingConeL1 = 0.5
        matchingConeHLT = 0.3

        # making MT variable
        MT = 0.
        myMuons = []
        for muon in muons:
             if muon.pt() < 10: continue
             if abs(muon.eta()) > 2.1: continue
             # apply vertex cuts which have been missing in muon selection before
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
        MT = (muon.pt() + met.pt()) * (muon.pt() + met.pt()) - (muon.px() + met.px()) * (muon.px() + met.px()) - (muon.py() + met.py()) * (muon.py() + met.py())
        MT = math.sqrt(MT)
        # Remove W+Jets
        mTH.Fill(MT)
        # normal MTCut
        if not invertMT:
            if MT > MTCut: continue
        else:
            # inverse MT cut
            if MT < MTCut: continue
                
        NEvtMuonMT += 1

   
        # select best tau
        cntTaus = 0
        tau = None
        bestZmass = 999999.0
        for iTau in offlineTaus:
            cntTaus += 1
            if iTau.pt() < 5: continue
            if abs(iTau.eta()) > 2.1: continue
            # check: select only taus in 2nd peak structure 
            # if ((iTau.pt() < 55) or  (iTau.pt() > 100)): continue
            # muon and tau with opposite charge
            if iTau.charge() != -1 * muon.charge(): continue
            # require seperation between muon and tau
            dR_mutau = deltaR(iTau.eta(), iTau.phi(), muon.eta(), muon.phi())
            if dR_mutau < 0.5 : continue
            # make sure tau has a leading track, which is needed later for the dZ cut
            if not iTau.leadPFChargedHadrCand().trackRef().isNonnull(): continue
            # offline taus: only 1prongs, only 3prongs or all?
            if not offlineProngs == 0:
                if not iTau.signalPFChargedHadrCands().size() == offlineProngs: continue
            NEvtMuonMTTau += 1
            Zmass = 2 * (iTau.mass() * muon.mass() + muon.energy() * iTau.energy() - muon.px() * iTau.px() - muon.py() * iTau.py() - muon.pz() * iTau.pz())
            Zmass = math.sqrt(Zmass)
            invM.Fill(Zmass)
            if tau == None:
                bestZmass = Zmass
                tau = iTau
                continue
            if iTau.pt() > tau.pt():
                bestZmass = Zmass
                tau = iTau
        
        # skip Z->mumu events, where 1 mu fakes a tau
        if ((maxInvM != 0) & (bestZmass > maxInvM)): continue
        
        nTaus.Fill(cntTaus)
        if tau == None: continue
        invMSel.Fill(bestZmass)
        
         # check which eta bins to run on
        doEtaBin = []
        for bin in etaBinList:
            if (abs(tau.eta()) >= bin[0]) & (abs(tau.eta()) <= bin[1]):
                doEtaBin.append(True)
            else: doEtaBin.append(False)
        
        NTauFilled += 1
        tauPt.Fill(tau.pt())
        tauEta.Fill(tau.eta())
       
        # loop on numerators
        myL1Taus = []
        myL1Jets = []
        myL2Taus = []
        myL2p5Taus = []
        myL3Taus = []
        
        for l1Tau in l1Taus:
             if l1Tau.pt() < L1TauPtCut: continue
             if abs(l1Tau.eta()) > 2.1: continue
             dr = deltaR(tau.eta(), tau.phi(), l1Tau.eta(), l1Tau.phi())
             myL1Taus.append((dr, l1Tau))            
        for l1Jet in l1Jets:
             if l1Jet.pt() < L1JetPtCut: continue
             if abs(l1Jet.eta()) > 2.1: continue
             dr = deltaR(tau.eta(), tau.phi(), l1Jet.eta(), l1Jet.phi())
             myL1Jets.append((dr, l1Jet))
        for l2Tau in l2Taus:
             if l2Tau.pt() < L2PtCut: continue
             if abs(l2Tau.eta()) > 2.1: continue
             dr = deltaR(tau.eta(), tau.phi(), l2Tau.eta(), l2Tau.phi())
             myL2Taus.append((dr, l2Tau))            
        for l2p5Tau in l2p5Taus:
             if l2p5Tau.pt() < L2PtCut: continue
             if abs(l2p5Tau.eta()) > 2.1: continue
             dr = deltaR(tau.eta(), tau.phi(), l2p5Tau.eta(), l2p5Tau.phi())
             myL2p5Taus.append((dr, l2p5Tau))
        for l3Tau in l3Taus:
             if l3Tau.pt() < L3PtCut: continue
             if abs(l3Tau.eta()) > 2.1: continue
             dr = deltaR(tau.eta(), tau.phi(), l3Tau.eta(), l3Tau.phi())
             myL3Taus.append((dr, l3Tau))    
        
        # store if object found
        foundL1 = False
        foundL2 = False
        foundL2p5 = False
        foundL3 = False
        
        foundEtaL1 = False
        foundEtaL2 = False
        foundEtaL2p5 = False
        foundEtaL3 = False
               
        # Filling L1 condition
        if ((L1TauPtCut == 0) & (L1JetPtCut == 0)):
            foundL1 = True
        else:
            if len(myL1Taus) > 0:
                 myL1Taus.sort()
                 if myL1Taus[0][0] < matchingConeL1:
                     foundL1 = True
                     if tau.pt() > offlinePtCut:
                         foundEtaL1 = True
                     for i in range(len(etaBinList)):
                        if doEtaBin[i] == True:
                            # fill resolution plots
                            resEtaL1[i].Fill(myL1Taus[0][1].eta() - tau.eta())
                            resPhiL1[i].Fill(deltaPhiWithSign(myL1Taus[0][1].phi(), tau.phi()))
                            resRL1[i].Fill(myL1Taus[0][0])
                            resPtL1[i].Fill(myL1Taus[0][1].pt() - tau.pt())
                            resRelptL1[i].Fill((myL1Taus[0][1].pt() - tau.pt()) / tau.pt())
                    

            if len(myL1Jets) > 0:
                 myL1Jets.sort()
                 if myL1Jets[0][0] < matchingConeL1:
                     foundL1 = True
                     if tau.pt() > offlinePtCut:
                         foundEtaL1 = True
                     for i in range(len(etaBinList)):
                        if doEtaBin[i] == True:
                            # fill resolution plots
                            resEtaL1[i].Fill(myL1Jets[0][1].eta() - tau.eta())
                            resPhiL1[i].Fill(deltaPhiWithSign(myL1Jets[0][1].phi(), tau.phi()))
                            resRL1[i].Fill(myL1Jets[0][0])
                            resPtL1[i].Fill(myL1Jets[0][1].pt() - tau.pt())
                            resRelptL1[i].Fill((myL1Jets[0][1].pt() - tau.pt()) / tau.pt())
                     

        # Filling L2 condition
        if len(myL2Taus) > 0:
             myL2Taus.sort()
             if myL2Taus[0][0] < matchingConeHLT:
                 foundL2 = True
                 if tau.pt() > offlinePtCut:
                     foundEtaL2 = True       
                 for i in range(len(etaBinList)):
                    if doEtaBin[i] == True:
                        # fill resolution plots
                        resEtaL2[i].Fill(myL2Taus[0][1].eta() - tau.eta())
                        resPhiL2[i].Fill(deltaPhiWithSign(myL2Taus[0][1].phi(), tau.phi()))
                        resRL2[i].Fill(myL2Taus[0][0])
                        resPtL2[i].Fill(myL2Taus[0][1].pt() - tau.pt())
                        resRelptL2[i].Fill((myL2Taus[0][1].pt() - tau.pt()) / tau.pt())
                 
        # Filling L2.5 condition
        if len(myL2p5Taus) > 0:
             myL2p5Taus.sort()
             if myL2p5Taus[0][0] < matchingConeHLT:
                 foundL2p5 = True
                 if tau.pt() > offlinePtCut:
                     foundEtaL2p5 = True
                 for i in range(len(etaBinList)):
                    if doEtaBin[i] == True:
                        # fill resolution plots
                        resEtaL2p5[i].Fill(myL2p5Taus[0][1].eta() - tau.eta())
                        resPhiL2p5[i].Fill(deltaPhiWithSign(myL2p5Taus[0][1].phi(), tau.phi()))
                        resRL2p5[i].Fill(myL2p5Taus[0][0])
                        resPtL2p5[i].Fill(myL2p5Taus[0][1].pt() - tau.pt())
                        resRelptL2p5[i].Fill((myL2p5Taus[0][1].pt() - tau.pt()) / tau.pt())
                 
         # Filling L3 condition         
        if len(myL3Taus) > 0:
             myL3Taus.sort()
             if myL3Taus[0][0] < matchingConeHLT:
                 if myL3Taus[0][1].leadPFChargedHadrCand().isNonnull():
                      if ((nProngsCut == 0) | (myL3Taus[0][1].signalPFChargedHadrCands().size() < nProngsCut)):
                           if myL3Taus[0][1].leadPFChargedHadrCand().trackRef().isNonnull():
                               # create histogram to check/mimic impact on L3 dZ cut
                               dZ = abs(tau.leadPFChargedHadrCand().trackRef().dz(offlinevertex) - myL3Taus[0][1].leadPFChargedHadrCand().trackRef().dz(offlinevertex))
                               deltaZTauTau.Fill(dZ)
                               if dZ < dZCut:
                                   # create histogram to check/mimic impact of matching between L3 and L1 objects at HLT
                                   dRl1l3 = matchL3toL1(myL3Taus[0], myL1Taus, myL1Jets)
                                   deltaRL1L3.Fill(dRl1l3)
                                   # if dRl1l3 < 0.5: # matching cut of 0.5 is hard-coded in HLT filter module
                                   foundL3 = True
                                   if tau.pt() > offlinePtCut:
                                       foundEtaL3 = True
                                   for i in range(len(etaBinList)):
                                       if doEtaBin[i] == True:
                                           # fill resolution plots
                                           resEtaL3[i].Fill(myL3Taus[0][1].eta() - tau.eta())
                                           resPhiL3[i].Fill(deltaPhiWithSign(myL3Taus[0][1].phi(), tau.phi()))
                                           resRL3[i].Fill(myL3Taus[0][0])
                                           resPtL3[i].Fill(myL3Taus[0][1].pt() - tau.pt())
                                           resRelptL3[i].Fill((myL3Taus[0][1].pt() - tau.pt()) / tau.pt())    
                                   # else:
                                   #    if (foundL1 & foundL2 & foundL2p5):
                                   #        print "*** L1L3 matching failed, min(dR)=", dRl1l3
        
        # fill efficiency graphs
        for i in range(len(etaBinList)):
            if doEtaBin[i] == True:
                effL1[i].SetPoint(NTauFilled,tau.pt(),  foundL1)
                effL2[i].SetPoint(NTauFilled,tau.pt(),  foundL2)
                effL2p5[i].SetPoint(NTauFilled,tau.pt(),foundL2p5)
                effL3[i].SetPoint(NTauFilled,tau.pt(),  foundL3)
                
                effL1L2[i].SetPoint(NTauFilled, tau.pt(),    (foundL1 and foundL2))
                effL1L2L2p5[i].SetPoint(NTauFilled, tau.pt(),(foundL1 and foundL2 and foundL2p5))
                effL2L2p5[i].SetPoint(NTauFilled, tau.pt(),  (foundL2 and foundL2p5))
                effL1_HLT[i].SetPoint(NTauFilled, tau.pt(),  (foundL1 and foundL2 and foundL2p5 and foundL3))
    
                effEtaL1[i].SetPoint(NTauFilled, tau.pt(),  foundEtaL1)  
                effEtaL2[i].SetPoint(NTauFilled, tau.pt(),  foundEtaL2)  
                effEtaL2p5[i].SetPoint(NTauFilled, tau.pt(),foundEtaL2p5)
                effEtaL3[i].SetPoint(NTauFilled, tau.pt(),  foundEtaL3)  
                
                effEtaL1L2[i].SetPoint(NTauFilled, tau.pt(),    (foundEtaL1 and foundEtaL2))                          
                effEtaL1L2L2p5[i].SetPoint(NTauFilled, tau.pt(),(foundEtaL1 and foundEtaL2 and foundEtaL2p5))            
                effEtaL2L2p5[i].SetPoint(NTauFilled, tau.pt(),  (foundEtaL2 and foundEtaL2p5))                        
                effEtaL1_HLT[i].SetPoint(NTauFilled, tau.pt(),  (foundEtaL1 and foundEtaL2 and foundEtaL2p5 and foundEtaL3))
                

  
    eventCount = NEvt
    
    print "Events: ", NEvt, ", with Muon: ", NEvtMuon, ", with good MT: ", NEvtMuonMT, ", with tau: ", NEvtMuonMTTau

    return  effL1, effL2, effL2p5, effL3, effL1L2, effL1L2L2p5, effL2L2p5, effL1_HLT, \
            effEtaL1, effEtaL2, effEtaL2p5, effEtaL3, effEtaL1L2, effEtaL1L2L2p5, effEtaL2L2p5, effEtaL1_HLT, \
            resEtaL1, resPhiL1, resRL1, resPtL1, resRelptL1, resEtaL2, resPhiL2, resRL2, resPtL2, resRelptL2, resEtaL2p5, resPhiL2p5, \
            resRL2p5, resPtL2p5, resRelptL2p5, resEtaL3, resPhiL3, resRL3, resPtL3, resRelptL3, \
            deltaZTauTau, mTH, invM, invMSel, nTaus, deltaRL1L3, tauPt, tauEta

def matchL3toL1(l3Tau, l1Taus, l1Jets):
    minDR = 100
    for l1Tau in l1Taus:
        dR = deltaR(l1Tau[1].eta(), l1Tau[1].phi(), l3Tau[1].eta(), l3Tau[1].phi())
        if dR < minDR: minDR = dR
    for l1Jet in l1Jets:
        dR = deltaR(l1Jet[1].eta(), l1Jet[1].phi(), l3Tau[1].eta(), l3Tau[1].phi())
        if dR < minDR: minDR = dR
    return minDR
    

def executeCode(cfgFileName):
    cfgFile = open(cfgFileName, 'r')
    cfg = imp.load_source('configuration', cfgFileName, cfgFile)
    cfgFile.close()
    
    print "++++ Process dataset ", cfg.dataset
    
    # settings
    invertMtCut = False
    if invertMtCut: print "Warning: mT cut is inversed"
    
    # parse info for file name
    offlProngsLabel = ""
    if cfg.offlProngs == 0: offlProngsLabel = ""
    else: offlProngsLabel = "_" + str(cfg.offlProngs) + "OfflineProngs"
    if cfg.invMCut == 0: invMLabel = ""
    else: invMLabel = "_maxInvM" + str(cfg.invMCut)
    if cfg.nProngs == 0: prongLabel = ""
    else: prongLabel = "_" + str(cfg.nProngs) + "Prongs"
    mtLabel = ""
    if invertMtCut: mtLabel = "Inversed"
    offlineLabel = cfg.offlineTauLabel.split("XX")[-1]
    ptThresholds = "ptThr"
    if ((cfg.l1tauPt != 0) | (cfg.l1jetPt != 0)): ptThresholds += ("-" + str(cfg.l1tauPt) + "-" + str(cfg.l1jetPt))
    ptThresholds += "-" + str(cfg.l2Pt) + "-" + str(cfg.l3Pt)
    outputFileName = "/user/kargoll/results/TurnOns/forMoriond2013/" + cfg.path + "/turnOnUnbinned_" + cfg.dataset + "_" + offlineLabel + "_L3" + cfg.l3Label[-10:] + "_" + ptThresholds + "_dZ" + str(cfg.dZ) + prongLabel + "_MT" + str(cfg.MT) + mtLabel + offlProngsLabel + invMLabel + "_smB.root"



    myhistos = histograms(cfg.input, cfg.offlineTauLabel, cfg.l1Label, cfg.l2Label, cfg.l2p5Label, cfg.l3Label, cfg.l1tauPt, cfg.l1jetPt, cfg.l2Pt, cfg.l3Pt, cfg.MT, invertMtCut, cfg.nProngs, cfg.dZ, cfg.offlProngs, cfg.invMCut, cfg.offlinePt, cfg.etaBins)
    newfile = ROOT.TFile(outputFileName, "RECREATE")
    nHistos = len(myhistos)
    for i in range(nHistos - 8):
        for j in range(len(cfg.etaBins)):
            myhistos[i][j].Write()
    for i in range(nHistos - 8, nHistos):
        myhistos[i].Write()
    newfile.Close()
    
    print nHistos, " histograms have been saved to ", outputFileName



# running the macro

# just to avoid opening windows
print "Let's go!"
ROOT.gROOT.SetBatch()

setTDRStyle(1)

# MSSM Higgs
# executeCode("config_DoubleTauProng1_2012A.py")
#executeCode("config_DoubleTauProng1_2012B.py")
#executeCode("config_DoubleTauProng1_2012Cv3.py")
#executeCode("config_DoubleTauProng1_2012Cv4.py")
#executeCode("config_DoubleTauProng1_2012Dpart1.py")
#executeCode("config_DoubleTauProng1_2012Dpart2.py")
# ## SM Higgs
executeCode("config_DoubleTauJet_2012A.py")
executeCode("config_DoubleTauJet_2012B.py")
executeCode("config_DoubleTauJet_2012Cv3.py")
executeCode("config_DoubleTauJet_2012Cv4.py")
executeCode("config_DoubleTauJet_2012Dpart1.py")
executeCode("config_DoubleTauJet_2012Dpart2.py")
