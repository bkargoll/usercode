import imp
import sys
import os

import math
from array import array
import ROOT

from tdrstyle import *

# tau turn ons
def eff2012IsoTau12fb(pt, eta):
    return (808.411*(0.764166*0.5*(ROOT.TMath.Erf((pt-33.2236)/2./0.97289/math.sqrt(pt))+1.))+ \
            4428.0*(0.802387*0.5*(ROOT.TMath.Erf((pt-38.0971)/2./0.82842/math.sqrt(pt))+1.))+ \
            1783.003*(0.818051*0.5*(ROOT.TMath.Erf((pt-37.3669)/2./0.74847/math.sqrt(pt))+1.))+ \
            5109.155*(0.796086*0.5*(ROOT.TMath.Erf((pt-37.3302)/2./0.757558/math.sqrt(pt))+1.)) \
            )/(808.411+4428.0+1783.003+5109.155)
            
def plotTauEff12fb(aVar, aPar):
    pt = aVar[0]
    eta = aPar[0]
    return eff2012IsoTau12fb(pt, eta)
            
def effDoubleTau1Prong2012ABC(pt, eta):
    return (808.411*(0.581803*0.5*(ROOT.TMath.Erf((pt-41.1051)/2./0.661146/math.sqrt(pt))+1.))+ \
            4428.0*(0.619611*0.5*(ROOT.TMath.Erf((pt-41.7388)/2./0.732977/math.sqrt(pt))+1.))+ \
            1783.003*(0.683073*0.5*(ROOT.TMath.Erf((pt-42.3529)/2./0.797681/math.sqrt(pt))+1.))+ \
            5109.155*(0.627328*0.5*(ROOT.TMath.Erf((pt-41.285)/2./0.694383/math.sqrt(pt))+1.)) \
            )/(808.411+4428.0+1783.003+5109.155)

def plotDoubleTau1ProngEff(aVar, aPar):
    pt = aVar[0]
    eta = aPar[0]
    return effDoubleTau1Prong2012ABC(pt, eta)

def effDoubleTau1Prong2012ABC1OfflProng(pt, eta):
    return (808.411*(0.73201*0.5*(ROOT.TMath.Erf((pt-41.0386)/2./0.663772/math.sqrt(pt))+1.))+ \
            4428.0*(0.764371*0.5*(ROOT.TMath.Erf((pt-41.2967)/2./0.71147/math.sqrt(pt))+1.))+ \
            1783.003*(0.840278*0.5*(ROOT.TMath.Erf((pt-41.7228)/2./0.754745/math.sqrt(pt))+1.))+ \
            5109.155*(0.7879*0.5*(ROOT.TMath.Erf((pt-40.9494)/2./0.66706/math.sqrt(pt))+1.)) \
            )/(808.411+4428.0+1783.003+5109.155)

def plotDoubleTau1ProngEff1OfflProng(aVar, aPar):
    pt = aVar[0]
    eta = aPar[0]
    return effDoubleTau1Prong2012ABC1OfflProng(pt, eta)

def eff2012IsoTau19fb(pt, eta):
    return (808.411*(0.764166*0.5*(ROOT.TMath.Erf((pt-33.2236)/2./0.97289/math.sqrt(pt))+1.))+\
            4428.0*(0.802387*0.5*(ROOT.TMath.Erf((pt-38.0971)/2./0.82842/math.sqrt(pt))+1.))+\
            1783.003*(0.818051*0.5*(ROOT.TMath.Erf((pt-37.3669)/2./0.74847/math.sqrt(pt))+1.))+\
            5109.155*(0.796086*0.5*(ROOT.TMath.Erf((pt-37.3302)/2./0.757558/math.sqrt(pt))+1.))+\
            4131*(0.828182*0.5*(ROOT.TMath.Erf((pt-37.6596)/2./0.830682/math.sqrt(pt))+1.))+\
            3143*(0.833004*0.5*(ROOT.TMath.Erf((pt-37.634)/2./0.777843/math.sqrt(pt))+1.))\
            )/(808.411+4428.0+1783.003+5109.155+4131+3143)

def plotTauEff19fb(aVar, aPar):
    pt = aVar[0]
    eta = aPar[0]
    return eff2012IsoTau19fb(pt, eta)

def eff2012IsoTau1prong19fb1OfflProng(pt, eta):
    return (808.411*(0.73201*0.5*(ROOT.TMath.Erf((pt-41.0386)/2./0.663772/math.sqrt(pt))+1.))+\
            4428.0*(0.764371*0.5*(ROOT.TMath.Erf((pt-41.2967)/2./0.71147/math.sqrt(pt))+1.))+\
            1783.003*(0.840278*0.5*(ROOT.TMath.Erf((pt-41.7228)/2./0.754745/math.sqrt(pt))+1.))+\
            5109.155*(0.7879*0.5*(ROOT.TMath.Erf((pt-40.9494)/2./0.66706/math.sqrt(pt))+1.))+\
            4131*(0.811053*0.5*(ROOT.TMath.Erf((pt-41.2314)/2./0.72215/math.sqrt(pt))+1.))+\
            3143*(0.802065*0.5*(ROOT.TMath.Erf((pt-41.0161)/2./0.654632/math.sqrt(pt))+1.))\
            )/(808.411+4428.0+1783.003+5109.155+4131+3143)
            
def plotTauEff2012IsoTau1prong19fb1OfflProng(aVar, aPar):
    pt = aVar[0]
    eta = aPar[0]
    return eff2012IsoTau1prong19fb1OfflProng(pt, eta)


def eff2012IsoTauD(pt, eta):
    return (4131*(0.828182*0.5*(ROOT.TMath.Erf((pt-37.6596)/2./0.830682/math.sqrt(pt))+1.))+\
            3143*(0.833004*0.5*(ROOT.TMath.Erf((pt-37.634)/2./0.777843/math.sqrt(pt))+1.))\
            )/(4131+3143)

def plotEff2012IsoTauD(aVar, aPar):
    pt = aVar[0]
    eta = aPar[0]
    return eff2012IsoTauD(pt, eta)



# jet turn ons
def eff2012Jet19fbABCD(pt, eta):
    if (abs(eta)<=2.1):
        return (808.411*(0.99212*0.5*(ROOT.TMath.Erf((pt-31.3706)/2./1.22821/math.sqrt(pt))+1.))+\
                4428.0*(0.99059*0.5*(ROOT.TMath.Erf((pt-32.1104)/2./1.23292/math.sqrt(pt))+1.))+\
                1783.003*(0.988256*0.5*(ROOT.TMath.Erf((pt-31.3103)/2./1.18766/math.sqrt(pt))+1.))+\
                5109.155*(0.988578*0.5*(ROOT.TMath.Erf((pt-31.6391)/2./1.22826/math.sqrt(pt))+1.))+\
                4131*(0.989049*0.5*(ROOT.TMath.Erf((pt-31.9836)/2./1.23871/math.sqrt(pt))+1.))+\
                3143*(0.988047*0.5*(ROOT.TMath.Erf((pt-31.6975)/2./1.25372/math.sqrt(pt))+1.))\
                )/(808.411+4428.0+1783.003+5109.155+4131+3143)
    elif (abs(eta)>2.1 and abs(eta)<=3.0):
        return (808.411*(0.969591*0.5*(ROOT.TMath.Erf((pt-36.8179)/2./0.904254/math.sqrt(pt))+1.))+\
                4428.0*(0.975932*0.5*(ROOT.TMath.Erf((pt-37.2121)/2./0.961693/math.sqrt(pt))+1.))+\
                1783.003*(0.990305*0.5*(ROOT.TMath.Erf((pt-36.3096)/2./0.979524/math.sqrt(pt))+1.))+\
                5109.155*(0.971612*0.5*(ROOT.TMath.Erf((pt-36.2294)/2./0.871726/math.sqrt(pt))+1.))+\
                4131*(0.977958*0.5*(ROOT.TMath.Erf((pt-37.131)/2./0.987523/math.sqrt(pt))+1.))+\
                3143*(0.968457*0.5*(ROOT.TMath.Erf((pt-36.3159)/2./0.895031/math.sqrt(pt))+1.))\
                )/(808.411+4428.0+1783.003+5109.155+4131+3143)
    else: return 0

def plotEff2012Jet19fbABCD(aVar, aPar):
    pt = aVar[0]
    eta = aPar[0]
    return eff2012Jet19fbABCD(pt, eta)

def eff2012JetABC(pt, eta):
    if (abs(eta)<=2.1):
        return (808.411*(0.99212*0.5*(ROOT.TMath.Erf((pt-31.3706)/2./1.22821/math.sqrt(pt))+1.))+\
                4428.0*(0.99059*0.5*(ROOT.TMath.Erf((pt-32.1104)/2./1.23292/math.sqrt(pt))+1.))+\
                1783.003*(0.988256*0.5*(ROOT.TMath.Erf((pt-31.3103)/2./1.18766/math.sqrt(pt))+1.))+\
                5109.155*(0.988578*0.5*(ROOT.TMath.Erf((pt-31.6391)/2./1.22826/math.sqrt(pt))+1.))\
                )/(808.411+4428.0+1783.003+5109.155)
    elif (abs(eta)>2.1 and abs(eta)<=3.0):
        return (808.411*(0.969591*0.5*(ROOT.TMath.Erf((pt-36.8179)/2./0.904254/math.sqrt(pt))+1.))+\
                4428.0*(0.975932*0.5*(ROOT.TMath.Erf((pt-37.2121)/2./0.961693/math.sqrt(pt))+1.))+\
                1783.003*(0.990305*0.5*(ROOT.TMath.Erf((pt-36.3096)/2./0.979524/math.sqrt(pt))+1.))+\
                5109.155*(0.971612*0.5*(ROOT.TMath.Erf((pt-36.2294)/2./0.871726/math.sqrt(pt))+1.))\
                )/(808.411+4428.0+1783.003+5109.155)
    else: return 0
    
def plotEff2012JetABC(aVar, aPar):
    pt = aVar[0]
    eta = aPar[0]
    return eff2012JetABC(pt, eta)

def eff2012JetD(pt, eta):
    if (abs(eta)<=2.1):
        return (4131*(0.989049*0.5*(ROOT.TMath.Erf((pt-31.9836)/2./1.23871/math.sqrt(pt))+1.))+\
                3143*(0.988047*0.5*(ROOT.TMath.Erf((pt-31.6975)/2./1.25372/math.sqrt(pt))+1.))\
                )/(4131+3143)
    elif (abs(eta)>2.1 and abs(eta)<=3.0):
        return (4131*(0.977958*0.5*(ROOT.TMath.Erf((pt-37.131)/2./0.987523/math.sqrt(pt))+1.))+\
                3143*(0.968457*0.5*(ROOT.TMath.Erf((pt-36.3159)/2./0.895031/math.sqrt(pt))+1.))\
                )/(4131+3143)
    else: return 0
    
def plotEff2012JetD(aVar, aPar):
    pt = aVar[0]
    eta = aPar[0]
    return eff2012JetD(pt, eta)

# general functions

def errorfunc(pt, norm, effTurnOn, width):
    #return norm * 0.5 * (ROOT.TMath.Erf((x - effTurnOn) / width / math.sqrt(2)) + 1.)
    return norm*0.5*(ROOT.TMath.Erf((pt-effTurnOn)/2./width/math.sqrt(pt))+1.)

def errorfit(ax,apar):
    return errorfunc(ax[0],apar[0],apar[1],apar[2])

def crystalballfunc(m, m0, sigma, alpha, n, norm):
    sqrtPiOver2 = 1.2533141373
    sqrt2 = 1.4142135624
    sig = abs(sigma)
    t = (m - m0) / sig
    if (alpha < 0):
        t = -t
    absAlpha = abs(alpha / sig)
    a = ROOT.TMath.Power(n / absAlpha, n) * math.exp(-0.5 * absAlpha * absAlpha)
    b = absAlpha - n / absAlpha
    arg = absAlpha / sqrt2
    if (arg > 5.):
        ApproxErf = 1
    elif (arg < -5.):
        ApproxErf = -1
    else:
        ApproxErf = ROOT.TMath.Erf(arg)
    leftArea = (1 + ApproxErf) * sqrtPiOver2
    rightArea = (a * 1 / ROOT.TMath.Power(absAlpha - b, n - 1)) / (n - 1)
    area = leftArea + rightArea
    if (t <= absAlpha):
        arg = t / sqrt2
        if (arg > 5.):
            ApproxErf = 1
        elif (arg < -5.):
            ApproxErf = -1
        else:
            ApproxErf = ROOT.TMath.Erf(arg)
        return norm * (1 + ApproxErf) * sqrtPiOver2 / area
    else:
        return norm * (leftArea + a * (1 / ROOT.TMath.Power(t - b, n - 1) - 1 / ROOT.TMath.Power(absAlpha - b, n - 1)) / (1 - n)) / area


def crystalballfit( ax, apar):
    return crystalballfunc(ax[0],apar[0],apar[1],apar[2],apar[3],apar[4])


def combineEffHistos(folder,lHistos,trigLevel):
    totalLumi = 0
    for histo in lHistos:
        totalLumi += histo[1]
        
    print "Total lumi is ", totalLumi, " /pb"
    
    #xbin = [5,10,15,17,18,19,20,25,30,35,40,45,50,55,60,120,200]
    if (trigLevel == "Jet"):
        xbin = [5,10,15,17,18,19,20,25,30,35,40,45,50,55,60,70,80,90,100,120,140,160,200]
        etaBin = "CentralForward"
    else:
        xbin = [5,10,15,17,18,19,20,22,24,26,28,30,32,34,36,38,40,42.5,45,47.5,50,52.5,55,60,70,80,100,120,200]
        etaBin = "Total"
    
    xbins = array.array('d',xbin)
    denominator = ROOT.TH1F("denominator","",len(xbin)-1,xbins)
    numerator = ROOT.TH1F("numerator","",len(xbin)-1,xbins)
        
    for histo in lHistos:
        file = ROOT.TFile(folder+histo[0],"READ")
        if (trigLevel == "Jet"):
            denomName = "denominatorJet_CentralForward"
        else:
            denomName = "denominator_Total"
        denom = ROOT.TH1F(file.Get(denomName))
        #denom.Sumw2()
        scale = (histo[1]/totalLumi)/denom.Integral()
        denom.Scale(scale)
        denominator.Add(denom)
        
        num = ROOT.TH1F(file.Get("numerator"+trigLevel+"_"+etaBin))
        #num.Sumw2()
        num.Scale(scale)
        numerator.Add(num)
        
    eff = ROOT.TGraphAsymmErrors(numerator,denominator,"cl=0.683 b(1,1) mode")
    return eff#, numerator, denominator

# global variables:
#efficiencyGraph = -1
#startFitFrom = 0

class unbinnedFit():
    
    def __init__(self,effGraph,startFitFrom):
        self.efficiencyGraph = effGraph
        self.startFitFrom = startFitFrom
    
    # fcn to be used by MINUIT
    def fcnUnbinned(self, npar, gin, f, par, iflag):
        logL = 0.
        
        # Read x and y values from graph
#        if efficiencyGraph == -1:
#            print "ERROR: efficiencyGraph not defined. Abort."
#            return
        x = self.efficiencyGraph.GetX()
        y = self.efficiencyGraph.GetY()
        # loop over entries in TGraph
        for i in range(self.efficiencyGraph.GetN()):
            # only use points above set threshold
            if x[i] < self.startFitFrom: continue
            # compute logL for given point
            if y[i] > 0.99: # event passed
                pEff = errorfunc(x[i], par[0], par[1], par[2])
                pEff = pEff if pEff > 0 else 0.00000001 
                logL = logL + math.log(pEff)
            elif y[i] < 0.01: # event failed
                pEff = 1 - errorfunc(x[i], par[0], par[1], par[2])
                pEff = pEff if pEff > 0 else 0.00000001
                logL = logL + math.log(pEff)
            else:
                print "Unexpected y value. Skip event..."
                continue
        # return the log likelihood by reference
        f[0] = -1 * logL
    
    def fitUnbinned(self):
        minuit = ROOT.TMinuit(5)    # initialize TMinuit with a maximum of 5 params
        minuit.SetFCN(self.fcnUnbinned)
        minuit.SetPrintLevel(-1)
        arglist = array( 'd', 10*[0.])
        arglist[0] = 1
        ierflg = ROOT.Long(1982)
        
        minuit.mnexcm("SET ERR", arglist, 1, ierflg)
        
        # Set starting values and step sizes for parameters
        vstart = array( 'd',[0.9, 30., 1.])
        step = array( 'd',[0.01, 0.1, 0.01])
        minuit.mnparm(0, "a1", vstart[0], step[0], 0.,1.0,ierflg)
        minuit.mnparm(1, "a2", vstart[1], step[1], 0.,0.,ierflg)
        minuit.mnparm(2, "a3", vstart[2], step[2], 0.,0.,ierflg)
        
        # Now ready for minimization step
        arglist[0] = 1000.
        arglist[1] = 10.
        
        minuit.mnexcm("MIGRAD", arglist, 2, ierflg)
        minuit.mnexcm("MINOS", arglist, 1, ierflg)
        x0  = ROOT.Double(0)
        dx0 = ROOT.Double(0)
        minuit.GetParameter(0,x0,dx0)
        y0  = ROOT.Double(0)
        dy0 = ROOT.Double(0)
        minuit.GetParameter(1,y0,dy0)
        z0  = ROOT.Double(0)
        dz0 = ROOT.Double(0)
        minuit.GetParameter(2,z0,dz0)
        
        # Print results
        amin = ROOT.Double(0.18)
        edm = ROOT.Double(0.19)
        errdef = ROOT.Double(0.20)
        nvpar = ROOT.Long(1983)
        nparx = ROOT.Long(1984)
        icstat = ROOT.Long(1985)
        minuit.mnstat( amin, edm, errdef, nvpar, nparx, icstat )
        minuit.mnprin( 3, amin )
        
#        func = ROOT.TF1("func",errorfunc,self.startFitFrom,1000.,3)
#        func.SetParameters(x0,y0,z0)
#        func.SetParNames("constant","coefficient","width")
#        func.SetParError(0,dx0)
#        func.SetParError(1,dy0)
#        func.SetParError(2,dz0)
        
        return [x0,y0,z0, dx0,dy0,dz0]