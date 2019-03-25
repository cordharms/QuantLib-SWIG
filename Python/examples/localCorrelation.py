
# Copyright (C) 2019, Cord Harms
#
# This file is part of QuantLib, a free-software/open-source library
# for financial quantitative analysts and developers - http://quantlib.org/
#
# QuantLib is free software: you can redistribute it and/or modify it under the
# terms of the QuantLib license.  You should have received a copy of the
# license along with this program; if not, please email
# <quantlib-dev@lists.sf.net>. The license is also available online at
# <http://quantlib.org/license.shtml>.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the license for more details.

from QuantLib import *
import numpy as np
import math
import matplotlib.pyplot as plt

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

class VolSurface:
    def __init__(self,todaysDate, spot, deltas, matAsDates, yieldCurve, divCurve):
        self.todaysDate = todaysDate
        self.spot = spot
        self.deltas = deltas
        self.matAsDates = matAsDates
        self.yieldCurve = yieldCurve
        self.divCurve = divCurve
        self.quotes = []
        self.strikes = []
        self.sabrStart = []
        self.smiledSurface = None
    def __calcStrikesFromDeltas(self):
        self.strikes = np.zeros((len(self.matAsDates),len(self.deltas)))
        for i in range(0,len(self.deltas)):
            for j in range(0,len(self.matAsDates)):
                if is_number(self.deltas[i]):
                    if self.deltas[i]<0:
                        bdc = BlackDeltaCalculator(Option.Put,
                            DeltaVolQuote.Spot,
                            self.spot.value(),
                            self.yieldCurve.discount(self.matAsDates[j]),
                            self.divCurve.discount(self.matAsDates[j]),
                            self.quotes[j][i]*math.sqrt(Actual365Fixed().yearFraction(self.todaysDate, self.matAsDates[j]))*0.01)
                        self.strikes[j][i] = bdc.strikeFromDelta(self.deltas[i])
                    else:
                        bdc = BlackDeltaCalculator(Option.Call,
                            DeltaVolQuote.Spot,
                            self.spot.value(),
                            self.yieldCurve.discount(self.matAsDates[j]),
                            self.divCurve.discount(self.matAsDates[j]),
                            self.quotes[j][i]*math.sqrt(Actual365Fixed().yearFraction(self.todaysDate, self.matAsDates[j]))*0.01)
                        self.strikes[j][i] = bdc.strikeFromDelta(self.deltas[i])
                else:
                    bdc = BlackDeltaCalculator(Option.Put,
                        DeltaVolQuote.Spot,
                        self.spot.value(),
                        self.yieldCurve.discount(self.matAsDates[j]),
                        self.divCurve.discount(self.matAsDates[j]),
                        self.quotes[j][i]*math.sqrt(Actual365Fixed().yearFraction(self.todaysDate, self.matAsDates[j]))*0.01)
                    self.strikes[j][i] = bdc.atmStrike(DeltaVolQuote.AtmFwd)
                    self.atmIndex = i
    def calcSABRSmiledSurface(self):
        self.__calcStrikesFromDeltas()
        self.sabrSections = VectorOfSABRInterpolatedSmileSection()
        for i in range(0,min(len(self.matAsDates),13)):
            sabr = SabrInterpolatedSmileSection(
                           self.matAsDates[i], self.spot.value(),self.strikes[i,:].tolist(),
                           False, self.quotes[i][self.atmIndex]*0.01, (self.quotes[i,:]*0.01).tolist(),
                           self.sabrStart[i][0]*1, self.sabrStart[i][1]*1, self.sabrStart[i][3]*1, 
                           self.sabrStart[i][2]*1,False,True,False,False)
            self.sabrSections.append(sabr)
            #print(sabr.alpha(), sabr.beta(), sabr.nu(), sabr.rho())
        self.smiledSurface = BlackVolTermStructureHandle(SmiledSurface(self.sabrSections,self.todaysDate,calendar, Following,Actual365Fixed()))
        self.interpolatedSurface = LocalVolTermStructureHandle(InterpolatedLocalVolSurface(self.smiledSurface,
            self.yieldCurve, self.divCurve, self.spot,100,100))
        
        
input("Press Enter to continue...")        

########################################################
# global data
########################################################

calendar = TARGET()
todaysDate = Date(10,October,2017);
Settings.instance().evaluationDate = todaysDate

useLocalCorr = False

########################################################
# market data
########################################################

underlying1 = QuoteHandle(SimpleQuote(1.0))
underlying2 = QuoteHandle(SimpleQuote(1.0))
underlyingCross = QuoteHandle(SimpleQuote(underlying1.value()/underlying2.value()))

yieldDom = YieldTermStructureHandle(FlatForward(0,calendar,0.0,Actual365Fixed()))
yield1 = YieldTermStructureHandle(FlatForward(0,calendar,0.0,Actual365Fixed()))
yield2 = YieldTermStructureHandle(FlatForward(0,calendar,0.0,Actual365Fixed()))

volDelta = [-0.05,-0.1,-0.15,-0.25,-0.35,"ATM",0.35,0.25,0.15,0.1,0.05]
volMat = ["1D","1W","2W","3W","1M","2M","3M","4M","5M","6M","9M","1Y","18M","2Y","3Y","4Y","5Y","6Y","7Y","10Y","15Y","20Y","25Y"]
#transform to dates
for i in range(0,len(volMat)): volMat[i] = calendar.advance(todaysDate,Period(volMat[i]))

vol1 = VolSurface(todaysDate, underlying1,volDelta, volMat,yieldDom,yield1)
vol2 = VolSurface(todaysDate, underlying2,volDelta, volMat,yieldDom,yield2)
volCross = VolSurface(todaysDate, underlyingCross,volDelta, volMat, yield2, yield1)

vol1.quotes = np.loadtxt("localCorrelation_vol1.txt", dtype='f', delimiter=',')
vol2.quotes = np.loadtxt("localCorrelation_vol2.txt", dtype='f', delimiter=',')
volCross.quotes = np.loadtxt("localCorrelation_volCross.txt", dtype='f', delimiter=',')
vol1.sabrStart = np.loadtxt("localCorrelation_sabrStartValues.txt", dtype='f', delimiter=',')
vol2.sabrStart = np.loadtxt("localCorrelation_sabrStartValues.txt", dtype='f', delimiter=',')
volCross.sabrStart = np.loadtxt("localCorrelation_sabrStartValues.txt", dtype='f', delimiter=',')

try:
    assert len(volDelta) == vol1.quotes.shape[1], "err: size of volDelta and columns of vol1 do not match"
    assert len(volDelta) == vol2.quotes.shape[1], "err: size of volDelta and columns of vol2 do not match"
    assert len(volDelta) == volCross.quotes.shape[1], "err: size of volDelta and columns of volCross do not match"
    assert len(volMat) == vol1.quotes.shape[0], "err: size of volMat and rows of vol1 do not match"
    assert len(volMat) == vol2.quotes.shape[0], "err: size of volMat and rows of vol2 do not match"
    assert len(volMat) == volCross.quotes.shape[0], "err: size of volMat and rows of volCross do not match"
    assert len(volMat) == vol1.sabrStart.shape[0], "err: size of volMat and rows of sabrStart do not match"
    assert len(volMat) == vol2.sabrStart.shape[0], "err: size of volMat and rows of sabrStart do not match"
    assert len(volMat) == volCross.sabrStart.shape[0], "err: size of volMat and rows of sabrStart do not match"
except AssertionError as e:
    print(e)
    quit()

########################################################
# derive blackVolSurface and stochastic processes
########################################################

vol1.calcSABRSmiledSurface()
vol2.calcSABRSmiledSurface()
volCross.calcSABRSmiledSurface()

gbsprocess1 = GeneralizedBlackScholesProcess(underlying1,yield1, yieldDom,vol1.smiledSurface,vol1.interpolatedSurface)
gbsprocess2 = GeneralizedBlackScholesProcess(underlying2,yield2, yieldDom,vol2.smiledSurface,vol2.interpolatedSurface)
gbsprocessCross = GeneralizedBlackScholesProcess(underlyingCross,yield2,yield1,volCross.smiledSurface,volCross.interpolatedSurface)

gbslist = GeneralizedBlackScholesProcessVector()
gbslist.append(gbsprocess1)
gbslist.append(gbsprocess2)

corr = DoubleVectorVector([[1,0.54],[0.54,1]])
aliases = ["EURUSD","GBPUSD"]
maprocess = None

if useLocalCorr:
    cts = CTSlocalInCrossCorrelationFX(gbslist,gbsprocessCross)
    calibrateFXlocalCTSWithParticleMethod(LocalCorrSurfaceABFFXHandle(cts),"QuarticKernel", 10000, 1,
			  1.0/80, 0.25, 3, 0.0848, -0.2, 0.01, 0.99, 50, 30)
    maprocess = LocalCorrelationBSModel(yieldDom,aliases, gbslist,LocalCorrTermStructureHandle(cts))
else:
    maprocess = MultiAssetBSModel(yieldDom,aliases, gbslist,corr)

schedule = Schedule(todaysDate,calendar.advance(todaysDate,Period("1y1m")),Period("4d"),calendar,Unadjusted,Unadjusted,DateGeneration.Forward,False)
times = DoubleVector([(schedule.dates()[i]-todaysDate)/365 for i in range(0,len(schedule.dates()))])

simulation = RealMCSimulation(maprocess,times,times,40000,1,True,True,False)
simulation.simulate();
simulation.calculateAssetAdjuster(times,aliases);

########################################################
# calculate asset smiles for a particular model
########################################################

#4:1M, 11:1Y
strikeRow = 11
frac = Actual365Fixed().yearFraction(todaysDate, vol1.matAsDates[strikeRow])

implVol = []
sabrVol = []
optionList = RealMCPayoffVector()   
    
for i in range(0,vol1.strikes.shape[1]):
    option = RealMCVanillaOption(frac,aliases[0],vol1.strikes[strikeRow][i],1)
    optionList.append(option)
price = RealMCPricer.NPVs(optionList,simulation)
for i in range(0,vol1.strikes.shape[1]):
    implVol.append(get_blackFormulaImpliedStdDev(Option.Call, vol1.strikes[strikeRow][i],1.0,price[i], 1.0,0.0)/math.sqrt(frac))
    sabrVol.append(vol1.smiledSurface.blackVol(vol1.matAsDates[strikeRow],vol1.strikes[strikeRow][i],True))
plt.subplot(2, 1, 1)
lineMC, = plt.plot(vol1.strikes[strikeRow],implVol,label='MC')
lineQ, = plt.plot(vol1.strikes[strikeRow],vol1.quotes[strikeRow]*0.01, label='quotes')
lineSABR, = plt.plot(vol1.strikes[strikeRow],sabrVol, label = 'sabr' )
plt.ylabel('vol in %')
plt.legend(handles=[lineMC, lineQ, lineSABR])


########################################################
# calculate cross smile for a particular model
########################################################

implVol = []
sabrVol = []
crossOptions = RealMCPayoffVector() 

asset = RealMCPayoffVector()
asset.append(RealMCAsset(frac,aliases[0]))
asset.append(RealMCAsset(frac,aliases[1])) 

for i in range(0,volCross.strikes.shape[1]):
    payoff = RealMCAxpy(-volCross.strikes[strikeRow][i],asset[1],asset[0])
    option =  RealMCPay(RealMCMax(payoff,RealMCFixedAmount(0)),frac)
    crossOptions.append(option)
price = RealMCPricer.NPVs(crossOptions,simulation)
for i in range(0,volCross.strikes.shape[1]):
    implVol.append(get_blackFormulaImpliedStdDev(Option.Call, volCross.strikes[strikeRow][i],1.0,price[i], 1.0,0.0)/frac)
    sabrVol.append(volCross.smiledSurface.blackVol(volCross.matAsDates[strikeRow],volCross.strikes[strikeRow][i],True))
plt.subplot(2, 1, 2)
lineMC, = plt.plot(volCross.strikes[strikeRow],implVol,label='MC')
lineQ, = plt.plot(volCross.strikes[strikeRow],volCross.quotes[strikeRow]*0.01, label='quotes')
lineSABR, = plt.plot(volCross.strikes[strikeRow],sabrVol, label = 'sabr' )
plt.ylabel('vol in %')
plt.legend(handles=[lineMC, lineQ, lineSABR])
plt.show()