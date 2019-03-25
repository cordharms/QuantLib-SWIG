
/*
 Copyright (C) 2019, Cord Harms

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

#ifndef quantlib_experimental_i
#define quantlib_experimental_i

%include templatemontecarlo.i

//RealStochasticProcesses:

%{
#include <ql/experimental/templatemodels/multiasset/multiassetbsmodel.hpp>

using QuantLib::MultiAssetBSModel;
typedef boost::shared_ptr<RealStochasticProcess> MultiAssetBSModelPtr;
%}

%rename(MultiAssetBSModel) MultiAssetBSModelPtr;
class MultiAssetBSModelPtr : public boost::shared_ptr<RealStochasticProcess> {
  public:
    %extend {
        MultiAssetBSModelPtr(const Handle<YieldTermStructure>&                                            termStructure,
			              const std::vector<std::string>&                                                 aliases,
			              const std::vector<GeneralizedBlackScholesProcessPtr>&                           processes,
			              const std::vector< std::vector<Real> >&                                         correlations) {
            std::vector<boost::shared_ptr<GeneralizedBlackScholesProcess>> gbslist(processes.size());
            for(size_t i=0;i<gbslist.size();i++){
                gbslist[i] = boost::dynamic_pointer_cast<GeneralizedBlackScholesProcess>(processes[i]);
            }
            return new MultiAssetBSModelPtr(new MultiAssetBSModel(termStructure,aliases,gbslist,correlations));
        };
    }
};

%{
#include <ql/experimental/termstructures/localCorrFX/localcorrsurfaceabffx.hpp>
#include <ql/experimental/termstructures/localcorrtermstructure.hpp>
#include <ql/experimental/termstructures/localCorrFX/CTSlocalInCrossCorrelationFX.hpp>
#include <ql/experimental/termstructures/localCorrFX/CTSlocalInCrossCovarianceFX.hpp>
#include <ql/experimental/termstructures/localCorrFX/CTSlocalInCrossNegSkewFX.hpp>
#include <ql/experimental/termstructures/localCorrFX/CTSlocalInCrossVolatilityFX.hpp>

using QuantLib::LocalCorrSurfaceABFFX;
using QuantLib::LocalCorrTermStructure;
using QuantLib::CTSlocalInCrossCorrelationFX;
using QuantLib::CTSlocalInCrossCovarianceFX;
using QuantLib::CTSlocalInCrossNegSkewFX;
using QuantLib::CTSlocalInCrossVolatilityFX;
%}

%ignore LocalCorrTermStructure;
class LocalCorrTermStructure : public Extrapolator {};
%ignore LocalCorrSurfaceABFFX;
class LocalCorrSurfaceABFFX : public LocalCorrTermStructure {};

%template(LocalCorrTermStructureBase) boost::shared_ptr<LocalCorrTermStructure>;
%template(LocalCorrSurfaceABFFXBase) boost::shared_ptr<LocalCorrSurfaceABFFX>;

%{
typedef boost::shared_ptr<LocalCorrSurfaceABFFX> CTSlocalInCrossCorrelationFXPtr;
typedef boost::shared_ptr<LocalCorrSurfaceABFFX> CTSlocalInCrossCovarianceFXPtr;
typedef boost::shared_ptr<LocalCorrSurfaceABFFX> CTSlocalInCrossNegSkewFXPtr;
typedef boost::shared_ptr<LocalCorrSurfaceABFFX> CTSlocalInCrossVolatilityFXPtr;
typedef boost::shared_ptr<LocalCorrSurfaceABFFX> LocalCorrSurfaceABFFXPtr;
%}

class LocalCorrSurfaceABFFXPtr : public boost::shared_ptr<LocalCorrSurfaceABFFX>, public boost::shared_ptr<LocalCorrTermStructure> {
  public:
    %extend {
        Matrix getLocalCorrelationSurface(Time t, std::vector<Real> assetGrid1, std::vector<Real> assetGrid2){
           return boost::dynamic_pointer_cast<LocalCorrSurfaceABFFX>(*self)->getLocalCorrelationSurface(t,assetGrid1,assetGrid2);
        };
    }
};

%template(LocalCorrSurfaceABFFXHandle) Handle<LocalCorrSurfaceABFFX>;
%template(LocalCorrTermStructureHandle) Handle<LocalCorrTermStructure>;

%rename(CTSlocalInCrossCorrelationFX) CTSlocalInCrossCorrelationFXPtr;
class CTSlocalInCrossCorrelationFXPtr : public LocalCorrSurfaceABFFXPtr {
  public:
    %extend {
        CTSlocalInCrossCorrelationFXPtr(const std::vector<GeneralizedBlackScholesProcessPtr>& processes,
								  const GeneralizedBlackScholesProcessPtr&					processToCal) {
            std::vector<boost::shared_ptr<GeneralizedBlackScholesProcess>> gbslist(processes.size());
            for(size_t i=0;i<gbslist.size();i++){
                gbslist[i] = boost::dynamic_pointer_cast<GeneralizedBlackScholesProcess>(processes[i]);
            }
            boost::shared_ptr<GeneralizedBlackScholesProcess> pToCal = boost::dynamic_pointer_cast<GeneralizedBlackScholesProcess>(processToCal);
            return new CTSlocalInCrossCorrelationFXPtr(new CTSlocalInCrossCorrelationFX(gbslist,pToCal));
        };
    }
};

%rename(CTSlocalInCrossCovarianceFX) CTSlocalInCrossCovarianceFXPtr;
class CTSlocalInCrossCovarianceFXPtr : public LocalCorrSurfaceABFFXPtr {
  public:
    %extend {
        CTSlocalInCrossCovarianceFXPtr(const std::vector<GeneralizedBlackScholesProcessPtr>& processes,
								  const GeneralizedBlackScholesProcessPtr&					processToCal) {
            std::vector<boost::shared_ptr<GeneralizedBlackScholesProcess>> gbslist(processes.size());
            for(size_t i=0;i<gbslist.size();i++){
                gbslist[i] = boost::dynamic_pointer_cast<GeneralizedBlackScholesProcess>(processes[i]);
            }
            boost::shared_ptr<GeneralizedBlackScholesProcess> pToCal = boost::dynamic_pointer_cast<GeneralizedBlackScholesProcess>(processToCal);
            return new CTSlocalInCrossCovarianceFXPtr(new CTSlocalInCrossCovarianceFX(gbslist,pToCal));
        };
    }
};

%rename(CTSlocalInCrossNegSkewFX) CTSlocalInCrossNegSkewFXPtr;
class CTSlocalInCrossNegSkewFXPtr : public LocalCorrSurfaceABFFXPtr {
  public:
    %extend {
        CTSlocalInCrossNegSkewFXPtr(const std::vector<GeneralizedBlackScholesProcessPtr>& processes,
								  const GeneralizedBlackScholesProcessPtr&					processToCal,
									double beta) {
            std::vector<boost::shared_ptr<GeneralizedBlackScholesProcess>> gbslist(processes.size());
            for(size_t i=0;i<gbslist.size();i++){
                gbslist[i] = boost::dynamic_pointer_cast<GeneralizedBlackScholesProcess>(processes[i]);
            }
            boost::shared_ptr<GeneralizedBlackScholesProcess> pToCal = boost::dynamic_pointer_cast<GeneralizedBlackScholesProcess>(processToCal);
            return new CTSlocalInCrossNegSkewFXPtr(new CTSlocalInCrossNegSkewFX(gbslist,pToCal,beta));
        };
    }
};

%rename(CTSlocalInCrossVolatilityFX) CTSlocalInCrossVolatilityFXPtr;
class CTSlocalInCrossVolatilityFXPtr : public LocalCorrSurfaceABFFXPtr {
  public:
    %extend {
        CTSlocalInCrossVolatilityFXPtr(const std::vector<GeneralizedBlackScholesProcessPtr>& processes,
								  const GeneralizedBlackScholesProcessPtr&					processToCal) {
            std::vector<boost::shared_ptr<GeneralizedBlackScholesProcess>> gbslist(processes.size());
            for(size_t i=0;i<gbslist.size();i++){
                gbslist[i] = boost::dynamic_pointer_cast<GeneralizedBlackScholesProcess>(processes[i]);
            }
            boost::shared_ptr<GeneralizedBlackScholesProcess> pToCal = boost::dynamic_pointer_cast<GeneralizedBlackScholesProcess>(processToCal);
            return new CTSlocalInCrossVolatilityFXPtr(new CTSlocalInCrossVolatilityFX(gbslist,pToCal));
        };
    }
};


%{
#include <ql/experimental/templatemodels/multiasset/localcorrelationbsmodel.hpp>

using QuantLib::LocalCorrelationBSModel;
using QuantLib::LocalCorrTermStructure;
typedef boost::shared_ptr<RealStochasticProcess> LocalCorrelationBSModelPtr;
%}

%rename(LocalCorrelationBSModel) LocalCorrelationBSModelPtr;
class LocalCorrelationBSModelPtr : public boost::shared_ptr<RealStochasticProcess> {
  public:
    %extend {
        LocalCorrelationBSModelPtr(const Handle<YieldTermStructure>&                                      termStructure,
			              const std::vector<std::string>&                                                 aliases,
			  			  const std::vector<GeneralizedBlackScholesProcessPtr>&                           processes,
						  const Handle<LocalCorrTermStructure>&											  localCorrTermStructure) {
            std::vector<boost::shared_ptr<GeneralizedBlackScholesProcess>> gbslist(processes.size());
            for(size_t i=0;i<gbslist.size();i++){
                gbslist[i] = boost::dynamic_pointer_cast<GeneralizedBlackScholesProcess>(processes[i]);
            }
            return new LocalCorrelationBSModelPtr(new LocalCorrelationBSModel(termStructure,aliases,gbslist,localCorrTermStructure));
        };
    }
};

%{
#include <ql/experimental/termstructures/Helper/ParticleMethodUtils.hpp>
#include <ql/experimental/termstructures/localCorrFX/localcorrsurfaceabfFX.hpp>

using QuantLib::blackFormulaImpliedStdDev;
using QuantLib::ParticleMethodUtils;
%}

    
%inline %{
    Real  get_blackFormulaImpliedStdDev(Option::Type optionType,
                                                Real strike,
                                                Real forward,
                                                Real blackPrice,
                                                Real discount,
                                                Real displacement=0.0,
                                                Real guess=0.1, Real accuracy=1.0e-6, Natural maxIter=100) {
        return blackFormulaImpliedStdDev(optionType,strike,forward,blackPrice,discount,displacement,guess,accuracy,maxIter);
    };
    void calibrateFXlocalCTSWithParticleMethod(Handle<LocalCorrSurfaceABFFX>& surface,const std::string& kernelIn, unsigned int numberOfPaths, Time maxTime,
			  Time deltaT, Time tMin, Real kappa, Real sigmaAVR, Real exponentN, Real gridMinQuantile,
			  Real gridMaxQuantile, unsigned int ns1, unsigned int ns2){
        ParticleMethodUtils::calibrateFX(surface,kernelIn,numberOfPaths,maxTime,deltaT,tMin,kappa,sigmaAVR,exponentN,gridMinQuantile,gridMaxQuantile,ns1,ns2);
    };
%}

#endif
