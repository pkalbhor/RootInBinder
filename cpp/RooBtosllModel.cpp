// Custom model in ROOT for B to sll process
#include "Riostream.h" 

#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
#include <math.h> 
#include "TMath.h" 
 
#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
 
class RooBtosllModel : public RooAbsPdf {
public:
    RooBtosllModel() {} ; 
    RooBtosllModel(const char *name, const char *title,
    RooAbsReal& _CosThetaL,
    RooAbsReal& _CosThetaK,
    RooAbsReal& _unboundAfb,
    RooAbsReal& _unboundFl);
    RooBtosllModel(const RooBtosllModel& other, const char* name=0) ;
    virtual TObject* clone(const char* newname) const { return new RooBtosllModel(*this,newname); }
    inline virtual ~RooBtosllModel() { }

protected:
    RooRealProxy CosThetaL ;
    RooRealProxy CosThetaK ;
    RooRealProxy unboundAfb ;
    RooRealProxy unboundFl ;
    Double_t evaluate() const ;

private:
    ClassDef(RooBtosllModel,1) // Your description goes here...
};

RooBtosllModel::RooBtosllModel(const char *name, const char *title, 
                       RooAbsReal& _CosThetaL,
                       RooAbsReal& _CosThetaK,
                       RooAbsReal& _unboundAfb,
                       RooAbsReal& _unboundFl) :
  RooAbsPdf(name,title), 
  CosThetaL("CosThetaL","CosThetaL",this,_CosThetaL),
  CosThetaK("CosThetaK","CosThetaK",this,_CosThetaK),
  unboundAfb("unboundAfb","unboundAfb",this,_unboundAfb),
  unboundFl("unboundFl","unboundFl",this,_unboundFl)
{ 
} 


RooBtosllModel::RooBtosllModel(const RooBtosllModel& other, const char* name) :  
  RooAbsPdf(other,name), 
  CosThetaL("CosThetaL",this,other.CosThetaL),
  CosThetaK("CosThetaK",this,other.CosThetaK),
  unboundAfb("unboundAfb",this,other.unboundAfb),
  unboundFl("unboundFl",this,other.unboundFl)
{ 
} 


Double_t RooBtosllModel::evaluate() const 
{ 
  // ENTER EXPRESSION IN TERMS OF VARIABLE ARGUMENTS HERE 
  Double_t fl = 0.5+TMath::ATan(unboundFl)/TMath::Pi();
  Double_t afb = 2.0*(1.-fl)*TMath::ATan(unboundAfb)/TMath::Pi();
  Double_t result = (9.0/16.0)*((0.5*(1.0-fl)*(1.0-CosThetaK*CosThetaK)*(1.0+CosThetaL*CosThetaL)) + (2.0*fl*CosThetaK*CosThetaK*(1.0-CosThetaL*CosThetaL)) + (afb*(1.0-CosThetaK*CosThetaK)*CosThetaL)) ;
  return result;
}