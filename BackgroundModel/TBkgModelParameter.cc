#include "TBkgModelParameter.hh"

ClassImp(TBkgModelParameter)

TBkgModelParameter::TBkgModelParameter()
{
	TBackgroundModel *f1 = new TBackgroundModel(500, 8000, 50, 0, false);
}

TBkgModelParameter::TBkgModelParameter(int fNParam)
{
	dNParam = fNParam
	
}


TBkgModelParameter::~TBkgModelParameter()
{}

