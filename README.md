Notes for the background model code (sometime in April, 2016)

### General
The fitter uses two classes:
TBkgModelParameter -- container for background model parameters. Primarily used to keep track of the histograms of each background model parameter.

TBackgroundModel -- Class with the actual fitter, fit is performed using TMinuit. If I were to re-write the code, I would probably try use TMinuit2 or RooFit.


### Running the code
I usually compile the two classes and then load the objects in a ROOT macro and run whatever functions I need

eg:
	gSystem->Load("TBkgModelParameter_cc.so");
	gSystem->Load("TBackgroundModel_cc.so");
	TBackgroundModel *f1 = new TBackgroundModel(500, 7000, 50, 0, false);
	f1->DoTheFit();
	f1->PrintParameters();

