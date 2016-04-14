Notes for the background model code (sometime in April, 2016)

### General
The fitter uses two classes:
TBkgModelParameter -- container for background model parameters. Primarily used to keep track of the histograms of each background model parameter, I loop through an array of these.

TBackgroundModel -- Class with the actual fitter, fit is performed using TMinuit.

### Running the code
I usually compile the two classes and then load the objects in a ROOT macro to run

eg:
	gSystem->Load("TBkgModelParameter_cc.so");
	gSystem->Load("TBackgroundModel_cc.so");

	TBackgroundModel *f1 = new TBackgroundModel(500, 7000, 50, 0, false);
	f1->DoTheFit();

