Notes for the background model code (sometime in April, 2016)

### General
The fitter uses two classes:

TBkgModelParameter -- container for background model parameters. Primarily used to keep track of the histograms of each background model parameter.

TBackgroundModel -- Class with the actual fitter, fit is performed using TMinuit. The header file describes most of the member functions in slight detail. 


### Required data/MC files
If you need the MC/data files, email me.

Model PDFs -- model PDFs are saved as normalized histograms (with 1 keV binning) in 4 files, the "Initialize" function loads all of the histograms from the files. Only Bulk gamma sources were generated using the double Gaussian detector response. 

	1) MCProduction_Bulk_1keV.root -- Bulk histograms, using single Gaussian smearing

	2) MCProduction_BulkCDR_1keV.root -- Bulk histograms, using double Gaussian (Calibration Detector Response)

	3) MCProduction_BulkCDRInternal_1keV.root -- Bulk histograms for internal shields using double Gaussian, this should have been combined with the previous file but I messed up producing the histograms and I didn't bother to remake the other file

	4) MCProduction_Surface_1keV.root -- Surface histograms, using single Gaussian smearing

Background Data -- I have taken the Reduced Unblinded ROOT files and reduced them further, storing only Energy, Multiplicity, and the PSA Cuts. The function "LoadData" should be edited if different data files/formats are used.

