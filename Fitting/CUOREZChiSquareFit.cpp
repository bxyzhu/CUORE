// Compile with eg: "g++ MySource.cc `root-config --libs --cflags` -o foo"
// CUORE-0 Montecarlo ChiSquared fit

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <stdio.h>

#include "TFile.h"
#include "TROOT.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"

#include "THStack.h"
#include "TMath.h"
#include "TH2F.h"
#include "TLegend.h"

#include "TRandom1.h"


#include <TMinuitMinimizer.h>

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#include <cmath>

using namespace std;

    // Probably want to create a function to load the file names here
    // void LoadData(const char* fileName)


    // Define chi square using simulation data
    Double_t FitChiSquare(const double* para)
    {
        // Load simulation data
        TFile      *Hist_FileName = new TFile(".root", "READ");
        
        // Define histograms for filling
        TH1F          *smoothed_hist_S1  = (TH1F*)Hist_FileName->Get("S1Spec_50_100");
        
        
        const Int_t         TotBinNbr = smoothed_hist_S1->GetSize();
        const Double_t      Max_xAxis = smoothed_hist_S1->GetXaxis()->GetXmax();
        const Double_t      Min_xAxis = smoothed_hist_S1->GetXaxis()->GetXmin();

        //std::cout<<"Max BIN: "<<Max_xAxis<<std::endl;
        
        
        Float_t    Data_yAxis_Content[TotBinNbr];
        Float_t  G4Data_yAxis_Content[TotBinNbr];
        
        
        TH1F      *hist_G4Data = new TH1F("hist_G4Data", "hist_G4Data", TotBinNbr, Min_xAxis, Max_xAxis);
        
        
     string rootFile = "SCENE_ND_G4_508keVNeutron_90deg_700-million_12052012_LAr-0-5keVTrigger_EJ301-100keVTrigger.root";
        
        
        //-------------------------------------------------------------------
        
        TFile     *f = new TFile(rootFile.c_str(), "READ");
        TTree *tData = (TTree*)f->Get("tData");

        
        Int_t            iInteractingVolume;
        
        vector<float>    *fDetector1_PeakStartTime      = new vector<float>;
        vector<float>    *fDetector1_PeakEndTime        = new vector<float>;
        vector<float>    *fDetector1_MaxPeakEnergy      = new vector<float>;
        vector<float>    *fDetector1_TotPeakEnergy      = new vector<float>;
        vector<float>    *fDetector1_MaxPeakTime        = new vector<float>;
        vector<float>    *fDetector1_MaxPeakPositionX   = new vector<float>;
        vector<float>    *fDetector1_MaxPeakPositionY   = new vector<float>;
        vector<float>    *fDetector1_MaxPeakPositionZ   = new vector<float>;
        vector<float>    *fDetector1_PeakDirectionX   = new vector<float>;
        vector<float>    *fDetector1_PeakDirectionY   = new vector<float>;
        vector<float>    *fDetector1_PeakDirectionZ   = new vector<float>;
        
        
        std::vector<std::string>    *sDetector1_PeakCreatorParticle = new std::vector<std::string>;
        std::vector<std::string>    *sDetector1_PeakCreatorProcess  = new std::vector<std::string>;
        
        Int_t           iDetector1_MultipleScattering;
        Int_t           iDetector1_PeakNbr;
        Float_t         fDetector1_TriggerEnergy;
        Float_t         fDetector1_TriggerTime;
        Float_t         fDetector1_TriggerPositionX;
        Float_t         fDetector1_TriggerPositionY;
        Float_t         fDetector1_TriggerPositionZ;
        Float_t         fDetector1_TotEneDeposit;
        char            sDetector1_TriggerParticle[30] = "";
        char            sDetector1_TriggerProcess[30] = "";
        char            sDetector1_TriggerParentParticle[30] = "";
        
        vector<float>    *fDetector2_PeakStartTime      = new vector<float>;  
        vector<float>    *fDetector2_PeakEndTime        = new vector<float>;
        vector<float>    *fDetector2_MaxPeakEnergy      = new vector<float>;
        vector<float>    *fDetector2_TotPeakEnergy      = new vector<float>;
        vector<float>    *fDetector2_MaxPeakTime        = new vector<float>;
        vector<float>    *fDetector2_MaxPeakPositionX   = new vector<float>;
        vector<float>    *fDetector2_MaxPeakPositionY   = new vector<float>;
        vector<float>    *fDetector2_MaxPeakPositionZ   = new vector<float>;
        
        std::vector<std::string>    *sDetector2_PeakCreatorParticle = new std::vector<std::string>;
        std::vector<std::string>    *sDetector2_PeakCreatorProcess  = new std::vector<std::string>;
        
        Int_t           iDetector2_PeakNbr;
        Float_t         fDetector2_TriggerEnergy;
        Float_t         fDetector2_TriggerTime;
        Float_t         fDetector2_TriggerPositionX;
        Float_t         fDetector2_TriggerPositionY;
        Float_t         fDetector2_TriggerPositionZ;
        Float_t         fDetector1_TriggerDirectionX;
        Float_t         fDetector1_TriggerDirectionY;
        Float_t         fDetector1_TriggerDirectionZ;
        Float_t         fDetector2_TotEneDeposit;
        char            sDetector2_TriggerParticle[30] = "";
        char            sDetector2_TriggerProcess[30] = "";
        
        
        tData->SetBranchAddress("iInteractingVolume",                   &iInteractingVolume);
        tData->SetBranchAddress("iDetector1_MultipleScattering",        &iDetector1_MultipleScattering);
        tData->SetBranchAddress("fDetector1_PeakStartTime",             &fDetector1_PeakStartTime);
        tData->SetBranchAddress("fDetector1_PeakEndTime",               &fDetector1_PeakEndTime);
        tData->SetBranchAddress("fDetector1_MaxPeakEnergy",             &fDetector1_MaxPeakEnergy);
        tData->SetBranchAddress("fDetector1_TotPeakEnergy",             &fDetector1_TotPeakEnergy);
        tData->SetBranchAddress("fDetector1_MaxPeakTime",               &fDetector1_MaxPeakTime) ;
        tData->SetBranchAddress("fDetector1_MaxPeakPositionX",          &fDetector1_MaxPeakPositionX);
        tData->SetBranchAddress("fDetector1_MaxPeakPositionY",          &fDetector1_MaxPeakPositionY);
        tData->SetBranchAddress("fDetector1_MaxPeakPositionZ",          &fDetector1_MaxPeakPositionZ);
        tData->SetBranchAddress("fDetector1_PeakDirectionX",            &fDetector1_PeakDirectionX);
        tData->SetBranchAddress("fDetector1_PeakDirectionY",            &fDetector1_PeakDirectionY);
        tData->SetBranchAddress("fDetector1_PeakDirectionZ",            &fDetector1_PeakDirectionZ);
        tData->SetBranchAddress("sDetector1_PeakCreatorParticle",       &sDetector1_PeakCreatorParticle);
        tData->SetBranchAddress("sDetector1_PeakCreatorProcess",        &sDetector1_PeakCreatorProcess);
        
        
        tData->SetBranchAddress("iDetector1_PeakNbr",                   &iDetector1_PeakNbr);
        tData->SetBranchAddress("fDetector1_TriggerEnergy",             &fDetector1_TriggerEnergy);
        tData->SetBranchAddress("fDetector1_TotEneDeposit",             &fDetector1_TotEneDeposit);
        tData->SetBranchAddress("fDetector1_TriggerTime",               &fDetector1_TriggerTime);
        tData->SetBranchAddress("fDetector1_TriggerPositionX",          &fDetector1_TriggerPositionX);
        tData->SetBranchAddress("fDetector1_TriggerPositionY",          &fDetector1_TriggerPositionY);
        tData->SetBranchAddress("fDetector1_TriggerPositionZ",          &fDetector1_TriggerPositionZ);
        tData->SetBranchAddress("fDetector1_TriggerDirectionX",         &fDetector1_TriggerDirectionX);
        tData->SetBranchAddress("fDetector1_TriggerDirectionY",         &fDetector1_TriggerDirectionY);
        tData->SetBranchAddress("fDetector1_TriggerDirectionZ",         &fDetector1_TriggerDirectionZ);
        tData->SetBranchAddress("sDetector1_TriggerParentParticle",     &sDetector1_TriggerParentParticle);
        tData->SetBranchAddress("sDetector1_TriggerParticle",           &sDetector1_TriggerParticle);
        tData->SetBranchAddress("sDetector1_TriggerProcess",            &sDetector1_TriggerProcess);
        
        
        
        
        
        tData->SetBranchAddress("fDetector2_PeakStartTime",             &fDetector2_PeakStartTime);
        tData->SetBranchAddress("fDetector2_PeakEndTime",               &fDetector2_PeakEndTime);
        tData->SetBranchAddress("fDetector2_MaxPeakEnergy",             &fDetector2_MaxPeakEnergy);
        tData->SetBranchAddress("fDetector2_TotPeakEnergy",             &fDetector2_TotPeakEnergy);
        tData->SetBranchAddress("fDetector2_MaxPeakTime",               &fDetector2_MaxPeakTime) ;
        tData->SetBranchAddress("fDetector2_MaxPeakPositionX",          &fDetector2_MaxPeakPositionX);
        tData->SetBranchAddress("fDetector2_MaxPeakPositionY",          &fDetector2_MaxPeakPositionY);
        tData->SetBranchAddress("fDetector2_MaxPeakPositionZ",          &fDetector2_MaxPeakPositionZ);
        tData->SetBranchAddress("sDetector2_PeakCreatorParticle",       &sDetector2_PeakCreatorParticle);
        tData->SetBranchAddress("sDetector2_PeakCreatorProcess",        &sDetector2_PeakCreatorProcess);
        
        
        tData->SetBranchAddress("iDetector2_PeakNbr",                   &iDetector2_PeakNbr);
        tData->SetBranchAddress("fDetector2_TriggerEnergy",             &fDetector2_TriggerEnergy);
        tData->SetBranchAddress("fDetector2_TotEneDeposit",             &fDetector2_TotEneDeposit);
        tData->SetBranchAddress("fDetector2_TriggerTime",               &fDetector2_TriggerTime);
        tData->SetBranchAddress("fDetector2_TriggerPositionX",          &fDetector2_TriggerPositionX);
        tData->SetBranchAddress("fDetector2_TriggerPositionY",          &fDetector2_TriggerPositionY);
        tData->SetBranchAddress("fDetector2_TriggerPositionZ",          &fDetector2_TriggerPositionZ);
        tData->SetBranchAddress("sDetector2_TriggerParticle",           &sDetector2_TriggerParticle);
        tData->SetBranchAddress("sDetector2_TriggerProcess",            &sDetector2_TriggerProcess);
        
        
        
        //---------------------------------------------------------------------- 
        
        
        
        //----------------------------------------------------------------------  
        
        /*--Time of Flight Calculation----
         Recoil Energy follows that:
         Er = 2En/(1+A)^2*(sin^2(theta) + A - cos(theta)*sqrt(A^2 - sin^2(theta)) ) 
         */
        
        
        const double ScatteringAngle = 90;          //--- Angle in degree --
        const double   NeutronEnergy = 530;         //--- Unit: keV ---
        const double     NeutronMass = 939565.38;   //--- Unit: keV ---
        const double  EJ301_Distance = 0.6;           //--- Unit: m ---
        //const double    LAr_Distance = 0.5;         //--- Unit: m ---
        const double    LAr_Distance = 0.7;           //--- Unit: m ---
        const double    Ar_AtomicNbr = 40;
        const double           Theta = ScatteringAngle/180*TMath::Pi();
        
        const double  RecoilEnergy = 2*NeutronEnergy/(1+Ar_AtomicNbr)/(1+Ar_AtomicNbr)*(TMath::Sin(Theta)*TMath::Sin(Theta) + Ar_AtomicNbr - TMath::Cos(Theta)*TMath::Sqrt(Ar_AtomicNbr*Ar_AtomicNbr - TMath::Sin(Theta)*TMath::Sin(Theta)));
        
        
        const double  ScatteredNeutronSpeed = 3*TMath::Sqrt(2*(NeutronEnergy-RecoilEnergy)/NeutronMass)/10; //--- Unit: m/ns ---
        const double           NeutronSpeed = 3*TMath::Sqrt(2*NeutronEnergy/NeutronMass)/10; //--- Unit: m/ns ---
        
        const double               Scattering_TOF = LAr_Distance/NeutronSpeed;
        const double             TOF_Window_Mid   = EJ301_Distance/ScatteredNeutronSpeed;
        const double             TOF_Window_Start = TOF_Window_Mid - 15;
        const double             TOF_Window_End   = TOF_Window_Mid + 15;
        const double      Scattering_Window_Start = Scattering_TOF - 15;
        const double      Scattering_Window_End   = Scattering_TOF + 15;
        
        
        
        const double    Detector1_TriggerThresholdLimit = 7;
        const double    Detector2_TriggerThresholdLimit = 20;
                
        
        const double    PE_Detector_Resolution = para[0];
        const double           Leff_LightYield = para[1];   //-- p.e./keVee --

        
         TRandom1  *RandomGenerator  = new TRandom1(0);
         TRandom1  *PoissonGenerator = new TRandom1(0);
         //----------------------------------------------------------------------    
        
        
        
        Int_t EntryNbr = (Int_t)tData->GetEntries();
        
              
        for(int i=0; i<EntryNbr; i++){
            
            tData->GetEntry(i);
            
            
            Float_t TOF = fDetector2_TriggerTime - fDetector1_TriggerTime;

            
            if( (fDetector2_TriggerEnergy>Detector2_TriggerThresholdLimit)&&(fDetector1_TriggerEnergy>Detector1_TriggerThresholdLimit)&&(TOF>=TOF_Window_Start)&&(TOF<TOF_Window_End)&&(fDetector1_TriggerTime>Scattering_Window_Start)&&(fDetector1_TriggerTime<Scattering_Window_End) )
            {
                Float_t DetectorSmeared_PE = RandomGenerator->Gaus(fDetector1_TotEneDeposit*Leff_LightYield, PE_Detector_Resolution*TMath::Sqrt(fDetector1_TotEneDeposit*Leff_LightYield));
                
                Float_t PMT_PE = PoissonGenerator->Poisson(DetectorSmeared_PE);
                
                //hist_G4Data->Fill(DetectorSmeared_PE);
                hist_G4Data->Fill(PMT_PE);
                
            }
        }
        
            
            Double_t ChiSquare = 0;

        //----------- Ratio uses the histogram area ratio -----------
        const double       Ratio = smoothed_hist_S1->Integral()/hist_G4Data->Integral();
        
        
            for(int i=Min_xAxis; i<=Max_xAxis; i++){
                
                Data_yAxis_Content[i]  = smoothed_hist_S1->GetBinContent(i);
              G4Data_yAxis_Content[i]  = hist_G4Data->GetBinContent(i);
                
            //ChiSquare = ChiSquare + (Data_yAxis_Content[i]- para[2]*G4Data_yAxis_Content[i])*(Data_yAxis_Content[i]- para[2]*G4Data_yAxis_Content[i])/(para[2]*G4Data_yAxis_Content[i]);
                
                if(G4Data_yAxis_Content[i]>0)
              ChiSquare = ChiSquare + (Data_yAxis_Content[i]- Ratio*G4Data_yAxis_Content[i])*(Data_yAxis_Content[i]- Ratio*G4Data_yAxis_Content[i])/(Ratio*G4Data_yAxis_Content[i]);
                
            }
                          
        Hist_FileName->Close();  
                    f->Close(); 
        
        return ChiSquare;
    }
  
        
int main(int argc, char **argv){
    
    const int TotLoopNbr = 100;
    
    // Input file name
 //   ifstream Input(argv[1]);

    // Output file name
    ofstream Output(argv[1]);
    
    Output<<"MinValue,  Detector_Resolution, Leff_LightYield"<<std::endl;
    
    for(int k=0; k<TotLoopNbr; k++){
 
    // Create Minimizer with various fitting algorithms       
    //ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit", "Migrad");
    ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit", "Simplex");
    //ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit", "Minuit2");

    // Set minimizer properties
    min->SetMaxFunctionCalls(500);
    min->SetMaxIterations(500);
    min->SetTolerance(0.1);
    
        const int ParaNbr = 2;    
        
    // Create function wrapper
    ROOT::Math::Functor f(&FitChiSquare, ParaNbr);
    
    double     step[ParaNbr] = {0.05, 0.05};
    double variable[ParaNbr] = {1.2, 1.9};
    
    min->SetFunction(f);
    
    // Set variable names and step size
    min->SetVariable(0,"Detector_Resolution",    variable[0], step[0]);
    min->SetVariable(1,"Leff_LightYield",        variable[1], step[1]);
    //min->SetVariable(2,"Amp",                    variable[2], step[2]);

        
        min->Minimize();
        
            const double *Result = min->X();
        
        Output<<min->MinValue()<<"  "<<Result[0]<<" "<<Result[1]<<"  "<<std::endl;
        
     /*   
    min->SetPrintLevel(2);
    min->PrintResults();
      */

    }
    
    return 0;
    
}
