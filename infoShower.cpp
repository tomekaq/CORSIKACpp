#include <crsRead/MCorsikaReader.h>

#include <crs/TSubBlock.h>
#include <crs/MRunHeader.h>
#include <crs/MEventHeader.h>
#include <crs/MEventEnd.h>
#include <crs/MParticleBlock.h>
#include <crs/MLongitudinalBlock.h>
#include <crs/MParticle.h>

#include <TFile.h>
#include <TH1D.h>
#include <TNtupleD.h>
#include <TH1.h>
#include <TH2.h>

#include <TRandom.h>
#include <TMath.h>
#include <iostream>
#include <sstream>
#include <map>
#include <vector>
#include <algorithm>

#include <stdlib.h>  

using namespace std;

int main (int argc, char **argv) 
{  
  if (argc<2) {
    cout << "  please specify Corsika file " << endl;
    return 1;
  }
    
  string fname(argv[1]);
  crsRead::MCorsikaReader cr(fname, 3);
  
  string inputFile;
  if (fname.rfind('/')==string::npos) {
    inputFile = fname;
  } else {
    inputFile = fname.substr(fname.rfind('/')+1, 
                             fname.size()-fname.rfind('/'));
  }
  

  crs::MRunHeader Run;
  while (cr.GetRun (Run)) {
    
    crs::MEventHeader Shower;
    while (cr.GetShower(Shower)) {
      
   
	ostringstream oFileName;

        cout << "---------------------------------\n"
           << " Shower info:\n";
           
 	cout << "Value of GetRunNumber: " << Shower.GetRunNumber ()<< endl; 	

 	cout << "Value of  DateStart:  " << Shower.GetDateStart()<< endl; 	
 	cout << "Value of  Version:  " << Shower.GetVersion () << endl; 
 	
 	cout << "Value of primary particle energy : " << Shower.GetEnergy () << endl;
        cout << "Type of primary particle : " <<	Shower.GetParticleId () << endl;
 	
 	cout << "Value of  FirstTarget: " << Shower.GetFirstTarget ()<< endl;
 	cout << "Value of  ZFirst: " << Shower.GetZFirst () << endl;
 	cout << "Value of  SpectralSlope: " << Shower.GetSpectralSlope () << endl;
 	
 	cout << "Value of  ArrayRotation: " << Shower.GetArrayRotation ()<< endl;

        cout << "Value of  Theta: "<<	Shower.GetTheta () << endl;
        cout << "Value of  ThetaMin: "<<	Shower.GetThetaMin () << endl;        
        cout << "Value of  ThetaMax: "<<	Shower.GetThetaMax () << endl;

 	cout << "Value of  Phi: " << Shower.GetPhi () << endl;
 	cout << "Value of  Phi Min: " << Shower.GetPhiMin () << endl; 	
 	cout << "Value of  Phi Max: " << Shower.GetPhiMax () << endl; 	

// wartosc seedów
        
        cout << "Type of HadronicLowEModell : "<<  	Shower.GetHadronicLowEModell ()  << endl;
	cout << "Type of HadronicHighEModell: "<< 	Shower.GetHadronicHighEModell () << endl;

        
	cout << "Value of GetCutoffHadrons: "<< Shower.GetCutoffHadrons ()<< endl;
	cout << "Value of GetCutoffMuons:  "<<   Shower.GetCutoffMuons () << endl;
	cout << "Value of GetCutoffElectrons "    <<   Shower.GetCutoffElectrons () << endl;
	cout << "Value of GetCutoffPhotons"    <<   Shower.GetCutoffPhotons () << endl;	
	
        
        cout << "Value of  GetMultipleScatteringStep: "<<  	Shower.GetMultipleScatteringStep () << endl;
        
	
	cout << "Value of  EFractionThinningH: "<<  	Shower.GetEFractionThinningH () << endl;
	cout << "Value of  EFractionThinningEM: "<< 	Shower.GetEFractionThinningEM () << endl;
	cout << "Value of  WMaxHadronic: "<<  	Shower.GetWMaxHadronic () << endl;
	cout << "Value of  WMaxEM: "  <<  	Shower.GetWMaxEM () << endl;
	cout << "Value of  RMaxThinning: "<<  	Shower.GetRMaxThinning () << endl;
	
    }
    
	crs::MEventEnd ShowerSummary;
	cr.GetShowerSummary(ShowerSummary);
      
	cout << "---------------------------------\n"
           << " Shower Summary info:\n";

	cout << "Value of  GetEventNumber: "<< ShowerSummary.GetEventNumber () << endl;
	cout << "Value of  GetPhotons: "<< 	ShowerSummary.GetPhotons () << endl;
	cout << "Value of  GetElectrons: "<<  ShowerSummary.GetElectrons () << endl;
	cout << "Value of  GetHadrons: "  <<  ShowerSummary.GetHadrons () << endl;
	cout << "Value of  GetMuons: "  <<  	ShowerSummary.GetMuons () << endl;
	cout << "Value of  GetParticles: "<<  ShowerSummary.GetParticles () << endl;    
  }
  return 0;
}
