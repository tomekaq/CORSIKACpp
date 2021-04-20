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

#include <iostream>
#include <sstream>
#include <map>
#include <vector>
#include <algorithm>

using namespace std;

struct ObsLevel {
  TNtupleD* data;
  double x;	
  double y;
  double x2;
  double y2;
  double w;
};



class Detector{
     public:
         double Xcenter;	
         double Ycenter;

         double Xmin;	
         double Ymin;

         double Xmax;	
         double Ymax;

         Detector (double x,double y,double d)
         {
            Xcenter = x;
            Ycenter = y;	
	    Xmin = Xcenter - d/2.0;
            Xmax = Xcenter + d/2.0;
            Ymin = Ycenter - d/2.0;
            Ymax = Ycenter + d/2.0;
         };
         void PrintSize()
         { 
            cout << "Lower border of X: "<<  Xmin << endl
                 << "Upper border of X: "<<  Xmax << endl             
                 << "Lower border of Y: "<<  Ymin << endl
                 << "Upper border of Y: "<<  Ymax << endl;
         };

         bool operator== (const Detector& d1) const 
         {
            if (d1.Xcenter > Xmin && d1.Xcenter < Xmax && d1.Ycenter > Ymin && d1.Ycenter < Ymax) 
                return true;
	    return false;
         };

         bool operator()(const Detector& d1) const
	 {
            if (d1.Xcenter > Xmin && d1.Xcenter < Xmax && d1.Ycenter > Ymin && d1.Ycenter < Ymax) 
               return true;
            return false;
         };
        
         bool operator<(const Detector& d1) const 
         { 
            if (d1.Xcenter > Xmin && d1.Xcenter < Xmax && d1.Ycenter > Ymin && d1.Ycenter < Ymax) 
		return true;
            return false;
         };
~Detector()
{}

};
 

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

//////////////////////// wczytywanie parametrów uzytkownika

   
  int nr = 100;
  int nc = 100;
  float m = 1000;
  float xmin = 0;
  float ymin = 0;
  
  float xmax = 1000;
  float ymax = 10000;
  float d= 100;
  
  cout <<"Enter the number of detectors in row:"<< endl;
  cin >> nr ;
  cout <<"Enter the number of detectors in column:"<< endl;
  cin >> nc ;
  cout <<"Enter the distance between detector [cm]:"<< endl;
  cin >> m ;

  cout <<"Enter your size of detector in cm:"<< endl;
  cin >> d ;
  
  xmin = -(nr*m)/2.0;  // center x point of first detector in a row
  ymin = -(nc*m)/2.0; // center y point of first detector
  
  xmax = (nr*m)/2.0;  // center x point of last detector in a row
  ymax = (nc*m)/2.0; // center y point of last detector

   cout << "---------------------------------\n";
 cout <<"Total number of detectors: "     << nr * nc<< endl;
 cout <<"Distance between detectors [cm]: " << m << endl;
 cout <<"Area of one detector [cm^2]: "     << d*d << endl;
 cout <<"Total area covered by detectors [cm^2]: "     <<  nr*nc*(d*d) << endl;
   cout << "---------------------------------\n";
  cout <<"Center of detector in first row, first column on x [cm]: "     << xmin<< endl;
  cout <<"Lower  range of this detector on x [cm] : "<< xmin - d<< endl;
  cout <<"Upper  range of this detector on x [cm] : "<< xmin + d<< endl;
   cout << "---------------------------------\n";
  cout <<"Center of detector in first row, first column on y [cm] :  "   << ymin << endl;
  cout <<"Lower  range of this detector on y : "<<  ymin - d<< endl;
  cout <<"Upper  range of this detector on y : "<< ymin + d<< endl;
   cout << "---------------------------------\n";


  cout <<"Center of detector in last row, last column on x [cm]: "<< xmax<< endl;
  cout <<"Lower  range of this detectors on x [cm]: "<< xmax - d<< endl;
  cout <<"Upper  range of this detectors on x [cm]: "<< xmax + d<< endl;
   cout << "---------------------------------\n";
  cout <<"Center detector in last row, last column on y [cm]: "<< ymax << endl;
  cout <<"Lower  range of this detectors on y [cm]: "<<  ymax - d<< endl;
  cout <<"Upper  range of this detectors on y [cm]: "<< ymax + d<< endl;
   cout << "---------------------------------\n";

  vector<Detector> vecDetectors;

  for (int i =0 ;i<nr;i++)
  {
     for (int j =0 ;j<nc;j++)
     {
	double k = xmin+i*m;
	double l = ymin+j*m;
        Detector detector(k,l,d);
        vecDetectors.push_back(detector);
     }
  }

  int numberofHits = 0;
  int ShowerCounter = 0;


  crs::MRunHeader Run;
  while (cr.GetRun (Run)) {
    
    crs::MEventHeader Shower;
    while (cr.GetShower(Shower)) {
      
      ++ShowerCounter;  

	ostringstream oFileName;
        oFileName << inputFile.c_str() << "_"
            << Shower.GetEventNumber () << ".root";

        cout << " Writing summary to output file: " << oFileName.str() << endl;
        const int nObsLevel = Shower.GetNObservationLevels();
 	
        map<int, ObsLevel> obsLevel;  
	TFile oFile (oFileName.str ().c_str (), "RECREATE");

      for (int iObsLevel=1; iObsLevel<=nObsLevel; ++iObsLevel) { 
        
        double height = Shower.GetObservationHeight(iObsLevel-1);            
        cout << " init obs-level " << iObsLevel << " at h=" << height << endl;
        ObsLevel emptyLevel;
        ostringstream tTitle, tName;
        tTitle << "Data at level " << iObsLevel;
        tName << "data_" << iObsLevel;
        emptyLevel.data = new TNtupleD(tName.str().c_str(), tTitle.str().c_str(), "id:e:x:y:w");
        obsLevel[iObsLevel] = emptyLevel;

      } // end loop observation levels
      


      crs::TSubBlock Data;
      while (cr.GetData (Data)) {
        
        switch (Data.GetBlockType ()) {
          
            case crs::TSubBlock::ePARTDATA:
            {
              const crs::MParticleBlock& ParticleData = Data;
              crs::MParticleBlock::ParticleListConstIterator iEntry;
              for (iEntry = ParticleData.FirstParticle();
                   iEntry != ParticleData.LastParticle();
                   ++iEntry) {
                
                   if (iEntry->IsParticle()) {
                  
                  crs::MParticle iPart(*iEntry);
                  
                  const int id    = iPart.GetParticleID();
                  const int level = iPart.GetObservationLevel();
                  const double w  = iPart.GetWeight();
                  const double e  = iPart.GetKinEnergy();
                  const double x  = iPart.GetX();
                  const double y  = iPart.GetY();
               

                  Detector particle(x,y,0);
 		  vector<Detector>::iterator it;
	          it = find(vecDetectors.begin(),vecDetectors.end(),particle);
                  if (it != vecDetectors.end())
                  {
			
			numberofHits++;	

                  obsLevel[level].data->Fill(id, e, x, y, w);
                  obsLevel[level].x  += x*w;
                  obsLevel[level].y  += y*w;
                  obsLevel[level].x2 += x*x*w;
                  obsLevel[level].y2 += y*y*w;
                  obsLevel[level].w  += w;
			
                  }


                }
              }
              
              break;
            }
            
            case crs::TSubBlock::eLONG:
              break;
              
            default:
              break;
        } // end data block    
       } // loop data
      oFile.cd();

      for (map<int, ObsLevel>::iterator iLevel = obsLevel.begin();
           iLevel != obsLevel.end();
           ++iLevel) {
        
	iLevel->second.data->Write();
        }

  
      crs::MEventEnd ShowerSummary;
      cr.GetShowerSummary(ShowerSummary);
      const double numberofParticles = ShowerSummary.GetParticles();		
      cout << "---------------------------------\n"
           << "  Number of particle = " << numberofParticles << endl;

   oFile.Close();
    }
  }
  cout << "Numbers of hit: " << numberofHits << endl;
  return 0;
}


		  
		
