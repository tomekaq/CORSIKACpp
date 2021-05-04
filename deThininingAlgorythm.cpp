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

struct ObsLevel {
  TNtupleD* data;
  double x;	
  double y;
  double x2;
  double y2;
  double w;
  double theta;
  double time;
};
//void clone(crs::MParticle particle);

class Interface
{
public:
  virtual Interface* clone() const = 0;
};



void PlotParticles (TH1D &hist,double variable);
//void PlotParticles2D (const crs::MParticleBlock &Data, TH2D &hist);
void PlotParticles2D (TH2D &hist,double varHorizontal, double varVertical,bool ScaleLOG);

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

 int noweczastki = 0;
int stareczastki=0;
  crs::MRunHeader Run;
  while (cr.GetRun (Run)) {
    
    crs::MEventHeader Shower;
    while (cr.GetShower(Shower)) {
      
   
	ostringstream oFileName;
        oFileName << inputFile.c_str() << "_dethinTest" << ".root";

        cout << " Writing summary to output file: " << oFileName.str() << endl;
        const int nObsLevel = Shower.GetNObservationLevels();
 	
        map<int, ObsLevel> obsLevel;  
	TFile oFile (oFileName.str ().c_str (), "RECREATE");


	TRandom *random       = new TRandom();

	TH1D hEnergy ("energy", "energy spectra of EM", 50, 0, 1);
	TH1D xHistogram ("x", "x distribution", 100, -10000, 10000);
        TH1D yHistogram ("y", "y distribution", 100, -10000, 10000);
	TH1D wHistogram ("weight", "weight distribution of particles", 500, -3, 10000 );
	//TH1D hEnergy ("energy", "energy spectra of EM", 50, -3, 1);
      //  TH2D* histogramWE = new TH2D("center", "" ,80, 0,10,80,0,6);		
	TH2D* histogramXY = new TH2D("y:x", "" ,100, -10000,10000,100,-10000,10000);		
	


      crs::TSubBlock Data;
      while (cr.GetData (Data)) {
        
        switch (Data.GetBlockType ()) {
          
            case crs::TSubBlock::ePARTDATA:
            { 
              const crs::MParticleBlock& ParticleData = Data;
              crs::MParticleBlock::ParticleListConstIterator iEntry;
              for (iEntry = ParticleData.FirstParticle();
                   iEntry != ParticleData.LastParticle();   ///tu ma być  ParticleData.LastParticle(); 
                   ++iEntry) {
                 
                   if (iEntry->IsParticle()) {
                  
                  crs::MParticle iPart(*iEntry);
                  
                  const int id    = iPart.GetParticleID();
                  const int level = iPart.GetObservationLevel();
                  const double w  = iPart.GetWeight();
                  const double e  = iPart.GetKinEnergy();
                  const double x  = iPart.GetX();
                  const double y  = iPart.GetY();
		  const double theta = iPart.GetTheta();
		  const double t  = iPart.GetTime();
               
		  if (obsLevel.count(level)==0) {
		    cout << " detected new obs-level " << level 
			 << ". Possibly INCLIN-Option " << endl;
		    ObsLevel emptyLevel;
		    ostringstream tTitle, tName;
		    tTitle << "Data at level " << level;
		    tName << "data_" << level;
		    obsLevel[level] = emptyLevel;
		    obsLevel[level].data = new TNtupleD(tName.str().c_str(), 
							tTitle.str().c_str(),
							"id:e:x:y:w:theta:time");
		  } 
		
		for(int i=1; i < w; i++)
		{
			double xNew;
			double yNew;
			double wNew;
			//cout << w<< " " << x + random->Gaus() << " " << y + random->Gaus() << endl;
			xNew = x + random->Gaus();				                
			yNew = y + random->Gaus();
			wNew = 1;
			obsLevel[level].data->Fill(id, e, xNew, yNew, wNew,theta,t);
                  	obsLevel[level].x  += xNew*wNew;
                  	obsLevel[level].y  += yNew*wNew;
                  	obsLevel[level].x2 += xNew*xNew*wNew;
                  	obsLevel[level].y2 += yNew*yNew*wNew;
                  	obsLevel[level].w  += wNew;
			
			PlotParticles (xHistogram, xNew);
			PlotParticles (yHistogram, yNew);
			PlotParticles (wHistogram, wNew);
			PlotParticles (hEnergy, e);
		}
		
		
		PlotParticles (xHistogram, x);
		PlotParticles (yHistogram, y);
		PlotParticles (wHistogram, w);
		PlotParticles (hEnergy, e);

		PlotParticles2D ( *histogramXY, x,y,0);

                obsLevel[level].data->Fill(id, e, x, y, w,theta,t);
                obsLevel[level].x  += x*w;
                obsLevel[level].y  += y*w;
                obsLevel[level].x2 += x*x*w;
                obsLevel[level].y2 += y*y*w;
                obsLevel[level].w  += w;

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
	hEnergy.Write();
	xHistogram.Write();
	yHistogram.Write();
	wHistogram.Write();
	histogramXY->Write();

      oFile.Close();
    }
  }

  return 0;
}

void PlotParticles (TH1D &hist,double variable) {

		hist.Fill (variable);    

}

void PlotParticles2D (TH2D &hist,double varHorizontal, double varVertical,bool ScaleLOG) {
	if (ScaleLOG) {
		hist.Fill (log10 (varHorizontal), varVertical);
	}
	else
	{
		hist.Fill (varHorizontal,varVertical);
        }
}
