/*
* panodisplay.C
* Displays images of simulated air showers taken by an array of PANOSETI telescopes 
* 
* This macro is a work of simulation. Any resemblance to analysis packages,
* living or dead, is purely coincidence.
*
* Author: Nik Korzoun
*/

#include "TCanvas.h"
#include "TColor.h"
#include "TEllipse.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH2.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TMultiGraph.h"
#include "TRandom3.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TTree.h"

#include "iostream"
#include "fstream"

/* dmod(A,B) - A modulo B (double) */
// from VASlamac.h
#define dmod(A,B) ((B)!=0.0?((A)*(B)>0.0?(A)-(B)*floor((A)/(B))\
							 :(A)+(B)*floor(-(A)/(B))):(A))

// Random seed
int seed = 200;
TRandom3 *r = new TRandom3(seed);

// Root file
TFile *f;
TTree *t;
int Ntel;

// Fixed telescope pointing (GrOptics definition)
// float fTel_Zenith = 20;
// float fTel_Azimuth = 180;

// Magnetic Declination in Degrees (unrotate CORSIKA i.e. ARRANG)
float fMag_Dec = 12.77;
// float fMag_Dec = 10.4; // FLWO

// Reconstructed params
float fShower_Xoffset = -99999.;
float fShower_Yoffset = -99999.;
float fShower_Az = -99999.;
float fShower_Ze = -99999.;
float fShower_Xcore = -99999.;
float fShower_Ycore = -99999.;
float fShower_stdP = -99999.;
float fShower_Chi2 = -99999.;

// unit conversion
const float petoadu = 16.0;

/*
* Read root file for displaying images
*/
void readFile(std::string rootfile){

    // load tree
    f = new TFile(rootfile.c_str());
    t = (TTree *)f->Get("tcors");

    // find number of telescopes
    t->Draw("telNumber","","goff");
    Ntel = t->GetV1()[0];
}

//! reduce large angle to intervall 0, 2*pi
// stolen from GM, corsikaIOreader
double redang( double iangle )
{
    if( iangle >= 0 )
    {
        iangle = iangle - int( iangle / ( 2. * M_PI ) ) * 2. * M_PI;
    }
    else
    {
        iangle = 2. * M_PI + iangle + int( iangle / ( 2. * M_PI ) ) * 2. * M_PI;
    }
    
    return iangle;
}

/*
 * Randomly spread arrival direction of photons to simulate PANOSETI PSF
*/
std::tuple<double,double> spread(double positionX, double positionY){

    // axial separation
    double R = TMath::Hypot(positionX,positionY);

    // apprxoimate PANOSETI PSF per axial separation
    double fwhm = 0.0064*(R*R) + 0.036*R + 0.38; //units of mm
    
    // convert from mm to degrees
    // 3.28mm per pixel, 9.9/32 degrees per pixel
    fwhm = 10*(fwhm/3.28)*(9.9/32);

    // https://en.wikipedia.org/wiki/Full_width_at_half_maximum
    double sigma = fwhm/2.355;

    return std::make_tuple(r->Gaus(positionX, sigma),r->Gaus(positionY, sigma));


}
/*
void testspread(double x, double y){
    // histogram
    //TH2I *test = new TH2I("test", "test", 32, -4.95, 4.95, 32, -4.95, 4.95 );
    TH2I *test = new TH2I("test", "test", 600, -4.95, 4.95, 600, -4.95, 4.95 );
    test->SetXTitle( "x" );
    test->SetYTitle( "y" );
    
    // fill image 10000 times
    for(int i=0; i<10000; i++){
        // scatter by PSF
        std::tuple<double,double> coords = spread(x,y);
        double xnew = std::get<0>(coords);
        double ynew = std::get<1>(coords);

        test->Fill(xnew,ynew,petoadu);
        
    }
    test->Draw("COLZ");
    test->ResetStats();
}
*/

/*
* Add night sky background roughly consistent with VERITAS, but scaled down
* to a PANOSETI telescope
*/
TH2I* addNSB(TH2I* image){

    // get telescope size to scale NSB
    t->Draw("telR","","goff");
    //double telrad = t->GetV1()[0];

    int Nbins = image->GetNcells();
    TH2I *newImage = (TH2I*)image->Clone();
    for(int i=1; i<= Nbins ;i++){
		newImage->AddBinContent(i, r->Poisson(1)); // NSB
	}
    image->Delete();
    return newImage;
}

/*
* Add electronic noise to each pixel
*/
TH2I* addElectronics(TH2I* image){

    int Nbins = image->GetNcells();
    TH2I *newImage = (TH2I*)image->Clone();
    for(int i=1; i<= Nbins ;i++){
        //newImage->AddBinContent(i, (int) r->Gaus(0,2.4));// electronics
        newImage->AddBinContent(i, (int)r->Gaus(0,10)); // mean of 8, but simulate pedestal subtraction to reduce to 0
	}
    image->Delete();
    return newImage;
}

/*
* Compute the fractional area of a square contained within a circle
* Assumes side of square has unit length
* Args:
*   R - radius of the circle
*   cx - x coordinate of the circle's origin
*   cy - y coordinate of the circle's origin
*   sx - x coordinate of the square's origin
*   sy - y coordinate of the square's origin    
*/
double intersectionalArea(int R, int cx, int cy, int sx, int sy){
    double xdiff = abs(sx-cx);
    double ydiff = abs(sy-cy);
    // check if square is fully contained in the circle
    if( pow(xdiff+0.5,2) + pow(ydiff+0.5,2) < R*R ){
        return 1.0;
    // check if square if fully outside the circle
    }else if( pow(xdiff-0.5,2) + pow(ydiff-0.5,2) > R*R ){
        return 0.0;
    // else integrate
    }else{
        //
        // ---- MANUAL INTEGRATION ---- SLOW ----
        //
        //
        /*

        // exploit symmetry to look at top right quartercircle
        if(sx < cx || sy < cy){
            sx = cx + xdiff;
            sy = cy + ydiff;
        }

        // if sx > sy, flip sx,sy so function is integrable
        if(xdiff > ydiff){
            sx = cx + ydiff;
            sy = cy + xdiff;
        }

        // integration bounds of square
        double xi = sx - 0.5;
        double xf = sx + 0.5;
        double yi = sy - 0.5;
        double yf = sy + 0.5;

        // integration bounds of circle
        TF1 circle = TF1("circle", "pow([0]*[0]-(x-[1])*(x-[1]),0.5)+[2]", xi, xf);
        circle.SetParameters(R,cx,cy);

        // draw
        circle.SetMinimum(yi);
        circle.SetMaximum(yf);
        circle.SetFillColor(kRed);
        circle.SetFillStyle(3004);
        //circle.Draw("FC");

        
        // integrate piecewise
        if(circle.Eval(xi) > yf){
            // find intersection point
            double xcrit = circle.GetX(yf);
            return yf*(xcrit-xi) + circle.Integral(xcrit, xf) - yi;
        }else{
            return circle.Integral(xi,xf) - yi;
        }
        */ 

        //
        // ---- LOOKUP INTEGRATION ---- FASTER ----
        // ---- VALID FOR 32x32 CAMERA SUBDIVIDING EACH PIXEL TO 5x5 AND APERTURE RADIUS 2 ----
        //
    
        // only five cases that are not 0,1
        if(xdiff==0 && ydiff==2){
            return 0.478967;
        }else if(xdiff==2 && ydiff==0){
            return 0.478967;
        }else if(xdiff==1 && ydiff==2){
            return 0.198797;
        }else if(xdiff==2 && ydiff==1){
            return 0.198797;
        }else if(xdiff==1 && ydiff==1){
            return 0.984969;
        }else{
            std::cout<<"WARNING: CANNOT FIND INTEGRATION"<<std::endl;
            std::cout<<"xdiff: "<<xdiff<<" ydiff: "<<ydiff<<std::endl;
            return 0;
        }
        
    }
    
}

/*
* Clean image according to ADU thresholds
*/
TH2I* clean(TH2I* image){
    
    int Nbins = image->GetNcells();
    TH2I *newImage = (TH2I*)image->Clone();
    int binsX = newImage->GetNbinsX(); // = 32
    int binsY = binsX; // square camera

    // aperture cleaning
    // https://arxiv.org/pdf/1506.07476.pdf, section 3.1

    // subdivide pixels into NxN subpixels where N = [angular pixel width]/[Aperture radius/2]
    // choose Aperture radius to be approximately the width of a gamma-ray shower - paper suggests 0.12 degrees
        // might want to increase this for panoseti - higher energy showers will be larger
    // therefore N = [0.31]/[0.06] ~=~ 5

    int N = 5;
    double apertureRadius = 0.12; // units of degrees
    int R = apertureRadius/0.06; // units of subpixels

    // get telescope size to scale NSB
    //t->Draw("telR","","goff");
    //double telrad = t->GetV1()[0];

    double NSB = petoadu*0.06; // mean value of NSB per pixel - ADU
    double readoutNoise = 10; // standard deviation of detector readout noise - ADU

    int SNR = 7; // signal to noise ratio required to keep pixel

    // create subdivided image
    TH2I* dividedImage = new TH2I("div", "div", binsX*N, -4.95, 4.95, binsY*N, -4.95, 4.95 );
    // loop over subpixels
    for(int i = 1; i<=binsX*N; i++){
        for(int j = 1; j<=binsY*N; j++){
            dividedImage->SetBinContent(i,j,image->GetBinContent(1+(i-1)/N,1+(j-1)/N));
        }
    }

    // remove pixels which are below image threshold
    std::vector<int> removeMe;

    // loop over pixels
    for(int i = 1; i<=binsX; i++){
        for(int j = 1; j<=binsY; j++){

            // check if one or more subpixels exceeds image threshold
            bool signal = false;
            int checkBin = newImage->GetBin(i,j);

            // loop over subpixels in pixel
            for(int k=N*(i-1)+1; k<= N*i; k++){
                for(int l=N*(j-1)+1; l<= N*j; l++){

                    double binSizeAvg = 0.;
                    double imageThreshold = 0.;

                    // equations 5, 6
                    for(int m=k-R; m<= k+R; m++){
                        for(int n=l-R; n<= l+R; n++){
                            // do not check if one pixel exceeded threshold, and do not select pixels outside the camera
                            if(!signal && m>=1 && m<=binsX*N && n>=1 && n<=binsY*N){
                                double w = intersectionalArea(R,k,l,m,n)/(N*N);
                                binSizeAvg += w * dividedImage->GetBinContent(m,n);
                                imageThreshold += w * (readoutNoise*readoutNoise + NSB);
                            }
                        }
                    }

                    imageThreshold = sqrt(imageThreshold);

                    // check if any subpixel exceeds image threshold
                    if(binSizeAvg > SNR * imageThreshold){
                        // if a subpixel already exceeds the threshold, we can move on to other pixels
                        signal = true;
                    }
                }
            }
            // if this point is reached, there is no signal in any subpixel of a pixel
            if(!signal){
                removeMe.push_back(checkBin);
            }
        }
    }

    // remove pixels which fail threshold check
    for(int i=0; i<(int)removeMe.size(); i++){
        newImage->SetBinContent(removeMe[i], 0);
    }
    
    // check for negative pixels
    removeMe.clear();
    for(int i=1; i<=binsX; i++){
		for(int j=1; j<=binsY; j++){
            int checkBin = newImage->GetBin(i,j);
            int binSize = newImage->GetBinContent(checkBin);
            
            if(binSize<0){
                removeMe.push_back(checkBin);
            }
        }
    }
    
    // remove pixels which are negative
    for(int i=0; i<(int)removeMe.size(); i++){
        newImage->SetBinContent(removeMe[i], 0);
    }

    // check for isolated pixels
    removeMe.clear();
    for(int i=1; i<=binsX; i++){
		for(int j=1; j<=binsY; j++){
            int checkBin = newImage->GetBin(i,j);
            int binSize = newImage->GetBinContent(checkBin);
            // make sure pixel has p.e. before checking to remove
            if(binSize!=0){
                bool remove = true;
                // get neighbors
                for (int p=i-1; p<=i+1; p++){
                    for (int q=j-1; q<=j+1; q++){
                        // do not add central pixel as neighbor
                        if(p!=i && q!=j){
                            // stay in bounds of image
                            if(p>=1 && p<=binsX && q>=1 && q<=binsY){
                                // find a neighbor with pixels in it
                                if(newImage->GetBinContent(newImage->GetBin(p,q)) != 0){
                                    remove = false;
                                }
                            }
                        }
                    }    
                }
                if (remove){
                    removeMe.push_back(checkBin);
                }
            }
        }
    }
    // remove pixels which are isolated
    for(int i=0; i<(int)removeMe.size(); i++){
        newImage->SetBinContent(removeMe[i], 0);
    }
    
    // discard image if there are fewer than 3 pixels
    int Nimagepix=0;
    for(int i = 1; i<=binsX; i++){
        for(int j = 1; j<=binsY; j++){
            int binSize = newImage->GetBinContent(i,j);
            if(binSize!=0){
                Nimagepix++;
            }
        }    
    }
    if(Nimagepix < 3){
        newImage->Reset();
    }

    image->Delete();
    dividedImage->Delete();
    return newImage;
    

    /*
    // this method is closer to how VERITAS works
    // threshold for image pixels

    // get telescope size to scale NSB
    t->Draw("telR","","goff");
    double telrad = t->GetV1()[0];

    int imageThreshold = (((telrad/0.25)*(telrad/0.25)*.06)+1)*5; // corresponds to ~5 sigma, where sigma is electronics noise (1 p.e.)
    //double imageThreshold = 5; // 5 corresponds to ~5 sigma, where sigma is electronics noise (1 p.e.)
    // pixel can have this many p.e. if it borders an image pixel
    int borderThreshold = (((telrad/0.25)*(telrad/0.25)*.06)+1)*3; // 6 corresponds to ~2.5 sigma, same as above
    //double borderThreshold = 2.5; // 2.5 corresponds to ~2.5 sigma, same as above

    int Nbins = image->GetNcells();
    TH2I *newImage = (TH2I*)image->Clone();
    int binsX = newImage->GetNbinsX();
    int binsY = binsX;

    // remove pixels which are below image threshold unless they neighbor a pixel above image threshold and are themselves above border threshold
    std::vector<int> removeMe;
	for(int i=1; i<=binsX; i++){
		for(int j=1; j<=binsY; j++){
            int checkBin = newImage->GetBin(i,j);
            int binSize = newImage->GetBinContent(checkBin);

            bool remove = true;
            // check if pixel is above image threshold
            if(binSize>=imageThreshold){
                remove = false;
            // check if pixel is above border threshold
            }else if(binSize>=borderThreshold){
                // check if a neighbor is above image threshold
                // get neighbors
                for (int p=i-1; p<=i+1; p++){
                    for (int q=j-1; q<=j+1; q++){
                        // do not add central pixel as neighbor
                        if(p!=i && q!=j){
                            // stay in bounds of image
                            if(p>=1 && p<=binsX && q>=1 && q<=binsY){
                                int neighbor = newImage->GetBinContent(newImage->GetBin(p,q));
                                // check if pixel borders a pixel above image threshold)
                                if (neighbor >= imageThreshold){
                                    remove = false;
                                } // else it gets removed
                            }
                        }
                    }    
                }
            }// else it gets removed

            // remove pixels
            if(remove){
                removeMe.push_back(checkBin);
            }
        }
    }
    // remove pixels which fail threshold check
    for(int i=0; i<(int)removeMe.size(); i++){
        newImage->SetBinContent(removeMe[i], 0);
    }
    
    // check for isolated pixels
    removeMe.clear();
    for(int i=1; i<=binsX; i++){
		for(int j=1; j<=binsY; j++){
            int checkBin = newImage->GetBin(i,j);
            int binSize = newImage->GetBinContent(checkBin);
            // make sure pixel has p.e. before checking to remove
            if(binSize!=0){
                bool remove = true;
                // get neighbors
                for (int p=i-1; p<=i+1; p++){
                    for (int q=j-1; q<=j+1; q++){
                        // do not add central pixel as neighbor
                        if(p!=i && q!=j){
                            // stay in bounds of image
                            if(p>=1 && p<=binsX && q>=1 && q<=binsY){
                                // find a neighbor with pixels in it
                                if(newImage->GetBinContent(newImage->GetBin(p,q)) != 0){
                                    remove = false;
                                }
                            }
                        }
                    }    
                }
                if (remove){
                    removeMe.push_back(checkBin);
                }
            }
        }
    }
    // remove pixels which are isolated
    for(int i=0; i<(int)removeMe.size(); i++){
        newImage->SetBinContent(removeMe[i], 0);
    }

    image->Delete();
    return newImage;
    */

    /*
    // This method flatly removes pixels below a certain threshold
    // This is primarily for testing and debugging
    
    int imageThreshold = 4;
    
    int Nbins = image->GetNcells();
    TH2I *newImage = (TH2I*)image->Clone();
    int binsX = newImage->GetNbinsX();
    int binsY = binsX;

    // remove pixels which are below image threshold 
    std::vector<int> removeMe;
	for(int i=1; i<=binsX; i++){
		for(int j=1; j<=binsY; j++){
            int checkBin = newImage->GetBin(i,j);
            int binSize = newImage->GetBinContent(checkBin);

            bool remove = true;
            // check if pixel is above image threshold
            if(binSize>=imageThreshold){
                remove = false;
            }else{// it gets removed
                removeMe.push_back(checkBin);
            }
        }
    }

    // remove pixels which fail threshold check
    for(int i=0; i<(int)removeMe.size(); i++){
        newImage->SetBinContent(removeMe[i], 0);
    }

    image->Delete();
    return newImage;

    */
}


/*
* Attempt image parameterization
* returns tuple which stores
* meanx, sigmax, meany, sigmay, angle, size, length, width
*/

std::tuple<double, double, double, double, double, double, double, double, double, double, double, double> parameterize(TH2I* image){
	//	Begin moment analysis

	double sumsig = 0;
	double sumxsig = 0;
	double sumysig = 0;
	double sumx2sig = 0;
	double sumy2sig = 0;
	double sumxysig = 0;
	double sumx3sig = 0;
	double sumy3sig = 0;
	double sumx2ysig = 0;
	double sumxy2sig = 0;

	double xmean = 0;
	double ymean = 0;
	double x2mean = 0;
	double y2mean = 0;
	double xymean = 0;
	double x3mean = 0;
	double y3mean = 0;
	double x2ymean = 0;
	double xy2mean = 0;

	// loop over all pixels
    int bins = image->GetNbinsX();
	for(int j=1; j<=bins; j++){
		for(int k=1; k<=bins; k++){

			double xi = image->GetXaxis()->GetBinCenter(j);
			double yi = image->GetYaxis()->GetBinCenter(k);

			const double si = image->GetBinContent(image->GetBin(j,k));
			sumsig+=si;

			const double sixi = si * xi;
			const double siyi = si * yi;

			sumxsig += sixi;
			sumysig += siyi;

			const double sixi2 = sixi * xi;
			const double siyi2 = siyi * yi;
			const double sixiyi = sixi * yi;

			sumx2sig += sixi2;
			sumy2sig += siyi2;
			sumxysig += sixiyi;

			sumx3sig += sixi2 * xi;
			sumy3sig += siyi2 * yi;
			sumx2ysig += sixi2 * yi;
			sumxy2sig += siyi2 * xi;
		}
	}

	// image parameter calculations
	if(sumsig > 0. ){
		xmean = sumxsig / sumsig;
		ymean = sumysig / sumsig;
		x2mean = sumx2sig / sumsig;
		y2mean = sumy2sig / sumsig;
		xymean = sumxysig / sumsig;
		x3mean = sumx3sig / sumsig;
		y3mean = sumy3sig / sumsig;
		x2ymean = sumx2ysig / sumsig;
		xy2mean = sumxy2sig / sumsig;
	}
	double xmean2 = xmean * xmean;
	double ymean2 = ymean * ymean;
	double meanxy = xmean * ymean;

	double sdevx2 = x2mean - xmean2;
	double sdevy2 = y2mean - ymean2;
	double sdevxy = xymean - meanxy;
	double sdevx3 = x3mean - 3.0*xmean*x2mean + 2.0*xmean*xmean2;
	double sdevy3 = y3mean - 3.0*ymean*y2mean + 2.0*ymean*ymean2;
	double sdevx2y = x2ymean - 2.0*xymean*xmean + 2.0*xmean2*ymean - x2mean*ymean;
	double sdevxy2 = xy2mean - 2.0*xymean*ymean + 2.0*xmean*ymean2 - xmean*y2mean;

	//Table 6 - Fegan, David J. (1997)
	double d = sdevy2 - sdevx2;
	double z = sqrt(d*d + 4.0*sdevxy*sdevxy);
	double u = 1.0 + d/z;
	double v = 2.0-u;
	double w = sqrt( (y2mean-x2mean)*(y2mean-x2mean) * 4.0*xymean*xymean );

    double dist = sqrt(xmean2 + ymean2); 
	double azwidth = sqrt( (xmean2*y2mean - 2.0*xmean*ymean*xymean + x2mean*ymean2) / (dist*dist) );
	double akwidth = sqrt( (x2mean + y2mean - w)/2.0 );

    // parameterize orientation
    double ac = (d+z)*ymean + 2.0*sdevxy*xmean;
	double bc = 2.0*sdevxy*ymean - (d-z)*xmean;
	double cc = sqrt(ac*ac + bc*bc);
	double cosphi = bc/cc;
	double sinphi = ac/cc;
    double tanphi = ((d+z)*ymean + 2.0*sdevxy*xmean) / (2.0*sdevxy*ymean - (d-z)*xmean);

	double phi = atan(tanphi);
    phi = redang(phi);

    double length = sqrt( (sdevx2 + sdevy2 + z)/2.0 );
	double width = sqrt( (sdevx2 + sdevy2 - z)/2.0 );
    double miss = fabs(-sinphi *xmean + cosphi*ymean);
    if(miss > dist){
        miss = dist; // weird rounding error
    }
	double sinalpha = miss/dist;
    double alpha = fabs(TMath::RadToDeg() * asin(sinalpha));

    phi = fabs(TMath::RadToDeg()*phi);
    return std::make_tuple(xmean, sqrt(sdevx2), ymean, sqrt(sdevy2), phi, sumsig, length, width, miss, dist, azwidth, alpha);
}

/*
* Print some useful information about the shower
*/
TString showerInfo(int eventNumber){

    if(!f){
        std::cout<< "error reading file, try readFile(\"rootfile.root\")" <<std::endl;
        return nullptr;
    }
    if(!t){
        std::cout<<"error reading tree"<<std::endl;
        return nullptr;
    }

    auto condition = Form("eventNumber==%d", eventNumber);
    t->Draw("energy:az:ze",condition,"goff");
    double energy = t->GetV1()[0];
    double az = t->GetV2()[0];
    //CORSIKA to GrOptics
    az=TMath::RadToDeg()*redang(M_PI - redang(TMath::DegToRad()*az - M_PI));
    az=az + fMag_Dec*TMath::DegToRad(); // unrotate array to correct for magnetic declination, ARRANG
    double ze = t->GetV3()[0];
    // double az = fTel_Azimuth;
    // double ze = fTel_Zenith;

    t->Draw("xCore:yCore",condition,"goff");
    // CORSIKA to GrOptics
    double xCore = -1*t->GetV2()[0];
    double yCore = t->GetV1()[0];


    TString info = Form(
        "=======================\n"
        "EVENT:\t\t%d\n"
        "ENERGY:\t\t%.2f GeV\n"
        "AZIMUTH:\t\t%.2f degrees\n"
        "ZENITH:\t\t%.2f degrees\n"
        "CORE LOCATION:\t%.2f m\t%.2f m\n",
        eventNumber, energy, az, ze, xCore, yCore);

    return info;
}

/*
* Draw a map of telescope and shower core positions
*/
TMultiGraph* eventMap(int eventNumber){

    if(!f){
        std::cout<< "error reading file, try readFile(\"rootfile.root\")" <<std::endl;
        return nullptr;
    }
    if(!t){
        std::cout<<"error reading tree"<<std::endl;
        return nullptr;
    }

    // read event data
    auto condition = Form("eventNumber==%d", eventNumber);
    t->Draw("telXpos:telYpos:xCore:yCore",condition,"goff");

    TMultiGraph* map = new TMultiGraph();
    TGraph* telescopes = new TGraph(Ntel);
    TGraph* shower = new TGraph(1);
    map->SetTitle("Event Map");

    // fill telescope positions
    TLatex *l1;
    TLatex *l2;
    telescopes->SetMarkerStyle(20);
    telescopes->SetMarkerSize(3);
    for(int i=0;i<Ntel;i++){

        // corsika X is North and Y is West
        // here we want to plot Y-North and X-East
        double x = -1*t->GetV2()[i];
        double y = t->GetV1()[i];
        telescopes->SetPoint(i,x,y);

        // label point
        l1 = new TLatex(x-30*cos(atan2(y,x)),y-30*sin(atan2(y,x)),Form("T%d",i+1));
        l1->SetTextSize(0.025);
        l1->SetTextFont(42);
        l1->SetTextAlign(21);
        /*l2 = new TLatex(x,y-75,Form("(%.2f,%.2f)",x,y));
        l2->SetTextSize(0.025);
        l2->SetTextFont(42);
        l2->SetTextAlign(21);*/

        telescopes->GetListOfFunctions()->Add(l1);
        //telescopes->GetListOfFunctions()->Add(l2);
        
    }
    // shower position
    //shower->SetMarkerStyle(5);
    //shower->SetMarkerSize(3);
    shower->SetMarkerStyle(47);
    shower->SetMarkerColor(kRed);
    // CORSIKA to GrOptics
    double x = -1*t->GetV4()[0];
    double y = t->GetV3()[0];
    shower->SetPoint(0,x,y);

    // label point
    /*
    l1= new TLatex(x,y-50,"shower core");
    l1->SetTextSize(0.025);
    l1->SetTextFont(42);
    l1->SetTextAlign(21);
    l2 = new TLatex(x,y-75,Form("(%.2f,%.2f)",x,y));
    l2->SetTextSize(0.025);
    l2->SetTextFont(42);
    l2->SetTextAlign(21);

    shower->GetListOfFunctions()->Add(l1);
    shower->GetListOfFunctions()->Add(l2);
    */
   
    // add to multigraph
    map->Add(telescopes);
    map->Add(shower);
    map->GetXaxis()->SetTitle( "E (m)" );
    map->GetYaxis()->SetTitle( "N (m)" );

    return map;
}
double slaDranrm( double angle )
/*
 **  - - - - - - - - - -
 **   s l a D r a n r m
 **  - - - - - - - - - -
 **
 **  Normalize angle into range 0-2 pi.
 **
 **  (double precision)
 **
 **  Given:
 **     angle     double      the angle in radians
 **
 **  The result is angle expressed in the range 0-2 pi (double).
 **
 **  Defined in slamac.h:  D2PI, dmod
 **
 **  Last revision:   19 March 1996
 **
 **  Copyright P.T.Wallace.  All rights reserved.
 */
{
	double w;
	
	w = dmod( angle, 6.2831853071795864769252867665590057683943387987502 ); // 2pi - D2PI - from VASlamac.h
	return ( w >= 0.0 ) ? w : w + 6.2831853071795864769252867665590057683943387987502; // 2pi - D2PI - from VASlamac.h
}

void slaDtp2s( double xi, double eta, double raz, double decz,
			   double* ra, double* dec )
/*
 **  - - - - - - - - -
 **   s l a D t p 2 s
 **  - - - - - - - - -
 **
 **  Transform tangent plane coordinates into spherical.
 **
 **  (double precision)
 **
 **  Given:
 **     xi,eta      double   tangent plane rectangular coordinates
 **                          (xi and eta are equivalent to VERITAS's
 **                          derotated camera coordinates Xderot
 **                          and Yderot, NOT the tangent plane RA/Dec,
 **                          (e.g. Xderot + wobbleWest + TargetRA))
 **     raz,decz    double   spherical coordinates of tangent point
 **
 **  Returned:
 **     *ra,*dec    double   spherical coordinates (0-2pi,+/-pi/2)
 **
 **  Called:  slaDranrm
 **
 **  Last revision:   3 June 1995
 **
 **  Copyright P.T.Wallace.  All rights reserved.
 */
{
	double sdecz, cdecz, denom;
	
	sdecz = sin( decz );
	cdecz = cos( decz );
	denom = cdecz - eta * sdecz;
	*ra = slaDranrm( atan2( xi, denom ) + raz );
	*dec = atan2( sdecz + eta * cdecz, sqrt( xi * xi + denom * denom ) );
}

/*
    "borrowed" from eventDisplay
    VSimpleStereoReconstructor.cpp
      
    reconstruction of shower direction
    Hofmann et al 1999, Method 1 (HEGRA method)
    shower direction by intersection of image axes
    shower core by intersection of lines connecting reconstruced shower
    direction and image centroids
    corresponds to rcs_method4 in VArrayAnalyzer
*/
bool reconstruct_direction( unsigned int i_ntel,
        double fTelElevation,
        double fTelAzimuth,
		double* img_size,
		double* img_cen_x,
		double* img_cen_y,
		double* img_phi,
		double* img_length,
		double* img_width)
{
	
	// make sure that all data arrays exist
	if( !img_size || !img_cen_x || !img_cen_y
			|| !img_phi || !img_width || !img_length)
	{
		//std::cout << "Missing data: cannot reconstruct event."<<std::endl;
		return false;
	}
	
	float xs = 0.;
	float ys = 0.;
	
	// fill data std::vectors for direction reconstruction
	std::vector< float > m;
	std::vector< float > x;
	std::vector< float > y;
	std::vector< float > s;
	std::vector< float > l;
	for( unsigned int i = 0; i < i_ntel; i++ )
	{
        // length == length and width == width protect against negative estimators of the variance - NK
		if( img_size[i] > 0. && img_length[i] == img_length[i] && img_width[i] == img_width[i]) 
		{
			s.push_back( img_size[i] );
			x.push_back( img_cen_x[i] );
			y.push_back( img_cen_y[i] );
			// in VArrayAnalyzer, we do a recalculatePhi. Is this needed (for LL)?
			// (not needed, but there will be a very small (<1.e-5) number of showers
			// with different phi values (missing accuracy in conversion from float
			// to double)
			if( cos(img_phi[i]) != 0. )
			{
				m.push_back( sin(img_phi[i]) / cos(img_phi[i]) );
			}
			else
			{
				m.push_back( 1.e9 );
			}
			if( img_length[i] > 0. )
			{
				l.push_back( img_width[i] / img_length[i] );
			}
			else
			{
				l.push_back( 1. );
			}
		}
	}
	// are there enough images the run an array analysis
	if( s.size() < 2 )
	{
		//std::cout << "Not enough images for reconstruction."<<std::endl;
		return false;
	}
	
	// don't do anything if angle between image axis is too small (for 2 images only)
    float fiangdiff;
	if( s.size() == 2 )
	{
		fiangdiff = -1.*fabs( atan( m[0] ) - atan( m[1] ) ) * TMath::RadToDeg();
	}
	else
	{
		fiangdiff = 0.;
	}
	
	///////////////////////////////
	// direction reconstruction
	////////////////////////////////////////////////
	// Hofmann et al 1999, Method 1 (HEGRA method)
	// (modified weights)
	
	float itotweight = 0.;
	float iweight = 1.;
	float ixs = 0.;
	float iys = 0.;
	float iangdiff = 0.;
	float b1 = 0.;
	float b2 = 0.;
	std::vector< float > v_xs;
	std::vector< float > v_ys;
	float fmean_iangdiff = 0.;
	float fmean_iangdiffN = 0.;
	
	for( unsigned int ii = 0; ii < m.size(); ii++ )
	{
		for( unsigned int jj = 1; jj < m.size(); jj++ )
		{
			if( ii >= jj )
			{
				continue;
			}
			
			// check minimum angle between image lines; ignore if too small

			iangdiff = fabs( atan( m[jj] ) - atan( m[ii] ) );
			if( iangdiff < 0 ||
					fabs( 180. * TMath::DegToRad() - iangdiff ) < 0 )
			{
				continue;
			}
			// mean angle between images
			if( iangdiff < 90. * TMath::DegToRad() )
			{
				fmean_iangdiff += iangdiff * TMath::RadToDeg();
			}
			else
			{
				fmean_iangdiff += ( 180. - iangdiff * TMath::RadToDeg() );
			}
			fmean_iangdiffN++;
			
			// weight is sin of angle between image lines
			iangdiff = fabs( sin( fabs( atan( m[jj] ) - atan( m[ii] ) ) ) );
			
			b1 = y[ii] - m[ii] * x[ii];
			b2 = y[jj] - m[jj] * x[jj];
			
			// line intersection
			if( m[ii] != m[jj] )
			{
				xs = ( b2 - b1 )  / ( m[ii] - m[jj] );
			}
			else
			{
				xs = 0.;
			}
			ys = m[ii] * xs + b1;

			iweight  = 1. / ( 1. / s[ii] + 1. / s[jj] ); // weight 1: size of images
			iweight *= ( 1. - l[ii] ) * ( 1. - l[jj] ); // weight 2: elongation of images (width/length)
			iweight *= iangdiff;                      // weight 3: angular differences between the two image axis
			iweight *= iweight;                       // use squared value

			ixs += xs * iweight;
			iys += ys * iweight;
			itotweight += iweight;
			
			v_xs.push_back( xs );
			v_ys.push_back( ys );
		}
	}
	// average difference between image pairs
	if( fmean_iangdiffN > 0. )
	{
		fmean_iangdiff /= fmean_iangdiffN;
	}
	else
	{
		fmean_iangdiff = 0.;
	}
	if( s.size() > 2 )
	{
		fiangdiff = fmean_iangdiff;
	}
	// check validity of weight
	if( itotweight > 0. )
	{
		ixs /= itotweight;
		iys /= itotweight;
		fShower_Xoffset = ixs;
		fShower_Yoffset = iys;
	}
	else
	{
		//std::cout << "Image weights invalid"<<std::endl;
        return false;
	}
	
    // (y sign flip!)
    fShower_Yoffset = -1.*fShower_Yoffset;

    double el = 0.;
	double az = 0.;
	slaDtp2s( -1.* fShower_Xoffset * TMath::DegToRad(),
							 fShower_Yoffset * TMath::DegToRad(),
							 fTelAzimuth * TMath::DegToRad(),
							 fTelElevation * TMath::DegToRad(),
							 &az, &el );
    
    fShower_Az = slaDranrm( az ) * TMath::RadToDeg();
    fShower_Ze = 90 - el* TMath::RadToDeg();

	return true;
}

/*
    "borrowed" from eventDisplay
    VGrIsuAnalyzer.cpp
      
    helper function for core reconstruction
*/

float rcs_perpendicular_dist( float xs, float ys, float xp, float yp, float m )
/* function to determine perpendicular distance from a point
   (xs,ys) to a line with slope m and passing through the point
   (xp,yp). Calculations in two dimensions.
*/
{
	float theta = 0.;
	float x = 0.;
	float y = 0.;
	float d = 0.;
	float dl = 0.;
	float dm = 0.;
	
	theta = atan( m );
	/* get direction cosines of the line from the slope of the line*/
	dl = cos( theta );
	dm = sin( theta );
	
	/* get x and y components of std::vector from (xp,yp) to (xs,ys) */
	x = xs - xp;
	y = ys - yp;
	
	/* get perpendicular distance */
	d = fabs( dl * y - dm * x );
	
	return d;
}

/*
    "borrowed" from eventDisplay
    VGrIsuAnalyzer.cpp
      
    helper function for core reconstruction
*/

int rcs_perpendicular_fit( std::vector<float> x, std::vector<float> y, std::vector<float> w, std::vector<float> m,
		unsigned int num_images, float* sx, float* sy, float* std )
/*
RETURN= 0 if no faults
ARGUMENT=x[10]     = x coor of point on line
         y[10]     = y coor of point on line
     w[10]     = weight of line
     m[10]     = slope of line
     num_images= number of lines
     sx        = x coor of point with minim. dist.
     sx        = y coor of point with minim. dist
     std       = rms distance from point to lines
    This procedure finds the point (sx,sy) that minimizes the square of the
perpendicular distances from the point to a set of lines.  The ith line
passes through the point (x[i],y[i]) and has slope m[i].
*/
{

	float totweight = 0.;
	float a1 = 0.;
	float a2 = 0.;
	float b1 = 0.;
	float b2 = 0.;
	float c1 = 0.;
	float c2 = 0.;
	float gamma = 0.;
	float D = 0.;
	float m2 = 0.;
	float d = 0.0;
	
	/* initialize variables */
	*sx = -999.;
	*sy = -999.;
	*std = 0.0;
	
	// check length of std::vectors
	
	if( x.size() == num_images && y.size() == num_images && w.size() == num_images && m.size() == num_images )
	{
		for( unsigned int i = 0; i < num_images; i++ )
		{
			totweight = totweight + w[i];
			
			m2 = m[i] * m[i];
			gamma  = 1.0 / ( 1. + m2 );
			
			/* set up constants for array  */
			D = y[i] - ( m[i] * x[i] );
			
			a1 = a1 + ( w[i] *  m2 * gamma );
			a2 = a2 + ( w[i] * ( -m[i] ) * gamma );
			b1 = a2;
			b2 = b2 + ( w[i] *  gamma );
			c1 = c1 + ( w[i] * D * m[i] * gamma );
			c2 = c2 + ( w[i] * ( -D ) * gamma );
			
		}
		/* do fit if have more than one telescope */
		if( ( num_images > 1 ) )
		{
			/* completed loop over images, now normalize weights */
			a1 = a1 / totweight;
			b1 = b1 / totweight;
			c1 = c1 / totweight;
			a2 = a2 / totweight;
			b2 = b2 / totweight;
			c2 = c2 / totweight;
			
			/*
			The source coordinates xs,ys should be solution
			of the equations system:
			a1*xs+b1*ys+c1=0.
			a2*xs+b2*ys+c2=0.
			*/
			
			*sx = -( c1 / b1 - c2 / b2 ) / ( a1 / b1 - a2 / b2 );
			*sy = -( c1 / a1 - c2 / a2 ) / ( b1 / a1 - b2 / a2 );
			
			/* std is average of square of distances to the line */
			for( unsigned int i = 0; i < num_images; i++ )
			{
				d = ( float )rcs_perpendicular_dist( ( float ) * sx, ( float ) * sy,
													 ( float )x[i], ( float )y[i], ( float )m[i] );
				*std = *std + d * d * w[i];
			}
			*std = *std / totweight;
		}
	}
	else
	{
		std::cout <<  "VGrIsuAnalyzer::rcs_perpendicular_fit error in std::vector length" << std::endl;
	}
	return 0;
}

/*
    "borrowed" from eventDisplay
    VGrIsuAnalyzer.cpp
      
    helper function for core reconstruction
*/

void setup_matrix( float matrix[3][3], float dl, float dm, float dn, bool bInvers )
{
	float sv = 0.;
	
	/* sv is the projection of the primary std::vector onto the xy plane */
	
	sv = sqrt( dl * dl + dm * dm );
	
	if( sv > 1.0E-09 )
	{
	
		/* rotation about z axis to place y axis in the plane
		   created by the vertical axis and the direction of the
		   incoming primary followed by a rotation about the new x
		   axis (still in the horizontal plane) until the new z axis
		   points in the direction of the primary.
		*/
		
		matrix[0][0] = -dm / sv;
		matrix[0][1] = dl / sv;
		matrix[0][2] = 0;
		
		matrix[1][0] = dn * dl / sv ;
		matrix[1][1] = dn * dm / sv;
		matrix[1][2] =  - sv;
		
		matrix[2][0] = -dl;
		matrix[2][1] = -dm;
		matrix[2][2] = -dn;
		
	}
	/* for verital incident showers, return identity matrix */
	else
	{
		matrix[0][0] = 1;
		matrix[0][1] = 0;
		matrix[0][2] = 0;
		
		matrix[1][0] = 0;
		matrix[1][1] = 1;
		matrix[1][2] = 0;
		
		matrix[2][0] = 0;
		matrix[2][1] = 0;
		matrix[2][2] = 1;
	}
	
	// invert matrix for rotations from shower coordinates into ground coordinates
	if( bInvers )
	{
		float temp = 0.;
		temp = matrix[0][1];
		matrix[0][1] = matrix[1][0];
		matrix[1][0] = temp;
		temp = matrix[0][2];
		matrix[0][2] = matrix[2][0];
		matrix[2][0] = temp;
		temp = matrix[1][2];
		matrix[1][2] = matrix[2][1];
		matrix[2][1] = temp;
	}
	
}

/*
    "borrowed" from eventDisplay
    VGrIsuAnalyzer.cpp
      
    helper function for core reconstruction
*/

void mtxmlt( float a[3][3], float b[3], float c[3] )
{
	for( int i = 0; i < 3; i++ )
	{
		c[i] = 0.0;
		for( int j = 0; j < 3; j++ )
		{
			c[i] += a[i][j] * b[j];
		}
	}
}

/*!

"borrowed" from eventDisplay
VSimpleStereoReconstructor.cpp
    
helper function for shower core reconstrction

RETURN=    None
ARGUMENT=  prim   =Simulated primary characteristics in the original
                   ground system.
          \par xfield  the X ground locations of a telescope
      \par yfield  the y ground locations of a telescope
      \par zfield  the z ground locations of a telescope
      \par xtelrot the telescope X location in the rotated reference frame.
      \par ytelrot the telescope Y location in the rotated reference frame.
      \par ztelrot the telescope Z location in the rotated reference frame.
      \par bInv do inverse rotation from shower coordinates into ground coordinates
Function to calculate the coor. of the primary and telescope in
the rotated frame.WHAT IS THE ROTATED FRAME? DOES THE ANALYSIS WORK EVEN
IF THERE IS NO SIMULATION SPECIFIC RECORD?
*/
void tel_impact( float xcos, float ycos, float xfield, float yfield, float zfield, float* xtelrot, float* ytelrot, float* ztelrot, bool bInv )
{
	float b[3] = { 0., 0., 0. };
	float c[3] = { 0., 0., 0. };
	float matrix[3][3] = { { 0., 0., 0. },
		{ 0., 0., 0. },
		{ 0., 0., 0. }
	};
	
	float dl = 0.;
	float dm = 0.;
	float dn = 0.;                               /*Direction cos of the primary in the ground frame*/
	
	/* determine the rotation matrix from setup_matrix */
	dl = xcos;
	dm = ycos;
	if( 1. - dl * dl - dm * dm < 0. )
	{
		dn = 0.;
	}
	else
	{
		dn = -sqrt( 1. - dl * dl - dm * dm );
	}
	setup_matrix( matrix, dl, dm, dn, bInv );
	for( unsigned int i = 0; i < 3; i++ )
	{
		c[i] = 0.;
	}
	
	/* determine the location of the telescope in the rotated frame */
	
	b[0] = xfield;
	b[1] = yfield;
	b[2] = zfield;
	if( c[0] )
	{
		dl = 0.;
	}
	if( c[1] )
	{
		dl = 0.;
	}
	if( c[2] )
	{
		dl = 0.;
	}
	
	mtxmlt( matrix, b, c );
	
	if( c[0] )
	{
		dl = 0.;
	}
	if( c[1] )
	{
		dl = 0.;
	}
	if( c[2] )
	{
		dl = 0.;
	}
	
	// (GM) small number check
	for( unsigned int i = 0; i < 3; i++ ) if( TMath::Abs( c[i] ) < 1.e-5 )
		{
			c[i] = 0.;
		}
		
	*xtelrot = c[0];
	*ytelrot = c[1];
	*ztelrot = c[2];
}

/*
* "borrowed" from eventDisplay
* VSimpleStereoReconstructor.cpp
*
* calculate shower core in ground coordinates and
* check validity of core reconstruction results
*/
bool fillShowerCore( float fTelElevation,
    float fTelAzimuth,
    float ximp,
    float yimp )
{
    // check validity
    if(!isnormal( ximp ) || !isnormal( yimp ) )
    {
        fShower_Xcore = -99999.;
        fShower_Ycore = -99999.;
        return false;
    }
    // reconstructed shower core in ground coordinates
    float i_xcos = 0.;
    float i_ycos = 0.;
    float zimp = 0.;
    float igz = 0.;
    // calculate z in shower coordinates (for z=0 in ground coordinates)
    if( fShower_Ze != 0. )
    {
        zimp = yimp / tan(( 90. - fShower_Ze ) * TMath::DegToRad() );
    }
    // calculate direction cosinii
    // taking telescope plane as reference plane.
    i_xcos = sin(( 90. - fTelElevation ) * TMath::DegToRad() )
             * sin(( fTelAzimuth - 180. ) * TMath::DegToRad() );
    if( fabs( i_xcos ) < 1.e-7 )
    {
        i_xcos = 0.;
    }
    i_ycos = sin(( 90. - fTelElevation ) * TMath::DegToRad() )
             * cos(( fTelAzimuth - 180. ) * TMath::DegToRad() );
    if( fabs( i_ycos ) < 1.e-7 )
    {
        i_ycos = 0.;
    }
    tel_impact( i_xcos, i_ycos, ximp, yimp, zimp, &fShower_Xcore, &fShower_Ycore, &igz, true );
    if( isinf( fShower_Xcore ) || isinf( fShower_Ycore )
            || TMath::IsNaN( fShower_Xcore ) || TMath::IsNaN( fShower_Ycore ) )
    {
        fShower_Xcore = -99999;
        fShower_Ycore = -99999;
        return false;
    }
    return true;
}

/*
    "borrowed" from eventDisplay
    VSimpleStereoReconstructor.cpp

    reconstruction of shower core
    Hofmann et al 1999, Method 1 (HEGRA method)
    shower core by intersection of lines connecting reconstruced shower
    direction and image centroids
    expected to be run after direction reconstruction
    corresponds to rcs_method4 in VArrayAnalyzer
*/
bool reconstruct_core( unsigned int i_ntel,
		double iShowerDir_xs,
		double iShowerDir_ys,
        double fTelElevation,
        double fTelAzimuth,
        double* iTelX,
		double* iTelY,
		double* iTelZ,
		double* img_size,
		double* img_cen_x,
		double* img_cen_y,
		double* img_width,
		double* img_length)
{
    // sign flip in reconstruction
	iShowerDir_ys *= -1.;
	
	// make sure that all data arrays exist
	if( !img_size || !img_cen_x || !img_cen_y
			|| !img_width || !img_length)
	{
        //std::cout << "Missing data: cannot reconstruct event."<<std::endl;
		return false;
	}
	
	////////////////////////////////////////////////
	// core reconstruction
	////////////////////////////////////////////////
	
	// calculated telescope positions in shower coordinates
	float i_xcos = sin( ( 90. - fTelElevation ) / TMath::RadToDeg() ) * sin( ( fTelAzimuth - 180. ) / TMath::RadToDeg() );
	float i_ycos = sin( ( 90. - fTelElevation ) / TMath::RadToDeg() ) * cos( ( fTelAzimuth - 180. ) / TMath::RadToDeg() );
	float i_xrot, i_yrot, i_zrot = 0.;
	
	float ximp = 0.;
	float yimp = 0.;
	float stdp = 0.;
	
	float i_cenx = 0.;
	float i_ceny = 0.;
	
	std::vector< float > m;
	std::vector< float > x;
	std::vector< float > y;
	std::vector< float > w;
	float iweight = 1.;
	
	for( unsigned int i = 0; i < i_ntel; i++ )
	{
        // length == length and width == width protect against negative estimators of the variance - NK
		if( img_size[i] > 0. && img_length[i] > 0. && img_length[i] == img_length[i] && img_width[i] == img_width[i])
		{
			// telescope coordinates
			// shower coordinates (telecope pointing)
			tel_impact( i_xcos, i_ycos, iTelX[i], iTelY[i], iTelZ[i], &i_xrot, &i_yrot, &i_zrot, false );
			x.push_back( i_xrot - iShowerDir_xs / TMath::RadToDeg() * i_zrot );
			y.push_back( i_yrot - iShowerDir_ys / TMath::RadToDeg() * i_zrot );
			
			// gradient of image
			i_cenx = img_cen_x[i] - iShowerDir_xs;
			i_ceny = img_cen_y[i] - iShowerDir_ys;
			if( i_cenx != 0. )
			{
				m.push_back( -1. * i_ceny / i_cenx );
			}
			else
			{
				m.push_back( 1.e9 );
			}
			// image weight
			iweight = img_size[i];
			iweight *= ( 1. - img_width[i] / img_length[i] );
			w.push_back( iweight * iweight );
		}
	}

    // check minimum angle between image lines; ignore if too small
    // Note difference to evndisp reconstruction: apply this here for 2-tel events only
    //
    double fAxesAngles_min = 0.; // NK - probably only matters if we ever implement disp. Setting to 0 for now.
    if( m.size() == 2 )
    {
        float iangdiff = fabs( atan( m[0] ) - atan( m[1] ) ) * TMath::RadToDeg();
        if( iangdiff < fAxesAngles_min || TMath::Abs( 180. - iangdiff ) < fAxesAngles_min )
        {
            fShower_Xcore = -99999.;
            fShower_Ycore = -99999.;
            fShower_Chi2 = 1.;
            return false;
        }
    }

	// Now call perpendicular_distance for the fit, returning ximp and yimp
	rcs_perpendicular_fit( x, y, w, m, ( int )w.size(), &ximp, &yimp, &stdp );
    
	// return to ground coordinates
    if( fillShowerCore( fTelElevation, fTelAzimuth, ximp, yimp ) )
    {
        fShower_Chi2 = 0.;
    }
    else
    {
        fShower_Chi2 = -1.;
    }
    fShower_stdP = stdp;
    return true;
}


/*
* Create an image in a single telescope for a given event number
* coordinate transformations done using GrOptics method GUtilityFuncts::sourceOnTelescopePlane
*/
TH2I* telEvent(int telNumber, int eventNumber){

    if(!f){
        std::cout<< "error reading file, try readFile(\"rootfile.root\")" <<std::endl;
        return nullptr;
    }
    if(!t){
        std::cout<<"error reading tree"<<std::endl;
        return nullptr;
    }

    // REL TO GrOptics COORDINATES
    
    // draw tree
    auto condition = Form("(telID==%d && eventNumber==%d )", telNumber, eventNumber);
    t->Draw("CX:CY:az:ze",condition,"goff");

    // read data from tree
    const int NCp = t->GetSelectedRows();
    auto cx = t->GetV1();
    auto cy = t->GetV2();
    
    auto prmAz = TMath::DegToRad()*t->GetV3()[0];
    prmAz = redang(M_PI - redang(prmAz - M_PI)); // CORSIKA to GrOptics
    auto prmZe = TMath::DegToRad()*t->GetV4()[0];
    // double telAz = prmAz;
    double telAz = prmAz + fMag_Dec*TMath::DegToRad(); // unrotate array to correct for magnetic declination, ARRANG
    double telZe = prmZe;

    // double telAz = fTel_Azimuth;
    // double telZe = fTel_Zenith;

    // sourceOnTelescopePlane
    double epsilon = numeric_limits<double>::epsilon();
    double xcos_t = sin(telZe)*sin(telAz);
    double ycos_t = sin(telZe)*cos(telAz);
    if (TMath::AreEqualAbs(xcos_t,0.0,epsilon)){xcos_t = 0.0;}
    if (TMath::AreEqualAbs(ycos_t,0.0,epsilon)){ycos_t = 0.0;}
    double zcos_t = sqrt(1-xcos_t*xcos_t-ycos_t*ycos_t);
    ROOT::Math::XYZVector nUnit_t(xcos_t,ycos_t,zcos_t);

    // and now rotation matrix to get to telescope coordinates
    ROOT::Math::Rotation3D rotM;
    ROOT::Math::RotationZ rz(telAz);
    ROOT::Math::RotationX rx(telZe);
    rotM = rx*rz;
    
    // fill image
    TH2I* image = new TH2I(Form("T%d",telNumber), Form("T%d",telNumber), 32, -4.95, 4.95, 32, -4.95, 4.95);
    for(int i=0; i<NCp; i++){
        double xcos_s = -1*cy[i]; // CORSIKA to GrOptics
        double ycos_s = cx[i]; // CORSIKA to GrOptics
        if (TMath::AreEqualAbs(xcos_s,0.0,epsilon)){xcos_s = 0.0;}
        if (TMath::AreEqualAbs(ycos_s,0.0,epsilon)){ycos_s = 0.0;}
        double zcos_s = sqrt(1-xcos_s*xcos_s-ycos_s*ycos_s);
        ROOT::Math::XYZVector nUnit_s(xcos_s,ycos_s,zcos_s);
        // some vector algebra to find the vector to the intersection point
        double dotP = nUnit_s.Dot(nUnit_t);
        ROOT::Math::XYZVector tTos_vecGC = (nUnit_s/dotP) - nUnit_t;
        ROOT::Math::XYZVector tTos_vecTC = rotM*tTos_vecGC;
        if (TMath::AreEqualAbs(tTos_vecTC.X(),0.0,numeric_limits<double>::epsilon())) {tTos_vecTC.SetX(0.0);}
        if (TMath::AreEqualAbs(tTos_vecTC.Y(),0.0,numeric_limits<double>::epsilon())) {tTos_vecTC.SetY(0.0);}
        if (TMath::AreEqualAbs(tTos_vecTC.Z(),0.0,numeric_limits<double>::epsilon())) {tTos_vecTC.SetZ(0.0);}
        double imgX = TMath::RadToDeg()*tTos_vecTC.X();
        double imgY = TMath::RadToDeg()*tTos_vecTC.Y();

        // scatter by PSF
        std::tuple<double,double> coords = spread(imgX,imgY);
        double x = std::get<0>(coords);
        double y = std::get<1>(coords);

        // sign flip telescope coordinates to camera coordinates: y-> -1*y
            // via GrOptics README:
                // "For historical reasons, the camera coordinate system's y axis
                //  is a reflection of the y-axis of the telescope coordinate"
    
        // okay but in GrOptics this turns into x-> -1*x, so is the reflection ON the y-axis, and not that the y coordinate is reflected? i.e. x-> -1*x???
        // https://github.com/groptics/GrOptics/blob/8fccd40fb8141f420c7e395626a623fd3baf7555/src/GArrayTel.cpp#L192

        image->Fill(-1*x,y,petoadu); 
    }

    image = addNSB(image);
    image = addElectronics(image);

    // trigger threshold
    if(image->GetMaximum()<6.5*petoadu){
        image->Reset();
    }else{
        image = clean(image);
        image->Draw("COLZ");
    }

    image->ResetStats();
    image->SetStats(0);

    image->GetXaxis()->SetLabelSize(0);
    image->GetYaxis()->SetLabelSize(0);
    image->GetXaxis()->SetTickLength(0);
    image->GetYaxis()->SetTickLength(0);
    
    return image;
    
}

/*
* get total signal in a pixel over all events in a single telescope
* check if cleaning is enabled and if pedestals are subtracted before running
*/
void paramPixel(){
    // check a file is loaded before trying to read data
    if(!f){
        std::cout << "No file loaded" << std::endl;
        return;
    }

    // openfile
    std::ofstream datafile;
    datafile.open("simpixel.csv", std::ios_base::app);

    // make all images in one telescope
    int N = t->GetEntries();

    const int tel = 1;
    for(int eventNumber=1; eventNumber<=N+1; eventNumber++){
        TH2I* image = telEvent(tel, eventNumber);
        int signal = image->GetBinContent(16,16); //central pixel
        // make sure image isnt empty
        if(image->GetSumOfWeights()!=0){
            datafile << signal << std::endl;
        }
        image->Delete();
        
    }
    datafile.close();
    // std::cout << "Parameterization completed " << std::endl;
}


/*
* Writes parameter distributions for each shower in a data file to CSV for making histograms like in Fegan 1997
*/
void paramCSV(bool reconstruct=false){

    // check a file is loaded before trying to read data
    if(!f){
        std::cout << "No file loaded" << std::endl;
        return;
    }

    // openfile
    std::ofstream datafile;
    std::string output = f->GetName();
    output = output.substr(0,output.size()-5)+".10xSig.csv";
    datafile.open(output);

    if(!reconstruct){
        datafile << "Event,Telescope,MeanX,StdX,MeanY,StdY,Phi,Size,Length,Width,Miss,Distance,Azwidth,Alpha,TrueAz,TrueZe,TrueXcore,TrueYcore,TrueEnergy" << std::endl;
    }else{
        datafile << "Event,Telescope,MeanX,StdX,MeanY,StdY,Phi,Size,Length,Width,Miss,Distance,Azwidth,Alpha,Az,Ze,Xcore,Ycore,stdP,TrueAz,TrueZe,TrueXcore,TrueYcore,TrueEnergy" << std::endl;
    }

    // make images and paramaterize every event in each telescope
    int N = t->GetEntries();
    // find event numbers
    t->Draw("eventNumber","","goff");
    int start = (int) t->GetV1()[0];
    int stop = (int) t->GetV1()[N-1];
    for(int eventNumber=start; eventNumber<=stop; eventNumber++){
        std::cout << "Parameterizing event "<< eventNumber << std::endl;

        double* meanx = new double[Ntel];
        double* stdx = new double[Ntel];
        double* meany = new double[Ntel];
        double* stdy = new double[Ntel];
        double* phi = new double[Ntel];
        double* phi_rad = new double[Ntel];
        double* size = new double[Ntel];
        double* length = new double[Ntel];
        double* width = new double[Ntel];
        double* miss = new double[Ntel];
        double* dist = new double[Ntel];
        double* azwidth = new double[Ntel];
        double* alpha = new double[Ntel];

        double* TelX = new double[Ntel];
        double* TelY = new double[Ntel];
        double* TelZ = new double[Ntel];

        for(int i=0; i<Ntel; i++){
            TH2I* image = telEvent(i+1, eventNumber);
            auto params = parameterize(image);
            image->Delete();

            meanx[i] = std::get<0>(params);
            stdx[i] = std::get<1>(params);
            meany[i] = std::get<2>(params);
            stdy[i] = std::get<3>(params);
            phi[i] = std::get<4>(params);
            phi_rad[i] = std::get<4>(params)*TMath::DegToRad();
            size[i] = std::get<5>(params);
            length[i] = std::get<6>(params);
            width[i] = std::get<7>(params);
            miss[i] = std::get<8>(params);
            dist[i] = std::get<9>(params);
            azwidth[i] = std::get<10>(params);
            alpha[i] = std::get<11>(params);

            // get telescope positions
            if(!f){
                std::cout<< "error reading file, try readFile(\"rootfile.root\")" <<std::endl;
            }
            if(!t){
                std::cout<<"error reading tree"<<std::endl;
            }

            // read event data
            auto condition = Form("eventNumber==%d", eventNumber);
            t->Draw("telXpos:telYpos:telZpos",condition,"goff");

            // convert from CORSIKA to GrOptics
            TelX[i]=-1*t->GetV2()[i];
            TelY[i]=t->GetV1()[i];
            TelZ[i]=t->GetV3()[i];
      
        }

        auto condition = Form("eventNumber==%d", eventNumber);
        t->Draw("energy:az:ze",condition,"goff");
        double energy = t->GetV1()[0];
        double az = t->GetV2()[0];
        // CORSIKA to GrOptics
        az=TMath::RadToDeg()*redang(M_PI - redang(TMath::DegToRad()*az - M_PI));
        az=az + fMag_Dec; // unrotate array to correct for magnetic declination, ARRANG
        double ze = t->GetV3()[0];
        // double az = fTel_Azimuth;
        // double ze = fTel_Zenith;
        
        t->Draw("xCore:yCore",condition,"goff");
        // CORSIKA to GrOptics
        double xCore = -1*t->GetV2()[0];
        double yCore = t->GetV1()[0];

        if(!reconstruct){
            // write data to file
            for(int i = 0; i<Ntel; i++){
                datafile << eventNumber << "," << i+1 << "," << meanx[i] << "," << stdx[i] << "," << meany[i] << "," << stdy[i] << "," << phi[i] <<","<< size[i] << "," << length[i] << "," << width[i] << "," << miss[i] 
                    << "," << dist[i] << "," << azwidth[i] << "," << alpha[i] << "," << az << "," << ze << "," << xCore 
                    << "," << yCore << "," << energy << std::endl;   
            }
        }else{
            // reconstruction
            if(reconstruct_direction(Ntel,90-ze,az,size,meanx,meany,phi_rad,length,width)){

                if(reconstruct_core(Ntel, fShower_Xoffset, fShower_Yoffset, 90-ze, az, TelX, TelY, TelZ, size, meanx, meany, width, length)){

                    // write data to file
                    for(int i = 0; i<Ntel; i++){
                        datafile << eventNumber << "," << i+1 << "," << meanx[i] << "," << stdx[i] << "," << meany[i] << "," << stdy[i] << "," << phi[i] <<","<< size[i] << "," << length[i] << "," << width[i] << "," << miss[i] 
                            << "," << dist[i] << "," << azwidth[i] << "," << alpha[i] << "," << fShower_Az << "," 
                            << fShower_Ze << "," << fShower_Xcore << "," << fShower_Ycore << "," << fShower_stdP << "," 
                            << az << "," << ze << "," << xCore << "," << yCore << "," << energy << std::endl;   
                    }
                }
            }else{
                // write data to file
                for(int i = 0; i<Ntel; i++){
                        datafile << eventNumber << "," << i+1 << "," << meanx[i] << "," << stdx[i] << "," << meany[i] << "," << stdy[i] << "," << phi[i] <<","<< size[i] << "," << length[i] << "," << width[i] << "," << miss[i] 
                            << "," << dist[i] << "," << azwidth[i] << "," << alpha[i] << "," << "nan" << "," 
                            << "nan" << "," << "nan" << "," << "nan" << "," << "nan" << "," 
                            << az << "," << ze << "," << xCore << "," << yCore << "," << energy << std::endl;   
                    }
            }
        }

    }

    datafile.close();
    std::cout << "Parameterization completed " << std::endl;
}

/*
* Show the process of drawing Cherenkov photons, simulating noise and cleaning the image, then paramaterizing
*/
void showClean(int telNumber, int eventNumber){

    if(!f){
        std::cout<< "error reading file, try readFile(\"rootfile.root\")" <<std::endl;
        return;
    }
    if(!t){
        std::cout<<"error reading tree"<<std::endl;
        return;
    }

    // plot image 
    TCanvas *c = new TCanvas("show cleaning","show cleaning",1600,360);
    gStyle->SetPalette(57); // reset to default palette (kBird)
	c->Divide(4,1,0.01,0.01);

    //
    // Cherenkov photons
    //

    c->cd(1);

    // draw tree
    auto condition = Form("(telID==%d && eventNumber==%d )", telNumber, eventNumber);
    t->Draw("CX:CY",condition,"goff");

    // read data from tree
    const int NCp = t->GetSelectedRows();
    auto imgX = t->GetV1();
    auto imgY = t->GetV2();

    TH2I* image = new TH2I("Only Cherenkov Photons", "Only Cherenkov Photons", 32, -4.95, 4.95, 32, -4.95, 4.95 );
    
    // fill image
    for(int i=0; i<NCp; i++){
        // scatter by PSF
        std::tuple<double,double> coords = spread(imgX[i],imgY[i]);
        double x = std::get<0>(coords);
        double y = std::get<1>(coords);

        // sign flip ground->sky
        image->Fill(x,-1*y,petoadu);
        
    }

    image->ResetStats();
    image->SetStats(0);

    image->GetXaxis()->SetLabelSize(0);
    image->GetYaxis()->SetLabelSize(0);
    image->GetXaxis()->SetTickLength(0);
    image->GetYaxis()->SetTickLength(0);

    image->DrawCopy("COLZ1","");

    //
    // Add noise
    //

    c->cd(2);
    image = addNSB(image);
    image = addElectronics(image);
    image->SetTitle("NSB Added");
    image->DrawCopy("COLZ1","");    

    //
    // Cleaning
    //

    c->cd(3);

    // trigger threshold
    if(image->GetMaximum()<6.5*petoadu){
        image->Reset();
    }else{
        image = clean(image);
    }

    image->SetTitle("After Aperture Cleaning");
    image->DrawCopy("COLZ1","");

    //
    // Parameterize
    //

    c->cd(4);
    image->SetTitle("Parameterized");
    image->DrawCopy("COLZ1","");

    auto params = parameterize(image);
    image->Delete();

    TEllipse *e = new TEllipse(std::get<0>(params), std::get<2>(params), std::get<6>(params), std::get<7>(params), 0, 360, std::get<4>(params));
    e->SetFillStyle(0);
    e->SetLineWidth(4);
    e->Draw("SAME");

}

/*
* Display timing gradient of images in all telescopes for a given event number, and plot
* arrival time as a function of Cherenkov photon position in camera
*/
void timegrad(int eventNumber){

    if(!f){
        std::cout<< "error reading file, try readFile(\"rootfile.root\")" <<std::endl;
    }
    if(!t){
        std::cout<<"error reading tree"<<std::endl;
    }

    // plot image from each telescope
    TCanvas *it = new TCanvas("Image Timing","Image Timing",800,720);
    it->Divide(ceil(sqrt(Ntel)),ceil(sqrt(Ntel)),0.01,0.01);

    gStyle->SetPalette(75); //kCherry

    // plot arrival times for each telescope
    TCanvas *h = new TCanvas("Arrival Times","Arrival Times",800,720);
    h->Divide(ceil(sqrt(Ntel)),ceil(sqrt(Ntel)),0.01,0.01);

    for(int i=0; i<Ntel; i++){
        gPad->SetTopMargin(0.1);
        gPad->SetBottomMargin(0.01);
        gPad->SetLeftMargin(0.01);

        // draw tree
        auto condition = Form("(telID==%d && eventNumber==%d )", i+1, eventNumber);
        t->Draw("CX:CY:CTime",condition,"goff");

        // read data from tree
        const int NCp = t->GetSelectedRows();
        auto imgX = t->GetV1();
        auto imgY = t->GetV2();
        auto imgT = t->GetV3();

        TH2D* time_grad = new TH2D(Form("T%d",i+1), Form("T%d",i+1), 32, -4.95, 4.95, 32, -4.95, 4.95 );
        TH2I* tmp = (TH2I*)time_grad->Clone();
        TH1D* arrival = new TH1D(Form("T%d",i+1), Form("T%d",i+1), 500, 0, 50 );

        // fill image
        for(int i=0; i<NCp; i++){
            // scatter by PSF
            std::tuple<double,double> coords = spread(imgX[i],imgY[i]);
            double x = std::get<0>(coords);
            double y = std::get<1>(coords);

            // sign flip ground->sky
            tmp->Fill(x,-1*y,petoadu);
            if(imgT[i]<1e-9){ // prevent empty bins
                time_grad->Fill(x,-1*y,1e-9);
            }else{
                time_grad->Fill(x,-1*y,imgT[i]);
            }
            
            arrival->Fill(imgT[i]);
        }

        // average arrival time in each pixel
        int bin=0;
        for(int i=1; i<=32; i++){
            for(int j=1; j<=32; j++){
                bin = time_grad->GetBin(i,j);
                if(tmp->GetBinContent(bin)){ //dont divide by zero
                    time_grad->SetBinContent(bin,time_grad->GetBinContent(bin)/tmp->GetBinContent(bin));
                }
            }
        }

        // apply cleaning 
        tmp = addNSB(tmp);
        tmp = addElectronics(tmp);

        // trigger threshold
        if(tmp->GetMaximum()<6.5*petoadu){
            tmp->Reset();
        }else{
            tmp = clean(tmp);
            //tmp->Draw("COLZ");
        }

        //figure out which pixels to remove
        for(int i=1; i<=32; i++){
            for(int j=1; j<=32; j++){
                bin = time_grad->GetBin(i,j);
                if(tmp->GetBinContent(bin)==0){ //dont divide by zero
                    time_grad->SetBinContent(bin,0);
                }
            }
        }

        arrival->GetXaxis()->SetTitle("Cherenkov photon arrival time (ns)");
        arrival->GetYaxis()->SetTitle("Cherenkov photon count");
        arrival->ResetStats();
        arrival->SetStats(0);

        time_grad->ResetStats();
        time_grad->SetStats(0);
        time_grad->GetZaxis()->SetTitle("average Cherenkov photon arrival time (ns)");
        int min = arrival->GetBinLowEdge(arrival->FindFirstBinAbove(0));
        //int max = 1 + arrival->GetBinLowEdge(arrival->FindLastBinAbove(0));
        //int max = time_grad->GetBinContent(arrival->GetBinLowEdge(arrival->FindLastBinAbove(0));
        //time_grad->GetZaxis()->SetRangeUser(min,max);
        time_grad->SetMinimum(min);

        time_grad->GetXaxis()->SetLabelSize(0);
        time_grad->GetYaxis()->SetLabelSize(0);
        time_grad->GetXaxis()->SetTickLength(0);
        time_grad->GetYaxis()->SetTickLength(0);

        tmp->Delete();

        it->cd(i+1);
        gPad->SetRightMargin(0.2);
        time_grad->DrawCopy("COLZ1","");
        h->cd(i+1);
        gPad->SetRightMargin(0.1);
        arrival->DrawCopy();

        time_grad->Delete();
        arrival->Delete();
    }
}
/*
* Display images from all telescope for a given event number, and the position of the
* shower core relative to the telescopes
*/
void panodisplay(int eventNumber){

    // debug
    //TStopwatch t;
    //t.Start();

    // plot image from each telescope
    TCanvas *c = new TCanvas("Array Event","Array Event",800,720);
    gStyle->SetPalette(57); // reset to default palette (kBird)
	c->Divide(ceil(sqrt(Ntel)),ceil(sqrt(Ntel)),0.01,0.01);
    
    double* meanx = new double[Ntel];
    double* stdx = new double[Ntel];
    double* meany = new double[Ntel];
    double* stdy = new double[Ntel];
    double* phi = new double[Ntel];
    double* phi_rad = new double[Ntel];
    double* size = new double[Ntel];
    double* length = new double[Ntel];
    double* width = new double[Ntel];
    double* miss = new double[Ntel];
    double* dist = new double[Ntel];
    double* alpha = new double[Ntel];

    double* TelX = new double[Ntel];
    double* TelY = new double[Ntel];
    double* TelZ = new double[Ntel];

    for(int i=0; i<Ntel; i++){
        c->cd(i+1);
        gPad->SetTopMargin(0.1);
        gPad->SetBottomMargin(0.01);
        gPad->SetLeftMargin(0.01);
        gPad->SetRightMargin(0.15);

        TH2I* image = telEvent(i+1, eventNumber);
        image->DrawCopy("COLZ1","");
        // parameterization
        auto params = parameterize(image);
        image->Delete();

        TEllipse *e = new TEllipse(std::get<0>(params), std::get<2>(params), std::get<6>(params), std::get<7>(params), 0, 360, std::get<4>(params));
	    e->SetFillStyle(0);
	    e->SetLineWidth(4);
	    e->Draw("SAME");

        meanx[i]=std::get<0>(params);
        stdx[i]=std::get<1>(params);
        meany[i]=std::get<2>(params);
        stdy[i]=std::get<3>(params);
        phi[i]=std::get<4>(params);
        phi_rad[i]=std::get<4>(params)*TMath::DegToRad();
        size[i]=std::get<5>(params);
        length[i]=std::get<6>(params);
        width[i]=std::get<7>(params);
        miss[i] = std::get<8>(params);
        dist[i] = std::get<9>(params);
        alpha[i]=std::get<11>(params);

        TString parameterInfo = Form(
        "=======================\n"
        "TELESCOPE:\t%d\n"
        "-----------------------\n"
        "MEAN-X:\t\t%f\n"
        "SIGMA-X:\t%f\n"
        "MEAN-Y:\t\t%f\n"
        "SIGMA-Y:\t%f\n"
        "PHI:\t\t%f\n"
        "SIZE:\t\t%f\n"
        "LENGTH:\t\t%f\n"
        "WIDTH:\t\t%f\n"
        "MISS:\t\t%f\n"
        "DIST:\t\t%f\n"
        "ALPHA:\t\t%f\n",
        i+1, meanx[i],stdx[i],meany[i],stdy[i],phi[i],size[i],length[i],width[i],miss[i],dist[i],alpha[i]);

        std::cout<<parameterInfo<<std::endl;

        // Visualize Alpha
        /*
        // image axis
        if( size[i] != 0. ){
            double m1 = ( sin(phi_rad[i]) / cos(phi_rad[i]) );
            TF1 *f = new TF1("f","[0]*(x-[1])+[2]",-5,5); 
            f->SetParameters(m1,meanx[i],meany[i]);
            f->SetLineColor(kBlack);
            f->SetLineWidth(2);
            f->Draw("SAME");

            // center to centroid
            double m2 = meany[i]/meanx[i];
            TF1 *c = new TF1("c","[0]*x",-5,5); 
            c->SetParameters(m2);
            c->SetLineColor(kRed);
            c->SetLineWidth(2);
            c->Draw("SAME");

        }
        */
        
        // get telescope positions
        if(!f){
            std::cout<< "error reading file, try readFile(\"rootfile.root\")" <<std::endl;
        }
        if(!t){
            std::cout<<"error reading tree"<<std::endl;
        }

        // read event data
        auto condition = Form("eventNumber==%d", eventNumber);
        t->Draw("telXpos:telYpos:telZpos",condition,"goff");

        // convert CORSIKA to GrOptics
        TelX[i]=-1*t->GetV2()[i];
        TelY[i]=t->GetV1()[i];
        TelZ[i]=t->GetV3()[i];
    }

    auto condition = Form("eventNumber==%d", eventNumber);
    t->Draw("az:ze",condition,"goff");
    double az = t->GetV1()[0];
    // CORSIKA to GrOptics
    az=TMath::RadToDeg()*redang(M_PI - redang(TMath::DegToRad()*az - M_PI));
    az=az + fMag_Dec; // unrotate array to correct for magnetic declination, ARRANG
    double ze = t->GetV2()[0];
    // double az = fTel_Azimuth;
    // double ze = fTel_Zenith;

    // showerInfo
    std::cout<< "Simulated Shower Params:" << std::endl << showerInfo(eventNumber) << std::endl ;

    // plot image of telescope and shower core positions
    TCanvas *m = new TCanvas("Event Map","Event Map",840,720);
    TMultiGraph *map = eventMap(eventNumber);

    gPad->SetTopMargin(0.1);
    gPad->SetBottomMargin(0.1);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);

    // reconstruction
    if(reconstruct_direction(Ntel,90-ze,az,size,meanx,meany,phi_rad,length,width)){
        std::cout<<"Reconstructed Direction: "<<fShower_Az <<", " << fShower_Ze <<std::endl;

        // draw telescopes
        if(reconstruct_core(Ntel, fShower_Xoffset, fShower_Yoffset, 90-ze ,az,TelX, TelY, TelZ, size, meanx, meany, width, length)){
            std::cout<<"Reconstructed Core: "<< fShower_Xcore << " m" << ", " << fShower_Ycore << " m" << ", +/- " << fShower_stdP << " m" << std::endl;

            // point
            // plot core
            TGraph* g = new TGraph(1);
            g->SetMarkerStyle(34);
            //g->SetMarkerSize(3);
            g->SetMarkerColor(kBlue);
            g->SetPoint(0,fShower_Xcore,fShower_Ycore);

            /*
            // label point
            TLatex *l1 = new TLatex(coreX,coreY-50,"reconstructed core");
            l1->SetTextSize(0.025);
            l1->SetTextFont(42);
            l1->SetTextAlign(21);
            TLatex *l2 = new TLatex(coreX,coreY-75,Form("(%.2f,%.2f)",coreX,coreY));
            l2->SetTextSize(0.025);
            l2->SetTextFont(42);
            l2->SetTextAlign(21);
            
            g->GetListOfFunctions()->Add(l1);
            g->GetListOfFunctions()->Add(l2);
            */
            g->Draw("AP SAME");

            map->Add(g);   
            
        }
    }

    map->Draw("AP");

    // debug
    //t.Stop();
    //t.Print();
}