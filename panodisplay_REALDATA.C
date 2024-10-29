/*
* panodisplay.C
* Displays images of simulated air showers taken by an array of PANOSETI telescopes 
* 
* This macro is a work of simulation. Any resemblance to analysis packages,
* living or dead, is purely coincidence.
*
* Author: Nik Korzoun
* Jamie is awesome
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
#include "string"

#include "Math/Vector3D.h"
#include "Math/Rotation3D.h"
#include "Math/EulerAngles.h"
#include "Math/RotationX.h"
#include "Math/RotationY.h"
#include "Math/RotationZ.h"
#include "Math/VectorUtil.h"
#include "TMatrixD.h"

using namespace ROOT::Math;

/* dmod(A,B) - A modulo B (double) */
// from VASlamac.h
#define dmod(A,B) ((B)!=0.0?((A)*(B)>0.0?(A)-(B)*floor((A)/(B))\
							 :(A)+(B)*floor(-(A)/(B))):(A))

// Random seed
int seed = 200;
TRandom3 *r = new TRandom3(seed);

// Root file
const char *prefix;
TFile *f;
TTree *t;
int Ntel=3; //hard-coded for now

Int_t           array_event_num;
Int_t           array_tel_event[3];
UShort_t        array_scope_id[3];
Short_t        array_pix_data[3][32][32];

// Reconstructed params
float fShower_Xoffset = -99999.;
float fShower_Yoffset = -99999.;
float fShower_Az = -99999.;
float fShower_Ze = -99999.;
float fShower_Xcore = -99999.;
float fShower_Ycore = -99999.;
float fShower_stdP = -99999.;

/*
* Read root file for displaying images
*/
void readFile(const char *infile_prefix){
    // load tree
    char array_infile_name[200];
    snprintf(array_infile_name,200,"%s.array",infile_prefix);
    f = new TFile(array_infile_name);
    t = (TTree*)f->Get("arraydata"); 
    prefix = infile_prefix;

    t->SetBranchAddress("array_event_num",&array_event_num);
    t->SetBranchAddress("array_tel_event",array_tel_event);
    t->SetBranchAddress("array_scope_id",array_scope_id);
    t->SetBranchAddress("array_pix_data",array_pix_data);
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
* Clean image according to p.e. thresholds
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

    double NSB = 1; // mean value of NSB per pixel - ADU
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
* Calculate offset based on approximated field drift velocity
* assumes drift velocity is linear with time
* times are obtained from file name
* times used in units of minutes
* arrays in function arguments are structured so elements are ordered like {before flip, after flip}
*/
//std::tuple<int,int> calcOffset(double pix_start[4], double pix_end[4], int time_start[2], int time_end[2]){
std::tuple<int,int> calcOffset(int initial_offset_x[2], int initial_offset_y[2], double drift_velocity_x[2], double drift_velocity_y[2]){

    // get time (approximate)
    // cout<<f->GetName()<<endl;
    int h = 60 * std::stoi( std::string(f->GetName()).substr(45,2));
    int m = std::stoi( std::string(f->GetName()).substr(47,2));
    int s = std::stoi( std::string(f->GetName()).substr(49,2)) / 60;
    int t = h+m+s;
    const int t0 = 310;

    // check meridian flip
    const int meridian = 441; // ~7:41 pm
    int x = 0;
    int y = 0;

    int x0 = initial_offset_x[0];
    int y0 = initial_offset_y[0];
    double vx = drift_velocity_x[0];
    double vy = drift_velocity_y[0];

    if(t>meridian){
        x0 = initial_offset_x[1];
        y0 = initial_offset_y[1];
        vx = drift_velocity_x[1];
        vy = drift_velocity_y[1];
    }

    x = x0+vx*(t-t0);
    y = y0+vx*(t-t0);

    return std::make_tuple(x,y);


    /*
    // get time (approximate)
    // cout<<f->GetName()<<endl;
    int h = 60 * std::stoi( std::string(f->GetName()).substr(45,2));
    int m = std::stoi( std::string(f->GetName()).substr(47,2));
    int s = std::stoi( std::string(f->GetName()).substr(49,2)) / 60;
    int t = h+m+s;

    // check meridian flip
    const int meridian = 441; // ~7:41 pm

    double vx = 0;
    double vy = 0;
    int x = 0;
    int y = 0;
    // calc velocity
    if(t<meridian){
        vx = (pix_end[0]-pix_start[0])/(time_end[0]-time_start[0]);
        vy = (pix_end[1]-pix_start[1])/(time_end[0]-time_start[0]);
        x = ((t-time_start[0])*vx)+pix_start[0];
        y = ((t-time_start[0])*vy)+pix_start[1];
    }else{
        vx = (pix_end[2]-pix_start[2])/(time_end[1]-time_start[1]);
        vy = (pix_end[3]-pix_start[3])/(time_end[1]-time_start[1]);
        x = ((t-time_start[1])*vx)+pix_start[2];
        y = ((t-time_start[1])*vy)+pix_start[3];
    }

    return std::make_tuple(x,y);
    */
}

/*
* Shift each pixel in an image by x,y
*/
TH2I* shift(TH2I* image, int x, int y){

    // clone image
    TH2I *newImage = (TH2I*)image->Clone();
    newImage->Reset();

    // loop over all bins
    int bins = image->GetNbinsX();
	for(int i=1; i<=bins; i++){
		for(int j=1; j<=bins; j++){
            if(i+x<bins && i-x>0 && j+y<bins && j-y>0){
                newImage->SetBinContent(i-x,j-y,image->GetBinContent(i,j));
            }
        }
    }
    image->Delete();
    return newImage;
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

    TMultiGraph* map = new TMultiGraph();
    TGraph* telescopes = new TGraph(Ntel);
    TGraph* shower = new TGraph(1);
    map->SetTitle("Event Map");

    // fill telescope positions
    TLatex *l1;
    TLatex *l2;
    telescopes->SetMarkerStyle(20);
    telescopes->SetMarkerSize(3);

    // Hard-wired for now, probably better to read in a .cfg file in the long term
    double* TelX = new double[Ntel]{80.48,-6.95,-73.52}; // K,D,C
    double* TelY = new double[Ntel]{53.59,-107.18,53.59};
    for(int i=0;i<Ntel;i++){
        double x=TelX[i];
        double y=TelY[i];
        telescopes->SetPoint(i,x,y);

        TString label="";
        switch(i){
            case 0:
                label="Kron";
                break;
            case 1:
                label="Dorm";
                break;
            case 2:
                label="Crocker";
                break;
        }

        // label point
        l1 = new TLatex(x-30*cos(atan2(y,x)),y-30*sin(atan2(y,x)),label);
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
	// telescope pointings
    // assume telescope pointing directly upwards and North
	//double fTelElevation = 60.; 
	//double fTelAzimuth   = 90.;
	
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
        double fTelElevation,
        double fTelAzimuth,
		double iShowerDir_xs,
		double iShowerDir_ys,
        double* iTelX,
		double* iTelY,
		double* iTelZ,
		double* img_size,
		double* img_cen_x,
		double* img_cen_y,
		double* img_width,
		double* img_length)
{
    // telescope pointings
    // assume telescope pointing directly upwards and North
	//double fTelElevation = 60.; 
	//double fTelAzimuth   = 90.;

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
	// Now call perpendicular_distance for the fit, returning ximp and yimp
	rcs_perpendicular_fit( x, y, w, m, ( int )w.size(), &ximp, &yimp, &stdp );
    fShower_Xcore = ximp;
    fShower_Ycore = yimp;
    fShower_stdP = stdp;
	
	// return to ground coordinates
    // NK - uh oh. idk what that means. Probably should have been using VSimpleStereoReconstructor::fillShowerCore

    return true;
}


/*
* Create an image in a single telescope for a given event number
* coordinate transformations done using GrOptics method GUtilityFuncts::sourceOnTelescopePlane
*/
// TH2I* telEvent(int telNumber, int eventNumber, int a, int b){
TH2I* telEvent(int telNumber, int eventNumber){
    if(!f){
        std::cout<< "error reading file, try readFile(\"rootfile.root\")" <<std::endl;
        return nullptr;
    }
    if(!t){
        std::cout<<"error reading tree"<<std::endl;
        return nullptr;
    }

    TString label="";
    switch(telNumber){
        case 1:
            label="Kron";
            break;
        case 2:
            label="Dorm";
            break;
        case 3:
            label="Crocker";
            break;
    }

    // from Jamie's arraydisplay.C
    t->GetEntry(eventNumber);
    TH2I* image = new TH2I(label,label, 32, -4.95, 4.95, 32, -4.95, 4.95);
    
    char pedvar_infile_name[200];
    TH2D *peds_2D_hist;
    TH2D *pedvars_2D_hist;
    if (array_scope_id[telNumber-1]>0) 
    {
        snprintf(pedvar_infile_name,200,"%s.T%d.pedvars",prefix,array_scope_id[telNumber-1]);
        TFile *pedvar_infile = TFile::Open(pedvar_infile_name, "read");
        // TFile *pedvar_infile=new TFile(pedvar_infile_name);
        peds_2D_hist=(TH2D*)pedvar_infile->Get("peds_2D_hist");
        pedvars_2D_hist=(TH2D*)pedvar_infile->Get("pedvars_2D_hist");
        peds_2D_hist->SetDirectory(nullptr);
        pedvars_2D_hist->SetDirectory(nullptr);
        delete pedvar_infile;
    }
    image->SetStats(0);
    for (int i=0;i<32;i++)
    {
        for (int j=0;j<32;j++)
        {
            double pixval=0;
            if (array_scope_id[telNumber-1]>0)
            {
                pixval=array_pix_data[telNumber-1][i][j]-peds_2D_hist->GetBinContent(i+1,j+1);
                //pixval=array_pix_data[telNumber][i][j];
                //pixval=(array_pix_data[telNumber][i][j]-peds_2D_hist->GetBinContent(i+1,j+1))/pedvars_2D_hist->GetBinContent(i+1,j+1);
                pixval=pixval;
                image->SetBinContent(i+1,j+1,pixval);
            }
        }
    }
    //image[jtel]->Draw("COLZ");

    image = clean(image);
    int x = 0;
    int y = 0;
    if(telNumber==1){// kron
        // no shift
    }else if(telNumber==2){ // dorm
        int x_initial[2] = {-3,9};
        int y_initial[2] = {-2,2};
        double vx_initial[2] = {0,0};
        double vy_initial[2] = {0,0};
        std::tuple<int,int> offset = calcOffset(x_initial,y_initial,vx_initial,vy_initial);
        x = std::get<0>(offset);
        y = std::get<1>(offset);
        // image = shift(image,x,y);
    }else if(telNumber==3){ //crocker
        // int x_initial[2] = {a,0};
        // int y_initial[2] = {b,0};
        // double vx_initial[2] = {0,0};
        // double vy_initial[2] = {0,0};
        int x_initial[2] = {-13,-16};
        int y_initial[2] = {-18,9};
        double vx_initial[2] = {-0.069,0.069};
        // double vx_initial[2] = {0,0};
        double vy_initial[2] = {0,-0.046};
        std::tuple<int,int> offset = calcOffset(x_initial,y_initial,vx_initial,vy_initial);
        x = std::get<0>(offset);
        y = std::get<1>(offset);
        // image = shift(image,x,y);
    }else{
        //error
        cout<<"telNumber "<<telNumber<< " is unknown!";
    }
    image->Draw("COLZ");
    

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
    datafile.open("pixels.csv", std::ios_base::app);

    // // make all images in one telescope
    // int N = t->GetEntries();

    // std::cout << "Processing file: "<< prefix << std::endl;
    // const int tel = 2; // Dorm
    // for(int eventNumber=1; eventNumber<=N+1; eventNumber++){
    //     TH2I* image = telEvent(tel, eventNumber);
    //     int signal = image->GetBinContent(20,5); //central pixel
    //     // make sure image isnt empty
    //     if(image->GetSumOfWeights()!=0){
    //         datafile << signal << std::endl;
    //     }
    //     image->Delete();
        
    // }
    // datafile.close();

    // make all images in one telescope
    int N = t->GetEntries();

    std::cout << "Processing file: "<< prefix << std::endl;
    const int tel = 2; // Dorm
    for(int eventNumber=1; eventNumber<=N+1; eventNumber++){
        TH2I* image = telEvent(tel, eventNumber);
        // make sure image isnt empty
        if(image->GetSumOfWeights()!=0){
            int binsX=image->GetNbinsX();
            int binsY=binsX;
            // loop over pixels
            for(int i = 1; i<=binsX; i++){
                for(int j = 1; j<=binsY; j++){
                    float signal = image->GetBinContent(i,j);
                    datafile<<signal;
                    if(i!=binsX || j!=binsY){
                        datafile<<',';
                    }else{
                        datafile<<std::endl;
                    }
                }
            }
        }
        
        image->Delete();
        
    }
    datafile.close();

    // std::cout << "Parameterization completed " << std::endl;
    std::cout << "Completed parameterizing file "<< prefix << std::endl;
}

/*
* writes information to csv needed to calculate circumcircle of three shower images
* tests shifting a,b
*/
// void paramCircumcircle(int a, int b){
void paramCircumcircle(){
    // check a file is loaded before trying to read data
    if(!f){
        std::cout << "No file loaded" << std::endl;
        return;
    }

    // openfile
    std::ofstream datafile;
    std::string output = prefix;
    datafile.open(output + ".circumcircle.all_corrections.csv");

    // std::ofstream datafile;
    // datafile.open("circumcircle.csv", std::ios_base::app); //append


    // make images and paramaterize every event in each telescope
    int N = t->GetEntries();

    std::cout << "Parameterizing file "<< prefix << std::endl;
    for(int eventNumber=1; eventNumber<=N+1; eventNumber++){
        // std::cout << "Parameterizing event "<< eventNumber << std::endl;

        double* meanx = new double[Ntel];
        double* meany = new double[Ntel];
        double* phi_rad = new double[Ntel];

        // make sure there are three images
        bool valid = true;

        for(int i=0; i<Ntel; i++){
            //TH2I* image = telEvent(i+1, eventNumber,a,b);
            TH2I* image = telEvent(i+1, eventNumber);
            auto params = parameterize(image);
            image->Delete();

            meanx[i] = std::get<0>(params);
            meany[i] = std::get<2>(params);
            phi_rad[i] = std::get<4>(params)*TMath::DegToRad();
            
            if(std::get<5>(params) < 100){ //size
                valid = false;
            }
        }
        // check there are three images
        if(valid){
            // datafile << eventNumber << ',' << a << ',' << b;
            datafile << eventNumber << ',' << "-13" << ',' << "-18";
            for(int i=0;i<Ntel;i++){
                datafile << ',' << meanx[i] << ',' << meany[i] << ',' << phi_rad[i];
            }
            datafile << std::endl;
        }
        
    }
    datafile.close();
    // std::cout << "Parameterization completed " << std::endl;
    std::cout << "Completed parameterizing file "<< prefix << std::endl;
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
    std::string output = prefix;
    datafile.open(output + ".csv");

    if(!reconstruct){
        datafile << "Event,Telescope,MeanX,StdX,MeanY,StdY,Phi,Size,Length,Width,Miss,Distance,Azwidth,Alpha,TrueAz,TrueZe,TrueXcore,TrueYcore,TrueEnergy" << std::endl;
    }else{
        datafile << "Event,Telescope,MeanX,StdX,MeanY,StdY,Phi,Size,Length,Width,Miss,Distance,Azwidth,Alpha,Az,Ze,Xcore,Ycore,stdP,TrueAz,TrueZe,TrueXcore,TrueYcore,TrueEnergy" << std::endl;
    }

    // make images and paramaterize every event in each telescope
    int N = t->GetEntries();

    for(int eventNumber=1; eventNumber<=N+1; eventNumber++){
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

        // Hard-wired for now, probably better to read in a .cfg file in the long term
        double* TelX = new double[Ntel]{80.48,-6.95,-73.52}; // K,D,C
        double* TelY = new double[Ntel]{53.59,-107.18,53.59};
        double* TelZ = new double[Ntel]{0.,0.,0.};

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
      
        }

        //auto condition = Form("eventNumber==%d", eventNumber);
        //t->Draw("energy:az:ze",condition,"goff");
        //double energy = t->GetV1()[0];
        //double az = t->GetV2()[0];
        //az = redang(M_PI - redang(az - M_PI)); // convert from CORSIKA
        //double ze = t->GetV3()[0];
        
        //t->Draw("xCore:yCore",condition,"goff");
        //double xCore = -1*t->GetV2()[0]; // convert from CORSIKA
        //double yCore = t->GetV1()[0]; // convert from CORSIKA

        if(!reconstruct){
            // write data to file
            for(int i = 0; i<Ntel; i++){
                datafile << eventNumber << "," << i+1 << "," << meanx[i] << "," << stdx[i] << "," << meany[i] << "," << stdy[i] << "," << phi[i] <<"," << size[i] << "," << length[i] << "," << width[i] << "," << miss[i] 
                    << "," << dist[i] << "," << azwidth[i] << "," << alpha[i] /*<< "," << az << "," << ze << "," << xCore 
                    << "," << yCore << "," << energy */<< std::endl;   
            }
        }else{
            // // reconstruction
            // if(reconstruct_direction(Ntel, 90-ze, az, size,meanx,meany,phi_rad,length,width)){

            //     if(reconstruct_core(Ntel, 90-ze, az, fShower_Xoffset, fShower_Yoffset, TelX, TelY, TelZ, size, meanx, meany, width, length)){

            //         // write data to file
            //         for(int i = 0; i<Ntel; i++){
            //             datafile << eventNumber << "," << i+1 << "," << meanx[i] << "," << stdx[i] << "," << meany[i] << "," stdy[i] << "," << phi[i] << "," << size[i] << "," << length[i] << "," << width[i] << "," << miss[i] 
            //                 << "," << dist[i] << "," << azwidth[i] << "," << alpha[i] << "," << fShower_Az << "," 
            //                 << fShower_Ze << "," << fShower_Xcore << "," << fShower_Ycore << "," << fShower_stdP /*<< "," 
            //                 << az << "," << ze << "," << xCore << "," << yCore << "," << energy */<< std::endl;   
            //         }
            //     }
            // }else{
            //     // write data to file
            //     for(int i = 0; i<Ntel; i++){
            //             datafile << eventNumber << "," << i+1 << "," << meanx[i] << "," << stdx[i] << "," << meany[i] << "," stdy[i] << "," << phi[i] << "," << size[i] << "," << length[i] << "," << width[i] << "," << miss[i] 
            //                 << "," << dist[i] << "," << azwidth[i] << "," << alpha[i] << "," << "nan" << "," 
            //                 << "nan" << "," << "nan" << "," << "nan" << "," << "nan" /*<< "," 
            //                 << az << "," << ze << "," << xCore << "," << yCore << "," << energy */<< std::endl;   
            //         }
            // }
        }

    }

    datafile.close();
    std::cout << "Parameterization completed " << std::endl;
}

/*
* Display images from all telescope for a given event number on a single image, and the position of the
* shower core relative to the telescopes
*/
void arraydisplay(int eventNumber){

    // debug
    //TStopwatch t;
    //t.Start();

    // plot image from each telescope
    gStyle->SetCanvasPreferGL(kTRUE);
    TCanvas *c = new TCanvas("Array Event","Array Event",800,720);
    gStyle->SetPalette(57, 0, 0.5); // reset to default palette (kBird), but lower opacity
    
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

    // Hard-wired for now, probably better to read in a .cfg file in the long term
    double* TelX = new double[Ntel]{80.48,-6.95,-73.52}; // K,D,C
    double* TelY = new double[Ntel]{53.59,-107.18,53.59};
    double* TelZ = new double[Ntel]{0.,0.,0.};


    int colors[3] = {kBlue+2,kCyan-7,kYellow-7};
    TEllipse *ellipses[3]={0};
    TF1 *axes[3]={0};
    TH2I* combined = new TH2I("combined","combined", 32, -4.95, 4.95, 32, -4.95, 4.95);
    for(int i=0; i<Ntel; i++){
        gPad->SetTopMargin(0.1);
        gPad->SetBottomMargin(0.1);
        gPad->SetLeftMargin(0.1);
        gPad->SetRightMargin(0.15);

        TH2I* image = telEvent(i+1, eventNumber);
        int bins = image->GetNbinsX();
        for(int j=1; j<=bins; j++){
            for(int k=1; k<=bins; k++){
                combined->AddBinContent(image->GetBin(j,k),image->GetBinContent(j,k));
            }
        }
        // parameterization
        auto params = parameterize(image);
        image->Delete();

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
        
        // image axis and ellipses
        if( size[i] != 0. ){
            TEllipse *e = new TEllipse(std::get<0>(params), std::get<2>(params), std::get<6>(params), std::get<7>(params), 0, 360, std::get<4>(params));
            e->SetFillStyle(0);
            e->SetLineWidth(4);
            e->SetLineColor(colors[i]);
            ellipses[i]=e;

            double m1 = ( sin(phi_rad[i]) / cos(phi_rad[i]) );
            TF1 *f = new TF1("f","[0]*(x-[1])+[2]",-5,5); 
            f->SetParameters(m1,meanx[i],meany[i]);
            f->SetLineColor(colors[i]);
            f->SetLineWidth(2);
            axes[i]=f;

            // center to centroid
            // double m2 = meany[i]/meanx[i];
            // TF1 *c = new TF1("c","[0]*x",-5,5); 
            // c->SetParameters(m2);
            // c->SetLineColor(kRed);
            // c->SetLineWidth(2);
            // c->Draw("SAME");

        }
        
        
    }
    combined->ResetStats();
    combined->SetStats(0);
    combined->Draw("COLZ1 SAME");

    // legend and parameterization
    auto legend = new TLegend(0.1,0.7,0.48,0.9);
    for(int i=0; i<Ntel; i++){
        TString label="";
        switch(i){
            case 0:
                label="Kron";
                break;
            case 1:
                label="Dorm";
                break;
            case 2:
                label="Crocker";
                break;
        }
        if(ellipses[i]){
            ellipses[i]->Draw("SAME");
            TText *t = new TText(ellipses[i]->GetX1(),ellipses[i]->GetY1(),label);
            t->SetTextAngle(ellipses[i]->GetTheta());
            t->SetTextColor(colors[i]);
            t->SetTextFont(42);
            t->SetTextSize(0.05);
            t->SetTextAlign(21);
            t->Draw("SAME");
        }
        if(axes[i]){
            axes[i]->Draw("SAME");
        }
    }

    // plot image of telescope and shower core positions
    TCanvas *m = new TCanvas("Event Map","Event Map",840,720);
    TMultiGraph *map = eventMap(eventNumber);

    gPad->SetTopMargin(0.1);
    gPad->SetBottomMargin(0.1);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);

    // get pointings
    //az = 
    //ze =

    /*
    // reconstruction
    if(reconstruct_direction(Ntel, 90-ze, az, size,meanx,meany,phi_rad,length,width)){
        std::cout<<"Reconstructed Direction: "<<fShower_Az <<", " << fShower_Ze <<std::endl;

        // draw telescopes
        if(reconstruct_core(Ntel, 90-ze, az, fShower_Xoffset, fShower_Yoffset, TelX, TelY, TelZ, size, meanx, meany, width, length)){
            std::cout<<"Reconstructed Core: "<< fShower_Xcore << " m" << ", " << fShower_Ycore << " m" << ", +/- " << fShower_stdP << " m" << std::endl;

            // point
            // plot core
            TGraph* g = new TGraph(1);
            g->SetMarkerStyle(34);
            //g->SetMarkerSize(3);
            g->SetMarkerColor(kBlue);
            g->SetPoint(0,fShower_Xcore,fShower_Ycore);

            g->Draw("AP SAME");

            map->Add(g);   
            
        }
    }
    */
    map->Draw("AP");

    // debug
    //t.Stop();
    //t.Print();
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
    TCanvas *c = new TCanvas("Array Event","Array Event",(int)1500/1.2,(int)500/1.2);//800,720);
    gStyle->SetPalette(57); // reset to default palette (kBird)
	//c->Divide(ceil(sqrt(Ntel)),ceil(sqrt(Ntel)),0.01,0.01);
    c->Divide(3,1,0.01,0.01);
    
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

    // Hard-wired for now, probably better to read in a .cfg file in the long term
    double* TelX = new double[Ntel]{80.48,-6.95,-73.52}; // K,D,C
    double* TelY = new double[Ntel]{53.59,-107.18,53.59};
    double* TelZ = new double[Ntel]{0.,0.,0.};

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
        
    }

    // plot image of telescope and shower core positions
    TCanvas *m = new TCanvas("Event Map","Event Map",840,720);
    TMultiGraph *map = eventMap(eventNumber);

    gPad->SetTopMargin(0.1);
    gPad->SetBottomMargin(0.1);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);

    // get pointings
    //az = 
    //ze =

    /*
    // reconstruction
    if(reconstruct_direction(Ntel, 90-ze, az, size,meanx,meany,phi_rad,length,width)){
        std::cout<<"Reconstructed Direction: "<<fShower_Az <<", " << fShower_Ze <<std::endl;

        // draw telescopes
        if(reconstruct_core(Ntel, 90-ze, az, fShower_Xoffset, fShower_Yoffset, TelX, TelY, TelZ, size, meanx, meany, width, length)){
            std::cout<<"Reconstructed Core: "<< fShower_Xcore << " m" << ", " << fShower_Ycore << " m" << ", +/- " << fShower_stdP << " m" << std::endl;

            // point
            // plot core
            TGraph* g = new TGraph(1);
            g->SetMarkerStyle(34);
            //g->SetMarkerSize(3);
            g->SetMarkerColor(kBlue);
            g->SetPoint(0,fShower_Xcore,fShower_Ycore);

            g->Draw("AP SAME");

            map->Add(g);   
            
        }
    }
    */
    map->Draw("AP");

    // debug
    //t.Stop();
    //t.Print();
}

