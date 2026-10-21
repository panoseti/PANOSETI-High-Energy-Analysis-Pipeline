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
#include "Math/Vector3D.h"
#include "Math/Rotation3D.h"
#include "Math/RotationX.h"
#include "Math/RotationZ.h"

#include <limits>
#include "iostream"
#include "fstream"

/* dmod(A,B) - A modulo B (double) */
// from VASlamac.h
#define dmod(A,B) ((B)!=0.0?((A)*(B)>0.0?(A)-(B)*floor((A)/(B))\
							 :(A)+(B)*floor(-(A)/(B))):(A))

// Random seed
int seed = 200;
TRandom3 *r = new TRandom3(seed);

// Set to true to use ray tracing with thin Fresnel lens model, false to use Gaussian approximation
const bool kUseRayTrace = true;

// Root file
TFile *f;
TTree *t;
int Ntel;

// Fixed telescope pointing (GrOptics definition)
// float fTel_Zenith = 20;
// float fTel_Azimuth = 180;

// Magnetic Declination in Degrees (unrotate CORSIKA i.e. ARRANG)
// float fMag_Dec = 12.77; // LICK
float fMag_Dec = 11.03;// PALOMAR
// float fMag_Dec = 10.4; // FLWO
// float fMag_Dec = 0;

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

// ============================================================================
// PANOSETI Thin Lens Optical Model Parameters & Tables
// ============================================================================

const double kOpticsF = 60.78;         // Focal length (cm)
const double kOpticsD = 46.09;         // Aperture diameter (cm)
const double kOpticsR = 23.045;        // Aperture radius (cm)
const double kOpticsRoughness = 0.0;   // Lens surface roughness (cm), 0 = disabled
const double kSigmaScattering = kOpticsRoughness/kOpticsF;

// Fresnel polynomial sag coefficients: y = sum(a_k * rho^(2k))
const double kPolyCoeffs[11] = {
    0.00000000000000000e+00, // a0
    1.70651353139298460e-02, // a1
    -3.53352632884415423e-06, // a2
    9.57163484153624398e-10, // a3
    3.55711029991929572e-13, // a4
    -1.12361085705273497e-15, // a5
    1.96912674836582728e-19, // a6
    9.79010455872597239e-22, // a7
    3.32575748290361466e-25, // a8
    -6.27878102919874693e-28, // a9
    -7.74304624855183072e-31  // a10
};

// Inverse CDF of photon energy (eV) for Cherenkov spectrum folded with full PDE
// (Atmospheric transmission at alt=10km, zn=30deg * PMMA transmission * SiPM PDE).
// Zero-tail trimmed to physical detection window [1.305 eV, 3.647 eV] (~340 to 950 nm, 201 points).
const int kNumUGrid = 201;
const double kInvCDFEnergy[201] = {
    1.305097, 1.405120, 1.449981, 1.486680, 1.519204,
    1.547630, 1.573380, 1.596913, 1.619006, 1.642698,
    1.661821, 1.680094, 1.697955, 1.715681, 1.732686,
    1.748250, 1.763344, 1.778084, 1.792426, 1.806646,
    1.819888, 1.832719, 1.845235, 1.857458, 1.869423,
    1.881172, 1.892775, 1.904193, 1.915422, 1.926377,
    1.937024, 1.947466, 1.957723, 1.967757, 1.977610,
    1.987181, 1.996568, 2.005831, 2.014978, 2.024016,
    2.032950, 2.041785, 2.050527, 2.059185, 2.067772,
    2.076315, 2.084811, 2.093248, 2.101612, 2.109908,
    2.117991, 2.125959, 2.133850, 2.141674, 2.149439,
    2.157156, 2.164818, 2.172418, 2.179967, 2.187447,
    2.194839, 2.202173, 2.209452, 2.216678, 2.223854,
    2.230979, 2.238057, 2.245086, 2.252068, 2.259010,
    2.265911, 2.272772, 2.279595, 2.286380, 2.293128,
    2.299846, 2.306532, 2.313189, 2.319814, 2.326411,
    2.332977, 2.339517, 2.346026, 2.352510, 2.358970,
    2.365408, 2.371824, 2.378220, 2.384596, 2.390950,
    2.397286, 2.403601, 2.409896, 2.416173, 2.422430,
    2.428669, 2.434889, 2.441093, 2.447280, 2.453458,
    2.459628, 2.465787, 2.471937, 2.478078, 2.484211,
    2.490338, 2.496462, 2.502582, 2.508699, 2.514813,
    2.520924, 2.527031, 2.533136, 2.539237, 2.545336,
    2.551432, 2.557525, 2.563617, 2.569710, 2.575803,
    2.581896, 2.587990, 2.594085, 2.600179, 2.606274,
    2.612370, 2.618466, 2.624564, 2.630663, 2.636764,
    2.642866, 2.648969, 2.655074, 2.661182, 2.667294,
    2.673412, 2.679535, 2.685663, 2.691795, 2.697935,
    2.704084, 2.710244, 2.716415, 2.722597, 2.728793,
    2.735005, 2.741233, 2.747478, 2.753740, 2.760019,
    2.766316, 2.772630, 2.778961, 2.785310, 2.791677,
    2.798063, 2.804472, 2.810905, 2.817363, 2.823845,
    2.830353, 2.836886, 2.843450, 2.850046, 2.856676,
    2.863339, 2.870040, 2.876786, 2.883577, 2.890412,
    2.897305, 2.904265, 2.911298, 2.918407, 2.925599,
    2.932890, 2.940294, 2.947809, 2.955475, 2.963291,
    2.971274, 2.979451, 2.987828, 2.996470, 3.005357,
    3.014575, 3.024096, 3.034036, 3.044433, 3.055327,
    3.066807, 3.078999, 3.092001, 3.106080, 3.121486,
    3.138554, 3.158061, 3.181192, 3.210602, 3.254194,
    3.646594
};

// Refractive index table for PMMA vs photon energy (1.0 to 6.0 eV, 51 points, step 0.1 eV)
const double kRefractiveIndexEMin = 1.00;
const double kRefractiveIndexEMax = 6.00;
const int kNumRefractiveIndexGrid = 51;
const double kRefractiveIndexGrid[51] = {
    1.466383, 1.467147, 1.467994, 1.468914, 1.469911, 1.470988, 1.472144, 1.473380,
    1.474700, 1.476102, 1.477590, 1.479164, 1.480827, 1.482580, 1.484426, 1.486365,
    1.488401, 1.490535, 1.492770, 1.495111, 1.497558, 1.500117, 1.502792, 1.505586,
    1.508467, 1.511482, 1.514626, 1.517903, 1.521315, 1.524868, 1.528565, 1.532412,
    1.536416, 1.540581, 1.544915, 1.549424, 1.554115, 1.558997, 1.564076, 1.569360,
    1.574865, 1.580594, 1.586560, 1.592774, 1.599251, 1.606005, 1.613040, 1.620392,
    1.628056, 1.636061, 1.644431
};

// Sample photon energy from precomputed inverse CDF
double sample_photon_energy() {
    double u = r->Rndm();
    double idx_d = u * (kNumUGrid - 1);
    int idx = (int)idx_d;
    if (idx >= kNumUGrid - 1) return kInvCDFEnergy[kNumUGrid - 1];
    if (idx < 0) return kInvCDFEnergy[0];
    double frac = idx_d - idx;
    return kInvCDFEnergy[idx] * (1.0 - frac) + kInvCDFEnergy[idx + 1] * frac;
}

// Get PMMA refractive index for a given energy (eV)
double get_refractive_index(double energy_eV) {
    if (energy_eV <= kRefractiveIndexEMin) return kRefractiveIndexGrid[0];
    if (energy_eV >= kRefractiveIndexEMax) return kRefractiveIndexGrid[kNumRefractiveIndexGrid - 1];
    double frac_idx = (energy_eV - kRefractiveIndexEMin) / (kRefractiveIndexEMax - kRefractiveIndexEMin) * (kNumRefractiveIndexGrid - 1);
    int idx = (int)frac_idx;
    double frac = frac_idx - idx;
    return kRefractiveIndexGrid[idx] * (1.0 - frac) + kRefractiveIndexGrid[idx + 1] * frac;
}

/*
 * Trace a photon through the PANOSETI Fresnel thin-lens optical model:
 * - Samples photon energy according to Cherenkov spectrum * atmosphere * PMMA * SiPM PDE
 * - Evaluates PMMA dispersion n(E)
 * - Uniformly samples entrance pupil impact position on circular aperture (diameter D = 46.09 cm)
 * - Refracts into lens at entry plane y = 0
 * - Refracts out of lens with normal defined by aspheric polynomial surface y = P(rho^2)
 * - Applies surface micro-roughness scattering if kOpticsRoughness > 0
 * - Propagates to focal plane at y = -F (F = 60.78 cm)
 *
 * Args:
 *   imgX_deg, imgY_deg: Incident photon direction in telescope frame (in degrees)
 * Returns:
 *   std::tuple<double, double>: Focal plane position converted back to angular degrees on sky
 */
std::tuple<double, double> sample_and_trace_photon(double imgX_deg, double imgY_deg) {
    const double nan = std::numeric_limits<double>::quiet_NaN();

    // Direction cosines in lens frame (optical axis is -Y, X is horizontal, Z is vertical)
    double tx = TMath::DegToRad() * imgX_deg;
    double tz = TMath::DegToRad() * imgY_deg;
    double dx = tx;
    double dz = tz;
    double dy = -1.0;
    double norm = sqrt(dx*dx + dy*dy + dz*dz);
    dx /= norm; dy /= norm; dz /= norm;

    // Sample photon energy and corresponding PMMA refractive index
    double energy = sample_photon_energy();
    double n_pmma = get_refractive_index(energy);

    // Uniformly sample entrance pupil disk (radius = D/2)
    double u1 = r->Rndm();
    double u2 = r->Rndm();
    double rad = kOpticsR * sqrt(u1);
    double phi = 2.0 * M_PI * u2;
    double x0 = rad * cos(phi);
    double z0 = rad * sin(phi);

    // 1. Refract into lens at y=0 (normal = (0, 1, 0))
    double cos_i = -dy; // since normal is (0, 1, 0)
    double n_ratio1 = 1.0 / n_pmma;
    double sin2_t1 = n_ratio1 * n_ratio1 * (1.0 - cos_i * cos_i);
    if (sin2_t1 > 1.0) {
        return std::make_tuple(nan, nan); // Total internal reflection fallback
    }
    double cos_t1 = sqrt(1.0 - sin2_t1);
    double lx = n_ratio1 * dx;
    double ly = n_ratio1 * dy + (n_ratio1 * cos_i - cos_t1);
    double lz = n_ratio1 * dz;
    double lnorm = sqrt(lx*lx + ly*ly + lz*lz);
    lx /= lnorm; ly /= lnorm; lz /= lnorm;

    // 2. Refract out of lens at aspheric polynomial surface y = P(rho^2)
    double rho2 = x0*x0 + z0*z0;
    double dp_du = 0.0;
    double rho2_pow = 1.0;
    for (int k = 1; k < 11; k++) {
        dp_du += k * kPolyCoeffs[k] * rho2_pow;
        rho2_pow *= rho2;
    }

    // Surface outward normal: N = (-2*x0*dp_du, 1, -2*z0*dp_du)
    double nx = -2.0 * dp_du * x0;
    double ny = 1.0;
    double nz = -2.0 * dp_du * z0;
    double nnorm = sqrt(nx*nx + ny*ny + nz*nz);
    nx /= nnorm; ny /= nnorm; nz /= nnorm;

    double cos_i2 = -(lx*nx + ly*ny + lz*nz);
    if (cos_i2 < 0) {
        nx = -nx; ny = -ny; nz = -nz;
        cos_i2 = -cos_i2;
    }
    double n_ratio2 = n_pmma;
    double sin2_t2 = n_ratio2 * n_ratio2 * (1.0 - cos_i2 * cos_i2);
    if (sin2_t2 > 1.0) {
        return std::make_tuple(nan, nan); // Total internal reflection fallback
    }
    double cos_t2 = sqrt(1.0 - sin2_t2);
    double ox = n_ratio2 * lx + (n_ratio2 * cos_i2 - cos_t2) * nx;
    double oy = n_ratio2 * ly + (n_ratio2 * cos_i2 - cos_t2) * ny;
    double oz = n_ratio2 * lz + (n_ratio2 * cos_i2 - cos_t2) * nz;
    double onorm = sqrt(ox*ox + oy*oy + oz*oz);
    ox /= onorm; oy /= onorm; oz /= onorm;

    // Surface micro-roughness scattering (Gaussian angular deviation)
    if (kSigmaScattering > 0.0) {
        double theta_scat = kSigmaScattering * sqrt(-2.0 * log(r->Rndm() + 1e-12));
        double phi_scat = 2.0 * M_PI * r->Rndm();

        // Scattered direction in a frame whose y axis is the outgoing ray
        double sin_th = sin(theta_scat);
        double wx = sin_th * cos(phi_scat), wy = cos(theta_scat), wz = sin_th * sin(phi_scat);

        // Rotate that frame's y axis onto (ox, oy, oz)
        double st = sqrt(ox*ox + oz*oz);
        if (st > 1e-12) {
            double cp = ox / st, sp = oz / st;
            double wr =  wx*cp + wz*sp;
            double wt = -wx*sp + wz*cp;
            double nx = wy*ox + wr*oy*cp - wt*sp;
            double ny = wy*oy - wr*st;
            double nz = wy*oz + wr*oy*sp + wt*cp;
            ox = nx; oy = ny; oz = nz;
        } else {                       // ray exactly along -y
            ox = wx; oy = -wy; oz = wz;
        }
    }

    // 3. Propagate to focal plane at y = -F
    if (oy >= 0.0) {
        return std::make_tuple(nan, nan); // Propagating away from focal plane
    }
    double s = -kOpticsF / oy;
    double x_fp = x0 + s * ox;
    double z_fp = z0 + s * oz;

    // Convert focal plane displacement (cm) to angular degrees on sky
    double out_x_deg = TMath::RadToDeg() * (x_fp / kOpticsF);
    double out_y_deg = TMath::RadToDeg() * (z_fp / kOpticsF);

    return std::make_tuple(out_x_deg, out_y_deg);
}

/*
 * Randomly spread arrival direction of photons to simulate PANOSETI PSF
*/
std::tuple<double,double> spread_gaussian(double positionX, double positionY){

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

// Retain spread() as an alias calling sample_and_trace_photon() or spread_gaussian()
std::tuple<double, double> spread(double positionX, double positionY) {
    if (kUseRayTrace) {
        return sample_and_trace_photon(positionX, positionY);
    } else {
        return spread_gaussian(positionX, positionY);
    }   
}

/*
void testspread(double x, double y){
    // histogram
    //TH2D *test = new TH2D("test", "test", 32, -4.95, 4.95, 32, -4.95, 4.95 );
    TH2D *test = new TH2D("test", "test", 600, -4.95, 4.95, 600, -4.95, 4.95 );
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
* Add noise consistent with typical pedvar to each pixel
*/
TH2D* addNoise(TH2D* image){

    int Nbins = image->GetNcells();
    TH2D *newImage = (TH2D*)image->Clone();
    for(int i=1; i<= Nbins ;i++){
        //newImage->AddBinContent(i, (int) r->Gaus(0,2.4));// electronics
        newImage->AddBinContent(i, (double)r->Gaus(0,10.0)); // mean of 8, but simulate pedestal subtraction to reduce to 0
	}
    image->Delete();
    return newImage;
}


/*
* Count the number of pixels with signal (non-zero content) in a TH2D image
*/
int countSignalPixels(TH2D* image) {
    int count = 0;
    for(int i = 1; i <= image->GetNbinsX(); i++) {
        for(int j = 1; j <= image->GetNbinsY(); j++) {
            if(image->GetBinContent(i, j) > 0) {
                count++;
            }
        }
    }
    return count;
}

/*
* Clean image according to image/border thresholds
*/
TH2D* clean(TH2D* image){
    // threshold cleaning for image pixels

    double imageThreshold = 4;
    double borderThreshold = 2; 

    double pedvarSize = 10.; // same as addNoise

    int Nbins = image->GetNcells();
    TH2D *newImage = (TH2D*)image->Clone();
    int binsX = newImage->GetNbinsX();
    int binsY = binsX;

    // remove pixels which are below image threshold unless they neighbor a pixel above image threshold and are themselves above border threshold
    std::vector<int> removeMe;
	for(int i=1; i<=binsX; i++){
		for(int j=1; j<=binsY; j++){
            int checkBin = newImage->GetBin(i,j);
            double binSize = newImage->GetBinContent(checkBin); // (pixdata-pedestal)
            double nsig = binSize/pedvarSize;

            bool remove = true;
            // check if pixel is above image threshold
            if(nsig>=imageThreshold){
                remove = false;
            // check if pixel is above border threshold
            }else if(nsig>=borderThreshold){
                // check if a neighbor is above image threshold
                // get neighbors
                for (int p=i-1; p<=i+1; p++){
                    for (int q=j-1; q<=j+1; q++){
                        // do not add central pixel as neighbor
                        if(!(p==i && q==j)){ 
                            // stay in bounds of image
                            if(p>=1 && p<=binsX && q>=1 && q<=binsY){
                                double neighbor = newImage->GetBinContent(newImage->GetBin(p,q));
                                double neighborSig = neighbor/pedvarSize;
                                // check if pixel borders a pixel above image threshold)
                                if (neighborSig >= imageThreshold){
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
    
    // check for isolated pixels and small islands with border pixels
    removeMe.clear();
    for(int i=1; i<=binsX; i++){
		for(int j=1; j<=binsY; j++){
            int checkBin = newImage->GetBin(i,j);
            double binSize = newImage->GetBinContent(checkBin);
            double nsig = binSize/pedvarSize;
            // make sure pixel has p.e. before checking to remove
            if(binSize!=0){
                int neighborCount = 0;
                int neighborBin = -1;
                // count neighbors
                for (int p=i-1; p<=i+1; p++){
                    for (int q=j-1; q<=j+1; q++){
                        // do not add central pixel as neighbor
                        if(!(p==i && q==j)){ 
                            // stay in bounds of image
                            if(p>=1 && p<=binsX && q>=1 && q<=binsY){
                                // find a neighbor with pixels in it
                                if(newImage->GetBinContent(newImage->GetBin(p,q)) != 0){
                                    neighborCount++;
                                    neighborBin = newImage->GetBin(p,q); // only removed if neighborCount == 1
                                }
                            }
                        }
                    }    
                }
                
                // remove isolated pixels
                if(neighborCount == 0){
                    removeMe.push_back(checkBin);
                }
                // remove 2-pixel islands if either pixel is below threshold
                else if(neighborCount == 1){

                    bool isBorderPixel = (nsig < imageThreshold);
                    if(isBorderPixel){
                        removeMe.push_back(checkBin);
                        removeMe.push_back(neighborBin);
                    }
                }
            }
        }
    }
    // remove pixels which are isolated or small islands with border pixels
    for(int i=0; i<(int)removeMe.size(); i++){
        newImage->SetBinContent(removeMe[i], 0);
    }

    image->Delete();
    return newImage;
}
/*
* Attempt image parameterization
* returns tuple which stores
* meanx, sigmax, meany, sigmay, angle, size, length, width
*/

std::tuple<double, double, double, double, double, double, double, double, double, double, double, double> parameterize(TH2D* image){
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
    az=az + fMag_Dec; // unrotate array to correct for magnetic declination, ARRANG
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
TH2D* telEvent(int telNumber, int eventNumber){

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
    double telAz = prmAz;
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

    // unrotate cx,cy from CORSIKA/ARRANG frame to true frame, magnetic declination
    double magDecRad = -fMag_Dec*TMath::DegToRad();

    // fill image
    TH2D* image = new TH2D(Form("T%d",telNumber), Form("T%d",telNumber), 32, -4.95, 4.95, 32, -4.95, 4.95);
    for(int i=0; i<NCp; i++){
        double xcos_s = -1*cy[i]*cos(magDecRad) + cx[i]*sin(magDecRad); // CORSIKA to GrOptics
        double ycos_s = cx[i]*cos(magDecRad) + cy[i]*sin(magDecRad); // CORSIKA to GrOptics
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

        // Apply PSF: raytraced Fresnel optical model or Gaussian approximation
        std::tuple<double,double> coords = spread(imgX,imgY);
        double x = std::get<0>(coords);
        double y = std::get<1>(coords);
        if (std::isnan(x) || std::isnan(y)) continue;

        // sign flip telescope coordinates to camera coordinates: y-> -1*y
            // via GrOptics README:
                // "For historical reasons, the camera coordinate system's y axis
                //  is a reflection of the y-axis of the telescope coordinate"
    
        // okay but in GrOptics this turns into x-> -1*x, so is the reflection ON the y-axis, and not that the y coordinate is reflected? i.e. x-> -1*x???
        // https://github.com/groptics/GrOptics/blob/8fccd40fb8141f420c7e395626a623fd3baf7555/src/GArrayTel.cpp#L192

        image->Fill(-1*x,y,petoadu); 
    }

    image = addNoise(image);

    // 2-pixel trigger in quadrants
    int npixtrigq1=0, npixtrigq2=0, npixtrigq3=0, npixtrigq4=0;
    for(int i=1; i<=16; i++){
        for(int j=1; j<=16; j++){
            if(image->GetBinContent(i,j)>6.5*petoadu) npixtrigq1++;
        }
    }
    for(int i=17; i<=32; i++){
        for(int j=1; j<=16; j++){
            if(image->GetBinContent(i,j)>6.5*petoadu) npixtrigq2++;
        }
    }
    for(int i=1; i<=16; i++){
        for(int j=17; j<=32; j++){
            if(image->GetBinContent(i,j)>6.5*petoadu) npixtrigq3++;
        }
    }
    for(int i=17; i<=32; i++){
        for(int j=17; j<=32; j++){
            if(image->GetBinContent(i,j)>6.5*petoadu) npixtrigq4++;
        }
    }

    bool trig=false;
    if(npixtrigq1>=2 || npixtrigq2>=2 || npixtrigq3>=2 || npixtrigq4>=2)
        trig=true;

    if(trig==false){
        image->Reset();
    }
    else{
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
    output = output.substr(0,output.size()-5)+".csv";
    datafile.open(output);

    if(!reconstruct){
        datafile << "Event,Telescope,MeanX,StdX,MeanY,StdY,Phi,Size,Npix,Length,Width,Miss,Distance,Azwidth,Alpha,TrueAz,TrueZe,TrueXcore,TrueYcore,TrueEnergy" << std::endl;
    }else{
        datafile << "Event,Telescope,MeanX,StdX,MeanY,StdY,Phi,Size,Npix,Length,Width,Miss,Distance,Azwidth,Alpha,Az,Ze,Xcore,Ycore,stdP,TrueAz,TrueZe,TrueXcore,TrueYcore,TrueEnergy" << std::endl;
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
        int* npix = new int[Ntel];
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
            TH2D* image = telEvent(i+1, eventNumber);
            auto params = parameterize(image);
            npix[i] = countSignalPixels(image);
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
                datafile << eventNumber << "," << i+1 << "," << meanx[i] << "," << stdx[i] << "," << meany[i] << "," << stdy[i] << "," << phi[i] <<","<< size[i] << "," << npix[i] << "," << length[i] << "," << width[i] << "," << miss[i]
                    << "," << dist[i] << "," << azwidth[i] << "," << alpha[i] << "," << az << "," << ze << "," << xCore
                    << "," << yCore << "," << energy << std::endl;
            }
        }else{
            // reconstruction
            if(reconstruct_direction(Ntel,90-ze,az,size,meanx,meany,phi_rad,length,width)){

                if(reconstruct_core(Ntel, fShower_Xoffset, fShower_Yoffset, 90-ze, az, TelX, TelY, TelZ, size, meanx, meany, width, length)){

                    // write data to file
                    for(int i = 0; i<Ntel; i++){
                        datafile << eventNumber << "," << i+1 << "," << meanx[i] << "," << stdx[i] << "," << meany[i] << "," << stdy[i] << "," << phi[i] <<","<< size[i] << "," << npix[i] << "," << length[i] << "," << width[i] << "," << miss[i]
                            << "," << dist[i] << "," << azwidth[i] << "," << alpha[i] << "," << fShower_Az << ","
                            << fShower_Ze << "," << fShower_Xcore << "," << fShower_Ycore << "," << fShower_stdP << ","
                            << az << "," << ze << "," << xCore << "," << yCore << "," << energy << std::endl;
                    }
                }
            }else{
                // write data to file
                for(int i = 0; i<Ntel; i++){
                        datafile << eventNumber << "," << i+1 << "," << meanx[i] << "," << stdx[i] << "," << meany[i] << "," << stdy[i] << "," << phi[i] <<","<< size[i] << "," << npix[i] << "," << length[i] << "," << width[i] << "," << miss[i]
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
    // Cherenkov photons - USE SAME COORDINATE TRANSFORMATION AS telEvent()
    //

    c->cd(1);

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
    double telAz = prmAz;
    double telZe = prmZe;

    // sourceOnTelescopePlane
    double epsilon = numeric_limits<double>::epsilon();
    double xcos_t = sin(telZe)*sin(telAz);
    double ycos_t = sin(telZe)*cos(telAz);
    if (TMath::AreEqualAbs(xcos_t,0.0,epsilon)){xcos_t = 0.0;}
    if (TMath::AreEqualAbs(ycos_t,0.0,epsilon)){ycos_t = 0.0;}
    double zcos_t = sqrt(1-xcos_t*xcos_t-ycos_t*ycos_t);
    ROOT::Math::XYZVector nUnit_t(xcos_t,ycos_t,zcos_t);

    // rotation matrix
    ROOT::Math::Rotation3D rotM;
    ROOT::Math::RotationZ rz(telAz);
    ROOT::Math::RotationX rx(telZe);
    rotM = rx*rz;

    // unrotate cx,cy from CORSIKA/ARRANG frame to true frame, magnetic declination
    double magDecRad = -fMag_Dec*TMath::DegToRad();

    TH2D* image = new TH2D("Only Cherenkov Photons", "Only Cherenkov Photons", 32, -4.95, 4.95, 32, -4.95, 4.95);

    // fill image
    for(int i=0; i<NCp; i++){
        double xcos_s = -1*cy[i]*cos(magDecRad) + cx[i]*sin(magDecRad); // CORSIKA to GrOptics
        double ycos_s = cx[i]*cos(magDecRad) + cy[i]*sin(magDecRad); // CORSIKA to GrOptics
        if (TMath::AreEqualAbs(xcos_s,0.0,epsilon)){xcos_s = 0.0;}
        if (TMath::AreEqualAbs(ycos_s,0.0,epsilon)){ycos_s = 0.0;}
        double zcos_s = sqrt(1-xcos_s*xcos_s-ycos_s*ycos_s);
        ROOT::Math::XYZVector nUnit_s(xcos_s,ycos_s,zcos_s);
        double dotP = nUnit_s.Dot(nUnit_t);
        ROOT::Math::XYZVector tTos_vecGC = (nUnit_s/dotP) - nUnit_t;
        ROOT::Math::XYZVector tTos_vecTC = rotM*tTos_vecGC;
        if (TMath::AreEqualAbs(tTos_vecTC.X(),0.0,numeric_limits<double>::epsilon())) {tTos_vecTC.SetX(0.0);}
        if (TMath::AreEqualAbs(tTos_vecTC.Y(),0.0,numeric_limits<double>::epsilon())) {tTos_vecTC.SetY(0.0);}
        if (TMath::AreEqualAbs(tTos_vecTC.Z(),0.0,numeric_limits<double>::epsilon())) {tTos_vecTC.SetZ(0.0);}
        double imgX = TMath::RadToDeg()*tTos_vecTC.X();
        double imgY = TMath::RadToDeg()*tTos_vecTC.Y();

        // Apply PSF: raytraced Fresnel optical model or Gaussian approximation
        std::tuple<double,double> coords = spread(imgX,imgY);
        double x = std::get<0>(coords);
        double y = std::get<1>(coords);
        if (std::isnan(x) || std::isnan(y)) continue;

        image->Fill(-1*x,y,petoadu); 
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
    image = addNoise(image);
    image->SetTitle("NSB Added");
    image->DrawCopy("COLZ1","");    

    //
    // Apply trigger
    //
    c->cd(3);

    // 2-pixel trigger in quadrants - SAME AS telEvent()
    int npixtrigq1=0, npixtrigq2=0, npixtrigq3=0, npixtrigq4=0;
    for(int i=1; i<=16; i++){
        for(int j=1; j<=16; j++){
            if(image->GetBinContent(i,j)>6.5*petoadu) npixtrigq1++;
        }
    }
    for(int i=17; i<=32; i++){
        for(int j=1; j<=16; j++){
            if(image->GetBinContent(i,j)>6.5*petoadu) npixtrigq2++;
        }
    }
    for(int i=1; i<=16; i++){
        for(int j=17; j<=32; j++){
            if(image->GetBinContent(i,j)>6.5*petoadu) npixtrigq3++;
        }
    }
    for(int i=17; i<=32; i++){
        for(int j=17; j<=32; j++){
            if(image->GetBinContent(i,j)>6.5*petoadu) npixtrigq4++;
        }
    }

    bool trig = (npixtrigq1>=2 || npixtrigq2>=2 || npixtrigq3>=2 || npixtrigq4>=2);

    if(!trig){
        image->Reset();
    } else {
        image = clean(image);
    }

    image->SetTitle("After Cleaning");
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
        TH2D* tmp = (TH2D*)time_grad->Clone();
        TH1D* arrival = new TH1D(Form("T%d",i+1), Form("T%d",i+1), 500, 0, 50 );

        // fill image
        for(int i=0; i<NCp; i++){
            // Apply PSF: raytraced Fresnel optical model or Gaussian approximation
            std::tuple<double,double> coords = spread(imgX[i],imgY[i]);
            double x = std::get<0>(coords);
            double y = std::get<1>(coords);
            if (std::isnan(x) || std::isnan(y)) continue;

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
        tmp = addNoise(tmp);

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
    int* npix = new int[Ntel];

    double* TelX = new double[Ntel];
    double* TelY = new double[Ntel];
    double* TelZ = new double[Ntel];

    for(int i=0; i<Ntel; i++){
        c->cd(i+1);
        gPad->SetTopMargin(0.1);
        gPad->SetBottomMargin(0.01);
        gPad->SetLeftMargin(0.01);
        gPad->SetRightMargin(0.15);

        TH2D* image = telEvent(i+1, eventNumber);
        image->DrawCopy("COLZ1","");
        // parameterization
        auto params = parameterize(image);
        npix[i] = countSignalPixels(image);
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
        "NPIX:\t\t%d\n"
        "LENGTH:\t\t%f\n"
        "WIDTH:\t\t%f\n"
        "MISS:\t\t%f\n"
        "DIST:\t\t%f\n"
        "ALPHA:\t\t%f\n",
        i+1, meanx[i],stdx[i],meany[i],stdy[i],phi[i],size[i],npix[i],length[i],width[i],miss[i],dist[i],alpha[i]);

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
