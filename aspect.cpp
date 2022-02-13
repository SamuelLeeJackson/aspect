/**
 * @file aspect.cpp
 * @author Samuel Lee Jackson (samuel.jackson@open.ac.uk), Ben Rozitis (ben.rozitis@open.ac.uk)
 * @brief Code to estimate the aspect related uncertainty in an individual asteroid phase curve
 * @version 0.1
 * @date 2021-10-25
 * 
 * @copyright Copyright (c) 2021
 * 
 */
#include <omp.h>
#include <time.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <random>
using namespace std;

// Constants

const double PI = 4.0 * atan(1.0);  // Value of PI
const double AU = 149.6e9;          // Astronomical unit in metres
const double FSUN = 1.896;          // Flux of Sun at 545 nm
const double Vzero = 3.631e-11;     // Zero point flux of V filter

// Shape Input Parameters
const int N = 24;
const int latspacing = 180.0/N;
const int ASTNUMVTX = (N*N)+2;
const int ASTNUMFACE = 2*N*N;

// Asteroid Input Parameters

const int NTAU = 18;                                                // Number of lightcurve timesteps
const int NUMTAX = 2; //4;                                                // Number of taxonomies to try
const string TAXONS[NUMTAX] = {"C", "S"}; //, "E", "M"};                  // Name of taxonomies
const double HAPKE_PARAMS[NUMTAX][5] = {
    {0.037, 1.03, 0.025, -0.47, (20.0 * (PI / 180.0))}, // (10) Hygiea
    {0.033, 1.40, 0.010, -0.25, (28.0 * (PI / 180.0))}, // (433) Eros
//    {0.660, 0.60, 0.027, -0.30, (28.0 * (PI / 180.0))}, // (44) Nysa
//    {0.226, 1.79, 0.041, -0.28, (28.0 * (PI / 180.0))}  // (216) Kleopatra
};

const int NMODEL = 1000;

// Multi-thread computing

const int num_threads = 2;  // Number of threads

// Cotangent function

double cot(double x) {
    return (1.0 / tan(x));
};

// NEA obliquity distribution PDF

double PDF(double x) {
    return (1.12*cos(x)*cos(x) - 0.32*cos(x) + 0.13)/(0.69*PI + sin(PI)*(0.56*cos(PI) - 0.32));
}

// Samples the obliquity from the PDF above, this is done by sampling a random obliquity from a uniform distribution, and sampling a random value of the PDF from a uniform distribution
// If the PDF at the randomly selected obliquity is greater than the random value of the PDF, then accept the sample
// This operates on the basis that where PDF(obl) is high, it is less likely that pdf_samp will be above this value, and this gives us more samples in this location
double sample_obliquity(std::uniform_real_distribution<double> obliquity_uniform, std::uniform_real_distribution<double> pdf_uniform, std::mt19937 mt_2) {
    bool sample_found = false;
    double obl;
    double pdf_samp;
    while (!sample_found) {
        obl = obliquity_uniform(mt_2);
        pdf_samp = pdf_uniform(mt_2);
        if (PDF(obl) > pdf_samp){
            sample_found = true;
        }
    }
    return obl;
};

// Class AstVertex defining position of a single vertex of the global asteroid shape model

class AstVertex {
   private:
    double _polarpos[2];
    double _equpos[3];
    double _heliopos[3];

   public:
    // Constructor
    AstVertex() {}

    double &rtnpolarpos(int n) { return _polarpos[n]; }
    double &rtnequpos(int n) { return _equpos[n]; }
    double &rtnheliopos(int n) { return _heliopos[n]; }

    // Destructor
    ~AstVertex() {}
};

// Class AstFacet defining surface properties of a single surface facet of the global asteroid shape model

class AstFacet {
   private:
    int _astvertices[3];
    double _heliomidpos[3];
    double _obsvector[3];
    double _illumvector[3];
    double _viewvector[3];
    double _normal[3];
    double _area;
    double _illumangle;
    double _viewangle;
    double _phaseangle;
    double _azangle;
    int _shadow;
    int _visible;

   public:
    // Constructor
    AstFacet() {}

    int &rtnastvertices(int n) { return _astvertices[n]; }
    double &rtnheliomidpos(int n) { return _heliomidpos[n]; }
    double &rtnobsvector(int n) { return _obsvector[n]; }
    double &rtnillumvector(int n) { return _illumvector[n]; }
    double &rtnviewvector(int n) { return _viewvector[n]; }
    double &rtnnormal(int n) { return _normal[n]; }
    double &rtnarea() { return _area; }
    double &rtnillumangle() { return _illumangle; }
    double &rtnviewangle() { return _viewangle; }
    double &rtnphaseangle() { return _phaseangle; }
    double &rtnazangle() { return _azangle; }
    int &rtnshadow() { return _shadow; }
    int &rtnvisible() { return _visible; }

    // Destructor
    ~AstFacet() {}
};

// Class Asteroid to control the asteroid global shape model

class Asteroid {
   private:
    int NV, NF;
    double _effective_diameter;
    double _tmatrix[3][3];
    double _polcentre[3];
    double _carcentre[3];
    double _obspolcentre[3];
    double _obscentre[3];
    double _rotvector[2];
    double _solar;
    double _Vmag[NTAU];
    double _phase[NTAU];
    double _avg_Vmag;
    double _ave_phase;
    double _min_Vmag;
    double _max_Vmag;

    AstVertex *V;
    AstFacet *F;

   public:
    AstVertex *beginv() const { return V; }
    AstVertex *endv() const { return V + NV; }

    AstFacet *beginf() const { return F; }
    AstFacet *endf() const { return F + NF; }

    // Asteroid constructor
    Asteroid(int n, int m) : NV(n), NF(m) {
        // Creates an array of AstVertex instances
        V = new AstVertex[n];
        // Creates and array of AstFacet instances
        F = new AstFacet[m];
    }

    // Get references to the private members of the class
    int &rtnnumvtx() { return NV; }
    int &rtnnumface() { return NF; }
    double &rtneffective_diameter() { return _effective_diameter; }
    double &rtntmatrix(int n, int m) { return _tmatrix[n][m]; }
    double &rtnpolcentre(int n) { return _polcentre[n]; }
    double &rtncarcentre(int n) { return _carcentre[n]; }
    double &rtnobspolcentre(int n) { return _obspolcentre[n]; }
    double &rtnobscentre(int n) { return _obscentre[n]; }
    double &rtnrotvector(int n) { return _rotvector[n]; }
    double &rtnsolar() { return _solar; }
    double &rtnVmag(int t) { return _Vmag[t]; }
    double &rtnphase(int t) { return _phase[t]; }
    double &rtnavg_Vmag() { return _avg_Vmag; }
    double &rtnavg_phase() { return _ave_phase; }
    double &rtnmin_Vmag() { return _min_Vmag; }
    double &rtnmax_Vmag() { return _max_Vmag; }
};

void createglobalmodel(Asteroid &A) {

    // Assign vertices
    AstVertex *V=A.beginv();

    int vertexnumber = 1;
    double lonspacing;

    // North pole vertex

    V->rtnpolarpos(0) = 90.0;
    V->rtnpolarpos(0) = 0.0;

    // Northern hemisphere vertices

    for (int i=1;i!=((N/2)+1);i++) {
        lonspacing=360.0/(4.0*i);
        for (int j=1;j!=(4*i)+1;j++) {
            (V+vertexnumber)->rtnpolarpos(0) = 90.0-(i*latspacing);
            (V+vertexnumber)->rtnpolarpos(1) = (j-1)*lonspacing;
            vertexnumber=vertexnumber+1;
        }
    }

    // Southern hemisphere vertices

    for (int i=((N/2)-1);i!=0;i--) {
        lonspacing=360.0/(4.0*i);
        for (int j=1;j!=(4*i)+1;j++) {
            (V+vertexnumber)->rtnpolarpos(0) = -90.0+(i*latspacing);
            (V+vertexnumber)->rtnpolarpos(1) = (j-1)*lonspacing;
            vertexnumber=vertexnumber+1;
        }
    }

    // South pole vertex

    (V+vertexnumber)->rtnpolarpos(0) = -90;
    (V+vertexnumber)->rtnpolarpos(1) = 0.0;

    // Assign facets

    AstFacet *F=A.beginf();

    int facenumber;
    facenumber=4;

    // North pole facets

    F->rtnastvertices(0)=1, F->rtnastvertices(1)=2, F->rtnastvertices(2)=0;
    (F+1)->rtnastvertices(0)=2, (F+1)->rtnastvertices(1)=3, (F+1)->rtnastvertices(2)=0;
    (F+2)->rtnastvertices(0)=3, (F+2)->rtnastvertices(1)=4, (F+2)->rtnastvertices(2)=0;
    (F+3)->rtnastvertices(0)=4, (F+3)->rtnastvertices(1)=1, (F+3)->rtnastvertices(2)=0;

    // Northern hemisphere facets

    for (int i=2;i!=((N/2)+1);i++) {
        for (int k=0;k!=3;++k) {
            for (int j=1;j!=i+1;j++) {
                (F+facenumber)->rtnastvertices(0)=(k*i)+(2*i*(i-1))+j;
                (F+facenumber)->rtnastvertices(1)=(k*i)+(2*i*(i-1))+j+1;
                (F+facenumber)->rtnastvertices(2)=(k*(i-1))+(2*(i-1)*(i-2))+j;
                facenumber=facenumber+1;
            }
        }
        for (int j=1;j!=i;j++) {
            (F+facenumber)->rtnastvertices(0)=(3*i)+(2*i*(i-1))+j;
            (F+facenumber)->rtnastvertices(1)=(3*i)+(2*i*(i-1))+j+1;
            (F+facenumber)->rtnastvertices(2)=(3*(i-1))+(2*(i-1)*(i-2))+j;
            facenumber=facenumber+1;
        }
        (F+facenumber)->rtnastvertices(0)=(2*i*(i+1));
        (F+facenumber)->rtnastvertices(1)=(2*i*(i-1))+1;
        (F+facenumber)->rtnastvertices(2)=(2*(i-1)*(i-2))+1;
        facenumber=facenumber+1;
    }

    for (int i=2;i!=((N/2)+1);i++) {
        for (int k=0;k!=3;k++) {
            for (int j=1;j!=i;j++) {
                (F+facenumber)->rtnastvertices(0)=(k*(i-1))+(2*(i-1)*(i-2))+j;
                (F+facenumber)->rtnastvertices(1)=(k*i)+(2*i*(i-1))+j+1;
                (F+facenumber)->rtnastvertices(2)=(k*(i-1))+(2*(i-1)*(i-2))+j+1;
                facenumber=facenumber+1;
            }
        }
        for (int j=1;j!=i-1;j++) {
            (F+facenumber)->rtnastvertices(0)=(3*(i-1))+(2*(i-1)*(i-2))+j;
            (F+facenumber)->rtnastvertices(1)=(3*i)+(2*i*(i-1))+j+1;
            (F+facenumber)->rtnastvertices(2)=(3*(i-1))+(2*(i-1)*(i-2))+j+1;
            facenumber=facenumber+1;
        }
        (F+facenumber)->rtnastvertices(0)=(2*i*(i-1));
        (F+facenumber)->rtnastvertices(1)=(2*i*(i+1));
        (F+facenumber)->rtnastvertices(2)=(2*(i-1)*(i-2))+1;
        facenumber=facenumber+1;
    }

    // Southern hemisphere facets

    for (int i=2;i!=((N/2)+1);i++) {
        for (int k=0;k!=3;k++) {
            for (int j=1;j!=i+1;j++) {
                (F+facenumber)->rtnastvertices(0)=(N*N)+(k*i)-(2*i*(i+1))+j+1;
                (F+facenumber)->rtnastvertices(1)=(N*N)+(k*i)-(2*i*(i+1))+j;
                (F+facenumber)->rtnastvertices(2)=(N*N)+(k*(i-1))-(2*i*(i-1))+j;
                facenumber=facenumber+1;
            }
        }
        for (int j=1;j!=i;j++) {
            (F+facenumber)->rtnastvertices(0)=(N*N)+(3*i)-(2*i*(i+1))+j+1;
            (F+facenumber)->rtnastvertices(1)=(N*N)+(3*i)-(2*i*(i+1))+j;
            (F+facenumber)->rtnastvertices(2)=(N*N)+(3*(i-1))-(2*i*(i-1))+j;
            facenumber=facenumber+1;
        }
        (F+facenumber)->rtnastvertices(0)=(N*N)-(2*i*(i+1))+1;
        (F+facenumber)->rtnastvertices(1)=(N*N)-(2*i*(i-1));
        (F+facenumber)->rtnastvertices(2)=(N*N)-(2*i*(i-1))+1;
        facenumber=facenumber+1;
    }

    for (int i=2;i!=((N/2)+1);i++) {
        for (int k=0;k!=3;k++) {
            for (int j=1;j!=i;j++) {
                (F+facenumber)->rtnastvertices(0)=(N*N)+(k*(i-1))-(2*i*(i-1))+j+1;
                (F+facenumber)->rtnastvertices(1)=(N*N)+(k*i)-(2*i*(i+1))+j+1;
                (F+facenumber)->rtnastvertices(2)=(N*N)+(k*(i-1))-(2*i*(i-1))+j;
                facenumber=facenumber+1;
            }
        }
        for (int j=1;j!=i-1;j++) {
            (F+facenumber)->rtnastvertices(0)=(N*N)+(3*(i-1))-(2*i*(i-1))+j+1;
            (F+facenumber)->rtnastvertices(1)=(N*N)+(3*i)-(2*i*(i+1))+j+1;
            (F+facenumber)->rtnastvertices(2)=(N*N)+(3*(i-1))-(2*i*(i-1))+j;
            facenumber=facenumber+1;
        }
        (F+facenumber)->rtnastvertices(0)=(N*N)-(2*i*(i-1))+1;
        (F+facenumber)->rtnastvertices(1)=(N*N)-(2*i*(i-1));
        (F+facenumber)->rtnastvertices(2)=(N*N)-(2*(i-1)*(i-2));
        facenumber=facenumber+1;
    }

    // South pole facets

    (F+facenumber)->rtnastvertices(0)=(N*N)-2, (F+facenumber)->rtnastvertices(1)=(N*N)-3, (F+facenumber)->rtnastvertices(2)=(N*N)+1;
    (F+facenumber+1)->rtnastvertices(0)=(N*N)-1, (F+facenumber+1)->rtnastvertices(1)=(N*N)-2, (F+facenumber+1)->rtnastvertices(2)=(N*N)+1;
    (F+facenumber+2)->rtnastvertices(0)=N*N, (F+facenumber+2)->rtnastvertices(1)=(N*N)-1, (F+facenumber+2)->rtnastvertices(2)=(N*N)+1;
    (F+facenumber+3)->rtnastvertices(0)=(N*N)-3, (F+facenumber+3)->rtnastvertices(1)=N*N, (F+facenumber+3)->rtnastvertices(2)=(N*N)+1;
}

// Determine asteroid pole geometry

void asteroidgeometry(Asteroid &A, double lambda_pole, double beta_pole, double a, double b, double c1, double c2) {
    double lambda, beta, theta;

    lambda = lambda_pole * (PI / 180.0);
    beta = beta_pole * (PI / 180.0);
    theta = (PI / 2.0) - beta;

    // Define rotation transformation matrix between body-centric equatorial and body-centric ecliptic coordinates
    A.rtntmatrix(0, 0) = cos(lambda) * cos(theta);
    A.rtntmatrix(0, 1) = -1.0 * sin(lambda);
    A.rtntmatrix(0, 2) = cos(lambda) * sin(theta);
    A.rtntmatrix(1, 0) = sin(lambda) * cos(theta);
    A.rtntmatrix(1, 1) = cos(lambda);
    A.rtntmatrix(1, 2) = sin(lambda) * sin(theta);
    A.rtntmatrix(2, 0) = -1.0 * sin(theta);
    A.rtntmatrix(2, 1) = 0.0;
    A.rtntmatrix(2, 2) = cos(theta);

    AstVertex *B;
    for (int v = 0; v < ASTNUMVTX; v++) {
        // For each vertex in model

        // Gets initial vertex and adds v to get pointer of required vertex
        B = A.beginv() + v;

        B->rtnequpos(0) = a * cos(B->rtnpolarpos(0) * (PI / 180.0)) * cos(B->rtnpolarpos(1) * (PI / 180.0));  // Body-centric cartesian x of vertex (scaled by a of ellipse)
        B->rtnequpos(1) = b * cos(B->rtnpolarpos(0) * (PI / 180.0)) * sin(B->rtnpolarpos(1) * (PI / 180.0));  // Body-centric cartesian y of vertex (scaled by b of ellipse)

        if (B->rtnpolarpos(0) >= 0.0) {
            // If latitude of vertex is in the northern hemisphere
            B->rtnequpos(2) = c1 * sin(B->rtnpolarpos(0) * (PI / 180.0));  // Body-centric cartesian z of vertex (scaled by c1)
        } else {
            B->rtnequpos(2) = c2 * sin(B->rtnpolarpos(0) * (PI / 180.0));  // Body-centric cartesian z of vertex (scaled by c2)
        }
    }
};

// Determine asteroid global shape model lightcurve geometry

void lightcurvegeometry(Asteroid &A, int step) {
    // Scale the solar flux by the heliocentric range squared to get flux at asteroid
    A.rtnsolar() = FSUN / (A.rtnpolcentre(0) * A.rtnpolcentre(0));

    A.rtncarcentre(0) = AU * A.rtnpolcentre(0) * cos(A.rtnpolcentre(2)) * cos(A.rtnpolcentre(1));  // heliocentric range (metres) * cos(heliocentric lat) * cos(heliocentric lon)  -> target cartesian heliocentric x
    A.rtncarcentre(1) = AU * A.rtnpolcentre(0) * cos(A.rtnpolcentre(2)) * sin(A.rtnpolcentre(1));  // heliocentric range (metres) * cos(heliocentric lat) * sin(heliocentric lon)  -> target cartesian heliocentric y
    A.rtncarcentre(2) = AU * A.rtnpolcentre(0) * sin(A.rtnpolcentre(2));                           // heliocentric range (metres) * sin(heliocentric lat)							-> target cartesian heliocentric z

    A.rtnobscentre(0) = A.rtncarcentre(0) - (AU * A.rtnobspolcentre(0) * cos(A.rtnobspolcentre(2)) * cos(A.rtnobspolcentre(1)));  // target cartesian heliocentric x - cartesian geocentric x -> observer cartesian heliocentric x
    A.rtnobscentre(1) = A.rtncarcentre(1) - (AU * A.rtnobspolcentre(0) * cos(A.rtnobspolcentre(2)) * sin(A.rtnobspolcentre(1)));  // target cartesian heliocentric y - cartesian geocentric y -> observer cartesian heliocentric x
    A.rtnobscentre(2) = A.rtncarcentre(2) - (AU * A.rtnobspolcentre(0) * sin(A.rtnobspolcentre(2)));                              // target cartesian heliocentric z - cartesian geocentric z -> observer cartesian heliocentric z

    double rotangle, rotequpos[3], eclipos[3];
    double normcalc[2][3], y[3], _n, _t, _x, _o;
    double imatrix[3][3];

    AstVertex *B;
    AstFacet *C;

    // Rotation phase of asteroid to compute light curve at
    rotangle = 2.0 * PI * step / NTAU;
    
    for (int v = 0; v < ASTNUMVTX; v++) {
        // For each vertex in model

        // Gets initial vertex and adds v to get pointer of required vertex
        B = A.beginv() + v;

        rotequpos[0] = (B->rtnequpos(0) * cos(rotangle)) - (B->rtnequpos(1) * sin(rotangle));  // Rotate body-centric cartesian x of vertex
        rotequpos[1] = (B->rtnequpos(0) * sin(rotangle)) + (B->rtnequpos(1) * cos(rotangle));  // Rotate body-centric cartesian y of vertex
        rotequpos[2] = B->rtnequpos(2);                                                        // Rotate body-centric cartesian z of vertex

        eclipos[0] = (A.rtntmatrix(0, 0) * rotequpos[0]) + (A.rtntmatrix(0, 1) * rotequpos[1]) + (A.rtntmatrix(0, 2) * rotequpos[2]);  // Convert body-centric equatorial x to body-centric ecliptic x (i.e. in reference frame of Sun & Earth)
        eclipos[1] = (A.rtntmatrix(1, 0) * rotequpos[0]) + (A.rtntmatrix(1, 1) * rotequpos[1]) + (A.rtntmatrix(1, 2) * rotequpos[2]);  // Convert body-centric equatorial x to body-centric ecliptic y
        eclipos[2] = (A.rtntmatrix(2, 0) * rotequpos[0]) + (A.rtntmatrix(2, 1) * rotequpos[1]) + (A.rtntmatrix(2, 2) * rotequpos[2]);  // Convert body-centric equatorial x to body-centric ecliptic z

        // Add ecliptic heliocentric position of body centre to ecliptic body-centric vertex positions
        B->rtnheliopos(0) = A.rtncarcentre(0) + eclipos[0];  // Heliocentric ecliptic x coordinate
        B->rtnheliopos(1) = A.rtncarcentre(1) + eclipos[1];  // Heliocentric ecliptic y coordinate
        B->rtnheliopos(2) = A.rtncarcentre(2) + eclipos[2];  // Heliocentric ecliptic z coordinate
    }
    
    for (int f = 0; f < ASTNUMFACE; f++) {
        // For each facet in model

        // Gets initial facet and adds f to get pointer of required facet
        C = A.beginf() + f;

        // C->rtnastvertices(n) gets the vertex number of the vertex
        // A.beginv() gets the pointer to the first vertex in memory
        // Since the vertices will be stored sequentially in memory, adding the vertex number to the starting position of the array in memory gives the pointer for the required vertex
        // A.beginv() + C->rtnastvertices(n) gives the pointer to the vertex
        // (A.beginv() + C->rtnastvertices(n))->rtnheliopos(n) gives the nth axis of the heliocentric position vector of the vertex
        C->rtnheliomidpos(0) = (((A.beginv() + C->rtnastvertices(0))->rtnheliopos(0)) + ((A.beginv() + C->rtnastvertices(1))->rtnheliopos(0)) + ((A.beginv() + C->rtnastvertices(2))->rtnheliopos(0))) / 3.0;  // Average heliocentric x position of facet using each vertex
        C->rtnheliomidpos(1) = (((A.beginv() + C->rtnastvertices(0))->rtnheliopos(1)) + ((A.beginv() + C->rtnastvertices(1))->rtnheliopos(1)) + ((A.beginv() + C->rtnastvertices(2))->rtnheliopos(1))) / 3.0;  // Average heliocentric y position of facet using each vertex
        C->rtnheliomidpos(2) = (((A.beginv() + C->rtnastvertices(0))->rtnheliopos(2)) + ((A.beginv() + C->rtnastvertices(1))->rtnheliopos(2)) + ((A.beginv() + C->rtnastvertices(2))->rtnheliopos(2))) / 3.0;  // Average heliocentric z position of facet using each vertex

        // Heliocentric facet position vector magnitude
        _t = sqrt((C->rtnheliomidpos(0) * C->rtnheliomidpos(0)) + (C->rtnheliomidpos(1) * C->rtnheliomidpos(1)) + (C->rtnheliomidpos(2) * C->rtnheliomidpos(2)));

        // Creates Earth position vector from asteroid reference frame
        C->rtnobsvector(0) = C->rtnheliomidpos(0) - A.rtnobscentre(0);  // Facet heliocentric x - observer heliocentric x - > Observer body-centric x
        C->rtnobsvector(1) = C->rtnheliomidpos(1) - A.rtnobscentre(1);  // Facet heliocentric y - observer heliocentric y - > Observer body-centric y
        C->rtnobsvector(2) = C->rtnheliomidpos(2) - A.rtnobscentre(2);  // Facet heliocentric z - observer heliocentric z - > Observer body-centric z

        // Observer position vector magnitude
        _o = sqrt((C->rtnobsvector(0) * C->rtnobsvector(0)) + (C->rtnobsvector(1) * C->rtnobsvector(1)) + (C->rtnobsvector(2) * C->rtnobsvector(2)));

        // Vertex 2 position - Vertex 1 position
        normcalc[0][0] = ((A.beginv() + C->rtnastvertices(1))->rtnheliopos(0)) - ((A.beginv() + C->rtnastvertices(0))->rtnheliopos(0));  // Heliocentric x of vertex 2 - heliocentric x of vertex 1
        normcalc[0][1] = ((A.beginv() + C->rtnastvertices(1))->rtnheliopos(1)) - ((A.beginv() + C->rtnastvertices(0))->rtnheliopos(1));  // Heliocentric y of vertex 2 - heliocentric y of vertex 1
        normcalc[0][2] = ((A.beginv() + C->rtnastvertices(1))->rtnheliopos(2)) - ((A.beginv() + C->rtnastvertices(0))->rtnheliopos(2));  // Heliocentric z of vertex 2 - heliocentric z of vertex 1

        // Vertex 3 position - Vertex 1 position
        normcalc[1][0] = ((A.beginv() + C->rtnastvertices(2))->rtnheliopos(0)) - ((A.beginv() + C->rtnastvertices(0))->rtnheliopos(0));  // Heliocentric x of vertex 3 - heliocentric x of vertex 2
        normcalc[1][1] = ((A.beginv() + C->rtnastvertices(2))->rtnheliopos(1)) - ((A.beginv() + C->rtnastvertices(0))->rtnheliopos(1));  // Heliocentric y of vertex 3 - heliocentric y of vertex 2
        normcalc[1][2] = ((A.beginv() + C->rtnastvertices(2))->rtnheliopos(2)) - ((A.beginv() + C->rtnastvertices(0))->rtnheliopos(2));  // Heliocentric z of vertex 3 - heliocentric z of vertex 2

        // Cross product of the two difference vectors to obtain the normal vector
        C->rtnnormal(0) = (normcalc[0][1] * normcalc[1][2]) - (normcalc[0][2] * normcalc[1][1]);
        C->rtnnormal(1) = (normcalc[0][2] * normcalc[1][0]) - (normcalc[0][0] * normcalc[1][2]);
        C->rtnnormal(2) = (normcalc[0][0] * normcalc[1][1]) - (normcalc[0][1] * normcalc[1][0]);

        // Magnitude of the facet normal position vector
        _n = sqrt((C->rtnnormal(0) * C->rtnnormal(0)) + (C->rtnnormal(1) * C->rtnnormal(1)) + (C->rtnnormal(2) * C->rtnnormal(2)));

        // Area of the facet is equal to half the magnitude of the facet normal (from cross product also equalling |a||b|sin(theta) = area of parallelogram from two vectors)
        C->rtnarea() = _n / 2.0;

        // Make facet normal into a unit vector
        C->rtnnormal(0) = C->rtnnormal(0) / _n;
        C->rtnnormal(1) = C->rtnnormal(1) / _n;
        C->rtnnormal(2) = C->rtnnormal(2) / _n;

        // Magnitude of the difference between vertex 2 and vertex 1 of the facet
        _x = sqrt((normcalc[0][0] * normcalc[0][0]) + (normcalc[0][1] * normcalc[0][1]) + (normcalc[0][2] * normcalc[0][2]));

        // Make difference between vertex 2 and vertex 1 of the facet into a unit vector
        normcalc[0][0] = normcalc[0][0] / _x;
        normcalc[0][1] = normcalc[0][1] / _x;
        normcalc[0][2] = normcalc[0][2] / _x;

        // Cross product between facet normal unit vector and unit vector describing
        y[0] = (C->rtnnormal(1) * normcalc[0][2]) - (C->rtnnormal(2) * normcalc[0][1]);
        y[1] = (C->rtnnormal(2) * normcalc[0][0]) - (C->rtnnormal(0) * normcalc[0][2]);
        y[2] = (C->rtnnormal(0) * normcalc[0][1]) - (C->rtnnormal(1) * normcalc[0][0]);

        // Looks like it creates a facet-centric coordinate system

        // A unit vector pointing from one vertex to another (x)
        imatrix[0][0] = normcalc[0][0];
        imatrix[0][1] = normcalc[0][1];
        imatrix[0][2] = normcalc[0][2];
        // The unit vector describing the cross product of the vertex difference vector and the facet normal unit vector (y)
        imatrix[1][0] = y[0];
        imatrix[1][1] = y[1];
        imatrix[1][2] = y[2];
        // The facet normal unit vector (z)
        imatrix[2][0] = C->rtnnormal(0);
        imatrix[2][1] = C->rtnnormal(1);
        imatrix[2][2] = C->rtnnormal(2);

        // This is a matrix-vector product
        // Does it make it a facet-centric position vector for the Sun
        // -1.0 multiplication is because rtnheliomidpos(0..2) gives a position vector from the Sun to the Target, we want Target to Sun
        C->rtnillumvector(0) = -1.0 * ((imatrix[0][0] * C->rtnheliomidpos(0)) + (imatrix[0][1] * C->rtnheliomidpos(1)) + (imatrix[0][2] * C->rtnheliomidpos(2)));
        C->rtnillumvector(1) = -1.0 * ((imatrix[1][0] * C->rtnheliomidpos(0)) + (imatrix[1][1] * C->rtnheliomidpos(1)) + (imatrix[1][2] * C->rtnheliomidpos(2)));
        C->rtnillumvector(2) = -1.0 * ((imatrix[2][0] * C->rtnheliomidpos(0)) + (imatrix[2][1] * C->rtnheliomidpos(1)) + (imatrix[2][2] * C->rtnheliomidpos(2)));

        // Calculates angle of illumination
        C->rtnillumangle() = acos(C->rtnillumvector(2) / _t);

        // If angle of illumination > 90 degrees then facet is not illuminated
        if (C->rtnillumangle() > (PI / 2.0)) {
            C->rtnshadow() = 0;
        } else {
            C->rtnshadow() = 1;
        }

        // This is a matrix-vector product
        // Does it make a facet-centric position vector for the Observer
        // -1.0 multiplication is because rtnobsvector(0..2) gives a position vector from the Observer to the Target, we want Target to Observer
        C->rtnviewvector(0) = -1.0 * ((imatrix[0][0] * C->rtnobsvector(0)) + (imatrix[0][1] * C->rtnobsvector(1)) + (imatrix[0][2] * C->rtnobsvector(2)));
        C->rtnviewvector(1) = -1.0 * ((imatrix[1][0] * C->rtnobsvector(0)) + (imatrix[1][1] * C->rtnobsvector(1)) + (imatrix[1][2] * C->rtnobsvector(2)));
        C->rtnviewvector(2) = -1.0 * ((imatrix[2][0] * C->rtnobsvector(0)) + (imatrix[2][1] * C->rtnobsvector(1)) + (imatrix[2][2] * C->rtnobsvector(2)));

        // Calculated facet viewing angle
        C->rtnviewangle() = acos(C->rtnviewvector(2) / _o);

        // if viewing angle > 90 degrees then facet is not visible
        if (C->rtnviewangle() > (PI / 2.0)) {
            C->rtnvisible() = 0;
        } else {
            C->rtnvisible() = 1;
        }

        // Returns phase angle of facet
        C->rtnphaseangle() = acos(((C->rtnillumvector(0) * C->rtnviewvector(0)) + (C->rtnillumvector(1) * C->rtnviewvector(1)) + (C->rtnillumvector(2) * C->rtnviewvector(2))) / (_t * _o));

        // Returns azimuth angle of facet (what is this?)
        C->rtnazangle() = acos((cos(C->rtnphaseangle()) - (cos(C->rtnillumangle()) * cos(C->rtnviewangle()))) / (sin(C->rtnillumangle()) * sin(C->rtnviewangle())));
        if (C->rtnazangle() != C->rtnazangle()) {
            C->rtnazangle() = 0.0;
        }
    }
};

// BEGIN HAPKE SCATTERING FUNCTIONS

// B function

double B_function(double B0, double h, double alpha) {
    return B0 / (1.0 + (tan(alpha / 2.0) / h));
};

// P function

double P_function(double G, double alpha) {
    return (1.0 - (G * G)) / pow((1.0 + (2.0 * G * cos(alpha)) + (G * G)), 1.5);
};

// Helper functions to make the mathematics more readable, corresponds directly to formulation in the paper and in Hapke (1993)

double chi_function(double theta) {
    return (1.0 / sqrt(1.0 + (PI * tan(theta) * tan(theta))));
}

double E1_function(double y, double theta) {
    return exp(-(2.0/PI) * cot(theta) * cot(y));
}

double E2_function(double y, double theta) {
    return exp(-(1.0/PI) * pow(cot(theta), 2) * pow(cot(y), 2));
}

double eta_function(double y, double theta) {
    return chi_function(theta) * (cos(y) + (sin(y) * tan(theta) * E2_function(y, theta) / (2.0 - E1_function(y, theta))));
}

// mu0dash function

double mu_0e_function(double i, double e, double az, double theta) {
    double mu0dash;

    if (i <= e) {
        mu0dash = (cos(az) * E2_function(e, theta)) + (pow(sin(az / 2.0), 2) * E2_function(i, theta));
        mu0dash = mu0dash / (2.0 - E1_function(e, theta) - ((az / PI) * E1_function(i, theta)));
        mu0dash = cos(i) + (sin(i) * tan(theta) * mu0dash);
        mu0dash = mu0dash * chi_function(theta);
        return mu0dash;
    } else {
        mu0dash = E2_function(i, theta) - (pow(sin(az / 2.0), 2) * E2_function(e, theta));
        mu0dash = mu0dash / (2.0 - E1_function(i, theta) - ((az / PI) * E1_function(e, theta)));
        mu0dash = cos(i) + (sin(i) * tan(theta) * mu0dash);
        mu0dash = mu0dash * chi_function(theta);
        return mu0dash;
    }
};

// mudash function

double mu_e_function(double i, double e, double az, double theta) {
    double mudash;

    if (i <= e) {
        mudash = E2_function(e, theta) - (pow(sin(az / 2.0), 2) * E2_function(i, theta));
        mudash = mudash / (2.0 - E1_function(e, theta) - ((az / PI) * E1_function(i, theta)));
        mudash = cos(e) + (sin(e) * tan(theta) * mudash);
        mudash = mudash * chi_function(theta);
        return mudash;
    } else {
        mudash = (cos(az) * E2_function(i, theta)) + (pow(sin(az / 2.0), 2) * E2_function(e, theta));
        mudash = mudash / (2.0 - E1_function(i, theta) - ((az / PI) * E1_function(e, theta)));
        mudash = cos(e) + (sin(e) * tan(theta) * mudash);
        mudash = mudash * chi_function(theta);
        return mudash;
    }
};

// F function

double F_function(double az) {
    return exp(-2.0 * tan(az / 2.0));
};

// S function

double S_function(double i, double e, double az, double theta) {
    double S;

    if (i <= e) {
        S = (mu_e_function(i, e, az, theta) / eta_function(e, theta)) * (cos(i) / eta_function(i, theta)) * chi_function(theta);
        S = S / (1.0 - F_function(az) + ((F_function(az) * chi_function(theta)) * (cos(i) / eta_function(i, theta))));
        return S;
    } else {
        S = (mu_e_function(i, e, az, theta) / eta_function(e, theta)) * (cos(i) / eta_function(i, theta)) * chi_function(theta);
        S = S / (1.0 - F_function(az) + ((F_function(az) * chi_function(theta)) * (cos(e) / eta_function(e, theta))));
        return S;
    }
};

// H function

double H_function(double w, double x) {
    return (1.0 + (2.0 * x)) / (1.0 + (2.0 * x * sqrt(1.0 - w)));
};

// rbi function

double rbi_function(double i, double e, double az, double alpha, double w, double B0, double h, double G, double theta) {
    double rbi;
    rbi = (w / (4.0 * PI)) * (mu_0e_function(i, e, az, theta) / (mu_0e_function(i, e, az, theta) + mu_e_function(i, e, az, theta)));
    rbi = rbi * (((1.0 + B_function(B0, h, alpha)) * P_function(G, alpha)) - 1.0 + (H_function(w, mu_0e_function(i, e, az, theta)) * H_function(w, mu_e_function(i, e, az, theta))));
    rbi = rbi * S_function(i, e, az, theta);
    return rbi;
};

// END HAPKE SCATTERING FUNCTIONS

// Code to observe asteroid reflected flux

void hapke_observe(Asteroid &A, int step, double ast_diameter, double ast_w, double ast_B0, double ast_h, double ast_G, double ast_theta) {
    double _d;
    double facet_flux;
    double astreflectedtotal;

    astreflectedtotal = 0.0;

    AstFacet *B = A.beginf();
    A.rtnphase(step) = B->rtnphaseangle();

    // For each facet in the model
    for (; B != A.endf(); ++B) {
        
        // Only calculate flux if both visible by observer and illuminated by the sun
        if (B->rtnshadow() != 0 && B->rtnvisible() != 0) {
            // Magnitude of view vector
            _d = sqrt((B->rtnviewvector(0) * B->rtnviewvector(0)) + (B->rtnviewvector(1) * B->rtnviewvector(1)) + (B->rtnviewvector(2) * B->rtnviewvector(2)));
            // Calculates total flux of facet, scaled by the magnitude of the view vector squared (inverse square law) - uses Hapke functions and passes hapke parameters in
            facet_flux = ((A.rtnsolar() * rbi_function(B->rtnillumangle(), B->rtnviewangle(), B->rtnazangle(), B->rtnphaseangle(), ast_w, ast_B0, ast_h, ast_G, ast_theta) * B->rtnarea() * 1.0e6 * cos(B->rtnviewangle())) / (_d * _d));
            if (facet_flux != facet_flux) {
                facet_flux = 0.0;
            }

            // Adds individual facet flux to the total
            astreflectedtotal = astreflectedtotal + facet_flux;
        }
    }

    // Scale reflected flux so that it is equivalent to object with diameter of 1km
    astreflectedtotal = astreflectedtotal * ((ast_diameter * ast_diameter) / (A.rtneffective_diameter() * A.rtneffective_diameter()));

    // Make it a reduced flux to 1AU from Sun and Earth
    astreflectedtotal = astreflectedtotal * A.rtnpolcentre(0) * A.rtnpolcentre(0) * A.rtnobspolcentre(0) * A.rtnobspolcentre(0);

    // Give Johnson V magnitude at given rotational phase
    A.rtnVmag(step) = -2.5 * log10(astreflectedtotal / Vzero);
};

// Main program loop

int main(int argc, char *argv[]) {

    if (argc != 3) {
        cerr << "Usage: " << argv[0] << " ASTEROID" << " OBS_LC_AMP" << endl;
        return -1;
    }

    const string asteroid = argv[1];
    const string min_lc_amp_str = argv[2];

    const double min_lc_amp = stod(min_lc_amp_str);
    const double min_ab_ratio = pow(10, 0.4 * min_lc_amp);
    // cout << min_ab_ratio << endl;

    // MAKE INPUT FILENAME FROM THE ABOVE STRING

    ifstream inclinations("numbered_asteroid_inclinations_20211014.txt");

    bool found = false;
    string delimiter = "\t";
    string inclination_str = "0.0";
    string line;
    int delimiter_pos;
    int end_line;
    string asteroid_name;
    while (std::getline(inclinations, line) && !found) {
        delimiter_pos = line.find(delimiter);
        end_line = line.find("\n");
        asteroid_name = line.substr(0, delimiter_pos);
        if (asteroid_name == asteroid) {
            inclination_str = line.substr(delimiter_pos + delimiter.length(), end_line);
            found = true;
        }
    }
    
    const double INCLINATION = stod(inclination_str) * PI/180.0;
    const double MAX_PDF = PDF(PI);

    std::random_device rd;
    std::mt19937 mt(rd());
    
    std::uniform_real_distribution<double> lambda(0.0, 2.0 * PI);
    std::uniform_real_distribution<double> obliquity(0, PI);
    std::uniform_real_distribution<double> pdf_dist(0, MAX_PDF);
    std::uniform_real_distribution<double> lambda_choose(0.0, 1.0);
    std::uniform_real_distribution<double> a(min_ab_ratio, 4.5);
    std::uniform_real_distribution<double> b(1.0, 2.0);
    std::uniform_real_distribution<double> c1(0.75, 1.25);
    std::uniform_real_distribution<double> c2(0.75, 1.25);

    double model_variants[NMODEL][6];
    double shape_ra;
    double shape_dec;
    double shape_obl;
    double shape_sine_beta;
    double shape_cos_lambda;
    double shape_lambda;
    double shape_beta;
    double shape_a;
    double shape_b;
    double shape_c1;
    double shape_c2;

    model_variants[0][0] = 0;
    model_variants[0][1] = 90;
    model_variants[0][2] = 1;
    model_variants[0][3] = 1;
    model_variants[0][4] = 1;
    model_variants[0][5] = 1;
    
    cout << 0 << "\t" << 0 << "\t" << 90 << "\t" << 1 << "\t" << 1 << "\t" << 1 << "\t" << 1 << endl;

    for (int i = 1; i < NMODEL; ++i) {
        
        shape_ra = lambda(mt);
        shape_obl = sample_obliquity(obliquity, pdf_dist, mt);
        shape_dec = (PI/2.0) - shape_obl;

        shape_sine_beta = sin(shape_dec) * cos(INCLINATION) - cos(shape_dec) * sin(INCLINATION) * sin(shape_ra);
        shape_beta = asin(shape_sine_beta);
        
        if (cos(shape_beta) == 0) {
            shape_lambda = 0.0;
        } else {
            shape_cos_lambda = (cos(shape_ra) * cos(shape_dec)) / cos(shape_beta);
            if (lambda_choose(mt) < 0.5) {
                shape_lambda = acos(shape_cos_lambda);
            } else {
                shape_lambda = (2 * PI) - acos(shape_cos_lambda);
            }
        }
        shape_beta = shape_beta * 180.0/PI;
        shape_lambda = shape_lambda * 180.0/PI;

        shape_a = a(mt);
        shape_b = b(mt);
        while ((shape_a/shape_b) < min_ab_ratio) {
            shape_b = b(mt);
        }
        shape_c1 = c1(mt);
        while (shape_c1 > shape_b) {
            shape_c1 = c1(mt);
        }
        shape_c2 = c2(mt);
        while (shape_c2 > shape_b) {
            shape_c2 = c2(mt);
        }

        model_variants[i][0] = shape_lambda;
        model_variants[i][1] = shape_beta;
        model_variants[i][2] = shape_a;
        model_variants[i][3] = shape_b;
        model_variants[i][4] = shape_c1;
        model_variants[i][5] = shape_c2;
        cout << i << "\t" << shape_lambda << "\t" << shape_beta << "\t" << shape_a << "\t" << shape_b << "\t" << shape_c1 << "\t" << shape_c2 << endl;
    }
    // cout << "GENERATED ALL SHAPES" << endl;

    const string readname = asteroid + "-geom.txt";

    const int NPOS = [&]() {
        int NPOS = 0;
        string line;
        ifstream myfile(readname);

        while (getline(myfile, line)) {
            ++NPOS;
        }
        return NPOS;
    }();

    const int num_input_cols = 7;

    double model_geometry[NPOS][num_input_cols];

    ifstream read(readname.c_str());
    // For each geometry
    for (int i = 0; i != NPOS; ++i) {
        // Read in each of the 6 parameters (heliocentric range, lon, lat, geocentric range, lon, lat, observed uncertainty)
        for (int j = 0; j != num_input_cols; ++j) {
            read >> model_geometry[i][j];
        }
    }

    ofstream Scurves(asteroid + "-S_type_models.txt");
    ofstream Ccurves(asteroid + "-C_type_models.txt");

omp_set_num_threads(num_threads);
#pragma	omp parallel for schedule(static) shared(model_variants, model_geometry)
    for (int tax = 0; tax < NUMTAX; ++tax) {

        // Initialise an asteroid instance with a number of vertices and faces
        Asteroid A(ASTNUMVTX, ASTNUMFACE);

        // Read in the model
        createglobalmodel(A);

        const double diameter = 1.0;                   // Asteroid diameter in km
        const double w = HAPKE_PARAMS[tax][0];         // Asteroid single scattering albedo
        const double B0 = HAPKE_PARAMS[tax][1];        // Asteroid opposition surge amplitude
        const double h = HAPKE_PARAMS[tax][2];         // Asteroid opposition surge width
        const double g = HAPKE_PARAMS[tax][3];         // Asteroid scattering asymmetry parameter
        const double theta = HAPKE_PARAMS[tax][4];     // Asteroid macroscopic roughness

        double H_array[NMODEL];
        double G_array[NMODEL];

        for (int i = 0; i != NMODEL; ++i) {
            // cout << i << endl;
            // Pass pole longitude and latitude to Asteroid instance to calculate pole geometry
            asteroidgeometry(A, model_variants[i][0], model_variants[i][1], model_variants[i][2], model_variants[i][3], model_variants[i][4], model_variants[i][5]);

            //Calculates the effective diameter of the ellipsoidal model using a, b, c1, & c2
            A.rtneffective_diameter() = pow((model_variants[i][2] * model_variants[i][3] * model_variants[i][4]), (1.0 / 3.0)) + pow((model_variants[i][2] * model_variants[i][3] * model_variants[i][5]), (1.0 / 3.0));

            // Loop over all viewing geometries to calculate average over
            for (int j = 0; j != NPOS; ++j) {
                A.rtnpolcentre(0) = model_geometry[j][0];                  // Heliocentric range
                A.rtnpolcentre(1) = model_geometry[j][1] * PI / 180.0;     // Heliocentric longitude
                A.rtnpolcentre(2) = model_geometry[j][2] * PI / 180.0;     // Heliocentric latitude
                A.rtnobspolcentre(0) = model_geometry[j][3];               // Geocentric range
                A.rtnobspolcentre(1) = model_geometry[j][4] * PI / 180.0;  // Geocentric longitude
                A.rtnobspolcentre(2) = model_geometry[j][5] * PI / 180.0;  //Geocentric latitude

                // Initialise average variable so it can be added to and divided later to calculate the average
                A.rtnavg_Vmag() = 0.0;
                A.rtnavg_phase() = 0.0;

                // For each step around one rotation
                for (int step = 0; step != NTAU; ++step) {
                    // Calculates illumination of each facet on the model (passing in a, b, c1, c2 for model scaling)
                    lightcurvegeometry(A, step);

                    // Calculate total asteroid flux at geometry at given rotation phase
                    hapke_observe(A, step, diameter, w, B0, h, g, theta);

                    // Sum all the Vmags at each step so the average can be found later
                    A.rtnavg_Vmag() = A.rtnavg_Vmag() + A.rtnVmag(step);
                    A.rtnavg_phase() = A.rtnavg_phase() + A.rtnphase(step);
                    
                    if (step == 0) {
                        // Min and max over lc are the same when it is the first step
                        A.rtnmin_Vmag() = A.rtnVmag(step);
                        A.rtnmax_Vmag() = A.rtnVmag(step);
                    } else {
                        // Check the vMag to see if it is smaller than the current minimum or larger than the current maximum, if so then make it the new minimum/maximum
                        if (A.rtnVmag(step) < A.rtnmin_Vmag()) {
                            A.rtnmin_Vmag() = A.rtnVmag(step);
                        }
                        if (A.rtnVmag(step) > A.rtnmax_Vmag()) {
                            A.rtnmax_Vmag() = A.rtnVmag(step);
                        }
                    }
                }

                // Divide the total Vmag over all light curve steps by the number of light curve steps
                // Gets the average Vmag over the lightcurve
                A.rtnavg_Vmag() = A.rtnavg_Vmag() / NTAU;
                A.rtnavg_phase() = A.rtnavg_phase() / NTAU;

                if (TAXONS[tax] == "C") {
                    Ccurves << i << "\t" << A.rtnavg_phase() << "\t" << A.rtnavg_Vmag() << endl;
                } else {
                    Scurves << i << "\t" << A.rtnavg_phase() << "\t" << A.rtnavg_Vmag() << endl;
                }
            }
        }
    }
};