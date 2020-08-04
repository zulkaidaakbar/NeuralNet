#include <vector>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <fstream>
#include <sstream>
#include "TLorentzVector.h"
#include "TLorentzRotation.h"

using namespace std;

typedef double Double_t;

double TProduct (TLorentzVector v1, TLorentzVector v2) {
	// Transverse product
  Double_t tv1v2;
  return tv1v2 = v1.Px() * v2.Px() + v1.Py() * v2.Py();
}

double model (
    
    const std::vector<double> inputs,
    const std::vector<double> params
)
{


const double k	= inputs[0];
const double QQ	= inputs[1];
const double xb	= inputs[2];
const double t	= inputs[3];
const double phi = inputs[4];
const double F1	= inputs[5];
const double F2	= inputs[6];
const double dvcs = inputs[7];


const double p0 = params[0];
const double p1 = params[1];
const double p2 = params[2];




Double_t ALP_INV = 137.0359998; 
Double_t PI = 3.1415926535;
Double_t RAD = PI / 180.;
Double_t M = 0.938272; 
Double_t GeV2nb = .389379*1000000; 

Double_t chisq = 0.0;
Double_t function = 0.0;
Double_t prob = 0.0;
Double_t nll = 0.0;


Double_t y, e, xi, tmin, kpr, gg, q, qp, po, pmag, cth, theta, sth, sthl, cthl, cthpr, sthpr, M2, tau;

  
TLorentzVector K, KP, Q, QP, D, p, P;
Double_t kkp, kq, kp, kpp;
Double_t kd, kpd, kP, kpP, kqp, kpqp, dd, Pq, Pqp, qd, qpd;
Double_t kk_T, kqp_T, kkp_T, kd_T, dd_T, kpqp_T, kP_T, kpP_T, qpP_T, kpd_T, qpd_T;

Double_t s;     
Double_t Gamma; 
Double_t jcob;  
Double_t AUUBH, BUUBH; 
Double_t AUUI, BUUI, CUUI; 
Double_t con_AUUBH, con_BUUBH, con_AUUI, con_BUUI, con_CUUI;  
Double_t bhAUU, bhBUU; 
Double_t iAUU, iBUU, iCUU; 
Double_t real_iAUU, real_iBUU, real_iCUU; 
Double_t xbhUU;
Double_t xIUU; 
Double_t real_xbhUU;
Double_t real_xIUU; 




M2 = M*M; 
y = QQ / ( 2. * M * k * xb ); 
  
gg = 4. * M2 * xb * xb / QQ; 
e = ( 1 - y - ( y * y * (gg / 4.) ) ) / ( 1. - y + (y * y / 2.) + ( y * y * (gg / 4.) ) ); 
xi = 1. * xb * ( ( 1. + t / ( 2. * QQ ) ) / ( 2. - xb + xb * t / QQ ) ); 
tmin = ( QQ * ( 1. - sqrt( 1. + gg ) + gg / 2. ) ) / ( xb * ( 1. - sqrt( 1. + gg ) + gg / ( 2.* xb ) ) ); 
kpr = k * ( 1. - y ); 

qp = t / 2. / M + k - kpr; 
po = M - t / 2. / M; 
pmag = sqrt( ( -t ) * ( 1. - t / 4. / M / M ) ); 
cth = -1. / sqrt( 1. + gg ) * ( 1. + gg / 2. * ( 1. + t / QQ ) / ( 1. + xb * t / QQ ) ); 
theta = acos(cth); 
sthl = sqrt( gg ) / sqrt( 1. + gg ) * ( sqrt ( 1. - y - y * y * gg / 4. ) ); 
cthl = -1. / sqrt( 1. + gg ) * ( 1. + y * gg / 2. ) ; 
tau = -0.25 * t / M2;
K.SetPxPyPzE( k * sthl, 0.0, k * cthl, k );
KP.SetPxPyPzE( K(0), 0.0, k * ( cthl + y * sqrt( 1. + gg ) ), kpr );
Q = K - KP;
p.SetPxPyPzE(0.0, 0.0, 0.0, M);

s = (p + K) * (p + K);
Gamma = 1. / ALP_INV / ALP_INV / ALP_INV / PI / PI / 16. / ( s - M2 ) / ( s - M2 ) / sqrt( 1. + gg ) / xb;
jcob = 1./ ( 2. * M * xb * K(3) ) * 2. * PI * 2.;

//
QP.SetPxPyPzE(qp * sin(theta) * cos( RAD*phi ), qp * sin(theta) * sin( RAD*phi ), qp * cos(theta), qp);
D = Q - QP; 
TLorentzVector pp = p + D; 
P = p + pp;
P.SetPxPyPzE(.5*P.Px(), .5*P.Py(), .5*P.Pz(), .5*P.E());

//
kkp  = K * KP;   
kq   = K * Q;    
kp   = K * p;    
kpp  = KP * p;   

kd   = K * D;  
kpd  = KP * D;   
kP   = K * P;    
kpP  = KP * P;   
kqp  = K * QP;   
kpqp = KP * QP;  
dd   = D * D;    
Pq   = P * Q;    
Pqp  = P * QP;   
qd   = Q * D;    
qpd  = QP * D;   

//

kk_T   = TProduct(K,K);
  kkp_T  = kk_T;
  kqp_T  = TProduct(K,QP);
  kd_T   = -1.* kqp_T;
  dd_T   = TProduct(D,D);
  kpqp_T = TProduct(KP,QP);
  kP_T   = TProduct(K,P);
  kpP_T  = TProduct(KP,P);
  qpP_T  = TProduct(QP,P);
  kpd_T  = TProduct(KP,D);
  qpd_T  = TProduct(QP,D);


//
AUUBH = ( (8. * M2) / (t * kqp * kpqp) ) * ( (4. * tau * (kP * kP + kpP * kpP) ) - ( (tau + 1.) * (kd * kd + kpd * kpd) ) );
BUUBH = ( (16. * M2) / (t* kqp * kpqp) ) * (kd * kd + kpd * kpd);

con_AUUBH = 2. * AUUBH * GeV2nb * jcob;
con_BUUBH = 2. * BUUBH * GeV2nb * jcob;

bhAUU = (Gamma/t) * con_AUUBH * ( F1 * F1 + tau * F2 * F2 );
bhBUU = (Gamma/t) * con_BUUBH * ( tau * ( F1 + F2 ) * ( F1 + F2 ) ) ;

xbhUU = bhAUU + bhBUU;
//


AUUI = -4.0 / (kqp * kpqp) * ( ( QQ + t ) * ( 2.0 * ( kP + kpP ) * kk_T   + ( Pq * kqp_T ) + 2.* ( kpP * kqp ) - 2.* ( kP * kpqp ) + kpqp * kP_T + kqp * kpP_T - 2.*kkp * kP_T ) +
                                                    ( QQ - t + 4.* kd ) * ( Pqp * ( kkp_T + kqp_T - 2.* kkp )  + 2.* kkp * qpP_T - kpqp * kP_T - kqp * kpP_T ) );
  // BUUI = 2.0 * xi / ( kqp * kpqp) * ( ( QQ + t ) * ( 2.* kk_T * ( kd + kpd ) + kqp_T * ( qd - kqp - kpqp + 2.*kkp ) + 2.* kqp * kpd - 2.* kpqp * kd ) +
  //                                                        ( QQ - t + 4.* kd ) * ( ( kk_T - 2.* kkp ) * qpd - kkp * dd_T - 2.* kd_T * kqp ) ) / tau;
  BUUI = 2.0 * xi / ( kqp * kpqp) * ( ( QQ + t ) * ( 2.* kk_T * ( kd + kpd ) + kqp_T * ( qd - kqp - kpqp + 2.*kkp ) + 2.* kqp * kpd - 2.* kpqp * kd ) +
                                                    ( QQ - t + 4.* kd ) * ( ( kk_T - 2.* kkp ) * qpd - kkp * dd_T - 2.* kd_T * kqp ) );
  CUUI = 2.0 / ( kqp * kpqp) * ( -1. * ( QQ + t ) * ( 2.* kkp - kpqp - kqp + 2.* xi * (2.* kkp * kP_T - kpqp * kP_T - kqp * kpP_T) ) * kd_T +
                                                  ( QQ - t + 4.* kd ) * ( ( kqp + kpqp ) * kd_T + dd_T * kkp + 2.* xi * ( kkp * qpP_T - kpqp * kP_T - kqp * kpP_T ) ) );


con_AUUI = AUUI * GeV2nb * jcob;
con_BUUI = BUUI * GeV2nb * jcob;
con_CUUI = CUUI * GeV2nb * jcob;



iAUU = (Gamma/(-t * QQ)) * cos( phi * RAD ) * con_AUUI * ( F1 * p1 + tau * F2 * p0 );
  
  iBUU = (Gamma/(-t * QQ)) * cos( phi * RAD ) * con_BUUI * ( F1 + F2 ) * ( p1 + p0 );
  iCUU = (Gamma/(-t * QQ)) * cos( phi * RAD ) * con_CUUI * ( F1 + F2 ) * p2;

xIUU = iAUU + iBUU + iCUU ;

function = xIUU + xbhUU + dvcs ;

return function;
    
}

// ----------------------------------------------------------------------------------------

double output_delta_derivative (
    
    const std::vector<double> inputs,
    const std::vector<double> params,
    double targetVals, 
    int n
)
{


const double k	= inputs[0];
const double QQ	= inputs[1];
const double xb	= inputs[2];
const double t	= inputs[3];
const double phi = inputs[4];
const double F1	= inputs[5];
const double F2	= inputs[6];
const double dvcs = inputs[7];


const double p0 = params[0];
const double p1 = params[1];
const double p2 = params[2];

double F = targetVals;
int par_int = n;




Double_t ALP_INV = 137.0359998; 
Double_t PI = 3.1415926535;
Double_t RAD = PI / 180.;
Double_t M = 0.938272; 
Double_t GeV2nb = .389379*1000000; 

Double_t chisq = 0.0;
Double_t function = 0.0;
Double_t prob = 0.0;
Double_t nll = 0.0;


Double_t y, e, xi, tmin, kpr, gg, q, qp, po, pmag, cth, theta, sth, sthl, cthl, cthpr, sthpr, M2, tau;

  
TLorentzVector K, KP, Q, QP, D, p, P;
Double_t kkp, kq, kp, kpp;
Double_t kd, kpd, kP, kpP, kqp, kpqp, dd, Pq, Pqp, qd, qpd;
//Double_t kk_T, kqp_T, kkp_T, kd_T, dd_T;
Double_t kk_T, kqp_T, kkp_T, kd_T, dd_T, kpqp_T, kP_T, kpP_T, qpP_T, kpd_T, qpd_T;

Double_t s;     
Double_t Gamma; 
Double_t jcob;  
Double_t AUUBH, BUUBH; 
Double_t AUUI, BUUI, CUUI; 
Double_t con_AUUBH, con_BUUBH, con_AUUI, con_BUUI, con_CUUI;  
Double_t bhAUU, bhBUU; 
Double_t iAUU, iBUU, iCUU; 
Double_t real_iAUU, real_iBUU, real_iCUU; 
Double_t xbhUU;
Double_t xIUU; 
Double_t real_xbhUU;
Double_t real_xIUU; 




M2 = M*M; 
y = QQ / ( 2. * M * k * xb ); 
  
gg = 4. * M2 * xb * xb / QQ; 
e = ( 1 - y - ( y * y * (gg / 4.) ) ) / ( 1. - y + (y * y / 2.) + ( y * y * (gg / 4.) ) ); 
xi = 1. * xb * ( ( 1. + t / ( 2. * QQ ) ) / ( 2. - xb + xb * t / QQ ) ); 
tmin = ( QQ * ( 1. - sqrt( 1. + gg ) + gg / 2. ) ) / ( xb * ( 1. - sqrt( 1. + gg ) + gg / ( 2.* xb ) ) ); 
kpr = k * ( 1. - y ); 

qp = t / 2. / M + k - kpr; 
po = M - t / 2. / M; 
pmag = sqrt( ( -t ) * ( 1. - t / 4. / M / M ) ); 
cth = -1. / sqrt( 1. + gg ) * ( 1. + gg / 2. * ( 1. + t / QQ ) / ( 1. + xb * t / QQ ) ); 
theta = acos(cth); 
sthl = sqrt( gg ) / sqrt( 1. + gg ) * ( sqrt ( 1. - y - y * y * gg / 4. ) ); 
cthl = -1. / sqrt( 1. + gg ) * ( 1. + y * gg / 2. ) ; 
tau = -0.25 * t / M2;
K.SetPxPyPzE( k * sthl, 0.0, k * cthl, k );
KP.SetPxPyPzE( K(0), 0.0, k * ( cthl + y * sqrt( 1. + gg ) ), kpr );
Q = K - KP;
p.SetPxPyPzE(0.0, 0.0, 0.0, M);

s = (p + K) * (p + K);
Gamma = 1. / ALP_INV / ALP_INV / ALP_INV / PI / PI / 16. / ( s - M2 ) / ( s - M2 ) / sqrt( 1. + gg ) / xb;
jcob = 1./ ( 2. * M * xb * K(3) ) * 2. * PI * 2.;

//
QP.SetPxPyPzE(qp * sin(theta) * cos( RAD*phi ), qp * sin(theta) * sin( RAD*phi ), qp * cos(theta), qp);
D = Q - QP; 
TLorentzVector pp = p + D; 
P = p + pp;
P.SetPxPyPzE(.5*P.Px(), .5*P.Py(), .5*P.Pz(), .5*P.E());

//
kkp  = K * KP;   
kq   = K * Q;    
kp   = K * p;    
kpp  = KP * p;   

kd   = K * D;  
kpd  = KP * D;   
kP   = K * P;    
kpP  = KP * P;   
kqp  = K * QP;   
kpqp = KP * QP;  
dd   = D * D;    
Pq   = P * Q;    
Pqp  = P * QP;   
qd   = Q * D;    
qpd  = QP * D;   

//
/*
kk_T = 0.5 * ( e / ( 1 - e ) ) * QQ;  
kkp_T = kk_T;  
kqp_T = ( QQ / ( sqrt( gg ) * sqrt( 1 + gg ) ) ) * sqrt ( (0.5 * e) / ( 1 - e ) ) * ( 1. + xb * t / QQ ) * sin(theta) * cos( RAD*phi );
kd_T = -1.* kqp_T;
dd_T = ( 1. - xi * xi ) * ( tmin - t );*/
kk_T   = TProduct(K,K);
  kkp_T  = kk_T;
  kqp_T  = TProduct(K,QP);
  kd_T   = -1.* kqp_T;
  dd_T   = TProduct(D,D);
  kpqp_T = TProduct(KP,QP);
  kP_T   = TProduct(K,P);
  kpP_T  = TProduct(KP,P);
  qpP_T  = TProduct(QP,P);
  kpd_T  = TProduct(KP,D);
  qpd_T  = TProduct(QP,D);


//
AUUBH = ( (8. * M2) / (t * kqp * kpqp) ) * ( (4. * tau * (kP * kP + kpP * kpP) ) - ( (tau + 1.) * (kd * kd + kpd * kpd) ) );
BUUBH = ( (16. * M2) / (t* kqp * kpqp) ) * (kd * kd + kpd * kpd);

con_AUUBH = 2. * AUUBH * GeV2nb * jcob;
con_BUUBH = 2. * BUUBH * GeV2nb * jcob;

bhAUU = (Gamma/t) * con_AUUBH * ( F1 * F1 + tau * F2 * F2 );
bhBUU = (Gamma/t) * con_BUUBH * ( tau * ( F1 + F2 ) * ( F1 + F2 ) ) ;

xbhUU = bhAUU + bhBUU;
//
/*
AUUI = -4.0 * cos( RAD*phi ) / (kqp * kpqp) * ( ( QQ + t ) * ( 2.0 * ( kP + kpP ) * kk_T   + ( Pq * kqp_T ) + 2.* ( kpP * kqp ) - 2.* ( kP * kpqp ) ) + ( QQ - t + 4.* kd ) * Pqp * ( kkp_T + kqp_T - 2.* kkp ) );
BUUI = 2.0 * xi * cos( RAD*phi ) / ( kqp * kpqp) * ( ( QQ + t ) * ( 2.* kk_T * ( kd + kpd ) + kqp_T * ( qd - kqp - kpqp + 2.*kkp ) + 2.* kqp * kpd - 2.* kpqp * kd ) + ( QQ - t + 4.* kd ) * ( ( kk_T - 2.* kkp ) * qpd - kkp * dd_T - 2.* kd_T * kqp ) ) / tau;
CUUI = 2.0 * cos( RAD*phi ) / ( kqp * kpqp) * ( -1. * ( QQ + t ) * ( 2.* kkp - kpqp - kqp ) * kd_T + ( QQ - t + 4.* kd ) * ( ( kqp + kpqp ) * kd_T + dd_T * kkp ) );
*/

AUUI = -4.0 / (kqp * kpqp) * ( ( QQ + t ) * ( 2.0 * ( kP + kpP ) * kk_T   + ( Pq * kqp_T ) + 2.* ( kpP * kqp ) - 2.* ( kP * kpqp ) + kpqp * kP_T + kqp * kpP_T - 2.*kkp * kP_T ) +
                                                    ( QQ - t + 4.* kd ) * ( Pqp * ( kkp_T + kqp_T - 2.* kkp )  + 2.* kkp * qpP_T - kpqp * kP_T - kqp * kpP_T ) );
  // BUUI = 2.0 * xi / ( kqp * kpqp) * ( ( QQ + t ) * ( 2.* kk_T * ( kd + kpd ) + kqp_T * ( qd - kqp - kpqp + 2.*kkp ) + 2.* kqp * kpd - 2.* kpqp * kd ) +
  //                                                        ( QQ - t + 4.* kd ) * ( ( kk_T - 2.* kkp ) * qpd - kkp * dd_T - 2.* kd_T * kqp ) ) / tau;
  BUUI = 2.0 * xi / ( kqp * kpqp) * ( ( QQ + t ) * ( 2.* kk_T * ( kd + kpd ) + kqp_T * ( qd - kqp - kpqp + 2.*kkp ) + 2.* kqp * kpd - 2.* kpqp * kd ) +
                                                    ( QQ - t + 4.* kd ) * ( ( kk_T - 2.* kkp ) * qpd - kkp * dd_T - 2.* kd_T * kqp ) );
  CUUI = 2.0 / ( kqp * kpqp) * ( -1. * ( QQ + t ) * ( 2.* kkp - kpqp - kqp + 2.* xi * (2.* kkp * kP_T - kpqp * kP_T - kqp * kpP_T) ) * kd_T +
                                                  ( QQ - t + 4.* kd ) * ( ( kqp + kpqp ) * kd_T + dd_T * kkp + 2.* xi * ( kkp * qpP_T - kpqp * kP_T - kqp * kpP_T ) ) );


con_AUUI = AUUI * GeV2nb * jcob;
con_BUUI = BUUI * GeV2nb * jcob;
con_CUUI = CUUI * GeV2nb * jcob;

/*iAUU = (Gamma/(-t * QQ)) * cos( RAD*phi ) * con_AUUI * ( F1 * p1 + tau * F2 * p0 );
iBUU = (Gamma/(-t * QQ)) * cos( RAD*phi ) * con_BUUI * tau * ( F1 + F2 ) * ( p1 + p0 );
iCUU = (Gamma/(-t * QQ)) * cos( RAD*phi ) * con_CUUI * ( F1 + F2 ) * p2;*/

iAUU = (Gamma/(-t * QQ)) * cos( phi * RAD ) * con_AUUI * ( F1 * p1 + tau * F2 * p0 );
  //iBUU = (Gamma/(-t * QQ)) * cos( phi * RAD ) * con_BUUI * tau * ( F1 + F2 ) * ( ReH + ReE );
  iBUU = (Gamma/(-t * QQ)) * cos( phi * RAD ) * con_BUUI * ( F1 + F2 ) * ( p1 + p0 );
  iCUU = (Gamma/(-t * QQ)) * cos( phi * RAD ) * con_CUUI * ( F1 + F2 ) * p2;

xIUU = iAUU + iBUU + iCUU ;

function = xIUU + xbhUU + dvcs ;

double deltaF = function - F;

double der_par0 = (Gamma/(-t * QQ)) * cos( phi * RAD ) * con_AUUI * ( tau * F2 ) + (Gamma/(-t * QQ)) * cos( phi * RAD ) * con_BUUI * ( F1 + F2 );
double der_par1 = (Gamma/(-t * QQ)) * cos( phi * RAD ) * con_AUUI * ( F1  ) + (Gamma/(-t * QQ)) * cos( phi * RAD ) * con_BUUI * ( F1 + F2 );
double der_par2 = (Gamma/(-t * QQ)) * cos( phi * RAD ) * con_CUUI * ( F1 + F2 ) ;

double derivative ;
if (n==0) {derivative = der_par0;}
if (n==1) {derivative = der_par1;}
if (n==2) {derivative = der_par2;}

double function_all;
function_all = (deltaF)*derivative;

return function_all;
    
}

// ----------------------------------------------------------------------------------------

class TrainingData
{
public:
	TrainingData(const string filename);
	bool isEof(void)
	{
		return m_trainingDataFile.eof();
	}
	void getTopology(vector<unsigned> &topology);

	// Returns the number of input values read from the file:
	unsigned getNextInputs(vector<double> &inputVals);
	unsigned getTargetOutputs(vector<double> &targetOutputVals);

private:
	ifstream m_trainingDataFile;
};

void TrainingData::getTopology(vector<unsigned> &topology)
{
	string line;
	string label;

	getline(m_trainingDataFile, line);
	stringstream ss(line);
	ss >> label;
	if(this->isEof() || label.compare("topology:") != 0)
	{
		abort();
	}

	while(!ss.eof())
	{
		unsigned n;
		ss >> n;
		topology.push_back(n);
	}
	return;
}

TrainingData::TrainingData(const string filename)
{
	m_trainingDataFile.open(filename.c_str());
}


unsigned TrainingData::getNextInputs(vector<double> &inputVals)
{
    inputVals.clear();

    string line;
    getline(m_trainingDataFile, line);
    stringstream ss(line);

    string label;
    ss >> label;
    if (label.compare("in:") == 0) {
        double oneValue;
        while (ss >> oneValue) {
            inputVals.push_back(oneValue);
        }
    }

    return inputVals.size();
}

unsigned TrainingData::getTargetOutputs(vector<double> &targetOutputVals)
{
    targetOutputVals.clear();

    string line;
    getline(m_trainingDataFile, line);
    stringstream ss(line);

    string label;
    ss>> label;
    if (label.compare("out:") == 0) {
        double oneValue;
        while (ss >> oneValue) {
            targetOutputVals.push_back(oneValue);
        }
    }

    return targetOutputVals.size();
}

struct Connection
{
	double weight;
	double deltaWeight;
};

class Neuron;

typedef vector<Neuron> Layer;

// ****************** class Neuron ******************

class Neuron
{
public:
	Neuron(unsigned numOutputs, unsigned myIndex);
	void setOutputVal(double val) { m_outputVal = val; }
	double getOutputVal(void) const { return m_outputVal; }
        //vector<double> getOutputValVec(void) const { return m_outputValVec;}
	void feedForward(const Layer &prevLayer);
	void feedForwardOutput(const Layer &prevLayer);
	void calcOutputGradients(const std::vector<double> inputs, const std::vector<double> params, double targetVals, int n);
	void calcHiddenGradients(const Layer &nextLayer);
	void updateInputWeights(Layer &prevLayer);
private:
	static double eta; // [0.0...1.0] overall net training rate
	static double alpha; // [0.0...n] multiplier of last weight change [momentum]
	static double transferFunction(double x);
	static double transferFunctionDerivative(double x);
	static double transferFunctionOutput(double x);
	static double transferFunctionDerivativeOutput(double x);
	// randomWeight: 0 - 1
	static double randomWeight(void) { return rand() / double(RAND_MAX); }
	//static double randomWeight(void) { return rand() % 31 - 15; }
	double sumDOW(const Layer &nextLayer) const;
	double m_outputVal;
	//vector<double> m_outputValVec;
	vector<Connection> m_outputWeights;
	unsigned m_myIndex;
	double m_gradient;
};

double Neuron::eta = 0.05; // overall net learning rate
double Neuron::alpha = 0.15; // momentum, multiplier of last deltaWeight, [0.0..n]

//***************************************************************************************************************


void Neuron::updateInputWeights(Layer &prevLayer)
{
	// The weights to be updated are in the Connection container
	// in the nuerons in the preceding layer

	for(unsigned n = 0; n < prevLayer.size(); ++n)
	{
		Neuron &neuron = prevLayer[n];
		double oldDeltaWeight = neuron.m_outputWeights[m_myIndex].deltaWeight;

		double newDeltaWeight = 
				// Individual input, magnified by the gradient and train rate:
				eta
				* neuron.getOutputVal()
				* m_gradient
				// Also add momentum = a fraction of the previous delta weight
				+ alpha
				* oldDeltaWeight;
		neuron.m_outputWeights[m_myIndex].deltaWeight = newDeltaWeight;
		//neuron.m_outputWeights[m_myIndex].weight += newDeltaWeight;
		neuron.m_outputWeights[m_myIndex].weight = neuron.m_outputWeights[m_myIndex].weight - newDeltaWeight;
	}
}

//***************************************************************************************************************
double Neuron::sumDOW(const Layer &nextLayer) const
{
	double sum = 0.0;

	// Sum our contributions of the errors at the nodes we feed

	for (unsigned n = 0; n < nextLayer.size() - 1; ++n)
	{
		sum += m_outputWeights[n].weight * nextLayer[n].m_gradient;
	}

	return sum;
}

void Neuron::calcHiddenGradients(const Layer &nextLayer)
{
	double dow = sumDOW(nextLayer);
	m_gradient = dow * Neuron::transferFunctionDerivative(m_outputVal);
}
void Neuron::calcOutputGradients(const std::vector<double> inputs, const std::vector<double> params, double targetVals, int n)
{
	
	m_gradient = output_delta_derivative(inputs, params, targetVals, n) * Neuron::transferFunctionDerivativeOutput(m_outputVal);
}

//--------------------------------------------------------
double Neuron::transferFunction(double x)
{
	// tanh - output range [-1.0..1.0]
	return tanh(x);
}

double Neuron::transferFunctionDerivative(double x)
{
	// tanh derivative
	return 1.0 - x * x;
}

double Neuron::transferFunctionOutput(double x)
{
	
	return x;
}

double Neuron::transferFunctionDerivativeOutput(double x)
{
	
	return 1.0 ;
}

//-----------------------------------------------------------

void Neuron::feedForward(const Layer &prevLayer)
{
	double sum = 0.0;

	// Sum the previous layer's outputs (which are our inputs)
    // Include the bias node from the previous layer.

	for(unsigned n = 0 ; n < prevLayer.size(); ++n)
	{
		sum += prevLayer[n].getOutputVal() * 
				 prevLayer[n].m_outputWeights[m_myIndex].weight;
	}

	m_outputVal = Neuron::transferFunction(sum);
}

void Neuron::feedForwardOutput(const Layer &prevLayer)
{
	double sum = 0.0;

	// Sum the previous layer's outputs (which are our inputs)
    // Include the bias node from the previous layer.

	for(unsigned n = 0 ; n < prevLayer.size(); ++n)
	{
		sum += prevLayer[n].getOutputVal() * 
				 prevLayer[n].m_outputWeights[m_myIndex].weight;
	}

	m_outputVal = Neuron::transferFunctionOutput(sum);
}


Neuron::Neuron(unsigned numOutputs, unsigned myIndex)
{
	for(unsigned c = 0; c < numOutputs; ++c){
		m_outputWeights.push_back(Connection());
		m_outputWeights.back().weight = randomWeight();
	}

	m_myIndex = myIndex;
}
// ****************** class Net ******************
class Net
{
public:
	Net(const vector<unsigned> &topology);
	void feedForward(const vector<double> &inputVals);
	void backProp(const vector<double> &inputVals, const vector<double> &targetVals, const vector<double> &resultVals);
	void getResults(vector<double> &resultVals) const;
	double getRecentAverageError(void) const { return m_recentAverageError; }

private:
	vector<Layer> m_layers; //m_layers[layerNum][neuronNum]
	double m_error;
	double m_recentAverageError;
	static double m_recentAverageSmoothingFactor;
};

double Net::m_recentAverageSmoothingFactor = 100.0; // Number of training samples to average over

void Net::getResults(vector<double> &resultVals) const
{
	resultVals.clear();

	for(unsigned n = 0; n < m_layers.back().size() - 1; ++n)
	{
		resultVals.push_back(m_layers.back()[n].getOutputVal());
	}
}

void Net::backProp(const std::vector<double> &inputVals, const std::vector<double> &targetVals, const std::vector<double> &resultVals)
{
	// Calculate overal net error (RMS of output neuron errors)
	vector<double> inputs = inputVals;
	vector<double> F_iuuxs = targetVals;
	vector<double> params = resultVals;
	Layer &outputLayer = m_layers.back();
	m_error = 0.0;

	for(unsigned n = 0; n < outputLayer.size() - 1; ++n)
	{
		double delta = F_iuuxs[n] - model(inputs,params);
		m_error += delta *delta;
	}
        
	m_error /= outputLayer.size() - 1; // get average error squared
	m_error = sqrt(m_error); // RMS

	// Implement a recent average measurement:

	m_recentAverageError = 
			(m_recentAverageError * m_recentAverageSmoothingFactor + m_error)
			/ (m_recentAverageSmoothingFactor + 1.0);
	// Calculate output layer gradients

	for(unsigned n = 0; n < outputLayer.size() - 1; ++n)
	{
		
		outputLayer[n].calcOutputGradients(inputs, params, targetVals[n], n);
	}
	// Calculate gradients on hidden layers

	for(unsigned layerNum = m_layers.size() - 2; layerNum > 0; --layerNum)
	{
		Layer &hiddenLayer = m_layers[layerNum];
		Layer &nextLayer = m_layers[layerNum + 1];

		for(unsigned n = 0; n < hiddenLayer.size(); ++n)
		{
			hiddenLayer[n].calcHiddenGradients(nextLayer);
		}
	}

	// For all layers from outputs to first hidden layer,
	// update connection weights

	for(unsigned layerNum = m_layers.size() - 1; layerNum > 0; --layerNum)
	{
		Layer &layer = m_layers[layerNum];
		Layer &prevLayer = m_layers[layerNum - 1];

		for(unsigned n = 0; n < layer.size() - 1; ++n)
		{
			layer[n].updateInputWeights(prevLayer);
		}
	}
}

void Net::feedForward(const vector<double> &inputVals)
{
	// Check the num of inputVals euqal to neuronnum expect bias
	assert(inputVals.size() == m_layers[0].size() - 1);

	// Assign {latch} the input values into the input neurons
	for(unsigned i = 0; i < inputVals.size(); ++i){
		m_layers[0][i].setOutputVal(inputVals[i]); 
	}

	// Forward propagate
	for(unsigned layerNum = 1; layerNum < m_layers.size(); ++layerNum){
		Layer &prevLayer = m_layers[layerNum - 1];
		for(unsigned n = 0; n < m_layers[layerNum].size() - 1; ++n){
			if (layerNum == (m_layers.size()-1))
			 {
			  m_layers[layerNum][n].feedForwardOutput(prevLayer);
			 }
			else
			 {
			  m_layers[layerNum][n].feedForward(prevLayer);
			 }
		}
	}
}
Net::Net(const vector<unsigned> &topology)
{
	unsigned numLayers = topology.size();
	for(unsigned layerNum = 0; layerNum < numLayers; ++layerNum){
		m_layers.push_back(Layer());
		// numOutputs of layer[i] is the numInputs of layer[i+1]
		// numOutputs of last layer is 0
		unsigned numOutputs = layerNum == topology.size() - 1 ? 0 :topology[layerNum + 1];

		// We have made a new Layer, now fill it ith neurons, and
		// add a bias neuron to the layer:
		for(unsigned neuronNum = 0; neuronNum <= topology[layerNum]; ++neuronNum){
			m_layers.back().push_back(Neuron(numOutputs, neuronNum));
			cout << "Made a Neuron!" << endl;
		}

		// Force the bias node's output value to 1.0. It's the last neuron created above
		m_layers.back().back().setOutputVal(1.0);
	}
}

void showVectorVals(string label, vector<double> &v)
{
	cout << label << " ";
	for(unsigned i = 0; i < v.size(); ++i)
	{
		cout << v[i] << " ";
	}
	cout << endl;
}
int main()
{
	TrainingData trainData("epoch_test7.txt");
	//e.g., {3, 2, 1 }
	vector<unsigned> topology;
	

	trainData.getTopology(topology);
	Net myNet(topology);

	vector<double> inputVals, targetVals, resultVals;
	int trainingPass = 0;
        
	while(!trainData.isEof())
	{

		++trainingPass;
		cout << endl << "Pass" << trainingPass;

		// Get new input data and feed it forward:
		if(trainData.getNextInputs(inputVals) != topology[0])
			break;
		showVectorVals(": Inputs :", inputVals);
//for (int i =0; i<1000; i++)
//         {
		myNet.feedForward(inputVals);

		// Collect the net's actual results:
		myNet.getResults(resultVals);
		showVectorVals("Outputs:", resultVals);

		// Train the net what the outputs should have been:
		trainData.getTargetOutputs(targetVals);
		//showVectorVals("Targets:", targetVals);
		assert(targetVals.size() == topology.back());

		myNet.backProp(inputVals, targetVals, resultVals);

		// Report how well the training is working, average over recnet
		//cout << "Net recent average error: "
		//     << myNet.getRecentAverageError() << endl;
		//cout <<resultVals[0]<<" "<<resultVals[1]<<" "<<resultVals[2]<<endl;
		cout<<model(inputVals, resultVals) - targetVals[0]<<" "<<model(inputVals, resultVals)<<endl;
                //cout<<output_delta_derivative (inputVals, resultVals,targetVals[0],0)<<endl; 
		//cout<<output_delta_derivative (inputVals, resultVals,targetVals[1],1)<<endl; 
		//cout<<output_delta_derivative (inputVals, resultVals,targetVals[2],2)<<endl; 
    
//	   }
	}

	cout << endl << "Done" << endl;

}
