#define ana_cxx
#include "ana.h"
ana::ana() {
}

ana::~ana() {
  // destructor implementation
}


Bool_t ana::PassEventSelection(float _vz, float _eta, float _pt, int _nclusters, bool _isconv, bool _truth_found_decay)
{
  bool ispass = true;
  if(fabs(_vz) > vzcut) ispass = false;
  if(_pt > highptcut || _pt < lowptcut) ispass = false;
  if(_nclusters == 0) ispass =false;
  if(_isconv) ispass = false;
  if(m_truth_found_decay && !_truth_found_decay) ispass = false;
  double etaminshifted = GetShiftedEta(_vz, etamin);
  double etamaxshifted = GetShiftedEta(_vz, etamax);
  if(_eta > etamaxshifted || _eta < etaminshifted) ispass = false;
  return ispass;
}

Double_t ana::GetShiftedEta(float _vz, float _eta)
{
  double theta = 2*atan(exp(-_eta));
  double z = radius / tan(theta);
  double zshifted = z - _vz;
  double thetashifted = atan2(radius,zshifted);
  double etashifted = -log(tan(thetashifted/2.0));
  return etashifted;
}

Float_t ana::deltaR(float eta1, float eta2, float phi1, float phi2)
{
  float dphi = (fabs(phi1-phi2) > M_PI) ? 2*M_PI - fabs(phi1-phi2) : fabs(phi1-phi2);
  float dR = sqrt((eta2-eta1) * (eta2-eta1) + dphi*dphi);
  return dR;
}

Int_t ana::findBin(double value)
{
  for (int i = 0; i < nPtBins; ++i) {
    if (value >= ptBins[i] && value < ptBins[i + 1]) {
      return  i;
    }
  }
  return -1;
}

