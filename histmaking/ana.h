#ifndef ana_h
#define ana_h

#include <TROOT.h>
#include "TClonesArray.h"
#include <TLorentzVector.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <vector>
#include "/sphenix/user/jpark4/Utility/Style_jaebeom.h"
#include "/sphenix/user/jpark4/Utility/commonUtility.h"

class ana {
  public :
    ana();
    ~ana();

    static Double_t crystalBall(Double_t *x, Double_t *par);
    virtual std::unordered_set<int> getSubsetIndices(int gridSize, int subsetSize);
    virtual Bool_t   PassEventSelection(float _vz, float _eta, float _pt, int _nclusters, bool _isconv, bool _truth_found_decay); 
    virtual void     SetTruthDecayFlag(bool truth_decay_flag){m_truth_found_decay = truth_decay_flag;}
    virtual Double_t GetShiftedEta(float _vz, float _eta);
    virtual Float_t  deltaR(float eta1, float eta2, float phi1, float phi2);

    float sPHENIX_posx = 0.6;
    float sPHENIX_posy = 0.85;
    float posy_diff = 0.05;

    static const int nDims=7;
    static const int NHIST = nDims*nDims;
    const int gridSize = 7;
    const double highptcut = 40; 
    const double lowptcut = 1; 
    const double vzcut = 30; 
    const double radius = 93; 
    static constexpr double etamin = -1;
    static constexpr double etamax = 1;

    const double dRcut = 0.05;
    const double erecotruthcut = 0.8;

    const double pi0mass = 0.135;
    const double etamass = 0.55;
    const double effpi0masslow = 0.05;
    const double effpi0masshigh = 0.2;
    const double effetamasslow = 0.4;
    const double effetamasshigh = 0.70;

    static const int nPtBins = 28;
    static constexpr double ptBins[nPtBins+1] = {2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 11, 12, 13, 14, 16, 18, 20, 22, 26, 30, 35, 40};
    virtual Int_t findBin(double value);

  private:
    bool m_truth_found_decay{true};
};

#endif

#ifdef ana_cxx

Double_t ana::crystalBall(Double_t *x, Double_t *par) {
  Double_t xx = x[0];
  Double_t mu = par[0];    // Mean
  Double_t sigma = par[1]; // Width (sigma)
  Double_t alpha = par[2]; // Transition point
  Double_t n = par[3];     // Power-law tail parameter
  Double_t norm1 = par[4];  // Normalization
  Double_t sigma2 = par[5];
  Double_t norm2 = par[6];

  Double_t t = (xx - mu) / sigma;
  Double_t cb;
  if (t > -alpha) {
    // Gaussian part
    cb = norm1 * exp(-0.5 *t*t);
  } else {
    // Power-law part
    Double_t a = pow(n / alpha, n) * exp(-0.5 * alpha * alpha);
    Double_t b = n / alpha - alpha;
    cb = norm1 * a * pow(b - t, -n);
  }

  Double_t t2 = (xx - mu) / sigma2;
  Double_t cb2;
  if (t2 > -alpha) {
    // Gaussian part
    cb2 = norm2 * exp(-0.5 *t2*t2);
  } else {
    // Power-law part
    Double_t a = pow(n / alpha, n) * exp(-0.5 * alpha * alpha);
    Double_t b = n / alpha - alpha;
    cb2 = norm2 * a * pow(b - t2, -n);
  }

  Double_t gauss = norm2 * exp(-0.5 * pow((xx - mu) / sigma2, 2));

  return cb + cb2;

}

std::unordered_set<int> ana::getSubsetIndices(int gridSize, int subsetSize) {
  std::unordered_set<int> indices;
  if (subsetSize > gridSize || subsetSize % 2 == 0) {
    std::cerr << "Invalid subset size!" << std::endl;
    return indices;
  }

  int center = (gridSize * gridSize - 1) / 2;
  int halfSubset = subsetSize / 2;

  for (int i = -halfSubset; i <= halfSubset; ++i) {
    for (int j = -halfSubset; j <= halfSubset; ++j) {
      int row = center / gridSize + i; 
      int col = center % gridSize + j; 
      if (row >= 0 && row < gridSize && col >= 0 && col < gridSize) {
        indices.insert(row * gridSize + col);
      }
    }
  }
  return indices;
}


#endif // #ifdef ana_cxx
