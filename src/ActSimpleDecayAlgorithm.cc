
#include "Activia/ActSimpleDecayAlgorithm.hh"
#include "Activia/ActNuclide.hh"
#include "Activia/ActProdNuclide.hh"
#include "Activia/ActTargetNuclide.hh"
#include "Activia/ActOutputSelection.hh"
#include "Activia/ActString.hh"
#include "Activia/ActAbsOutput.hh"   // pour outputLineOfText / outputTable
#include "Activia/ActTarget.hh"      // (nécessaire pour ActTarget*)
#include "Activia/ActOutputTable.hh" // (nécessaire pour ActOutputTable)
#include <iostream>
#include <cmath>                     // pour std::exp

using std::cout;
using std::endl;

ActSimpleDecayAlgorithm::ActSimpleDecayAlgorithm(ActTarget* target,
                                                 ActProdNuclideList* prodList,
                                                 ActTime* times)
  : ActAbsDecayAlgorithm(target, prodList, times),
    _chains() {}

ActSimpleDecayAlgorithm::~ActSimpleDecayAlgorithm() {}

void ActSimpleDecayAlgorithm::writeOutputPreamble() {
  ActString line("Using Bateman 2-member decay algorithm (k->j) with short-lived side-branches folded in.");
  if (_outputData) _outputData->outputLineOfText(line);
}

std::map<ActNuclide*, double> ActSimpleDecayAlgorithm::buildProductionRates() {
  std::map<ActNuclide*, double> R;

  // Somme des taux de production sur tous les isotopes cibles
  std::vector<ActProdXSecData*> allXSec = _theTarget->getXSections();
  for (ActProdXSecData* px : allXSec) {
    if (!px) continue;
    auto xmap = px->getXSecData();   // map<ActNuclide*, ActXSecGraph>
    for (auto &kv : xmap) {
      double r = kv.second.getTotalProdRate();
      R[kv.first] += r;
    }
  }
  return R;
}

void ActSimpleDecayAlgorithm::addShortLivedSideBranchesTo(std::map<ActNuclide*, double>& R) {
  const int nProd = _theProdList->getNProdNuclides();
  for (int i = 0; i < nProd; ++i) {
    ActProdNuclide* p = _theProdList->getProdNuclide(i);
    if (!p) continue;
    ActNuclide* j = p->getProduct();
    if (!j) continue;

    double Rj = R[j];
    const int nsb = p->getNSideBranches();
    for (int s = 0; s < nsb; ++s) {
      ActNuclide* sb = p->getSideBranch(s);
      if (!sb) continue;
      Rj += R[sb];
    }
    R[j] = Rj;
  }
}

void ActSimpleDecayAlgorithm::calculateDecays(ActAbsOutput* output) {
  _outputData = output;

  // Charger les paires parent->fils depuis decayChains.dat
  _chains.load("decayChains.dat", _theProdList);
  _decayEnergies.load("decayEnergy.dat", _theProdList);

  // Temps
  const double texp = _times->getExposureTime();
  const double tdec = _times->getDecayTime();

  // Taux de production et pliage des side-branches
  std::map<ActNuclide*, double> R = buildProductionRates();
  addShortLivedSideBranchesTo(R);

  // Table de sortie 
  std::vector<ActString> cols;
  cols.emplace_back("Z");
  cols.emplace_back("A");
  cols.emplace_back("R_j");
  cols.emplace_back("R_k"); 
  cols.emplace_back("dNdt"); 
  cols.emplace_back("Nj_exp");
  cols.emplace_back("Nk_exp");
  cols.emplace_back("Q_EC_keV");
  cols.emplace_back("Q_beta_keV");
  ActOutputTable table("Bateman2_Yields", cols);

  // Calcul isotope par isotope
  const int nProd = _theProdList->getNProdNuclides();
  for (int i = 0; i < nProd; ++i) {
    ActProdNuclide* p = _theProdList->getProdNuclide(i);
    if (!p) continue;
    ActNuclide* j = p->getProduct();
    if (!j) continue;

    const double lj = lam(j->getHalfLife());
    const double Rj = R[j];

    // Direct : N_j(t_exp)
    const double Nj_dir_exp = (lj > 0.0 ? (Rj / lj) * (1.0 - std::exp(-lj * texp)) : 0.0);

    // Chaîné 
    double Nj_chain_exp_sum = 0.0;
    double Nk_cool_sum      = 0.0;
    double Rk_sum           = 0.0;
    double Nk_texp_sum      = 0.0;

    const auto& parents = _chains.parentsOf(j);
    for (const auto& link : parents) {
      ActNuclide* k = link.parent;
      const double b  = link.branching;
      const double lk = lam(k->getHalfLife());
      const double Rk = R[k];

      Rk_sum += Rk;

      // Expo 
      Nj_chain_exp_sum += Nj_from_k_exposure(Rk, lj, lk, b, texp);

      // N_k(t_exp)
      const double Nk_texp = (lk > 0.0 ? (Rk / lk) * (1.0 - std::exp(-lk * texp)) : 0.0);
      Nk_texp_sum += Nk_texp;

      // Refroidissement 
      Nk_cool_sum += Nj_chain_cooling_add(Nk_texp, lj, lk, b, tdec);
    }

    const double Nj_texp = Nj_dir_exp + Nj_chain_exp_sum;

    // N_j(t_dec)
    const double Nj_tdec = (Nj_dir_exp + Nj_chain_exp_sum) * std::exp(-lj * tdec) + Nk_cool_sum;

    // dNdt = lambda_j * N_j(t_dec)
    const double dNdt = lj * Nj_tdec; 

    std::vector<double> row;
    row.push_back(static_cast<double>(j->getZ()));
    row.push_back(j->getA());
    row.push_back(Rj);
    row.push_back(Rk_sum); 
    row.push_back(dNdt);
    row.push_back(Nj_texp);
    row.push_back(Nk_texp_sum); 
    auto qv = _decayEnergies.getQValues(j);
    row.push_back(qv.QEC);
    row.push_back(qv.Qbeta);
    table.addRow(row);
  }

  // Impression
  if (_outputData) {
    writeOutputPreamble();
    _outputData->outputTable(table);
  }
}
