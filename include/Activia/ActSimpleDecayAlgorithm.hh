
#ifndef ACT_SIMPLE_DECAY_ALGORITHM_HH
#define ACT_SIMPLE_DECAY_ALGORITHM_HH

#include "Activia/ActAbsDecayAlgorithm.hh"
#include "Activia/ActProdNuclideList.hh"
#include "Activia/ActTarget.hh"
#include "Activia/ActTime.hh"
#include "Activia/ActConstants.hh"
#include "Activia/ActXSecGraph.hh"
#include "Activia/ActProdXSecData.hh"
#include "Activia/ActOutputTable.hh"
#include "Activia/ActDecayChains.hh"
#include "Activia/ActDecayEnergies.hh"
#include <map>
#include <vector>
#include <string>
#include <cmath>


// Algorithme Bateman (chaîne 2 membres k->j) + side-branches short-lived (déversement immédiat).
// Sortie: table ROOT/ASCII avec colonnes "Z", "A", "dNdt" (per kg per day).

class ActSimpleDecayAlgorithm : public ActAbsDecayAlgorithm {
public:
  ActSimpleDecayAlgorithm(ActTarget* target,
                          ActProdNuclideList* prodList,
                          ActTime* times);
  virtual ~ActSimpleDecayAlgorithm();

  virtual void calculateDecays(ActAbsOutput* output);
  void writeOutputPreamble();

private:
  ActDecayChains _chains;

  ActDecayEnergies _decayEnergies;


  std::map<ActNuclide*, double> buildProductionRates();
  void addShortLivedSideBranchesTo(std::map<ActNuclide*, double>& R);

  inline double lam(double halfLifeDays) const {
    return (halfLifeDays > 0.0 ? ActConstants::ln2 / halfLifeDays : 0.0);
  }

  // Contribution parent k -> j accumulée pendant l’exposition
  inline double Nj_from_k_exposure(double Rk, double lj, double lk, double b, double texp) const {
    if (lk <= 0.0) return 0.0;
    if (std::fabs(lk - lj) < 1e-12) {
      return b * Rk * ((1.0 - std::exp(-lj*texp))/lj - texp * std::exp(-lj*texp));
    }
    return b * Rk * ((1.0 - std::exp(-lj*texp))/lj
                   - (std::exp(-lj*texp) - std::exp(-lk*texp)) / (lk - lj));
  }

  // Ajout pendant le refroidissement
  inline double Nj_chain_cooling_add(double Nk_texp, double lj, double lk, double b, double tdec) const {
    if (lk <= 0.0) return 0.0;
    if (std::fabs(lk - lj) < 1e-12) {
      return b * lk * Nk_texp * (tdec * std::exp(-lj * tdec));
    }
    return b * lk * Nk_texp * ((std::exp(-lk*tdec) - std::exp(-lj*tdec)) / (lj - lk));
  }
};

#endif // ACT_SIMPLE_DECAY_ALGORITHM_HH
