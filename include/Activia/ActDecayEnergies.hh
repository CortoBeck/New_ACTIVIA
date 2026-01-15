#ifndef ACT_DECAY_ENERGIES_HH
#define ACT_DECAY_ENERGIES_HH

#include "Activia/ActNuclide.hh"
#include "Activia/ActProdNuclideList.hh"
#include <map>
#include <string>

class ActDecayEnergies {
public:
  struct QValues {
    double QEC;    // keV
    double Qbeta;  // keV
  };

  ActDecayEnergies() {}
  ~ActDecayEnergies() {}

  void load(const std::string& filename, ActProdNuclideList* prodList);
  QValues getQValues(ActNuclide* nuclide) const;

private:
  ActNuclide* findProductNuclide(ActProdNuclideList* prodList, int Z, double A) const;

  std::map<ActNuclide*, QValues> _qmap;
  static const QValues _zero;
};

#endif
