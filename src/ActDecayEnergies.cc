#include "Activia/ActDecayEnergies.hh"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>

using std::cout;
using std::endl;

const ActDecayEnergies::QValues ActDecayEnergies::_zero = {0.0, 0.0};

ActNuclide* ActDecayEnergies::findProductNuclide(ActProdNuclideList* prodList,
                                                 int Z, double A) const {
  if (!prodList) return nullptr;
  const int nProd = prodList->getNProdNuclides();
  for (int i = 0; i < nProd; ++i) {
    ActProdNuclide* p = prodList->getProdNuclide(i);
    if (!p) continue;
    ActNuclide* nu = p->getProduct();
    if (nu->getZ() == Z && std::fabs(nu->getA() - A) < 1e-6) return nu;
  }
  return nullptr;
}

void ActDecayEnergies::load(const std::string& filename,
                            ActProdNuclideList* prodList) {
  _qmap.clear();

  std::ifstream in(filename.c_str());
  if (!in.is_open()) {
    cout << "ActDecayEnergies::load: cannot open " << filename << endl;
    return;
  }

  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;

    // Format:
    // Z A QEC Qbeta
    std::istringstream iss(line);
    int Z;
    double A, QEC, Qbeta;
    if (!(iss >> Z >> A >> QEC >> Qbeta)) continue;

    ActNuclide* nu = findProductNuclide(prodList, Z, A);
    if (!nu) continue;

    _qmap[nu] = {QEC, Qbeta};
  }
}

ActDecayEnergies::QValues
ActDecayEnergies::getQValues(ActNuclide* nuclide) const {
  auto it = _qmap.find(nuclide);
  if (it == _qmap.end()) return _zero;
  return it->second;
}
