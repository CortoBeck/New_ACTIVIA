
#include "Activia/ActDecayChains.hh"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>

using std::cout;
using std::endl;

const std::vector<ActDecayChains::ParentLink> ActDecayChains::_empty = {};

ActNuclide* ActDecayChains::findProductNuclide(ActProdNuclideList* prodList, int Z, double A) const {
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

void ActDecayChains::load(const std::string& filename, ActProdNuclideList* prodList) {
  _parentsOf.clear();

  std::ifstream in(filename.c_str());
  if (!in.is_open()) {
    cout << "ActDecayChains::load: cannot open " << filename << endl;
    return;
  }

  std::string line;
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    if (line[0] == '#') continue;

    // Format: Z_k A_k  Z_j A_j  branching
    std::istringstream iss(line);
    int Zk, Zj;
    double Ak, Aj, b;
    if (!(iss >> Zk >> Ak >> Zj >> Aj >> b)) continue;

    ActNuclide* k = findProductNuclide(prodList, Zk, Ak);
    ActNuclide* j = findProductNuclide(prodList, Zj, Aj);
    if (!k || !j) {
      // Parent ou fils non listé comme "produit" dans decayData.dat → ignorer la ligne
      continue;
    }

    _parentsOf[j].push_back(ParentLink{ k, b });
  }
}

const std::vector<ActDecayChains::ParentLink>& ActDecayChains::parentsOf(ActNuclide* child) const {
  auto it = _parentsOf.find(child);
  if (it == _parentsOf.end()) return _empty;
  return it->second;
}
