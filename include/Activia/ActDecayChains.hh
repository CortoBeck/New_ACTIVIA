
#ifndef ACT_DECAY_CHAINS_HH
#define ACT_DECAY_CHAINS_HH

#include "Activia/ActProdNuclideList.hh"
#include "Activia/ActNuclide.hh"
#include <map>
#include <vector>
#include <string>


class ActDecayChains {
public:
  struct ParentLink {
    ActNuclide* parent;
    double      branching;
  };

  ActDecayChains() : _parentsOf() {}
  ~ActDecayChains() {}

  void load(const std::string& filename, ActProdNuclideList* prodList);
  const std::vector<ParentLink>& parentsOf(ActNuclide* child) const;

private:
  ActNuclide* findProductNuclide(ActProdNuclideList* prodList, int Z, double A) const;

  std::map<ActNuclide*, std::vector<ParentLink>> _parentsOf;
  static const std::vector<ParentLink> _empty;
};

#endif // ACT_DECAY_CHAINS_HH
