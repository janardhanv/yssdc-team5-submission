#include "pair.h"
#include <string>
#include <vector>


bool CheckPairsForBounds(const Pairs& pairs, const PairsLimit& pairsLimit, const size_t locusLength);
bool CheckPairsForOrder(const Pairs& pairs);
bool CheckPairsForComplementarity(const std::string& locus, const Pairs& pairs);
bool CheckPairsForPseudoknots(const Pairs& pairs);
bool CheckPairsForUniqueness(const Pairs& pairs);
bool CheckPairsForHairpins(const Pairs& pairs);
bool CheckPairsForHelixLength(const Pairs& pairs);
bool CheckPairsForFAT(const std::string& locus, const PairsLimit& pairsLimit);


class BoundsChecker
{
public:
    bool CheckAndAdd(const PairsLimit& segment);

private:
    std::vector<PairsLimit> m_segments;
};
