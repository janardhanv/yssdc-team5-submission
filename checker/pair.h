#ifndef PAIR_H
#define PAIR_H

#include <vector>
#include <cstdio>


struct Pair
{
    unsigned int  firstNucleotide;
    unsigned int secondNucleotide;
};

typedef std::vector<Pair> Pairs;


struct PairsLimit
{
    unsigned int firstPosition;
    unsigned int lastPosition;
};


inline void PrintPair(const Pair& pair, const unsigned int index)
{
    printf("%u: (%u, %u)\n",
        index + 1,
        pair.firstNucleotide + 1,
        pair.secondNucleotide + 1);
}


inline unsigned int GetLength(const PairsLimit& pairsLimit)
{
    return pairsLimit.lastPosition + 1 - pairsLimit.firstPosition;
}

#endif // PAIR_H
