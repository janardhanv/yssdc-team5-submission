#include "checking.h"

#include <stack>
#include <set>

const unsigned int HAIRPIN_LENGTH_MIN = 5;
const unsigned int HAIRPIN_LENGTH_MAX = 21;
const unsigned int HELIX_LENGTH_MIN = 4;
const unsigned int GENE_LENGTH_MIN = 100;
const double FAT_MIN = 0.35;
const double FAT_MAX = 0.65;


class DoublesFinder
{
public:
    bool Contains(const unsigned int index) const
    {
        return m_foundNumbers.find(index) != m_foundNumbers.end();
    }

    void Add(const unsigned int index)
    {
        m_foundNumbers.insert(index);
    }

    bool AddAndCheck(const unsigned int index)
    {
        if (Contains(index))
        {
            return false;
        }
        Add(index);
        return true;
    }

private:
    std::set<unsigned int> m_foundNumbers;
};


bool CheckPairsForBounds(const Pairs& pairs, const PairsLimit& pairsLimit, const size_t locusLength)
{
    if ( !(pairsLimit.firstPosition < pairsLimit.lastPosition &&
           pairsLimit.lastPosition < locusLength) )
    {
        printf("Error. Invalid bound in block.\n");
        return false;
    }
    unsigned int maxNucleotidePosition = 0;
    for (size_t i = 0; i < pairs.size(); ++i)
    {
        const Pair& pair = pairs[i];
        if ( !(pairsLimit.firstPosition <= pair.firstNucleotide &&
               pair.secondNucleotide <= pairsLimit.lastPosition) )
        {
            printf("Error. Out of bounds in pair:\n");
            PrintPair(pair, i);
            return false;
        }
        if (pair.secondNucleotide > maxNucleotidePosition)
        {
            maxNucleotidePosition = pair.secondNucleotide;
        }
    }

    if (pairsLimit.firstPosition != pairs.front().firstNucleotide)
    {
        printf("Error. There is unpaired nucleotide on the beginning of block.\n");
        return false;
    }
    if (pairsLimit.lastPosition != maxNucleotidePosition)
    {
        printf("Error. There is unpaired nucleotide on the end of block.\n");
        return false;
    }
    const unsigned int blockLength = pairsLimit.lastPosition + 1 - pairsLimit.firstPosition;
    if (blockLength < GENE_LENGTH_MIN)
    {
        printf("Error. Block is too short, length = %u.\n", blockLength);
        return false;
    }
    return true;
}


bool CheckPairsForOrder(const Pairs& pairs)
{
    for (size_t i = 0; i < pairs.size(); ++i)
    {
        const Pair& pair = pairs[i];
        if (pair.firstNucleotide >= pair.secondNucleotide)
        {
            printf("Error. Invalid pair:\n");
            PrintPair(pair, i);
            return false;
        }
    }
    for (size_t i = 1; i < pairs.size(); ++i)
    {
        const Pair& pairPrev = pairs[i-1];
        const Pair& pairNext = pairs[i];
        if (pairNext.firstNucleotide <= pairPrev.firstNucleotide)
        {
            printf("Error. Unsorted pairs:\n");
            PrintPair(pairPrev, i-1);
            PrintPair(pairNext, i);
            return false;
        }
    }
    return true;
}


bool CheckPairsForUniqueness(const Pairs& pairs)
{
    DoublesFinder doublesFinder;
    for (size_t i = 0; i < pairs.size(); ++i)
    {
        const Pair& pair = pairs[i];
        if (!doublesFinder.AddAndCheck(pair.firstNucleotide))
        {
            printf("Error. Nucleotide index %u occurs at least twice.\n", pair.firstNucleotide + 1);
            return false;
        }
        if (!doublesFinder.AddAndCheck(pair.secondNucleotide))
        {
            printf("Error. Nucleotide index %u occurs at least twice.\n", pair.secondNucleotide + 1);
            return false;
        }
    }
    return true;
}


bool CheckPairsForComplementarity(const std::string& locus,
                                  const Pairs& pairs)
{
    for (size_t i = 0; i < pairs.size(); ++i)
    {
        const Pair& pair = pairs[i];
        char c1 = locus[pair.firstNucleotide];
        char c2 = locus[pair.secondNucleotide];
        if (c1 > c2)
        {
            std::swap(c1, c2);
        }
        if (('a' == c1 && 't' != c2) ||
            ('c' == c1 && 'g' != c2) )
        {
            printf("Error. Non-complement symbols in pair:\n");
            printf("(%c-%c) --- ", locus[pair.firstNucleotide], locus[pair.secondNucleotide]);
            PrintPair(pair, i);
            return false;
        }
    }
    return true;
}


bool FindPairByNucleotide(const Pairs pairs, const unsigned int requiredNucleotide,
                          size_t& resultIndex, Pair* resultPair)
{
    for (size_t index = 0; index < pairs.size(); ++index)
    {
        const Pair& pair = pairs[index];
        if (pair.firstNucleotide  == requiredNucleotide ||
            pair.secondNucleotide == requiredNucleotide)
        {
            resultIndex = index;
            *resultPair = pair;
            return true;
        }
    }
    return false;
}


bool CheckPairsForPseudoknots(const Pairs& pairs)
{
    std::stack<unsigned int> stack;
    for (size_t i = 0; i < pairs.size(); ++i)
    {
        const Pair& pair = pairs[i];
        const unsigned int firstNucl  = pair.firstNucleotide;
        const unsigned int secondNucl = pair.secondNucleotide;
        while (!stack.empty() && stack.top() < firstNucl)
        {
            stack.pop();
        }
        if (!stack.empty() && secondNucl >= stack.top())
        {
            printf("Error. Pseudoknot in the following lines:\n");
            PrintPair(pair, i);
            size_t secondIndex;
            Pair secondPair;
            FindPairByNucleotide(pairs, stack.top(), secondIndex, &secondPair);
            PrintPair(secondPair, secondIndex);
            return false;
        }
        stack.push(secondNucl);
    }
    return true;
}


bool CheckHairpin(const Pair& pair, const size_t index)
{
    const unsigned int hairpinLength = pair.secondNucleotide - pair.firstNucleotide;
    if ( !(HAIRPIN_LENGTH_MIN <= hairpinLength && hairpinLength <= HAIRPIN_LENGTH_MAX) )
    {
        printf("Error. Inadmissible hairpin length = %u; see pair ", hairpinLength);
        PrintPair(pair, index);
        return false;
    }
    return true;
}


bool CheckPairsForHairpins(const Pairs& pairs)
{
    for (size_t i = 0; i + 1 < pairs.size(); ++i)
    {
        const Pair& currPair = pairs[i];
        const Pair& nextPair = pairs[i+1];
        if (nextPair.firstNucleotide > currPair.secondNucleotide)
        {
            if (!CheckHairpin(currPair, i))
            {
                return false;
            }
        }
    }
    const size_t maxIndex = pairs.size() - 1;
    return CheckHairpin(pairs[maxIndex], maxIndex);
}


bool CheckPairsForHelixLength(const Pairs& pairs)
{
    unsigned int helixLength = 1;
    for (size_t i = 1; i < pairs.size(); ++i)
    {
        const Pair& prevPair = pairs[i-1];
        const Pair& currPair = pairs[i];
        if (prevPair.firstNucleotide + 1 == currPair.firstNucleotide)
        {
            ++helixLength;
        }
        else
        {
            if (helixLength < HELIX_LENGTH_MIN)
            {
                printf("Error. Helix is too short (length = %u); see pair ", helixLength);
                PrintPair(currPair, i);
                return false;
            }
            helixLength = 1;
        }
    }
    return true;
}


bool CheckPairsForFAT(const std::string& locus, const PairsLimit& pairsLimit)
{
    unsigned int countAandT = 0;
    for (unsigned int i = pairsLimit.firstPosition;
         i <= pairsLimit.lastPosition;
         ++i)
    {
        const char nucleotide = locus[i];
        if ('a' == nucleotide || 't' == nucleotide)
        {
            ++countAandT;
        }
    }
    const double fat = countAandT * 1.0 / GetLength(pairsLimit);
    if ( !(FAT_MIN <= fat && fat <= FAT_MAX) )
    {
        printf("Error. Inadmissible F_AT = %.3f.\n", fat);
        return false;
    }
    printf("* F_AT = %.3f\n", fat);
    return true;
}


bool BoundsChecker::CheckAndAdd(const PairsLimit& segment)
{
    bool result = true;
    for (std::vector<PairsLimit>::const_iterator it = m_segments.begin();
            it != m_segments.end();
            ++it)
    {
        if (segment.firstPosition <= it->firstPosition &&
            segment.lastPosition  >= it->firstPosition)
        {
            result = false;
            break;
        }
        if (segment.firstPosition <= it->lastPosition &&
            segment.lastPosition  >= it->lastPosition)
        {
            result = false;
            break;
        }
        if (segment.firstPosition >= it->firstPosition &&
            segment.lastPosition  <= it->lastPosition)
        {
            result =  false;
            break;
        }
    }
    if (!result)
    {
        printf("Error. Current block intersects another block.\n");
        return false;
    }
    m_segments.push_back(segment);
    return true;
}
