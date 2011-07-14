#include "gbk_reader.h"
#include "answer_reader.h"
#include "checking.h"

#include <map>

const unsigned int F_PAIRED_MINS_COUNT = 4;
const unsigned int F_PAIRED_MINS[] = {55, 60, 65, 70};


void ShowHelp()
{
    printf("Command-line parameters:\n");
    printf("\trnachecker GBK_FILENAME ANSWER_FILENAME\n");
}


bool ReadGbkFile(const char* const filename, std::string* locus)
{
    printf("Reading %s...\n", filename);
    GbkReader gr(filename);
    if (!gr.IsReady())
    {
        printf("Can't open GBK file!\n");
        return false;
    }
    if (!gr.ReadLocus(locus))
    {
        printf("Can't read GBK file!\n");
        return false;
    }
    printf("Done.\n");
    return true;
}


bool CheckAll(const AnswerReader& answerReader, const std::string& locus)
{
    if (!answerReader.ReadHeader())
    {
        return false;
    }

    std::map<unsigned int, unsigned int> totalPairsCount;
    std::map<unsigned int, unsigned int> totalBlocksLength;
    for (size_t i = 0; i < F_PAIRED_MINS_COUNT; ++i)
    {
        const unsigned int category = F_PAIRED_MINS[i];
        totalPairsCount[category] = 0;
        totalBlocksLength[category] = 0;
    }

    BoundsChecker boundsChecker;
    for (size_t count = 1; ; ++count)
    {
        printf("Reading block #%lu...\n", count);
        Pairs pairs;
        PairsLimit pairsLimit;
        const ReadResult readResult = answerReader.ReadBlock(&pairs, &pairsLimit);
        if (readResult == END_OF_FILE)
        {
            printf("Finished. All blocks read.\n");
            break;
        }
        else if (readResult == ERROR)
        {
            return false;
        }
        if (pairs.size() == 0)
        {
            continue;
        }

        // printf("%s\n", locus.substr(pairsLimit.firstPosition, GetLength(pairsLimit)).c_str());

        printf("Checking block #%lu...\n", count);
        if (!boundsChecker.CheckAndAdd(pairsLimit) ||
            !CheckPairsForBounds(pairs, pairsLimit, locus.length()) ||
            !CheckPairsForOrder(pairs) ||
            !CheckPairsForUniqueness(pairs) ||
            !CheckPairsForComplementarity(locus, pairs) ||
            !CheckPairsForPseudoknots(pairs) ||
            !CheckPairsForHairpins(pairs) ||
            !CheckPairsForHelixLength(pairs) ||
            !CheckPairsForFAT(locus, pairsLimit) )
        {
            return false;
        }

        const unsigned int blockLength = pairsLimit.lastPosition + 1 - pairsLimit.firstPosition;
        const double fPaired = pairs.size() * 2.0 / blockLength;
        printf("* fPaired = %.3f\n", fPaired);
        for (size_t i = 0; i < F_PAIRED_MINS_COUNT; ++i)
        {
            const unsigned int category = F_PAIRED_MINS[i];
            const double fPairedMin = category * 0.01;
            if (fPaired >= fPairedMin)
            {
                totalPairsCount[category] += pairs.size();
                totalBlocksLength[category] += blockLength;
            }
        }
    }

    printf("OK\n");
    printf("-------------------------------------------------------------------\n");
    printf("fPairedMin ; Total number of base pairings ; Total length of blocks\n");
    for (size_t i = 0; i < F_PAIRED_MINS_COUNT; ++i)
    {
        const unsigned int category = F_PAIRED_MINS[i];
        const double fPairedMin = category * 0.01;
        printf("%.2f ; %u ; %u\n", fPairedMin, totalPairsCount[category], totalBlocksLength[category]);
    }
    return true;
}


int main(int argc, char* argv[])
{
    printf("Rna-Checker ver. 1.14.16.05\n");
    if (argc != 3)
    {
        ShowHelp();
        return 2;
    }
    
    std::string locus;
    if (!ReadGbkFile(argv[1], &locus))
    {
        return 1;
    }

    printf("Reading %s...\n", argv[2]);
    AnswerReader ar(argv[2]);
    if (!ar.IsReady())
    {
        printf("Can't open %s!\n", argv[2]);
        return 1;
    }

    return CheckAll(ar, locus) ? 0 : 1;
}
