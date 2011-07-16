#include "answer_reader.h"
#include <cstring>
#include <sstream>

const unsigned int BUF_LENGTH = 4096;
const unsigned int MAX_INDEX = 1000000000;


bool IsComment(const char* const str)
{
    if (strlen(str) == 0)
    {
        return true;
    }
    return (str[0] == '#');
}


bool RemoveEndOfLine(char* const str)
{
    const unsigned int length = strlen(str);
    if (0 == length)
    {
        return false;
    }
    char& lastChar = str[length - 1];
    if (13 == lastChar || 10 == lastChar) //  "\r\n" in Windows, "\n" in Unix
    {
        lastChar = '\0';
        return true;
    }
    return false;
}


AnswerReader::AnswerReader(const char* const inputName)
{
    m_inFile = fopen(inputName, "r");
}


AnswerReader::~AnswerReader()
{
    if (m_inFile)
    {
        fclose(m_inFile);
    }
}


bool AnswerReader::IsReady() const
{
    return m_inFile != NULL;
}


bool AnswerReader::ReadString(char* const outStr) const
{
    do
    {
        if (!fgets(outStr, BUF_LENGTH, m_inFile))
        {
            return false;
        }
        while (RemoveEndOfLine(outStr))
        {}
    } while (IsComment(outStr));
    return true;
}


bool AnswerReader::ReadHeader() const
{
    char currentBuf[BUF_LENGTH];
    if (!ReadString(currentBuf))
    {
        printf("IO error. Invalid header of answer file.\n");
        return false;
    }
    if (strlen(currentBuf) != 0 && '>' == currentBuf[0])
    {
        printf("%s\n", currentBuf);
        return true;
    }
    else
    {
        printf("IO error. Invalid header of answer file.\n");
        return false;
    }
}


ReadResult AnswerReader::ReadBlock(Pairs* pairs, PairsLimit* pairsLimit) const
{
    static char currentBuf[BUF_LENGTH];
        
    if (!ReadString(currentBuf))
    {
        // printf("Can't read next block.\n");
        return END_OF_FILE;
    }
    std::stringstream parser(currentBuf);
    unsigned int count;
    parser >> pairsLimit->firstPosition >> pairsLimit->lastPosition >> count;
    printf("* %u pairs in block.\n", count);
    if (count > MAX_INDEX)
    {
        printf("IO error. Bad count of pairs.\n");
        return ERROR;
    }
    --pairsLimit->firstPosition;
    --pairsLimit->lastPosition;
    if (pairsLimit->firstPosition > MAX_INDEX || 
        pairsLimit->lastPosition  > MAX_INDEX )
    {
        printf("IO error. Out of bounds in block declaration.\n");
        return ERROR;
    }

    pairs->resize(count);
    for (unsigned int index = 0; index < count; ++index)
    {
        Pair pair;
        if (!ReadString(currentBuf))
        {
            printf("IO error. Can't read pair #%u", index + 1);
            return ERROR;
        }
        parser.clear();
        parser.str(currentBuf);
        parser >> pair.firstNucleotide >> pair.secondNucleotide;
        --pair.firstNucleotide;
        --pair.secondNucleotide;
        if (pair.firstNucleotide  > MAX_INDEX || 
            pair.secondNucleotide > MAX_INDEX)
        {
            printf("IO error. Out of bounds in pair ");
            PrintPair(pair, index);
            return ERROR;
        }
        (*pairs)[index] = pair;
    }
    return SUCCESS;
}
