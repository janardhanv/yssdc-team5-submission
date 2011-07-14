#include <vector>
#include <cstdio>

#include "pair.h"


enum ReadResult
{
    SUCCESS,
    END_OF_FILE,
    ERROR
};


class AnswerReader
{
public:
    AnswerReader(const char* const inputName);
    ~AnswerReader();
    bool IsReady() const;
    bool ReadHeader() const;
    ReadResult ReadBlock(Pairs* pairs, PairsLimit* pairsLimit) const;

private:
    bool ReadString(char* const outStr) const;

    FILE* m_inFile;
};
