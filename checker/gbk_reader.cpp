#include "gbk_reader.h"

const size_t BUF_LENGTH = 4096;

// forward declaration
bool RemoveEndOfLine(char* const str);


GbkReader::GbkReader(const char* const inputName)
{
    m_inFile = fopen(inputName, "r");
}


GbkReader::~GbkReader()
{
    if (m_inFile != NULL)
    {
        fclose(m_inFile);
    }
}


bool GbkReader::IsReady() const
{
    return (m_inFile != NULL);
}


bool GbkReader::ReadLocus(std::string* locus) const
{
    locus->clear();
    size_t length;
    if (fscanf(m_inFile, "%lu\n", &length) != 1)
    {
        printf("IO error. Invalid header of GBK file.\n");
        return false;
    }
    printf("Locus length = %lu\n", length);
    char buf[BUF_LENGTH];
    do
    {
        if (!fgets(buf, BUF_LENGTH, m_inFile))
        {
            printf("IO error. Length of locus read = %lu\n", locus->length());
            return false;
        }
        RemoveEndOfLine(buf);
        *locus += buf;
    } while (locus->length() != length);
    return true;
}
