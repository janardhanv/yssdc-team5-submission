#include <string>
#include <cstdio>


class GbkReader
{
public:
    GbkReader(const char* const inputName);
    ~GbkReader();
    bool IsReady() const;
    bool ReadLocus(std::string* locus) const;

private:
    FILE* m_inFile;
};
