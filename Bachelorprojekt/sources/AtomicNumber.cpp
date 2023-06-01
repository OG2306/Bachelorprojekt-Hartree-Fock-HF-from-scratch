#include "AtomicNumber.h"
#include <cstring>

int findatomicnumber(const char* const symbol)
{
    for (int i = 0; i < 118; i++) {
        if (strcmp(elements[i], symbol) == 0) return i + 1;
    }
    return 0; //Returnerer 0 hvis symbolet ikke genkendes
}