#include <iostream>

using namespace std;

extern "C" 
{
    const char* convseq (const char* file)
    {
        ifstream f (file);
        string s;
        f >>s >>s;

        f.close();
    }
}