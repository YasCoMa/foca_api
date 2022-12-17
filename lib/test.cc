#include <iostream>

#include <fstream>
#include <sstream>
#include <boost/serialization/map.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <cctype>
#include <dlfcn.h>

using namespace std;



int main ()
{
    string file = "/home/tavinbio/public_html/foca_backend/data/variants_mutations_unique_all.ser";

    multimap<string,vector<string>> mmvtab;

    ifstream f (file, ios::binary);

    boost::archive::binary_iarchive bia (f);

    bia >> mmvtab;

    f.close();

    return 0;
}








