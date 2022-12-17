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

//Convert vector to string
inline string vjson (vector<string> &v)
{
    return "\"ORF\":\"" + v[0] + "\",\"lastdate\":\"" + v[1] + "\",\"Nucl_Modified\":\"" + v[2] + "\",\"AA_Modified\":\"" + v[3] + "\",\"Count\":\"" + v[4] + "\",\"Proportion\":\"" + v[5] + "\"";
}

//get utr position
inline int posutr (string u)
{
    string num = "";
    for (auto c: u)
        if (isdigit(c)) 
            num += c;
        else
            return atoi (num.c_str());

    return 0;
}

extern "C" {
//Unique mutation by variants
    const char* get_unique_mutations (const char* pais, const char* lins, int ldate, double prop)
    {
   	    multimap<string,vector<string>> mmvtab;

        string ppais = pais;
        string file;

        if (ppais=="world")
            file = "/home/tavinbio/public_html/foca_backend/data/variants_mutations_unique_all.ser";
        else if (ppais=="brasil")
            file = "";


        ifstream f (file, ios::binary);
        boost::archive::binary_iarchive bia (f);
        bia >> mmvtab;

        f.close();

        set<string> sli;
	    string s;
    	stringstream ssl (lins);
	    while (ssl >> s) sli.insert (s);

        string stab = "";
	    pair<multimap<string,vector<string>>::iterator,multimap<string,vector<string>>::iterator> ret;
	    multimap<string,vector<string>>::iterator it;

        int pos;

        string str;

	    for (auto l: sli){ 
            ret = mmvtab.equal_range (l);
            for (it=ret.first; it!=ret.second; ++it){
                pos = posutr(it->second[2]);
                if ( (atof(it->second[5].c_str()) > prop) && (atoi(it->second[1].c_str()) > ldate) )
                    stab += "{\"lineage\":\"" + it->first + "\"," + vjson(it->second) + "},";   
            }
        }

        if (stab != "") stab.pop_back();

        stab = "{\"data\":[" + stab + "]}";

	    char* ctab = new char [stab.size()+1];
        strcpy (ctab, stab.c_str());

	    return ctab;
    }


}
