// Read a file and parse each line according to the format file.
// Output integer matrices with lookup maps in either ascii text 
// or matlab binary form. 

#include <time.h>
#include <string.h>
#include <iostream>
#include <vector>

using namespace std;

#ifdef __GNUC__
#include <stdint.h>
typedef uint64_t uint64;
typedef int64_t int64;
typedef int32_t int32;
#else
typedef unsigned __int64 uint64;
typedef __int64 int64;
typedef __int32 int32;
#endif

struct eqstr {
  bool operator()(const char* s1, const char* s2) const {
    return strcmp(s1, s2) == 0;
  }
};

struct ltstr {
  bool operator()(const char* s1, const char* s2) const {
    return strcmp(s1, s2) == 0;
  }
};

#ifdef __GNUC__
//#include <unordered_map>
//typedef unordered_map<string, int> strhash;
#include <hash_map>
typedef __gnu_cxx::hash_map<const char*, int,  __gnu_cxx::hash<const char*>, eqstr > strhash;
#else
#include <hash_map>
typedef stdext::hash_map<const char*, int,  stdext::hash_compare<const char*, std::less<const char *> > > strhash;
#endif

typedef vector<char *> unhash;
typedef vector<int> ivector;
typedef vector<int64> divector;
