#include <assert.h>
#include <math.h>
#include <string>
#include <iostream> 
#include <iomanip>
#include <fstream> 
#include <sstream> 
#include <vector> 
#include <algorithm>
#include <stdexcept>
#include <time.h>
#ifdef WITHGZS
#include "gzstream.h"
#endif

#define BUFSIZE 1048576

#ifdef __GNUC
#include <stdint.h>
typedef uint64_t uint64;
typedef int64_t int64;
typedef int32_t int32;
#else
typedef unsigned __int64 uint64;
typedef __int64 int64;
typedef __int32 int32;
#endif


template <class T>
std::string to_string(T t)
{
  std::ostringstream oss;
  oss << t;
  return oss.str();
}

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

#ifdef __GNUC
#include <unordered_map>
typedef std::unordered_map<const char*, int,  __gnu_cxx::hash<const char*>, eqstr> strhash;
#else
#include <hash_map>
typedef std::hash_map<const char*, int,  std::hash<const char*>, eqstr> strhash;
#endif

using std::cerr;
using std::cout;
using std::endl;
using std::setw;
using std::string;
using std::min;
using std::max;
using std::ifstream;
using std::ofstream;
using std::istream;
using std::ostream;
using std::ios_base;
using std::istringstream;
using std::ostringstream;
using std::sort;
using std::stable_sort;

#ifdef DEBUG
template <class T>
class vector : public std::vector<T> {
public:
  vector(unsigned size) : std::vector<T>(size) {}
  vector() : std::vector<T>() {}
  inline T & operator[](unsigned i) {
    if (i >= std::vector<T>::size()) {
      cout << "index out of range " << i << "," << std::vector<T>::size() << endl;
      throw std::out_of_range("");
    }
    return std::vector<T>::at(i);}
  inline const T & operator[](unsigned i) const {
    if (i >= std::vector<T>::size()) {
      cout << "index out of range " << i << "," << std::vector<T>::size() << endl;
      throw std::out_of_range("");
    }
    return std::vector<T>::at(i);}
};
#else
using std::vector;
#endif

typedef struct quadint {
  uint64 top;
  uint64 bottom;
} qint;


template <class T>
class vpair {
public:
  int ind;
  T val;
  vpair(int i, T v): ind(i), val(v) {};
  vpair(): ind(0), val(0) {};
  ~vpair() {};
};

template <class T>
bool indxlt(const vpair<T> & v1, const vpair<T> & v2) {
  return (v1.ind < v2.ind);
}

template <class T>
bool valult(const vpair<T> & v1, const vpair<T> & v2) {
  return (v1.val < v2.val);
}

typedef vector<char *> unhash;
typedef vpair<float> fpair;
typedef vpair<int> ipair;
typedef vector<fpair> sfvector;
typedef vector<ipair> ipvector;
typedef vector<int> ivector;
typedef vector<int64> divector;
typedef vector<qint> qvector;
typedef vector<ivector> imatrix;
typedef vector<divector> dimatrix;
typedef vector<float> fvector;
typedef vector<double> dvector;
typedef vector<string> svector;
typedef vector<fvector> fmatrix;

class fieldtype {
 public:
  ivector iv;
  divector div;
  qvector qv;
  fvector fv;
  dvector dv;
  imatrix im;
  dimatrix dim;
 fieldtype() : iv(0), div(0), fv(0), dv(0), im(0), dim(0) {};
  int writeInts(string);
  int writeDInts(string);
  int writeQInts(string);
  int writeFloats(string);
  int writeDoubles(string);
  int writeIntVecs(string);
  int writeIntVecsTxt(string);
  int writeDIntVecs(string);
};

typedef vector<fieldtype> ftvector;

class stringIndexer {
 public:
  strhash htab;
  unhash unh;
  ivector count;
  int size;
  char * linebuf;
  stringIndexer(int sz) : htab(sz), unh(sz), count(sz), size(sz) {
    linebuf = new char[BUFSIZE];
  };
  stringIndexer() : htab(0), unh(0), count(0), size(0) {
    linebuf = new char[BUFSIZE];
  };
  stringIndexer(const stringIndexer & si);
  ~stringIndexer();
  stringIndexer(char *);
  int writeMap(string);
  int checkword(char *);
  ivector checkstring(char *, const char *);
  ivector checkstrings(char **, const char *, int);
  ivector checkgroup(char ** here, const char * delim2, int len);
  divector checkdgroup(char ** here, const char * delim2, int len);
  stringIndexer shrink(int);
  stringIndexer shrink_to_size(int);  
  ivector indexMap(stringIndexer &);
};

typedef vector<stringIndexer> srivector;

int parseLine(char *, int, const char *, srivector &, ivector &);
int parseLineX(char *, int, const char *, ivector &, svector &,
	       srivector &, ftvector &, int);
string getfield(char * line, const char * delim1, int k);

enum ftypes {
  ftype_int = 1,
  ftype_dint,
  ftype_qhex,
  ftype_float,
  ftype_double,
  ftype_word,
  ftype_string,
  ftype_date,
  ftype_mdate,
  ftype_dt,
  ftype_mdt,
  ftype_group,
  ftype_igroup,
  ftype_digroup
};

#ifdef WITHGZS
istream * open_in(string);
ostream * open_out(string);
istream * open_in_buf(string, char *, int);
ostream * open_out_buf(string, char *, int);
#endif
