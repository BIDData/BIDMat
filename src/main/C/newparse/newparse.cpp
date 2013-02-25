// Read a file and parse each line according to the format file.
// Output integer matrices with lookup maps in either ascii text 
// or matlab binary form. 

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "newparse.h"
#include "gzstream.h"

extern int yylex(void);
extern FILE*   yyin;

ivector wcount;
ivector tokens;
unhash unh;
strhash htab;
int numlines=0;

void strtoupper(char * str) {
  if (str)
    while (*str) {
      *str = toupper(*str);
      str++;
    }
}

void strtolower(char * str) {
  if (str)
    while (*str) {
      *str = tolower(*str);
      str++;
    }
}

int checkword(char * str) {
  strtolower(str);
  int userno;
  if (htab.count(str)) {
    userno = htab[str];
    wcount[userno-1]++;
  } else {
    try {
      char * newstr = new char[strlen(str)+1];
      strcpy(newstr, str);
      wcount.push_back(1);
      unh.push_back(newstr);
      userno = unh.size();
      htab[newstr] = userno;
    } catch (bad_alloc) {
      cerr << "stringIndexer:checkstr: allocation error" << endl;
      throw;
    }
  }
  //  fprintf(stderr, "token %s (%d)\n", str, userno);
  tokens.push_back(userno);
  return userno;
}

void addtok(int tok) {
  tokens.push_back(tok | 0x80000000);
}

// For dates of form YYYY/MM/DD HH:MM:SS
int parsedate(char * str) {
  struct tm tnow;
  int i, fields[6];
  float secs;
  char * next;
  time_t tt;
  const char * delims = "T-+/: ";
  fields[0]=1900;
  fields[1]=1;
  fields[2]=0;
  fields[3]=0;
  fields[4]=0;
  fields[5]=0;
  next = strpbrk(str, delims);
  for (i = 0; i < 5 && next && *next; i++) {
    *next++ = 0;
    next += strspn(next, delims);
    sscanf(str, "%d", &fields[i]);
    str = next;
    next = strpbrk(str, delims);
  }
  if (str)
    sscanf(str, "%f", &secs);
  tnow.tm_year = fields[0]-1900;
  tnow.tm_mon = fields[1]-1;
  tnow.tm_mday = fields[2];
  tnow.tm_hour = fields[3];
  tnow.tm_min = fields[4];
  tnow.tm_sec = (int)secs;
  tnow.tm_isdst = 0;
  tt = mktime(&tnow);
  tnow.tm_year = 100;   // 
  tnow.tm_mon = 0;
  tnow.tm_mday = 1;
  tnow.tm_hour = 0;
  tnow.tm_min = 0;
  tnow.tm_sec = 0;
  tnow.tm_isdst = 0;
  tt = tt - mktime(&tnow) + 24*3600;
  return tt;
}


istream * open_in(string ifname) {
  istream * ifstr = NULL;
  int opened = 0;
  if (ifname.rfind(".gz") == ifname.length() - 3) {
    ifstr = new igzstream(ifname.c_str(), ios_base::in);
    opened = ((igzstream *)ifstr)->rdbuf()->is_open(); 
  }
  else {
    ifstr = new ifstream(ifname.c_str(), ios_base::in);
    opened = ((ifstream *)ifstr)->is_open(); 
  }
  if (!*ifstr || !opened) {
    cerr << "Couldnt open input file " << ifname << endl;
    throw;
  }
  return ifstr;
}

ostream * open_out(string ofname, int level) {
  ostream * ofstr = NULL;
  int opened = 0;
  if (ofname.rfind(".gz") == ofname.length() - 3) {
    ofstr = new ogzstream(ofname.c_str(), ios_base::out, level);
    opened = ((ogzstream *)ofstr)->rdbuf()->is_open(); 
  } else {
    ofstr = new ofstream(ofname.c_str(), ios_base::out);
    opened = ((ofstream *)ofstr)->is_open(); 
  }
  if (!*ofstr || !opened) {
    cerr << "Couldnt open output file " << ofname << endl;
    throw;
  }
  return ofstr;
}

istream * open_in_buf(string ifname, char * buffer, int buffsize) {
  istream * ifstr = NULL;
  int opened = 0;
  if (ifname.rfind(".gz") == ifname.length() - 3) {
    ifstr = new igzstream(ifname.c_str(), ios_base::in);
    opened = ((igzstream *)ifstr)->rdbuf()->is_open(); 
    ((igzstream *)ifstr)->rdbuf()->pubsetbuf(buffer, buffsize);
  }
  else {
    ifstr = new ifstream(ifname.c_str(), ios_base::in);
    opened = ((ifstream *)ifstr)->is_open(); 
    ((ifstream *)ifstr)->rdbuf()->pubsetbuf(buffer, buffsize);
  }
  if (!*ifstr || !opened) {
    cerr << "Couldnt open input file " << ifname << endl;
    throw;
  }
  return ifstr;
}

ostream * open_out_buf(string ofname, int buffsize, int level) {
  char *buffer = new char[buffsize];
  ostream * ofstr = NULL;
  int opened = 0;
  if (ofname.rfind(".gz") == ofname.length() - 3) {
    ofstr = new ogzstream(ofname.c_str(), ios_base::out, level);
    opened = ((ogzstream *)ofstr)->rdbuf()->is_open(); 
    ((ogzstream *)ofstr)->rdbuf()->pubsetbuf(buffer, buffsize);
  } else {
    ofstr = new ofstream(ofname.c_str(), ios_base::out);
    opened = ((ofstream *)ofstr)->is_open(); 
    ((ofstream *)ofstr)->rdbuf()->pubsetbuf(buffer, buffsize);
  }
  if (!*ofstr || !opened) {
    cerr << "Couldnt open output file " << ofname << endl;
    throw;
  }
  delete [] buffer;
  return ofstr;
}

void closeos(ostream *ofs) {
  ofstream * of1 = dynamic_cast<ofstream *>(ofs);
  if (of1) of1->close();
  ogzstream * of2 = dynamic_cast<ogzstream *>(ofs);
  if (of2) of2->close();
}

int writeIntVec(ivector & im, string fname, int buffsize, int level) {
  int fmt, nrows, ncols, nnz;
  ostream *ofstr = open_out_buf(fname.c_str(), buffsize, level);
  fmt = 110;
  nrows = im.size();
  ncols = 1;
  nnz = nrows;
  ofstr->write((const char *)&fmt, 4);
  ofstr->write((const char *)&nrows, 4);
  ofstr->write((const char *)&ncols, 4);
  ofstr->write((const char *)&nnz, 4);
  ofstr->write((const char *)im.data(), 4 * nrows);
  closeos(ofstr);
  return 0;
}

int writeBVec(unhash & unh, string fname, int buffsize, int level) {
  int i, s, fmt, nrows, ncols, nnz;
  ostream *ofstr = open_out_buf(fname.c_str(), buffsize, level);
  fmt = 301; // 3=sparse(no rows), 0=byte, 1=int
  ncols = unh.size();
  ivector cols;
  cols.push_back(0);
  for (i=0, nrows=0, nnz=0; i<ncols; i++) {
    s = strlen(unh[i]);
    nrows = max(nrows, s);
    nnz = nnz + s;
    cols.push_back(nnz);
  }
  ofstr->write((const char *)&fmt, 4);
  ofstr->write((const char *)&nrows, 4);
  ofstr->write((const char *)&ncols, 4);
  ofstr->write((const char *)&nnz, 4);
  ofstr->write((const char *)cols.data(), 4 * (ncols+1));
  for (i=0; i<ncols; i++) {
    ofstr->write(unh[i], cols[i+1] - cols[i]);
  }
  closeos(ofstr);
  return 0;
}


int main(int argc, char ** argv) {
  int iarg=1, membuf=1048576, complevel=3;
  char *here;
  char *ifname = NULL;
  string odname="", dictname = "";
  while (iarg < argc) {
    if (strncmp(argv[iarg], "-i", 2) == 0) {
      ifname = argv[++iarg];
    } else if (strncmp(argv[iarg], "-o", 2) == 0) {
      odname = argv[++iarg];
    } else if (strncmp(argv[iarg], "-d", 2) == 0) {
      dictname = argv[++iarg];
    } else if (strncmp(argv[iarg], "-c", 2) == 0) {
      complevel = strtol(argv[++iarg],NULL,10);
    } else if (strncmp(argv[iarg], "-s", 2) == 0) {
      membuf = strtol(argv[++iarg],NULL,10);
    } else {
      cout << "Unknown option " << argv[iarg] << endl;
      exit(1);
    }
    iarg++;
  }
  if (dictname.size() == 0) dictname = odname;
  here = strtok(ifname, " ,");
  while (here != NULL) {
    if (strstr(here, ".gz") - here == strlen(here) - 3) {
      yyin = popen( (string("gunzip -c ")+here).c_str(), "r" );
    } else {
      yyin = fopen( here, "r" );
    }
    fprintf(stderr, "\nScanning %s\n", here);
    yylex();
    if (strstr(here, ".gz") - here == strlen(here) - 3) {
      pclose(yyin);
    } else {
      fclose(yyin);
    }
    fprintf(stderr, "\r%05d lines", numlines);
    writeIntVec(tokens, odname+here, membuf, complevel);
    tokens.clear();
    numlines = 0;
    here = strtok(NULL, " ,");
  }
  fprintf(stderr, "\nWriting Dictionary\n");
  writeIntVec(wcount, dictname+"wcount.gz", membuf, complevel);
  writeBVec(unh, dictname+"dict.gz", membuf, complevel);
}
