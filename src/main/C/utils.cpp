#include "utils.h"

void indxsort(sfvector & sfv) {
  sort(sfv.begin(), sfv.end(), indxlt<float>);
}

void valusort(sfvector & sfv) {
  sort(sfv.begin(), sfv.end(), valult<float>);
}

void indxsort(sfvector::iterator p1, sfvector::iterator p2) {
  sort(p1, p2, indxlt<float>);
}

void valusort(sfvector::iterator p1, sfvector::iterator p2) {
  sort(p1, p2, valult<float>);
}

void indxsorts(sfvector & sfv) {
  stable_sort(sfv.begin(), sfv.end(), indxlt<float>);
}

void valusorts(sfvector & sfv) {
  stable_sort(sfv.begin(), sfv.end(), valult<float>);
}

void indxsorts(sfvector::iterator p1, sfvector::iterator p2) {
  stable_sort(p1, p2, indxlt<float>);
}

void valusorts(sfvector::iterator p1, sfvector::iterator p2) {
  stable_sort(p1, p2, valult<float>);
}

const int stroffset = 'A' - 'a';

void strtoupper(char * str) {
  if (str)
    while (*str) {
      if (*str & 128)
	*str = ' ';
      if (*str <= 'z' && *str >= 'a')
	*str += stroffset;
      str++;
    }
}

int stringIndexer::
checkword(char * str) {
  int userno;
  char * newstr;
  if (htab.count(str)) {
    userno = htab[str];
    count[userno-1]++;
  } else {
    try {
      newstr = new char[strlen(str)+1];
      strcpy(newstr, str);
      userno = ++size;
      htab[newstr] = userno;
      count.push_back(1);
      unh.push_back(newstr);
    } catch (std::bad_alloc) {
      cerr << "stringIndexer:checkstr: allocation error" << endl;
      throw;
    }
  }
  return userno;
}

ivector stringIndexer::
checkstring(char * here, const char * delim2) {
  char * next;
  ivector out(0);
  here += strspn(here, delim2);
  //  strtoupper(here);
  while (here && *here) {
    next = strpbrk(here, delim2);
    if (next) {
      *(next++) = 0;
      next += strspn(next, delim2);
    }
    out.push_back(checkword(here));
    here = next;
  }    
  return out;
}

ivector stringIndexer::
checkstrings(char ** here, const char * delim2, int len) {
  int i;
  char * next;
  ivector out(0);
  //  strtoupper(*here);
  for (i = 0; i < len; i++) {
    next = strpbrk(*here, delim2);
    if (next) {
      *(next++) = 0;
    }
    if (strlen(*here) > 0)
      out.push_back(checkword(*here));
    else
      out.push_back(0);
    *here = next;
  }    
  return out;
}

ivector stringIndexer::
checkgroup(char ** here, const char * delim2, int len) {
  int i;
  char * next;
  ivector out(0);
  int64 val;
  next = strpbrk(*here, delim2);
  if (next == *here) (*here)++;
  for (i = 0; i < len && *here; i++) {
    next = strpbrk(*here, delim2);
    if (next) {
      *(next++) = 0;
    }
    if (strlen(*here) > 0) {
      sscanf(*here, "%d", &val);
      out.push_back(val);
    } else
      out.push_back(0);
    *here = next;
  }    
  return out;
}


divector stringIndexer::
checkdgroup(char ** here, const char * delim2, int len) {
  int i;
  char * next;
  divector out(0);
  int64 val;
  next = strpbrk(*here, delim2);
  if (next == *here) (*here)++;
  for (i = 0; i < len && *here; i++) {
    next = strpbrk(*here, delim2);
    if (next) {
      *(next++) = 0;
    }
    if (strlen(*here) > 0) {
      sscanf(*here, "%lld", &val);
      out.push_back(val);
    } else
      out.push_back(0);
    *here = next;
  }    
  return out;
}

// Need a deep copy constructor for stringIndexer for this to work
stringIndexer stringIndexer::
shrink(int threshold) {
  int i, newc;
  char * newstr;
  for (i = 0, newc = 0; i < size; i++) {
    if (count[i] >= threshold) {
      newc++;
    }
  }
  stringIndexer retval(newc);
  try {
    for (i = 0, newc = 0; i < size; i++) {
      if (count[i] >= threshold) {
	newstr = new char[strlen(unh[i])+1];
	strcpy(newstr, unh[i]);
	retval.count[newc] = count[i];
	retval.unh[newc] = newstr;
	retval.htab[newstr] = ++newc;
      }
    }
  } catch (std::bad_alloc) {
    cerr << "stringIndexer:shrink: allocation error" << endl;
    throw;
  }
  return retval;
}

stringIndexer stringIndexer::
shrink_to_size(int newsize) {
  int i, newc, minth=1, maxth=0, midth;
  for (i = 0, newc = 0; i < size; i++) {
    if (count[i] >= maxth) {
      maxth = count[i];
    }
  }
  while (maxth > (minth+1)) {
    midth = (maxth + minth) >> 1;
    for (i = 0, newc = 0; i < size; i++) {
      if (count[i] >= midth) {
	newc++;
      }
    }
    if (newc > newsize) {
      minth = midth;
    } else {
      maxth = midth;
    }
  }
  return shrink(maxth);
}

ivector stringIndexer::
indexMap(stringIndexer & A) {
  int i;
  ivector remap(size);
  for (i = 0; i < size; i++) {
    if (A.htab.count(unh[i])) {
      remap[i] = A.htab[unh[i]];
    } else {
      remap[i] = 0;
    }
  }
  return remap;
}

int parseLine(char *line, int lineno, char *delim, srivector &srv, ivector &out) {
  int i;
  char * here, * next;
  here = line;
  for (i = 0; i < srv.size(); i++) {
    next = strpbrk(here, delim);
    if (!next && i < (srv.size() - 1)) {
      cerr << "parseLine: format error line " << lineno << endl;
      cerr << "  contents: " << line << " ... " << here << endl;
      throw;
    }
    if (next) {
      *(next++) = 0;
    }
    out[i] = srv[i].checkword(here);
    here = next;
  }
  return 0;
}

// Return a unix local time for dates of form MM-DD-YYYY
int parsemdt(char * str) {
  struct tm tnow;
  int day=0, month=1, year=1900;
  char * next;
  time_t tt;
  const char * delims = "-/: ";
  next = strpbrk(str, delims);
  if (next) {
    *next++ = 0;
    sscanf(str, "%d", &month);
    str = next;
    next = strpbrk(str, delims);
    if (next) {
      *next++ = 0;
      sscanf(str, "%d", &day);
      sscanf(next, "%d", &year);
    }
  }
  tnow.tm_year = year-1900;
  tnow.tm_mon = month-1;
  tnow.tm_mday = day;
  tnow.tm_isdst = 0;
  tnow.tm_sec = 0;
  tnow.tm_min = 0;
  tnow.tm_hour = 0;
  tnow.tm_isdst = 0;
  tt = mktime(&tnow);
  tnow.tm_year = 70;   // Jan 2, 1970 in GMT
  tnow.tm_mon = 0;
  tnow.tm_mday = 2;
  tnow.tm_isdst = 0;
  tt = tt - mktime(&tnow) + 24*3600;
  return tt;
}

// For dates of form YYYY-MM-DD
int parsedt(char * str) {
  struct tm tnow;
  int day=0, month=1, year=1900;
  char * next;
  time_t tt;
  const char * delims = "-/: ";
  next = strpbrk(str, delims);
  if (next) {
    *next++ = 0;
    sscanf(str, "%d", &year);
    str = next;
    next = strpbrk(str, delims);
    if (next) {
      *next++ = 0;
      sscanf(str, "%d", &month);
      sscanf(next, "%d", &day);
    }
  }
  tnow.tm_year = year-1900;
  tnow.tm_mon = month-1;
  tnow.tm_mday = day;
  tnow.tm_isdst = 0;
  tnow.tm_sec = 0;
  tnow.tm_min = 0;
  tnow.tm_hour = 0;
  tnow.tm_isdst = 0;
  tt = mktime(&tnow);
  tnow.tm_year = 70;   // Jan 2, 1970 in GMT
  tnow.tm_mon = 0;
  tnow.tm_mday = 2;
  tnow.tm_isdst = 0;
  tt = tt - mktime(&tnow) + 24*3600;
  return tt;
}

// For dates of form YYYY/MM/DD HH:MM:SS
int parsedate(char * str) {
  struct tm tnow;
  int i, fields[6];
  float secs;
  char * next;
  time_t tt;
  const char * delims = "-/: ";
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
    //    printf("%d ", fields[i]);
    str = next;
    next = strpbrk(str, delims);
  }
  if (str)
    sscanf(str, "%f", &secs);
  //  printf("%s %f %d\n", str, secs, (int)secs);
  tnow.tm_year = fields[0]-1900;
  tnow.tm_mon = fields[1]-1;
  tnow.tm_mday = fields[2];
  tnow.tm_hour = fields[3];
  tnow.tm_min = fields[4];
  tnow.tm_sec = (int)secs;
  tnow.tm_isdst = 0;
  tt = mktime(&tnow);
  tnow.tm_year = 70;   // Jan 2, 1970 in GMT
  tnow.tm_mon = 0;
  tnow.tm_mday = 2;
  tnow.tm_hour = 0;
  tnow.tm_min = 0;
  tnow.tm_sec = 0;
  tnow.tm_isdst = 0;
  tt = tt - mktime(&tnow) + 24*3600;
  return tt;
}

// For dates of form MM/DD/YYYY HH:MM:SS
int parsemdate(char * str) {
  struct tm tnow;
  int i, fields[6];
  float secs;
  char * next;
  time_t tt;
  const char * delims = "-/: ";
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
    //    printf("%d ", fields[i]);
    str = next;
    next = strpbrk(str, delims);
  }
  if (str)
    sscanf(str, "%f", &secs);
  //  printf("%s %f %d\n", str, secs, (int)secs);
  tnow.tm_year = fields[2]-1900;
  tnow.tm_mon = fields[0]-1;
  tnow.tm_mday = fields[1];
  tnow.tm_hour = fields[3];
  tnow.tm_min = fields[4];
  tnow.tm_sec = (int)secs;
  tnow.tm_isdst = 0;
  tt = mktime(&tnow);
  tnow.tm_year = 70;   // Jan 2, 1970 in GMT
  tnow.tm_mon = 0;
  tnow.tm_mday = 2;
  tnow.tm_hour = 0;
  tnow.tm_min = 0;
  tnow.tm_sec = 0;
  tnow.tm_isdst = 0;
  tt = tt - mktime(&tnow) + 24*3600;
  return tt;
}

int delimcount(char * line, const char * delims) {
  int cnt=0;
  while (line = strpbrk(line, delims)) {
    cnt++;
    line++;
  }
  return cnt;
}

string getfield(char * line, const char * delim1, int k) {
  char * here, * next;
  int i;
  next = line - 1;
  here = line;
  for (i = 0; i < k; i++) {
    if (next) 
      here = next+1;
    else
      return "";
    next = strpbrk(here, delim1);
  }
  if (next)
    return string(here, (int)(next - here));
  else
    return string(here);
}    


int parseLineX(char * line, int lineno, const char * delim1, ivector & tvec, 
	       svector & delims, srivector & srv, ftvector & out, int grpsize) {
  int i, ival;
  int64 dival;
  qint qval;
  float fval;
  double dval;
  char * here, * next;

  here = line;
  for (i = 0; i < tvec.size(); i++) {
    if (i < tvec.size()-1) {
      next = strpbrk(here, delim1);
      if (!next) {
	cerr << "parseLineX: format error line " << lineno << endl;
	cerr << "  contents: " << line << " ... " << here << endl;
	throw 10;
      }
      *(next++) = 0;
    }
    switch (tvec[i]) {
    case ftype_int:
      sscanf(here, "%d", &ival);
      out[i].iv.push_back(ival);
      break;
    case ftype_dint:
      sscanf(here, "%lld", &dival);
      out[i].div.push_back(dival);
      break;
    case ftype_qhex:
      sscanf(here, "%16Lx%16Lx", &qval.top, &qval.bottom);
      out[i].qv.push_back(qval);
      break;
    case ftype_float:
      sscanf(here, "%f", &fval);
      out[i].fv.push_back(fval);
      break;
    case ftype_double:
      sscanf(here, "%lf", &dval);
      out[i].dv.push_back(dval);
      break;
    case ftype_word:
      here += strspn(here, " ");
      out[i].iv.push_back(srv[i].checkword(here));
      break;
    case ftype_string:
      out[i].im.push_back(srv[i].checkstring(here, delims[i].c_str()));
      break;
    case ftype_dt:  
      ival = parsedt(here);
      if (ival < 0)
	printf("\nWarning: bad dt on line %d\n", lineno);
      out[i].iv.push_back(ival);
      break;
    case ftype_mdt:
      ival = parsemdt(here);
      if (ival < 0)
	printf("\nWarning: bad mdt on line %d\n", lineno);
      out[i].iv.push_back(ival);
      break;
    case ftype_date:
      ival = parsedate(here);
      if (ival < 0)
	printf("\nWarning: bad date on line %d\n", lineno);
      out[i].iv.push_back(ival);
      break;
    case ftype_mdate:
      ival = parsemdate(here);
      if (ival < 0)
	printf("\nWarning: bad date on line %d\n", lineno);
      out[i].iv.push_back(ival);
      break;
    case ftype_group:
      *(next-1) = *delim1;
      out[i].im.push_back(srv[i].checkstrings(&here, delim1, grpsize));
      next = here;
      break;
    case ftype_igroup:
      //      *(next-1) = *delim1;
      out[i].im.push_back(srv[i].checkgroup(&here, delims[i].c_str(), grpsize));
      next = here;
      break;
    case ftype_digroup:
      //      *(next-1) = *delim1;
      out[i].dim.push_back(srv[i].checkdgroup(&here, delims[i].c_str(), grpsize));
      next = here;
      break;
    default:
      break;
    }
    here = next;
  }
  return 0;
}


stringIndexer::
stringIndexer(char * fname) : htab(0), unh(0), count(0), size(0) {
  int idx, cnt, n = 0;
  char * here;
  ifstream ifstr(fname, ios_base::in);
  linebuf = new char[BUFSIZE];
  while (!ifstr.bad() && !ifstr.eof()) {
    n++;
    ifstr.getline(linebuf, BUFSIZE-1);
    linebuf[BUFSIZE-1] = 0;
    if (ifstr.fail()) {
      ifstr.clear();
      ifstr.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    }
    if (strlen(linebuf) > 2) {
      sscanf(linebuf, "%d %d", &idx, &cnt);
      here = strchr(linebuf, ' ');
      here = strchr(++here, ' ');
      here++;
      if (idx != n) {
	cerr << "stringIndexer load: file format error" << endl;
	throw;
      }
      checkword(here);
      count[n-1] = cnt;
    }
  }  
}

// Deep copy constructor
stringIndexer::
stringIndexer(const stringIndexer & si) : 
  htab(0), unh(si.size), count(si.size), size(si.size) 
{
  int i;
  char *str, *newstr;
  for (i = 0; i < size; i++) {
    str = si.unh[i];
    newstr = new char[strlen(str)+1];
    strcpy(newstr, str);
    unh[i] = newstr;
    htab[newstr] = i+1;
    count[i] = si.count[i];
  }
  linebuf = new char[BUFSIZE];
}

stringIndexer::
~stringIndexer() {
  int i;
  for (i=0; i<unh.size(); i++) {
    if (unh[i]) {
      delete [] unh[i];
      unh[i] = NULL;
    }
  }
  if (linebuf) {
    delete [] linebuf;
    linebuf = NULL;
  }
}

int stringIndexer::
writeMap(string fname) {
  int i;
  char * here;
  ofstream ofstr(fname.c_str(), ios_base::out);
  for (i = 0; i < size; i++) {
    sprintf(linebuf, "%d %d %s\n", i+1, count[i], unh[i]);
    ofstr << linebuf;
  }
  return 0;
}

void writeStrVector(char * linebuf, ivector & iv) {
  char * here = linebuf;
  int j;
  for (j = 0; j < (iv.size()-1); j++) {
    sprintf(here, "%d ", iv[j]);
    here += strlen(here);
  }
  sprintf(here, "%d\n", iv[j]);
}
  
int fieldtype::
writeInts(string fname) {
  int i;
  char linebuf[20];
  ofstream ofstr(fname.c_str(), ios_base::out);
  for (i = 0; i < iv.size(); i++) {
    sprintf(linebuf, "%d\n", iv[i]);
    ofstr << linebuf;
  }
  ofstr.close();
  return 0;
}

int fieldtype::
writeDInts(string fname) {
  int i;
  char linebuf[20];
  ofstream ofstr(fname.c_str(), ios_base::out);
  for (i = 0; i < div.size(); i++) {
    sprintf(linebuf, "%lld\n", div[i]);
    ofstr << linebuf;
  }
  ofstr.close();
  return 0;
}


int fieldtype::
writeQInts(string fname) {
  int i;
  char linebuf[40];
  ofstream ofstr(fname.c_str(), ios_base::out);
  for (i = 0; i < qv.size(); i++) {
    sprintf(linebuf, "%llu%llu\n", qv[i].top, qv[i].bottom);
    ofstr << linebuf;
  }
  ofstr.close();
  return 0;
}

int fieldtype::
writeFloats(string fname) {
  int i;
  char linebuf[20];
  ofstream ofstr(fname.c_str(), ios_base::out);
  for (i = 0; i < fv.size(); i++) {
    sprintf(linebuf, "%f\n", fv[i]);
    ofstr << linebuf;
  }
  return 0;
}

int fieldtype::
writeDoubles(string fname) {
  int i;
  char linebuf[20];
  ofstream ofstr(fname.c_str(), ios_base::out);
  for (i = 0; i < dv.size(); i++) {
    sprintf(linebuf, "%f\n", dv[i]);
    ofstr << linebuf;
  }
  return 0;
}

int fieldtype::
writeIntVecsTxt(string fname) {
  int i, j;
  char * linebuf, * here;
  linebuf = new char[BUFSIZE];
  ofstream ofstr(fname.c_str(), ios_base::out);
  for (i = 0; i < im.size(); i++) {
    ivector & ivv = im[i];
    here = linebuf;
    for (j = 0; j < (ivv.size()-1); j++) {
      sprintf(here, "%d ", ivv[j]);
      here += strlen(here);
    }
    sprintf(here, "%d\n", ivv[j]);
    ofstr << linebuf;
  }
  return 0;
}

int fieldtype::
writeIntVecs(string fname) {
  int i, j, nrows, fmt, ncols, nnz, v;
  ofstream ofstr(fname.c_str(), ios_base::out | ios_base::binary);
  ncols = 2;
  nnz = 0;
  for (i = 0, nrows = 0; i < im.size(); i++) nrows += im[i].size();
  fmt = 110;
  ofstr.write((const char *)&fmt, 4);
  ofstr.write((const char *)&nrows, 4);
  ofstr.write((const char *)&ncols, 4);
  ofstr.write((const char *)&nnz, 4);
  for (i = 0; i < im.size(); i++) {
    for (j = 0; j < im[i].size(); j++) {
      ofstr.write((const char *)&i, 4);
    }
  }
  for (i = 0; i < im.size(); i++) {
    ivector & ivv = im[i];
    for (j = 0; j < ivv.size(); j++) {
      v = ivv[j];
      ofstr.write((const char *)&v, 4);
    }
  }
  ofstr.close();
  return 0;
}

int fieldtype::
writeDIntVecs(string fname) {
  int i, j;
  char * linebuf, * here;
  linebuf = new char[BUFSIZE];
  ofstream ofstr(fname.c_str(), ios_base::out);
  for (i = 0; i < dim.size(); i++) {
    divector & divv = dim[i];
    here = linebuf;
    for (j = 0; j < (divv.size()-1); j++) {
      sprintf(here, "%lld ", divv[j]);
      here += strlen(here);
    }
    sprintf(here, "%lld\n", divv[j]);
    ofstr << linebuf;
  }
  return 0;
}


#ifdef WITHGZS
    

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

ostream * open_out(string ofname) {
  ostream * ofstr = NULL;
  int opened = 0;
  if (ofname.rfind(".gz") == ofname.length() - 3) {
    ofstr = new ogzstream(ofname.c_str(), ios_base::out);
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

ostream * open_out_buf(string ofname, char * buffer, int buffsize) {
  ostream * ofstr = NULL;
  int opened = 0;
  if (ofname.rfind(".gz") == ofname.length() - 3) {
    ofstr = new ogzstream(ofname.c_str(), ios_base::out);
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
  return ofstr;
}
  
#endif
