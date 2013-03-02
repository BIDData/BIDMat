/* Scanner for Twitter XML with emoticons */

%{
  extern int checkword(char *);
  extern void addtok(int tok);
  extern int parsedate(char * str);
  extern int numlines;
  
%}

%option never-interactive
%option noyywrap

LETTER	   [a-zA-Z_]
DIGF	   [0-9][0-9][0-9][0-9]
DIGT	   [0-9][0-9]
DIGIT	   [0-9]
PUNCT	   [;:,.?!]

%% 

-?{DIGIT}+    {
  long long iv = strtoll(yytext, NULL, 10);
  addtok(iv);
  iv = iv >> 31;
  if (iv > 0 || iv < -1) {
    addtok(iv);
  }
}	     
     
-?{DIGIT}+"."{DIGIT}*   {
  float f = strtof(yytext, NULL);
  int iv = *((int *)(&f));
  addtok(iv >> 1);
}

{DIGF}"-"{DIGT}"-"{DIGT}"T"{DIGT}":"{DIGT}":"{DIGT}("-"|"+"){DIGT}":"{DIGT}       {
  int tt = parsedate(yytext);
  addtok(tt);
}

{LETTER}+    {
  int iv = checkword(yytext);
	}

"<"{LETTER}+">"    {
  int iv = checkword(yytext);
	}

"</"{LETTER}+">"    {
  int iv = checkword(yytext);
	}

[:;]-[>)] {
  int iv = checkword(yytext);
	  }

[:;]-[<(] {
  int iv = checkword(yytext);
	  }

{PUNCT}	  {
  int iv = checkword(yytext);
	  }

"..""."*  {
  char ell[] = "...";
  int iv  = checkword(ell);
	  }

[\n]	  {
	  numlines++;
	  if (numlines % 1000000 == 0) {
	  fprintf(stderr, "\r%05d lines", numlines);
	  }	  
	  }

.         {}

%%

