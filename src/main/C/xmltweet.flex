/* Scanner for Twitter XML with emoticons */

%{
  extern int checkword(char *);
  extern int numlines;
%}

%option never-interactive
%option noyywrap

LETTER	   [a-zA-Z_]
DIGIT	   [0-9]
PUNCT	   [;:,.?!]

%% 

-?{DIGIT}+    {
  int iv = checkword(yytext);
       }	     
     
-?{DIGIT}+"."{DIGIT}*   {
  int iv = checkword(yytext);
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
	  if (numlines % 1000 == 0) {
	  fprintf(stderr, "\r%05d lines", numlines);
	  }	  
	  }


%%

