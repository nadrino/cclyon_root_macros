#ifndef KEYREAD_C 
#define KEYREAD_C

#include "keyread.h"

//class to read inputs from card files

using namespace std;

keyread::keyread(const char* afile){
  fname = afile;
}

double keyread::getKeyD(string key){
  return dmap[key];
}

float keyread::getKeyF(string key){
  return fmap[key];
}

int keyread::getKeyI(string key){
  return imap[key];
}

TString keyread::getKeyS(string key){
  return smap[key];
}

void keyread::readFile(){
  ifstream file(fname);
  if (!file.is_open()){
    cout<<fname<<" is not a good file "<<endl;
    return;
  }
  char line[LINESIZEMAX];
  int  nread;
  TString sline;
  int nline=0;
  while(!file.eof()){
    nline++;
    file.getline(line,LINESIZEMAX);
    for (int i=0;i<file.gcount();i++){
      sline.Append(line[i]);
    }
    processLine(sline);
    sline.Clear();
  }
}

void keyread::processLine(TString sline){
  TString stype;
  TString skey;
  TString sval;
  int ival;
  float fval;
  double dval;
  ///get key type
  if (sline(0)=='i') stype="i";
  else if (sline(0)=='f') stype="f";
  else if (sline(0)=='d') stype="d";
  else if (sline(0)=='s') stype="s";
  else if (sline(0)=='$'){
    return;
  }
  else if (sline(0)==' '){
    return;
  }
  else {
    return;
  }
  ///get key name
  int ichar=2;
  while (sline(ichar)!='='){
    // ignore spaces
    if(sline(ichar)!=' ') skey.Append(sline(ichar));
    ichar++;
  }
  ///get key value
  int jchar=ichar+1;
  while (sline(jchar)!=';'){
    // ignore spaces
    if(sline(ichar)!=' ') sval.Append(sline(jchar));
    jchar++;
  }
  //place key value in appropriate map
  if (stype(0)=='i'){
    ival = sval.Atoi();
    imap[skey.Data()] =  ival;
  }
  if (stype(0)=='f'){
    fval = sval.Atof(); 
    fmap[skey.Data()]=fval;
  }
  if (stype(0)=='d'){
    dval=sval.Atof();
    dmap[skey.Data()]=dval;
  }
  if (stype(0)=='s'){
    smap[skey.Data()]=sval;
  }
   
}




#endif 
