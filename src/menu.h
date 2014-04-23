/*******************************************************************************
    * SMATR 2.O :  Standardised Makjor Axis Tests and Routines
    * http://www.bio.mq.edu.au/ecology/SMATR/
    * Defnition of menu class 
    * Date   : 24/10/06                                 
    * Copyright (C) 2006 Daniel Falster
    *
    * This program is free software; you can redistribute it and/or
    * modify it under the terms of the GNU General Public License
    * as published by the Free Software Foundation, version 2, as 
    * stated at http://www.gnu.org/licenses/gpl.html
    * 
    * This program is distributed in the hope that it will be useful,
    * but WITHOUT ANY WARRANTY; without even the implied warranty of
    * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    * GNU General Public License for more details.
    * 
    * Contact details: Daniel Falster, dfalster [at] bio.mq.edu.au
    * Dept Biological Sciences, Macquarie University 2109, Australia
*****************************************************************************/

#include <string>
#include <vector>
#include <map>
#include <list>
#include <iostream>
#include <fstream>
#include <iomanip> // for formatting
#include <math.h>
#include <stdio.h>

#define allow_ftest 0      //Allow Ftest instead of Walds for comparison of fitted and resiudal values
#ifdef _WXOSX
      #define NEWLINE '\r'     //MAC OSX newline character
#else
     #define NEWLINE '\n'      //PC newline character
#endif
using namespace std;

//USEFUL FUNCTIONS 
void Print_String(string name, ofstream &file, int width);
double StrTodouble(string one);
bool StrCompare(string number1, string number2);
bool IsInteger(string number);
bool isDecimal(string number);
bool isBetween(string number, int low, int upp);
bool isBetween_dec(string number, float low, float upp);
string read_file(ifstream& File, char stop);

//data structures used in menu calss
typedef vector<string> Variables;
typedef map<string,int> Map;
typedef vector<Map> Filters;
struct profile{vector<Variables> data; string fname; int no_lines; int transformation; int stat_choice; int WALD;  int Intercept; int ME; int iterations; int CI; double p; int p_method; int exclusions; Filters filter; double b; double a;};

class Menu {
  public:
     Menu();                          // constructor
     void set();
     void display();
     void display_input();
     bool isExit();
     bool inputData();
     void reEnter();
     void display_data();
     void clear_profile();
     profile& getProfile();
     string getOutfilenme();
  private:
     profile sim_profile;
     string inFilename, outFilename, in2;
     string string_in;
     char action;
     Filters::iterator FilI;
     Map::iterator MI;
     char redo_filters, redo_files;
};


