/*******************************************************************************
    * SMATR 2.O :  Standardised Makjor Axis Tests and Routines
    * http://www.bio.mq.edu.au/ecology/SMATR/
    * main.cpp                                  *
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

#include "SMAstats.h"

int main(void)
{
bool hetero;
ofstream Outfile; 
string out_name;
//Menu data structure
Menu menu;               

while(!menu.isExit())
    {menu.set();
    if(menu.isExit()) exit(1);
    SMAstats stats;
   //output file operations
    if(Outfile.is_open() && out_name!=menu.getOutfilenme()) // if name of output file changed - open new file
         Outfile.close();
    if(!Outfile.is_open()) //if no file open open one
       {out_name=menu.getOutfilenme();
        Outfile.open(out_name.c_str()); if(!Outfile) {cerr<<"File "<<out_name<<" not opened"<<endl; system("pause"); exit(EXIT_FAILURE);}}
    //input data from file fit SMA slopes to groups
    stats.input(menu.getProfile(),Outfile);
    //fit lines to individual groups
    stats.fit_individual(menu.getProfile()); 
    //print data to file
    stats.print_individual_details(menu.getProfile(),Outfile);  Outfile<<endl;
    //test for common slope 
    hetero = stats.test_common(menu.getProfile(),Outfile);
    if(!hetero)     //run post-hoc comparions
                stats.slope_postHoc(menu.getProfile(),Outfile);
    else           //run ANCOVA tests for common slope
                stats.ANCOVA(menu.getProfile(),Outfile);
    //print line break
    for(int i=0; i<100; i++) Outfile<<"-"; Outfile<<"\n"<<endl;
    menu.reEnter();
    }
Outfile.close();
return 0;
}
