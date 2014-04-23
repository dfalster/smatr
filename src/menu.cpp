/*******************************************************************************
    * SMATR 2.O :  Standardised Makjor Axis Tests and Routines
    * http://www.bio.mq.edu.au/ecology/SMATR/
    * Implementation of menu class 
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
#include "menu.h"




//USEFUL FUNCTIONS 
//Print strings of a given width, where reaminaing length is blank space. 
void Print_String(string name, ofstream &file, int width)
    {file << name; for(int i=name.length(); i< width; i++) file<<' ';  return;}
//Converts a string to double
double StrTodouble(string one)
    {char temp[12]; int i; for(i=0; i<one.length(); i++) temp[i]=one[i]; temp[i]='\0'; return atof(temp);}
//Compares two strings as doubleing point numbers
bool StrCompare(string number1, string number2)    
       {return bool(StrTodouble(number1)> StrTodouble(number2));}
       
bool isInteger(string number)
     {for(int i=0; i<number.length(); i++)
             if(!isdigit(number[i])) return 0;
     return 1;}

bool isBetween(string number, int low, int upp)
     {if(!isInteger(number)) return 0;
     int i = atoi(number.c_str());
     return bool(i>= low && i<=upp);}

bool isBetween_dec(string number, float low, float upp)
     {if(!isDecimal(number)) return 0;
     float f = atof(number.c_str());
     return bool(f>= low && f<=upp);}
     
bool isDecimal(string number)
     {for(int i=0; i<number.length(); i++)
             {switch(i){
                     case 0:   if(!isdigit(number[i]) && number[i]!='.' && number[i]!='-' ) return 0; break;
                     default:  if(!isdigit(number[i]) && number[i]!='.') return 0;                    break;}
             }
     return 1;}
     
//scans a sequence from file until stop character or end of line is reached       
string read_file(ifstream& File, char stop)
    {string temp;    char ch;
    File.get(ch); 
    while(!File.eof() && ch!=NEWLINE && ch!=stop ) {
                      temp+=ch; File.get(ch);
                   //   if(ch==NEWLINE) cout<<"R "<<endl;
                      } return temp;}

Menu::Menu()
       {
       in2=inFilename="input.txt";
       outFilename="output.txt";
       sim_profile.no_lines=0;
       sim_profile.stat_choice=1;
       sim_profile.CI =95;
       sim_profile.p=0.05;
       sim_profile.WALD=1;
       sim_profile.ME=0;
       sim_profile.Intercept=1;
       sim_profile.exclusions=6;
       sim_profile.p_method =1;
       sim_profile.iterations=1000;
       sim_profile.transformation=2;
       sim_profile.b='x'; //hypothesised slope
       sim_profile.a='x'; //hypothesised intercept
       redo_filters='y';
       redo_files='n';
       Variables temp;
       temp.insert(temp.end(),"none");
       for(int i=0; i<4; i++)sim_profile.data.insert(sim_profile.data.end(),temp);
       action='0';
       }

void Menu::set(void)
       {
       while(action!='l'&&!isExit())
           {
           display();
           cin >> string_in; 
           action = tolower(string_in[0]);
           switch(action)
                {
                case'i':
                        cout<<"\nEnter name of input filename-> ";
                        in2=inFilename;
                        cin >>inFilename;
                        if(inFilename==in2) redo_files='n'; else redo_files='y';
                        break;
                case'o':
                        cout<<"\nEnter name of output filename-> ";
                        cin >>outFilename;
                        break;
                case'm':
                        cout<<"Choose the type of line you wish to fit: \n\t0. OLS - Ordinary Least Squares Regression\n\t1. SMA - standardised major axis (aka RMA)\n\t2. MA - major axis"<<endl;
                        cin >>string_in;
                        while(!isBetween(string_in, 0, 2))
                                   {cout <<"Bad entry, try again:\t"; cin >>string_in;}
                        sim_profile.stat_choice = atoi(string_in.c_str());
                        break;
                #if allow_ftest
                case'w':
                        cout<<"Choose routines for comparing Fitted & Residuals: \n\t0. F statistic (ANOVA)\n\t1. WALD Statistic (recommended)"<<endl;
                        cin >>string_in;
                        while(!isBetween(string_in, 0, 1))
                                   {cout <<"Bad entry, try again:\t"; cin >>string_in;}
                        sim_profile.WALD = atoi(string_in.c_str());
                        break; 
                #endif
                case't':
                        cout<<"Include intercept: \n\t0. No (force through origin)\n\t1. Yes (default)"<<endl;
                        cin >>string_in;
                        while(!isBetween(string_in, 0, 1))
                                   {cout <<"Bad entry, try again:\t"; cin >>string_in;}
                        sim_profile.Intercept= atoi(string_in.c_str());
                        break; 
                case'e':
                        cout<<"Include measurement error: \n\t0. No \n\t1. Yes"<<endl;
                         cin >>string_in;
                        while(!isBetween(string_in, 0, 1))
                                   {cout <<"Bad entry, try again:\t"; cin >>string_in;}
                        sim_profile.ME = atoi(string_in.c_str());
                         break; 
                case'c':
                        cout<<"\nEnter confidence level (1-100) -> ";
                        cin >>string_in;
                        while(!isBetween(string_in, 0, 100))
                                   {cout <<"Bad entry, try again:\t"; cin >>string_in;}
                        sim_profile.CI= atoi(string_in.c_str());
                        break;
                case'p':
                        cout<<"\nEnter critical p value -> ";
                        cin >>string_in;
                        while(!isBetween_dec(string_in, 0.0, 1.0))
                                   {cout <<"Bad entry, try again:\t"; cin >>string_in;}
                        sim_profile.p= atof(string_in.c_str());
                        break;
                case'r':
                        cout<<"p-value for common slope test calculated with resampling or chi-squared distn: \n\t0. Chi-squared\n\t1. Resampling (default)"<<endl;
                        cin >>string_in;
                        while(!isBetween(string_in, 0, 1))
                                   {cout <<"Bad entry, try again:\t"; cin >>string_in;}
                        sim_profile.p_method = atoi(string_in.c_str());
                        break; 
                case'g':
                        cout<<"\nEnter minimum group size (1-50) -> ";
                        cin >>string_in;
                        while(!isBetween(string_in, 1, 50))
                                   {cout <<"Bad entry, try again:\t"; cin >>string_in;}
                        sim_profile.exclusions= atoi(string_in.c_str());
                        break;
                case'n':
                        cout<<"\nEnter no iterations -> ";
                        cin >>string_in;
                        while(!isInteger(string_in))
                                   {cout <<"Bad entry, try again:\t"; cin >>string_in;}
                        sim_profile.iterations = atoi(string_in.c_str());
                        break;
                case 'x':
                        cout<<"\nGoodBye"<<endl;
                        break;
                case 'l':
                         if(inputData()==0)
                             cout<<"Data failure!\n";
                         break;
                default:break;
                }
           }
       return;
       }

void Menu::display(void)
     {
      #ifdef _WIN32
      system("cls");
      #endif
      cout<<"\nSMATR: Standardised Major Axis Tests & Routines, v2.0\nby D Falster, I Wright & D Warton 2006\n\nCURRENT SETTINGS"<<endl;
      cout<<"I - Input File:       "<<inFilename<<endl;
      cout<<"O - Output File:      "<<outFilename<<endl;
      cout<<"M - Method:           "; if(sim_profile.stat_choice==0) cout<<"OLS"; else if(sim_profile.stat_choice==1) cout<<"SMA"; else cout << "MA"; cout <<endl;
      cout<<"T - inTercept:        "; if(sim_profile.Intercept==1) cout<<"Yes"; else cout<<"No"; cout <<endl;
      cout<<"E - measuremet Error: "; if(sim_profile.ME==1) cout<<"Yes"; else cout<<"No"; cout <<endl;
      #if allow_ftest
      cout<<"W - Comparions:       "; if(sim_profile.WALD==1) cout<<"WALD"; else cout<<"ANOVA"; cout <<endl;
      #endif
      cout<<"C - Confidence int:   "<<sim_profile.CI<<"%"<<endl;
      cout<<"G - min Group size:   "<<sim_profile.exclusions<<endl;
      cout<<"R - Resample:         "; if(sim_profile.p_method==1) cout<<"Yes"; else cout<<"No"; cout <<endl;    
      cout<<"N - No. iterations:   "; if(sim_profile.p_method==1) cout<<sim_profile.iterations; else cout<<"na"; cout<<endl;
      cout<<"P - critical P-value: "<<sim_profile.p<<endl;
      cout<<"\nChoose a letter to change current setting. Press L to load data or x to exit"<<endl;
      return;
      }

void Menu::display_input(void)
      {int n=3;
      #ifdef _WIN32
      system("cls");
      #endif
      cout<<"\nDATA INPUT - settings"<<endl;
      cout<<"G - Grouping variable:  "; cout<< sim_profile.data[0][0]<<endl;
      cout<<"Y - y variable:         "; cout<< sim_profile.data[1][0]<< endl;
      cout<<"X - x variable:         "; cout<< sim_profile.data[2][0]<<endl;
      cout<<"F - Filter variable:    "; while(n < sim_profile.data.size()){cout<<sim_profile.data[n][0]<<' '; n++;} cout<<endl;
      cout<<"T - Transformation:     "; if(sim_profile.transformation==0) cout<<"log Y vs log X"; else if (sim_profile.transformation==1) cout << "ln Y vs ln X"; else cout<< "linear Y vs linear X"; cout <<endl;
      cout<<"B - H0 slope:           "; if(sim_profile.b!='x') cout <<sim_profile.b<<endl; else cout<<"none"<<endl;
      if(sim_profile.Intercept){cout<<"A - H0 intercept:       "; if(sim_profile.a!='x') cout <<sim_profile.a<<endl; else cout<<"none"<<endl;}
      cout<<"D - Display data        ";
      cout<<"\nChoose a letter to change current setting. Press L to load data"<<endl;
      return;
      }

bool Menu::isExit(void)
       {return bool(action=='x');}

void Menu::reEnter(void)
       {action=1;
       return;}

profile& Menu::getProfile(void)
       {return sim_profile;}

void Menu::clear_profile(void)
     {
     sim_profile.no_lines=0;
     redo_filters='y';
     redo_files='n';
     Variables temp;
     temp.insert(temp.end(),"none");
     sim_profile.data.erase(sim_profile.data.begin(),sim_profile.data.end());
     for(int i=0; i<4; i++)sim_profile.data.insert(sim_profile.data.end(),temp);
     return;}
     

bool Menu::inputData()
       {
       int n=0;
       string line;
       action ='0';
       char ch;
       Variables variables_head;  //variable names at top of file
       vector<Variables>::iterator K;
       
       if(redo_files=='y') clear_profile();
       //open file
       ifstream inFile(inFilename.c_str(), ios::in);
       if(!inFile) {cerr<<"File "<<inFilename<<" not opened"<<endl; 
                   #ifdef _WIN32
                   system("pause");
                   #endif
                   exit(EXIT_FAILURE);}
       sim_profile.fname = inFilename;
       //get variables names
       Variables::iterator I =variables_head.begin();
       while(inFile.peek()!=NEWLINE)
            {ch=inFile.get();
            if(ch=='\t')     {I=variables_head.insert(variables_head.end(),line);line.erase();}
            else line=line+ch;}
       ch=inFile.get(); n++;;
       I=variables_head.insert(variables_head.end(),line);line.erase();
       //count number lines
       if(sim_profile.no_lines==0)  
            {
         //   while(getline(inFile,line)) n++; 
            while(!inFile.eof())
                {n++;
                line = read_file(inFile, NEWLINE);   }
            sim_profile.no_lines=n;
            }
      
       Variables data_line;
       //input data
       while(action!='l')
           {
          // string temp;
           int x;
           display_input();
            
           cin >> string_in; 
           action = tolower(string_in[0]);
           //action=tolower(cin.get());
           if(action!='l' && action!='a' && action!='b' && action!='d'  ) //display variables
                {cout<< "Variable list"<<endl;for(int i=0; i<variables_head.size(); i++) cout << i <<": "<< variables_head[i]<<"\t"; cout<<"\n"<<endl; }
           switch(action)
                {case'd':
                        display_data();
                        break;
                case'g':
                        cout <<"Choose grouping variables (press c to clear)"<<endl;
                        cin>> string_in;
                        data_line.erase(data_line.begin(), data_line.end());
                        if(string_in[0]=='c')
                              {K=sim_profile.data.begin();K=sim_profile.data.erase(K);
                               data_line.insert(data_line.end(),"none");
                               K=sim_profile.data.insert(K, data_line);}
                        else
                           {
                           while(!isBetween(string_in, 0, variables_head.size()-1))
                                   {cout <<"Bad entry, try again:\t"; cin >>string_in;}
                           x= atoi(string_in.c_str());
                           //get data from file
                           data_line.erase(data_line.begin(), data_line.end());
                           //return to start of file
                           inFile.clear(); inFile.seekg(0, ios::beg);                   
                           for(int i=0; i< sim_profile.no_lines; i++)                         //do for each row
                                   {for(int j=0; j< x; j++)       read_file(inFile, '\t'); //scan to appropriate column
                                    data_line.push_back(read_file(inFile, '\t'));          //read data
                                    if(x<variables_head.size()-1) read_file(inFile, NEWLINE); //read to end of line
                                   }
                           
                           if(sim_profile.data[0][0]=="none")
                               {K=sim_profile.data.begin();K=sim_profile.data.erase(K);
                                K=sim_profile.data.insert(K, data_line);
                                }
                           else
                               for(int i=0; i< sim_profile.no_lines; i++)
                                   {sim_profile.data[0][i]+=' ';
                                    sim_profile.data[0][i]+= data_line[i];}
                           }
                        break;
                case'y':
                        cout <<"Choose Y variable"<< endl;
                        cin >>string_in;
                           while(!isBetween(string_in, 0, variables_head.size()-1))
                                   {cout <<"Bad entry, try again:\t"; cin >>string_in;}
                        x= atoi(string_in.c_str());
                        data_line.erase(data_line.begin(), data_line.end());
                        //return to start of file
                        inFile.clear(); inFile.seekg(0, ios::beg);                     
                        for(int i=0; i< sim_profile.no_lines; i++)                         //do for each row
                                   {
                                   for(int j=0; j< x; j++)       read_file(inFile, '\t'); //scan to appropriate column
                                   data_line.push_back(read_file(inFile, '\t'));          //read data
                                   if(x<variables_head.size()-1) read_file(inFile, NEWLINE); //read to end of line
                                   }
                        K=sim_profile.data.begin();K++;
                        K=sim_profile.data.erase(K);
                        sim_profile.data.insert(K, data_line);
                        break;
                case'x':
                        cout <<"Choose X variable"<< endl;
                        cin >>string_in;
                        while(!isBetween(string_in, 0, variables_head.size()-1))
                                   {cout <<"Bad entry, try again:\t"; cin >>string_in;}
                        x= atoi(string_in.c_str());
                        //get data from file
                        data_line.erase(data_line.begin(), data_line.end());
                        //return to start of file
                        inFile.clear(); inFile.seekg(0, ios::beg);                      
                        for(int i=0; i< sim_profile.no_lines; i++)                         //do for each row
                                   {for(int j=0; j< x; j++)       read_file(inFile, '\t'); //scan to appropriate column
                                    data_line.push_back(read_file(inFile, '\t'));          //read data
                                    if(x<variables_head.size()-1) read_file(inFile, NEWLINE); //read to end of line
                                   }
                        K=sim_profile.data.begin();K++; K++;
                        K=sim_profile.data.erase(K);
                        sim_profile.data.insert(K, data_line);
                        break;
                        
                case'f':
                        redo_filters='y';
                        while(sim_profile.filter.size()!=0)            //erase existing filters
                                   sim_profile.filter.erase(sim_profile.filter.begin());
                        cout <<"Choose filter variables (press c to clear)"<< endl;
                        cin>> string_in;
                        data_line.erase(data_line.begin(), data_line.end());
                        if(tolower(string_in[0])=='c')
                              {while(sim_profile.data.size()>3)
                                   {K=sim_profile.data.end(); K--;K=sim_profile.data.erase(K);}
                               data_line.insert(data_line.end(),"none");
                               K=sim_profile.data.insert(K, data_line);
                               }
                        else
                           {
                           while(!isBetween(string_in, 0, variables_head.size()-1))
                                   {cout <<"Bad entry, try again:\t"; cin >>string_in;}
                           x= atoi(string_in.c_str());
                          //get data from file
                           data_line.erase(data_line.begin(), data_line.end());
                           //return to start of file
                          inFile.clear(); inFile.seekg(0, ios::beg);                      
                          for(int i=0; i< sim_profile.no_lines; i++)                         //do for each row
                                    {for(int j=0; j< x; j++)       read_file(inFile, '\t'); //scan to appropriate column
                                    data_line.push_back(read_file(inFile, '\t'));          //read data
                                    if(x<variables_head.size()-1) read_file(inFile, NEWLINE); //read to end of line
                                   }
                                
                           if(sim_profile.data[3][0]=="none")
                               {K=sim_profile.data.end(); K--;
                                K=sim_profile.data.erase(K);
                                K=sim_profile.data.insert(K, data_line);}
                           else
                               {K=sim_profile.data.end(); K--;
                               K=sim_profile.data.insert(K, data_line);}
                           }
                        break;
                 case't':
                        cout<<"Choose a transformation:\n\t0. log Y vs log X\n\t1. ln Y vs ln X\n\t2. linear Y vs linear X"<<endl;
                        cin >>string_in;
                        while(!isBetween(string_in, 0, 2))
                                   {cout <<"Bad entry, try again:\t"; cin >>string_in;}
                        sim_profile.transformation= atoi(string_in.c_str());
                        break;
                 case'b':
                        cout<<"\n Enter hypothesised slope value HO: b = ";
                        cin >>string_in;
                        while(!isDecimal(string_in))
                                   {cout <<"Bad entry, try again:\t"; cin >>string_in;}
                        sim_profile.b= atof(string_in.c_str());
                        break;
                 case'a':
                        cout<<"\n Enter hypothesised intercept value HO: a = ";
                        cin >>string_in;
                        while(!isDecimal(string_in))
                                   {cout <<"Bad entry, try again:\t"; cin >>string_in;}
                        sim_profile.a= atof(string_in.c_str());
                        break;
                 case 'l':
                         if(sim_profile.data[1][0]=="none" || sim_profile.data[2][0]=="none")
                             {cout<< "\nCan't load data because missing X or Y variable"<<endl;
                             action = 'b'; 
                             #ifdef _WIN32
                             system("pause");
                             #endif
                             break;}
                         #ifdef _WIN32
                         system("cls");
                         #endif
                          //if no group variables fill group column with "none"
                         if(sim_profile.data[0][0]=="none")
                            while(sim_profile.data[0].size() < sim_profile.no_lines)
                                    sim_profile.data[0].insert(sim_profile.data[0].begin(),"none");
                         char action2;
                         //if needed create filter maps
                         if(sim_profile.data[3][0]!="none" && redo_filters=='y')//are there any filter variables?
                             {
                             for(int n =3; n < sim_profile.data.size(); n++) //do for each variable
                                 {
                                 cout<<"\nFiltering variable " << sim_profile.data[n][0]<<"\nValues:\t";
                                 for(int i=1;i < min(sim_profile.no_lines,30); i++) cout<<sim_profile.data[n][i]<<" ";
                                 cout<<"\n\n Is it:\n\tC - categorical";
                                 cout<<"\n\tS - scalar"<<endl;
                                 cin>> action2;
                                 if(action2=='s') //for scaler variables tranform to categorical
                                     {string min, max, Temp;
                                      cout<<"Enter minimum: "; cin >> min;
                                      cout<<"Enter maximum: "; cin >> max;
                                      for(int i=1;i < sim_profile.no_lines; i++)
                                           {if(StrCompare(sim_profile.data[n][i],max))
                                                 sim_profile.data[n][i]= "> " +max;
                                           else if(StrCompare(min,sim_profile.data[n][i]))
                                                 sim_profile.data[n][i]= "< " +min;
                                           else
                                               {Temp.erase();Temp= "> " + min; Temp+= " & < "; Temp+=max;
                                               sim_profile.data[n][i]= Temp;}
                                            }
                                     }
                                 Map temp;
                                 cout<<"Choose y or n:\n";
                                 for(int i=1;i < sim_profile.no_lines; i++)
                                      {
                                      MI =temp.find(sim_profile.data[n][i]);
                                      if(MI==temp.end())
                                          {cout << "remove " <<sim_profile.data[n][i]<<"? ";
                                          cin >> action2;
                                          if (action2== 'y') temp.insert(make_pair(sim_profile.data[n][i], 0));
                                          else temp.insert(make_pair(sim_profile.data[n][i],1));}
                                      }
                                 sim_profile.filter.insert(sim_profile.filter.end(),temp);
                                 }
                             }
                         redo_filters='n';
                         break;
                 default:break;
                }
        }
     inFile.close();   
     return 1;}

void Menu::display_data()
    {
     int i=0, f=0;
     char c;
     cin.get(c);
     while(i < sim_profile.no_lines)
             {cout<<"Line\tGroups\tY\tX\tFilters\n";
             for(i=f; i < min(f+21, sim_profile.no_lines); i++)
                  {
                  cout<<i<<"\t";
                  for(int j=0; j<sim_profile.data.size(); j++)
                          if(sim_profile.data[j][0]=="none") cout << "-\t";
                          else cout <<sim_profile.data[j][i]<<"\t";
                  cout <<endl;
                  }
             f+=21;
             cout<<"\nPress enter to continue"; cin.get(c);
             #ifdef _WIN32
             system("cls");
             #endif
             }
   }
   
string Menu::getOutfilenme()
    {return outFilename;}




