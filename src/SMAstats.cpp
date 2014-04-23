/*******************************************************************************
    * SMATR 2.O :  Standardised Makjor Axis Tests and Routines
    * http://www.bio.mq.edu.au/ecology/SMATR/
    * Implementation of class SMAstats                                   *
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

//Requires GNU GSL v1.8 
// Download from http://sourceforge.net/projects/gnuwin32/

//requirements for linear algebra using GNU GSL 
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#define using_GNU 1   // 1 = using GNU GSL, 0 = Rmath lib
#if using_GNU
    //requirements if using GNU GSL for p-values
    #include <gsl/gsl_randist.h>
    #include <gsl/gsl_cdf.h>
    #define pinv_df_max 200         //
#else        
   //requirements for R statistical library used for probability distributions
   #define MATHLIB_STANDALONE 1
   #include <R/Rmath.h>
   #define HAVE_EXPM1 1
   #define HAVE_LOG1P 1
   // Functions required by R stats library for library to function
   double expm1(double x) {return(exp(x)-1);} //calculates exponential of small numbers. 
   double log1p(double x) {return(log(1+x));} //calculates log of small numbers. Required by R stats library for library to function
#endif             

//USEFUL FUNCTIONS 
int signDF(double value) {if (value<0.0) return -1;  else return 1;}//tests sign of a decimal, returns 1 for +ve, -1 for -ve

//Data structures used in POST-HOC routines to make an ordered list of groups. 
struct PH_data {double slope; int groupNo; string group; string letters; int n;};
typedef list<PH_data> Grouplist;

//IMPLEMENATION OF CLASS-------------------------------------------------------------

SMAstats::SMAstats() {max_namelength=7;}

void SMAstats::input(profile& data, ofstream &file)
       {
       #ifdef _WIN32
       system("cls");
       #endif
       int j, in, count=0;
       Group::iterator GrpI;
       Map::iterator FilI;
       
       summary_stats temp;
       //add single element
       temp.rawX.insert(temp.rawX.begin(), 0);  temp.rawY.insert(temp.rawY.begin(), 0);            

       //put data into maps
       groups.clear(); //to remove any previous data
       for(int i=1;i < data.no_lines; i++)
           {in =bool(data.data[1][i].size()!= 0 && data.data[2][i].size()!=0 && data.data[0][i].size()!=0); //check if empty cell??
           if(data.data[0][i].size()==1)   in=bool(data.data[1][i]!=" " && in);//remove data if grouping variable missing
           if(data.data[1][i].size()==1)   in=bool(data.data[1][i]!=" " && in); //remove data if y variable missing
           if(data.data[2][i].size()==1)   in=bool(data.data[2][i]!=" " && in); //remove data if x yariable missing
           for(j=0; j<data.filter.size(); j++)   in= bool((data.filter[j].find(data.data[3+j][i]))->second && in); //check if data point has been filtered
                    
           //Add data to groups
           if(in)//only use if not filtered
              {GrpI=groups.find(data.data[0][i]);    //search fro group name
              if(GrpI==groups.end()) //no previous entry
                    {
                    temp.rawX[0] = StrTodouble(data.data[2][i]);  temp.rawY[0] = StrTodouble(data.data[1][i]); 
                    groups.insert(make_pair(data.data[0][i], temp));                    
                    if(data.data[0][i].size()> max_namelength) max_namelength= data.data[0][i].size();}
              else //previous entry- add to existsing list
                 {GrpI->second.rawX.insert(GrpI->second.rawX.begin(),StrTodouble(data.data[2][i]));
                  GrpI->second.rawY.insert(GrpI->second.rawY.begin(),StrTodouble(data.data[1][i]));}              
              }
           else count++; //counts number of filtered data pairs
           }
       if(max_namelength%6==0) max_namelength++; //used for formatting. Don't want column width to be divisible by 6 - stuffs up tab formatting
             
       //Print summary infomation to file
       stat_choice = data.stat_choice; 
       intercept_flag = data.Intercept;
       ME_flag = data.ME;
      
       if(stat_choice== 0) file << "OLS analysis";  else if(stat_choice== 1) file << "SMA analysis"; else file << "MA analysis";
       if(!intercept_flag) file << "\t fitted through origin";
       if(ME_flag) file << "\t fitted using Methods of Moments";
       file<<"\nFilename:\t"<<data.fname;
       file<<"\nGrp:\t"<<data.data[0][0];
       if(data.transformation==0) /*log log*/ file<<"\tY: log "<< data.data[1][0]<<"\tX: log "<<data.data[2][0]<<endl;
            else if(data.transformation==1) /*ln ln */ file<<"\tY: ln "<< data.data[1][0]<<"\tX: ln "<<data.data[2][0]<<endl;
            else /*lin lin*/ file<<"\tY: "<< data.data[1][0]<<"\tX: "<<data.data[2][0]<<endl;
       file <<"No. iterations:\t"<<data.iterations<<"\tCritical p:\t"<<data.p<<"\tConf.interval:\t"<<data.CI<<"%"<<endl;
       file << "Filters:\n";
       for(j=0; j<data.filter.size(); j++)
                 {file<<data.data[j+3][0]<<"\t";
                 for(FilI=data.filter[j].begin(); FilI!=data.filter[j].end(); FilI++)
                      {file<<FilI->first<<" ("<<FilI->second<<")\t";}
                 file << endl;}
       file << "Filtered or empty data points:\t"<<count<<"\nGroups with less than "<< data.exclusions<<" data points (removed):\t";

       //remove sites with too few data points
       GrpI=groups.begin(); 
       while(GrpI!=groups.end())
           {GrpI->second.n = GrpI->second.rawX.size();   
            if(GrpI->second.n < data.exclusions)  //remove groups wiht too few data point
                {cout<< "Group "<<GrpI->first<<" erased with "<<GrpI->second.n<<" elements"<<endl;
                file<<GrpI->first<<"\t"; groups.erase(GrpI); 
                GrpI=groups.begin();}
           else GrpI++;}
       count =0;
       //allocate memory and run transformations
       for(GrpI=groups.begin(); GrpI!=groups.end(); GrpI++)
          {  //allocate memroy to  vectors for fitted values
            GrpI->second.Fitted.insert(GrpI->second.Fitted.begin(), GrpI->second.n, 0);   
            GrpI->second.Resid.insert(GrpI->second.Resid.begin(), GrpI->second.n, 0);   
            //run transformations
            for(int i=0;i < GrpI->second.n; i++)//apply transformation if any
                   {switch(data.transformation)
                      {case 0:    //log-log
                              GrpI->second.rawX[i]=log10(GrpI->second.rawX[i]); GrpI->second.rawY[i]=log10(GrpI->second.rawY[i]); break;
                      case 1:    //ln-ln
                             GrpI->second.rawX[i]=log(GrpI->second.rawX[i]); GrpI->second.rawY[i]=log(GrpI->second.rawY[i]); break;
                      case 2:default:break;} }//lin-lin                              
           count +=GrpI->second.n;}//counts total number datapoints

      file<<"\nNo. groups:\t"<<groups.size()<<"\nNo. datapoints:\t" <<count<< endl;
      //enter measurment error
      file << "Measurement Error: ";
      if(ME_flag) 
          {cout<< "Enter error variance estimate for ";  
                  if(data.transformation==0) cout<<"log "; else if(data.transformation==1) cout<<"ln ";
                  cout << data.data[2][0]<<": "; cin>>ME.second; //X
          cout<< "Enter error variance estimate for ";
                 if(data.transformation==0) cout<<"log "; else if(data.transformation==1) cout<<"ln ";
                  cout << data.data[1][0]<<": "; cin>>ME.first; //Y
          file<< data.data[2][0]<< " = "<<ME.second<<"\t"<<data.data[1][0]<< " = "<<ME.first<< endl;}
      else
          {file <<"Not included"<< endl; ME.first=ME.second=0.0; }      
      return; }

void SMAstats::displayData()
       {cout << "Data groups"<<endl;
       for(Group::iterator GrpI=groups.begin(); GrpI!=groups.end(); GrpI++)
           {cout<<GrpI->first<<"\t";
           for(int i=0; i<GrpI->second.rawX.size(); i++)   cout <<"("<<GrpI->second.rawX[i]<<","<<GrpI->second.rawY[i]<<")  ";
           cout<<endl;
           cout<<"Size: "<<GrpI->second.rawX.size()<<"  SXX: "<<GrpI->second.s2_X<<"  Slope: "<<GrpI->second.slope<<"  R2: "<<GrpI->second.R2<<endl;
           }
       }
     
void SMAstats::fit_individual(profile& data) {
     fit_lines(groups, 1);
     fit_CIs(groups, data.CI);}


void SMAstats::print_individual_details(profile& data, ofstream &file)
      {
       double R2res, F,T,df, p; //for slope comparison against hypothesised slope
       W_B = max(7, Dec_B+4); W_A = max(7, Dec_A+4);  //FORMATTING 
     //print summaries to file
       if(data.stat_choice== 0) file << "\nOLS "; else if(data.stat_choice== 1) file << "\nSMA "; else file << "\nMA ";
       file<<"results\n"<<setiosflags(ios::left)<<setw(max_namelength)<<"Group"<<"\t"<<setw(7)<<"n"<<"\t"<<setw(7)<<"R2"<<"\t"<<setw(7)<<"p"<<"\t"<<setw(W_B)<<"Slope"<<"\t"<<setw(W_B)<<"LowCI"<<"\t"<<setw(W_B)<<"UppCI";
       file<<"\t"<<setiosflags(ios::left)<<setw(W_A)<<"Interc"<<"\t"<<setw(W_A)<<"LowCI"<<"\t"<<setw(W_A)<<"UppCI"<<"\t"<<setw(7)<<"Ymean"<<"\t"<<setw(7)<<"Xmean"<<"\t";
       if(data.b!='x') //if hypothesised slope is not 0
           file<<setw(W_B)<<"H0_b"<<"\t"<<setw(7)<<"F"<<"\t"<<setw(7)<<"p"<<"\t";
       if(data.a!='x' && intercept_flag) //if hypothesised slope is not 0
           file<<setw(W_A)<<"H0_a"<<"\t"<<setw(7)<<"T"<<"\t"<<setw(7)<<"p";
       file<<endl;
       
       for (Group::iterator GrpI=groups.begin(); GrpI!=groups.end(); GrpI++) 
           {//print summary stats to file
           Print_String(GrpI->first, file, max_namelength); file<<"\t"<<setprecision(0)<<GrpI->second.n<<"\t"<<setiosflags(ios::fixed)<<setprecision(3)<<setw(7)<<GrpI->second.R2<<"\t"<<setprecision(3)<<setw(7)<<GrpI->second.p<<"\t"<<setprecision(Dec_B)<<setw(W_B)<<GrpI->second.slope<<"\t";
           file<<setw(W_B)<<GrpI->second.slopeCI.first<<"\t"<<setw(W_B)<< GrpI->second.slopeCI.second<<"\t";   
           file<<setw(7)<<setprecision(Dec_A)<<setw(W_A)<<GrpI->second.intercept<<"\t"<<setw(W_A)<<GrpI->second.interceptCI.first<<"\t"<<setw(W_A)<<GrpI->second.interceptCI.second;  
           file<<"\t"<<setprecision(3)<<setw(7)<<GrpI->second.Ymean<<"\t"<<setw(7)<<GrpI->second.Xmean<<"\t";      
           //test if slope diff to hypothesised slope - BUT only do this if hypothesised slope is not 0
           if(data.b!='x')
               {R2res = pow(covar_rf(GrpI->second.s2_Y, GrpI->second.s2_X, GrpI->second.s_XY, data.b, GrpI->second.n, intercept_flag),2)
                        / var_fitted(GrpI->second.s2_Y, GrpI->second.s2_X, GrpI->second.s_XY, data.b, GrpI->second.n, intercept_flag)
                        / var_resid(GrpI->second.s2_Y, GrpI->second.s2_X, GrpI->second.s_XY, data.b, GrpI->second.n, intercept_flag);
                F= R2res/(1-R2res)*(GrpI->second.n-2);
                if(intercept_flag) df = GrpI->second.n-2; else df = GrpI->second.n-1; 
                #if using_GNU      //Using GSL library
                p = gsl_cdf_fdist_Q(F ,1,df); 
                #else        //Using R math library
                p = pf(F ,1,df,0,0);
                #endif
                file <<setw(W_B)<<data.b<<"\t"<<setw(7)<<F<<"\t"<<setprecision(3)<<setw(7)<<p<<"\t";}//Get probability using inverse function for f distribution
           //test if interecpt diff to hypothesised value BUT only do this if hypothesised value is not 0
           if(intercept_flag && data.a!='x')
               {T = (GrpI->second.intercept - data.a)/ sqrt(GrpI->second.s2_A);
                #if using_GNU
                p = 2.0*gsl_cdf_tdist_P(-fabs(T),GrpI->second.n-2); //Using GSL library
                #else        
                p = 2.0*pt(-fabs(T),GrpI->second.n-2,1,0);
                #endif
                file <<setw(W_A)<<data.a<<"\t"<<setw(7)<<T<<"\t"<<setprecision(3)<<setw(7)<<p;}//Get probability using inverse function for t distribution     
           file<<endl;
           }
        }

bool SMAstats::test_common(profile& data, ofstream &file)
     {
     double slope,p, TS, p2, TS2;
     int n=0, precis= (int)ceil(log10((double)data.iterations)); //decimal precision for p-values
     pair<double, double> CIs;
     
     if(groups.size()==1) 
         {cout <<"\nSingle group: no test for common slope.\nSee file ouput for line fitted to all data\n"<<endl;
          file <<"Single group: no test for common slope"<<endl;
          #ifdef _WIN32
          system("pause");
          #endif
          return 1;}
          
     file <<"TEST FOR COMMON SLOPE ACROSS GROUPS"<<endl;
     cout <<"\nTEST FOR COMMON SLOPE ACROSS GROUPS"<<endl;
     //fit common slope and claulate signifcicance
     fit_common(groups, data.iterations, slope);
     TS =common_TestStat(groups, slope,0);
     p=common_pval(groups, data.iterations, slope, data.p_method);
     //if not signifcantly heterogenous, calculate CI for common slope
     if(p>data.p) CIs = commonCI(groups, data.CI, slope, data.iterations);
     
     //count total datapoints
     for(Group::iterator GrpI=groups.begin(); GrpI!=groups.end(); GrpI++) n+=GrpI->second.n;
     //output to file & screen
     file <<setiosflags(ios::left)<<setw(7)<<"Grps"<<"\t"<<setw(7)<<"N"<<"\t"<<setw(W_B)<<"Slope"<<"\t"<<setw(W_B)<<"LowCI"<<"\t"<<setw(W_B)<<"UppCI"<<"\t";
     cout <<setiosflags(ios::left)<<setw(7)<<"Grps"<<"\t"<<setw(7)<<"N"<<"\t"<<setw(W_B)<<"Slope"<<"\t"<<setw(W_B)<<"LowCI"<<"\t"<<setw(W_B)<<"UppCI"<<"\t";
     if(data.b!=0 && p>data.p) //if hypothesised slope is not 0
           {file<<setw(W_B)<<"H0_b"<<"\t"<<setw(7)<<"X2"<<"\t"<<setw(7)<<"p";
            cout<<setw(W_B)<<"H0_b"<<"\t"<<setw(7)<<"X2"<<"\t"<<setw(7)<<"p";}
     file<<"\n"<<setprecision(0)<<setw(7)<<groups.size()<<"\t"<<setw(7)<<n<<"\t"<<setprecision(Dec_B)<<setw(W_B)<<slope<<"\t";
     cout<<"\n"<<setprecision(0)<<setw(7)<<groups.size()<<"\t"<<setw(7)<<n<<"\t"<<setprecision(Dec_B+1)<<setw(W_B)<<slope<<"\t";
     if(p>data.p)
           {file<<setw(W_B)<<CIs.first<<"\t"<<setw(W_B)<<CIs.second<<"\t";
           cout<<setw(W_B)<<CIs.first<<"\t"<<setw(W_B)<<CIs.second<<"\t";
           //test if slope diff to hypothesised slope - BUT only do this if hypothesised slope is not 0
           if(data.b!=0)
                        {TS2= common_TestStat(groups, slope, data.b);
                         #if using_GNU
                         p2 = gsl_cdf_chisq_Q(TS2,1);  //Using GSL library
                         #else        
                         p2 = pchisq(TS2, 1, 0,0); //Using R math library
                         #endif
                        file <<setw(W_B)<<data.b<<"\t"<<setprecision(3)<<setw(7)<<TS2<<"\t"<<setw(7)<<p2;
                        cout <<setw(W_B)<<data.b<<"\t"<<setprecision(3)<<setw(7)<<TS2<<"\t"<<setw(7)<<p2;}//Get probability using inverse function for f distribution
           }
     file<<endl; cout<<endl;   
        
     file <<setiosflags(ios::left)<<setprecision(3)<<"\nTest statistic:\t"<<TS<<"\tp =\t"<<setprecision(precis)<<p<<setprecision(3)<<endl;
     cout <<setiosflags(ios::left)<<setprecision(4)<<"\nTest statistic:\t"<<TS<<"\tp =\t"<<setprecision(precis+1)<<p<<setprecision(4)<<endl;
     common=slope;  
     return bool(p>data.p);}

void SMAstats::slope_postHoc(profile& data, ofstream &file)
     {
     char resp, letter='1';
     struct testdata {int df; double X; double slope; double p;};
     testdata testmatrix[groups.size()][groups.size()];
     int precis=(int)ceil(log10((double) data.iterations));  //decimal precision for p-values
     Grouplist grouplist;
     Grouplist::iterator Mstart, MSig, Mwork;
     Group Pair;
     Group::iterator It1, It2;
     int i,j, count=0;
     double SLOPE;
     PH_data temp;
     
     cout<<"\nRun post-hoc multiple slope comparisons between groups? (y or n)"<<endl;
     cin>> resp; resp=tolower(resp);
     switch(resp)
                 {
                 case'y':
                         file << "\nPOST-HOC MULTIPLE COMPARISON OF SLOPES AMONG GROUPS"<<endl;
                         cout << "\nRunning permutation tests for "<<setprecision(0)<<(int)((pow((double) groups.size(),2)-groups.size())/2)<<setprecision(3)<<" comparisons\nCount - ";
                         i=0;
                         for(It1=groups.begin(); It1!=groups.end(); It1++)
                              {j=0;
                              for(It2=groups.begin(); It2!=It1; It2++)
                                  {Pair.insert(make_pair(It1->first, It1->second));
                                   Pair.insert(make_pair(It2->first, It2->second));
                                  if(!fit_common(Pair,data.iterations,SLOPE))
                                       {cout<<"slope didn't converge for groups: "<<It1->first<<" & "<<It2->first<<", check data."<<endl; 
                                       #ifdef _WIN32
                                       system("pause");
                                       #endif
                                       testmatrix[j][i].slope=testmatrix[i][j].slope=0.0;}
                                  else {testmatrix[j][i].slope=testmatrix[i][j].slope=SLOPE;}
                                  testmatrix[j][i].X=testmatrix[i][j].X=common_TestStat(Pair , SLOPE,0);
                                  testmatrix[j][i].df=testmatrix[i][j].df=1;
                                  testmatrix[j][i].p=testmatrix[i][j].p=common_pval(Pair,data.iterations,testmatrix[i][j].slope, data.p_method);
                                  Pair.clear();
                                  j++; count++;
                                  cout<<count<<" ";}
                              i++;}
                         cout <<endl;
                         //print matix to file
                         //full matrix
                         file<<"Common slope, test stat, p-value\n"<<setw(max_namelength)<<setiosflags(ios::left)<<"Group"<<"\t";
                         for(It1=groups.begin(); It1!=groups.end(); It1++) {Print_String(It1->first, file, 18+precis); file<<"\t";} file<<endl;
                         It1=groups.begin();
                         for(i=0; i<groups.size(); i++)
                                  {Print_String(It1->first, file, max_namelength); file<<"\t"; It1++;
                                  for(j=0; j<groups.size(); j++)
                                        {if(j==i) file<<"("<<setw(16+precis)<<setprecision(1)<<1.0<<")\t";
                                        else file<<setiosflags(ios::fixed)<<setprecision(3)<<"("<<setw(6)<<testmatrix[i][j].slope<<","<<setw(6)<<testmatrix[i][j].X<<","<<setprecision(precis)<<testmatrix[i][j].p<<")\t";}
                                  file<<endl;
                                  }
                         //p-values matrix
                         file<<"\np-values only\n"<<setw(max_namelength)<<"Group"<<"\t";
                         for(It1=groups.begin(); It1!=groups.end(); It1++) {Print_String(It1->first, file, max_namelength); file<<"\t";} file<<endl;
                         It1=groups.begin();
                         for(i=0; i<groups.size(); i++)
                                  {Print_String(It1->first, file,max_namelength); file<<"\t"; It1++;
                                  for(j=0; j<groups.size(); j++)
                                        {if(j==i) file<<setiosflags(ios::fixed)<<setprecision(1)<<setw(max_namelength)<<1.0<<"\t";
                                        else file<<setiosflags(ios::fixed)<<setprecision(precis)<<setw(max_namelength)<<testmatrix[i][j].p<<"\t";}
                                  file<<endl;
                                  }
                         //make significance lists - prints matrix with columns and rows ordered by slope
                        It1=groups.begin(); i=0;
                        while(It1!=groups.end())   //put groups into list sorted by slopes
                              {
                              temp.groupNo=i; temp.group=It1->first; temp.slope=It1->second.slope; temp.n=It1->second.n;
                              Mstart=grouplist.begin();
                              while(Mstart!=grouplist.end()&& It1->second.slope > Mstart->slope)
                                  {Mstart++;}
                              Mstart=grouplist.insert(Mstart, temp);
                              i++; It1++;}
                        //print matrix
                         file <<"\nSorted pair-wise significance matrix for critical p = "<<setprecision(3)<<data.p<<"; x is indicated where p > "<<data.p<<"\nCode\t"<<setw(max_namelength)<<setiosflags(ios::left)<<"Group"<<"\t"<<setw(8)<<"slope"<<"\t"<<"n"<<"\t|";
                         for(Mwork=grouplist.begin(); Mwork!=grouplist.end(); Mwork++)
                                     {file<<letter; letter++; } file<<"|";
                         letter='1';
                         for(Mwork=grouplist.begin(); Mwork!=grouplist.end(); Mwork++)
                             {file<<"\n"<<letter<<"\t"; Print_String(Mwork->group, file, max_namelength);/*<<setw(max_namelength)<<Mwork->group*/ file<<"\t"<<setiosflags(ios::fixed)<<setprecision(3)<<setw(8)<<Mwork->slope<<"\t"<<Mwork->n<<"\t|";
                             letter++;
                             for(Mstart=grouplist.begin(); Mstart!=grouplist.end(); Mstart++)
                                {if(Mstart==Mwork) file<<'0';
                                 else {if(testmatrix[Mstart->groupNo][Mwork->groupNo].p >= data.p) file <<'x';
                                       else file <<' ';}
                                }
                             file<<"|";
                             }
                         file<<endl;
                         break;
                 default:case'n':
                         file<< "\nNo post hoc comparisons between groups"<< endl;
                         cout<< "\nNo post hoc comparisons between groups"<< endl;
                         break;
                 }
     return;}

void SMAstats::ANCOVA(profile& data, ofstream &file)
     {
     char resp,resp2;
     string filename;
     ofstream FILE;
     Group::iterator GrpI;
     double grandX=0.0, grandY=0.0, grandR=0.0, grandF=0.0, alpha,p, MSd;
     double n=0, theta,y;
     
     if(groups.size()==1) return;
          
     cout<<"\nRun comparisons for lines with common slope? (y or n)"<<endl;
     cin>> resp; resp=tolower(resp);
     switch(resp)
           {
           default:case'n':
                  file<< "\nComparisons of lines with common slope not run"<< endl;
                  cout<< "\nComparisons of lines with common slope not run"<< endl;
                  break;
           case'y':
                  file << "\nCOMPARISON OF LINES WITH COMMON SLOPE"<<endl;
                 //calculate residual and fitted scores using common slope
                 for(GrpI=groups.begin(); GrpI!=groups.end(); GrpI++)
                        {for(int i=0;i < GrpI->second.n; i++)
                              {GrpI->second.Resid[i] = residual_axis(GrpI->second.rawY[i], GrpI->second.rawX[i], common);
                               GrpI->second.Fitted[i] = fitted_axis(GrpI->second.rawY[i], GrpI->second.rawX[i], common); }
                        
                        GrpI->second.Fmean=sum(GrpI->second.Fitted)/GrpI->second.n;  //average distance along slope
                        GrpI->second.Rmean=sum(GrpI->second.Resid)/GrpI->second.n;   // average elevation of slope
                        GrpI->second.s2_F = var_fitted(GrpI->second.s2_Y, GrpI->second.s2_X, GrpI->second.s_XY, common, GrpI->second.n, intercept_flag);          // var.variance(GrpI->second.Fitted, GrpI->second.Fitted);
                        GrpI->second.s2_R = var_resid(GrpI->second.s2_Y, GrpI->second.s2_X, GrpI->second.s_XY, common, GrpI->second.n, intercept_flag);          // var.variance(GrpI->second.Resid, GrpI->second.Resid);
                        GrpI->second.s_RF = covar_rf(GrpI->second.s2_Y, GrpI->second.s2_X, GrpI->second.s_XY, common, GrpI->second.n, intercept_flag);          // var.variance(GrpI->second.Resid, GrpI->second.Resid);           
                        GrpI->second.s2_B = var_slope(GrpI->second.s2_Y, GrpI->second.s2_X, GrpI->second.s_XY,  GrpI->second.s2_F, GrpI->second.s2_R, common, GrpI->second.n); //slope variance with common slope
                        }
                //calucalte grandmean X,Y, F, R
                 for(GrpI=groups.begin(); GrpI!=groups.end(); GrpI++)
                        {n+=GrpI->second.n;
                         grandX+=GrpI->second.Xmean*GrpI->second.n; grandY+=GrpI->second.Ymean*GrpI->second.n; grandR+=sum(GrpI->second.Resid); grandF+=sum(GrpI->second.Resid);}
                 grandX=grandX/n; grandY=grandY/n; grandR=grandR/n; grandF=grandF/n;                
                 //print summary to file
                 file<<"Common slope:\t"<<setprecision(Dec_B)<<common<<"\t"<<"Grand mean X:"<<"\t"<<setprecision(3)<<grandX<<"\t"<<"Grand mean Y:"<<"\t"<<grandY<<"\t"<<"Grand mean F:"<<"\t"<<setprecision(3)<<grandF<<"\t"<<"Grand mean R:"<<"\t"<<grandR<<endl;
                 file<<setiosflags(ios::left)<<setw(max_namelength)<<"\nGroup"<<"\tn\tR2\t"<<setw(7)<<"Interc"<<"\t"<<setw(7)<<"YgrandX"<<"\t"<<setw(7)<<"XgrandY"<<"\t"<<setw(7)<<"Ymean"<<"\t"<<setw(7)<<"Xmean"<<"\t"<<setw(7)<<"Fmean"<<"\t"<<setw(7)<<"Rmean"<<"\t"<<endl;               
               //output x at grand Y, Yat grand X, distance along common etc
                for(GrpI=groups.begin(); GrpI!=groups.end(); GrpI++)
                       {alpha=GrpI->second.Ymean-common*GrpI->second.Xmean;
                        Print_String(GrpI->first, file, max_namelength); file<<"\t"<<GrpI->second.n<<"\t"<<setiosflags(ios::fixed)<<setprecision(3)<<GrpI->second.R2<<"\t"<<setw(7)<<alpha<<"\t"<<setw(7)<<alpha+common*grandX<<"\t"<<setw(7)<<(grandY-alpha)/common<<"\t"<<setw(7)<<GrpI->second.Ymean<<"\t"<<setw(7)<<GrpI->second.Xmean<<"\t"<<setw(7)<<GrpI->second.Fmean<<"\t"<<setw(7)<<GrpI->second.Rmean<<endl;}
                 n=0;
                //run comparisons on Fitted and Residual scores tests for elevation and shift along slope
                 resp2='d';
                 while(resp2!='x')
                      {switch(resp2)
                             {case's':
                                     cout<<"\n Enter name of output filename-> ";
                                     cin >>filename;
                                     FILE.open(filename.c_str()); if(!FILE) {cerr<<"File "<<filename<<" not opened"<<endl; 
                                     #ifdef _WIN32
                                     system("pause"); 
                                     #endif
                                     exit(EXIT_FAILURE);}
                                     FILE<< "Group\tY\tX\tR\tF"<<endl;
                                     for(GrpI=groups.begin(); GrpI!=groups.end(); GrpI++)
                                          for(int i=GrpI->second.n-1; i>=0; i--)
                                                  FILE<<GrpI->first<<"\t"<<GrpI->second.rawY[i]<<"\t"<<GrpI->second.rawX[i]<<"\t"<<GrpI->second.Resid[i]<<"\t"<<GrpI->second.Fitted[i]<<endl;
                                     FILE.close();
                                     file<< "\nData transformed";
                                     switch(stat_choice){
                                             case 0: file << "as R=Y-Bx, F= X for B= "<<common;  break; //OLS
                                             case 1: file << "as R=Y-Bx, F= Y+BX for B= "<<common; break;  //SMA
                                             case 2: file << "as R=Y-Bx, F= BY+X for B= "<<common; break;}  //MA
                                     file<<"\nTransformed data saved to file "<<filename<<endl;
                                     resp2='d';break;
                             case'f': file<< "\nTesting for shifts along the common slope";
                                      cout<< "\nTesting for shifts along the common slope";
                                     #if !allow_ftest
                                         p=walds(groups,2,MSd,1,file, data.p);
                                     #else
                                     if(data.WALD) p=walds(groups,2,MSd,1,file, data.p);
                                     else p=ANOVA(groups,2, MSd,1, file, data.p);
                                     #endif
                                     if(p<data.p)
                                           {file <<"Significant shift along common slope!"<<endl;
                                            cout <<"Significant shift along common slope!\n\nRun post-hoc multiple comparisons (y or n)?"<<endl;
                                            cin >>resp2;
                                            if(resp2=='y')
                                                {file<<"Running post-hoc multiple comparisons"<<endl;
                                                ANCOVA_post(groups, data.WALD, 2, MSd, data.p, file);}
                                            else file<<"No post-hoc multiple comparisons"<<endl;
                                            }
                                     else  {file <<"No shift along common slope"<<endl; cout <<"No shift along common slope"<<endl;
                                           #ifdef _WIN32
                                           system("pause");
                                           #endif
                                           }
                                     resp2='d';break;
                             case'r':
                                     if(!intercept_flag) {cout<<"Not possible for lines fitted through the origin";}
                                     else   {file<< "\nTesting for shifts in elevation between groups";
                                            cout<< "\nTesting for shifts in elevation between groups";
                                             #if !allow_ftest
                                                 p=walds(groups,3,MSd,1,file, data.p);
                                             #else
                                            if(data.WALD) p=walds(groups,3,MSd,1,file, data.p);
                                            else p=ANOVA(groups,3, MSd,1, file, data.p);
                                            #endif
                                            if(p<data.p)
                                                        {file <<"Significant elevation shift between groups!"<<endl;
                                                        cout <<"Significant elevation shift between groups!\n\nRun post-hoc multiple comparisons (y or n)?"<<endl;
                                                        cin >>resp2;
                                                        if(resp2=='y')
                                                           {file<<"\nRunning post-hoc multiple comparisons"<<endl;
                                                           ANCOVA_post(groups, data.WALD, 3, MSd, data.p, file);}
                                                        else  file<<"No post-hoc multiple comparisons"<<endl;}
                                            else  {file <<"No elevation shift between groups"<<endl; cout <<"No elevation shift between groups"<<endl;
                                                   #ifdef _WIN32
                                                   system("pause");
                                                   #endif
                                                  }
                                            }
                                     resp2='d';break;
                             case'x':break;
                             default:case'd':
                                     #ifdef _WIN32
                                     system("cls");
                                     #endif 
                                     cout<<"ANCOVA OPTIONS\n\n"<<endl;
                                     cout<<"S - Save transformed data to disk"<<endl;
                                     if(intercept_flag) cout<<"R - compare Residuals for shifts in elevation"<<endl;
                                     cout<<"F - compare Fitted values for shifts along the common slope"<<endl;
                                     cout<<"X - eXit"<<endl;
                                     resp2=tolower(cin.get());
                                     break;
                             }
                      }
           break;}
    return;
    }


void SMAstats::fit_lines(Group& datagroup, bool resid_flag)
     {
     for (Group::iterator GrpI=datagroup.begin(); GrpI!=datagroup.end(); GrpI++)
           {//calculate useful group summary data
           GrpI->second.Xmean=sum(GrpI->second.rawX)/GrpI->second.n;
           GrpI->second.Ymean=sum(GrpI->second.rawY)/GrpI->second.n;
           //caluclate X,Y variance
           GrpI->second.s2_X =variance(GrpI->second.rawX, GrpI->second.rawX, intercept_flag)-ME.second;
           GrpI->second.s2_Y =variance(GrpI->second.rawY, GrpI->second.rawY, intercept_flag)-ME.first;
           GrpI->second.s_XY =variance(GrpI->second.rawX, GrpI->second.rawY, intercept_flag);         
           //calculate R2 
           GrpI->second.R2   =pow(GrpI->second.s_XY,2)/(GrpI->second.s2_Y*GrpI->second.s2_X);
           //calculate slope, 
           GrpI->second.slope  = B_est(GrpI->second.s2_Y, GrpI->second.s2_X, GrpI->second.s_XY);
           if(resid_flag)
                   {//intercept
                   if(intercept_flag)   GrpI->second.intercept = GrpI->second.Ymean-GrpI->second.slope*GrpI->second.Xmean;
                   else                 GrpI->second.intercept = 0;
                   //calculate fitted and resiudal scores
                   for(int i=0;i < GrpI->second.n; i++)
                           {GrpI->second.Resid[i] = residual_axis(GrpI->second.rawY[i], GrpI->second.rawX[i], GrpI->second.slope);
                           GrpI->second.Fitted[i] = fitted_axis(GrpI->second.rawY[i], GrpI->second.rawX[i], GrpI->second.slope); }}
           //Variance of Fitted and residual scores 
           GrpI->second.s2_F = var_fitted(GrpI->second.s2_Y, GrpI->second.s2_X, GrpI->second.s_XY, GrpI->second.slope, GrpI->second.n, intercept_flag);          // var.variance(GrpI->second.Fitted, GrpI->second.Fitted);
           GrpI->second.s2_R = var_resid(GrpI->second.s2_Y, GrpI->second.s2_X, GrpI->second.s_XY, GrpI->second.slope, GrpI->second.n, intercept_flag);          // var.variance(GrpI->second.Resid, GrpI->second.Resid);
           GrpI->second.s_RF = covar_rf(GrpI->second.s2_Y, GrpI->second.s2_X, GrpI->second.s_XY, GrpI->second.slope, GrpI->second.n, intercept_flag);          // var.variance(GrpI->second.Resid, GrpI->second.Resid);
           GrpI->second.s2_B = var_slope(GrpI->second.s2_Y, GrpI->second.s2_X, GrpI->second.s_XY,  GrpI->second.s2_F, GrpI->second.s2_R, GrpI->second.slope, GrpI->second.n);
           }
     return; }     


void SMAstats::fit_CIs(Group& datagroup, double CI)
     {
     Dec_A=Dec_B=0;      //used in formatting decimal points
     double alpha = 1-CI/100.0;  
     for (Group::iterator GrpI=datagroup.begin(); GrpI!=datagroup.end(); GrpI++)
           {
           //caluculate the p-val  - does an ANOVA on predicted Y values (treatment), compared observed Y (control)
           double SSreg=0.0, SStot=0.0;
           for(int i=0;i < GrpI->second.n; i++)
                    {SStot+= pow(GrpI->second.rawY[i]- GrpI->second.Ymean,2);
                     SSreg += pow( (GrpI->second.Ymean-GrpI->second.s_XY/GrpI->second.s2_X*GrpI->second.Xmean + GrpI->second.s_XY/GrpI->second.s2_X *GrpI->second.rawX[i]) - GrpI->second.Ymean,2);}
           //Get probability using inverse function for f distribution
           #if using_GNU
                        GrpI->second.p=gsl_cdf_fdist_Q( SSreg /((SStot-SSreg)/(GrpI->second.n-2)),1,GrpI->second.n-2);//using GSL library
           #else        
                        GrpI->second.p=pf( SSreg /((SStot-SSreg)/(GrpI->second.n-2)),1,GrpI->second.n-2,0,0);//using R mathlibrary
           #endif
           //calculate slope CIs, Intercept CIs etc
           double B, df;
           if(intercept_flag) df =GrpI->second.n-2;    else  df =GrpI->second.n-1;
           switch(stat_choice)    //slope CI
               {case 0: //OLS
                 B = sqrt(GrpI->second.s2_B);
                 #if using_GNU  //Using GSL library
                 B *= gsl_cdf_tdist_Pinv(1-alpha/2,df);  
                 #else //Using R math library       
                 B *= qt(1-alpha/2,df,1,0);  
                 #endif
                 GrpI->second.slopeCI.first = GrpI->second.slope - B;
                 GrpI->second.slopeCI.second = GrpI->second.slope + B;
                 break;
               case 1:    //SMA
                 B= (1-GrpI->second.R2)/df;
                 #if using_GNU //Using GSL library
                 B*= gsl_cdf_fdist_Pinv(1-alpha,1,df); 
                 #else        //Using R math library
                 B*= qf(1-alpha,1,df,1,0); 
                 #endif
                 GrpI->second.slopeCI.first = min(GrpI->second.slope*(sqrt(B+1)-sqrt(B)),GrpI->second.slope*(sqrt(B+1)+sqrt(B))); 
                 GrpI->second.slopeCI.second = max(GrpI->second.slope*(sqrt(B+1)-sqrt(B)),GrpI->second.slope*(sqrt(B+1)+sqrt(B))); 
                 break;
               case 2:    //MA
                 B= 1/df *(GrpI->second.s2_X * GrpI->second.s2_Y - pow(GrpI->second.s_XY,2));
                 #if using_GNU  //Using GNU GSL library
                 B *= gsl_cdf_fdist_Pinv(1-alpha,1,df); 
                 #else    //Using R math library
                 B*= qf(1-alpha,1,df,1,0); 
                 #endif
                 double brackets = (GrpI->second.s2_Y-GrpI->second.s2_X) +sqrt(pow(GrpI->second.s2_Y-GrpI->second.s2_X,2)+4*pow(GrpI->second.s_XY,2) - 4*B);
                 GrpI->second.slopeCI.first = 1/(2*(GrpI->second.s_XY +sqrt(B)))*brackets;
                 GrpI->second.slopeCI.second = 1/(2*(GrpI->second.s_XY -sqrt(B)))*brackets;
                 break;
               }
           if(intercept_flag)   {//intercept sample varaince
                                GrpI->second.s2_A = GrpI->second.s2_R/GrpI->second.n + pow(GrpI->second.Xmean,2)*GrpI->second.s2_B;
                                //intercept CI
                                B= sqrt(GrpI->second.s2_A);
                                #if using_GNU //Using GSL library
                                B*= gsl_cdf_tdist_Pinv(1-alpha/2,GrpI->second.n-2); 
                                #else        //Using R library
                                B *= qt(1-alpha/2,GrpI->second.n-2,1,0);
                                #endif
                                GrpI->second.interceptCI.first = GrpI->second.intercept - B;
                                GrpI->second.interceptCI.second = GrpI->second.intercept + B;}
           else                 {GrpI->second.s2_A = 0; GrpI->second.interceptCI.first = 0; GrpI->second.interceptCI.second = 0;}                  
          //check for formatting - if slope is very shallow need more decimals
          Dec_A = max(Dec_A,(int)(4.0-log10(fabs(GrpI->second.intercept)))); 
          Dec_B = max(Dec_B,(int)(4.0-log10(fabs(GrpI->second.slope)))); 
          }
      return; }


bool SMAstats::fit_common(Group& new_groups,int iterations, double &SLOPE)
     {
     double slope=0.0, slope2=0.0, error= 1.0;
     int i, count=0, flag=1;
     double s2y_P=0.0, s2x_P=0.0, sxy_P=0.0, n=0.0, s2r, s2f, weight;
     Group::iterator Gr;         
     //inital slope estimate -use pooled sum of squares to estimate slope without group struture
     for(Gr =new_groups.begin(); Gr != new_groups.end(); Gr++)
        {s2y_P+=Gr->second.s2_Y*(Gr->second.n-1); s2x_P+=Gr->second.s2_X*(Gr->second.n-1); sxy_P+=Gr->second.s_XY*(Gr->second.n-1);}                 
        slope=B_est(s2y_P, s2x_P, sxy_P);
     //iteratively recalculate slope
     while(flag==1)
        {count++;
         s2y_P=s2x_P=sxy_P=0.0;
         //recalculate varaiance & weighting        
        for(Gr=new_groups.begin(); Gr!=new_groups.end(); Gr++)
             {s2r= var_resid(Gr->second.s2_Y, Gr->second.s2_X, Gr->second.s_XY , slope, 0, 0);
              s2f= var_fitted(Gr->second.s2_Y, Gr->second.s2_X, Gr->second.s_XY , slope, 0, 0);
              if(stat_choice==2) weight=Gr->second.n*(1/s2r - 1/s2f);
              else               weight=Gr->second.n*(1/s2r + 1/s2f);
              //pooled variances
              s2y_P+=(weight*Gr->second.s2_Y);  s2x_P+=(weight*Gr->second.s2_X); sxy_P+=(weight*Gr->second.s_XY);}
        slope2= B_est(s2y_P, s2x_P, sxy_P); //returns new slope estimate    
        error = fabs(slope - slope2);
        if(error < 0.0001) //slopes have converged
             {flag=0; SLOPE=slope2;}
        else if(count == 200) //hasn't converged
            {cout<<"SLOPES didn't converge"<<slope<<" "<<slope2<<endl; 
             #ifdef _WIN32
            system("pause");
            #endif
            flag=0;}
         slope=slope2;     
         }
     return bool (count < iterations); //returns true if slope stablised before exhausting all starting points
}


//Calculates CI for common slope using algoritm from Warton et al 2006, Appendix D
pair<double, double> SMAstats::commonCI(Group& new_groups, double CI, double slope, int interations){
      double lowCI = 0,uppCI=0, error = 0.01;
      Group::iterator GrpI; 
      //Initail estimate for common slope
      double sxx=0,sxy=0,syy=0, B; int n=0;
      for(GrpI=groups.begin(); GrpI != groups.end(); GrpI++)
        {sxx+=GrpI->second.s2_X; syy+=GrpI->second.s2_Y; sxy+=GrpI->second.s_XY; n+=GrpI->second.n;}    
      B= (1-(pow(sxy,2)/(sxx*syy)))/(n-groups.size()-1); 
      #if using_GNU //Using GSL library
      B *= gsl_cdf_fdist_Qinv(1-CI/100.0,1,n-groups.size()-1); 
      #else        //Using R math library
      B *= qf(1-CI/100.0,1,n-groups.size()-1,0,0); 
      #endif
      uppCI= slope + (max(sqrt(B+1)-sqrt(B),sqrt(B+1)+sqrt(B))-1)*signDF(sxy)*sqrt(syy/sxx);
      lowCI= slope - (1-min(sqrt(B+1)-sqrt(B),sqrt(B+1)+sqrt(B)))*signDF(sxy)*sqrt(syy/sxx);
    
     //refine using optimisation
     double TS_inner, TS_mid, TS_outer, B_inner, B_outer, B_mid;
     double TS_crit;
     #if using_GNU //Using GSL library
     TS_crit = gsl_cdf_chisq_Qinv(1-CI/100.0,1); 
     #else        //Using R math library
     TS_crit = qchisq(1- CI/100, 1,0,0);
     #endif
     bool flag, low, mid, top;  //flags to indicate which needs recalculating  
     //solve for upper
     flag=low=mid=top=1; n=0;
     B_inner = uppCI * (1-(uppCI/slope-1)/2); B_outer = uppCI * (1 +(uppCI/slope-1)*2); //estimate bounding range based on inital estiamtes
     while(flag) //Finds root via bisection method 
        {if(low) {TS_inner = common_TestStat(new_groups, slope, B_inner); low=0;}
        if(top) {TS_outer = common_TestStat(new_groups, slope, B_outer); top=0;}
        if(mid) {B_mid =  B_inner+ (B_outer-B_inner)/2;
                TS_mid = common_TestStat(new_groups, slope, B_mid); mid=0;}
        if(fabs(TS_mid-TS_crit)< error)   {flag=0; uppCI=B_mid;}     //Convergence             
        else if(TS_inner < TS_crit && TS_mid > TS_crit)      {B_outer=B_mid; TS_outer = TS_mid; mid=1;}         //Root lies between low & middle
        else if(TS_mid < TS_crit && TS_outer > TS_crit)     {B_inner=B_mid; TS_inner = TS_mid; mid=1;}            //Root lies between middle & top
        else if(TS_inner < TS_crit && TS_outer < TS_crit) {B_inner=B_mid; TS_inner=TS_mid;     B_outer*=2; top=mid=1;}  //Root is > outer  -->try larger interval
        else if(TS_inner > TS_crit && TS_outer > TS_crit)  {B_outer=B_mid; TS_outer=TS_mid; B_inner/=2; low=mid=1;}  //Root is < inner --> try smaller interval
        if(n>100) flag=0;    
        n++;}
    //solve for lower CI
     flag=top=mid=low=1; n=0;
     B_inner = slope ; B_outer = lowCI + (lowCI-slope)/2;     //estimate bounding range based on inital estiamtes
     while(flag){//Finds root via bisection method 
        if(low) {TS_inner = common_TestStat(new_groups, slope, B_inner); low=0;}
        if(top) {TS_outer = common_TestStat(new_groups, slope, B_outer); top=0;}
        if(mid) {B_mid =  B_inner+ (B_outer-B_inner)/2;
                TS_mid = common_TestStat(new_groups, slope, B_mid); mid=0;}
        if(fabs(TS_mid-TS_crit)< error)   {flag=0; lowCI=B_mid;}     //Convergence             
        else if(TS_inner < TS_crit && TS_mid > TS_crit)      {B_outer=B_mid; TS_outer = TS_mid; mid=1;}         //Root lies between low & middle
        else if(TS_mid < TS_crit && TS_outer > TS_crit)     {B_inner=B_mid; TS_inner = TS_mid; mid=1;}            //Root lies between middle & top
        else if(TS_inner < TS_crit && TS_outer < TS_crit) {B_inner=B_mid; TS_inner=TS_mid;     B_outer/=2; top=mid=1;}  //Root is > top  -->try larger interval
        else if(TS_inner > TS_crit && TS_outer > TS_crit)  {B_outer=B_mid; TS_outer=TS_mid; B_inner*=2; low=mid=1;}  //Root is > top  -->try larger interval
        if(n>100) flag=0;    
        n++;}
    return(make_pair(lowCI, uppCI));}

double SMAstats::common_TestStat(Group& new_groups, double Slope, double H0_b)
     {
     double teststat=0.0, r2rf, df;
     
     Group::iterator GrpI;
     if(H0_b==0)    //Likelihood ratio statistic for SMA & MA , Warton et al Addpenix  D eqn 42
         {for(GrpI=new_groups.begin(); GrpI!=new_groups.end(); GrpI++)
                 {r2rf= pow(covar_rf(GrpI->second.s2_Y, GrpI->second.s2_X, GrpI->second.s_XY , Slope, 0, 0),2)
                   /var_resid(GrpI->second.s2_Y, GrpI->second.s2_X, GrpI->second.s_XY , Slope, 0, 0)
                   /var_fitted(GrpI->second.s2_Y, GrpI->second.s2_X, GrpI->second.s_XY , Slope, 0, 0);      
                   if(intercept_flag) df =GrpI->second.n-2;    else  df =GrpI->second.n-1;
                   teststat+=((df-0.5)*log( 1 - r2rf));}
          return -teststat;}          
     else    //Likelihood ratio statistic for comparing group slope to H0_b, Warton et al Addpenix  D eqn 48
         {for(GrpI=new_groups.begin(); GrpI!=new_groups.end(); GrpI++)
                 {
                 if(intercept_flag) df =GrpI->second.n-2;    else  df =GrpI->second.n-1;
 
                 teststat+=(df-0.5)*log(
                             (var_fitted(GrpI->second.s2_Y, GrpI->second.s2_X, GrpI->second.s_XY , H0_b, 0, 0)
                             *var_resid(GrpI->second.s2_Y, GrpI->second.s2_X, GrpI->second.s_XY , H0_b, 0, 0)
                             /pow(kb(H0_b),2))
                             /(
                              var_fitted(GrpI->second.s2_Y, GrpI->second.s2_X, GrpI->second.s_XY , Slope, 0, 0)
                              *var_resid(GrpI->second.s2_Y, GrpI->second.s2_X, GrpI->second.s_XY , Slope, 0, 0)
                              /pow(kb(Slope),2)));}
         return teststat;}
     }

double SMAstats::common_pval(Group& in_group, int iterations, double common_slope, int resample)
     {
     vector<double>  Res;
     double p=1.0, slope,r;
     int i,j,k;
     Group new_groups = in_group; //copy groups structure
     Group::iterator GrpI2;
     double teststat = common_TestStat(in_group, common_slope, 0);
     if(!resample)// chi-squared distn
         {
         #if using_GNU      //Using GSL library
         return  gsl_cdf_chisq_Q(teststat, new_groups.size()-1);                 
         #else        //Using R math library
         return  pchisq(teststat, new_groups.size()-1,0,0);                 
         #endif
         }
     else {//use resampling
         //Calculate position along line using common slope and residuals
          //cout <<"\nRunning "<<iterations<<" permutations to calculate p-value........."<<endl;
          for(GrpI2=new_groups.begin(); GrpI2!= new_groups.end(); GrpI2++)
               {
               if(intercept_flag) GrpI2->second.intercept = GrpI2->second.Ymean - common_slope*GrpI2->second.Xmean;
               else GrpI2->second.intercept = 0;
               for(i=0;i < GrpI2->second.n; i++)//calculate fitted and resiudal scores
                   { GrpI2->second.Fitted[i] = fitted_axis(GrpI2->second.rawY[i], GrpI2->second.rawX[i], common_slope); 
                     GrpI2->second.Resid[i] = residual_axis(GrpI2->second.rawY[i], GrpI2->second.rawX[i], common_slope);                        
                    Res.push_back(GrpI2->second.Resid[i] -  GrpI2->second.intercept);} //note for resampling, residuals are resampled rather than residual scores  --> intercept subtracted from residual score 
               }
         //resmapling procedure
         for(int j=0; j<iterations-1; j++)
                  {random_shuffle(Res.begin(),Res.end()); // shufles residuals      
                  k=0;
                  for(GrpI2=new_groups.begin(); GrpI2 != new_groups.end(); GrpI2++)
                       {
                       for(i=0; i<GrpI2->second.n; i++)
                            {//Calculate new values for x & y
                            switch(stat_choice){
                                 case 0: GrpI2->second.rawY[i] = Res[k] + GrpI2->second.Fitted[i] * common_slope  + GrpI2->second.intercept;
                                      GrpI2->second.rawX[i] = GrpI2->second.rawX[i]; break;
                                 case 1:  GrpI2->second.rawY[i] = 0.5*(GrpI2->second.Fitted[i]+Res[k]+GrpI2->second.intercept);
                                      GrpI2->second.rawX[i] = 0.5/common_slope*(GrpI2->second.Fitted[i] - Res[k] -GrpI2->second.intercept); break;
                                 case 2: GrpI2->second.rawY[i] = 1/(1+pow(common_slope,2))*(common_slope * GrpI2->second.Fitted[i]+Res[k]-GrpI2->second.intercept);
                                      GrpI2->second.rawX[i] = 1/(1+pow(common_slope,2))*(GrpI2->second.Fitted[i] - common_slope *( Res[k] -GrpI2->second.intercept)); break;
                                      }          
                            k++;}
                       }
                  //fit slopes to data
                  fit_lines(new_groups, 0);
                  if(fit_common(new_groups,iterations,slope)) //only counts data if slopes converge, otherwise skips to next permutation
                         {if(common_TestStat(new_groups, slope, 0) > teststat) p+=1.0;}
                  else {j--; }
                  }
         p=p/iterations;
         Res.clear(); 
         return p;
         }
     }


double SMAstats::walds(Group& new_groups, int column, double& MSd, bool print, ofstream &file, double p_crit)
     {
     if(print) {file<< " using WALD statistic"<<endl; cout<< " using WALD statistic"<<endl;}
     Group::iterator GrpI;    
     int i,j,k,g= new_groups.size();
     double b, teststat, p, acom=0, temp=0, s2B_com =0, s2Ri;
     gsl_vector *A = gsl_vector_calloc(g);
     gsl_vector *s2R = gsl_vector_calloc(g);
     gsl_vector *X = gsl_vector_calloc(g);
     gsl_matrix *s2A = gsl_matrix_alloc (g, g);
     gsl_matrix *L = gsl_matrix_calloc (g-1, g);   
     gsl_matrix *I = gsl_matrix_alloc (g-1, g-1); gsl_matrix_set_identity (I);      
    //calculate vectors A, s2R, X 
     i=0; 
     for (GrpI=new_groups.begin(); GrpI!=new_groups.end(); GrpI++)
        {s2B_com+=1/GrpI->second.s2_B;
        switch(column)
           {case 2:  //Compare Fitted vals for testing shift along slope
                //gsl_vector_set (A, i,  sum(GrpI->second.Fitted)/GrpI->second.n); 
                gsl_vector_set (A, i,  fitted_axis(GrpI->second.Ymean, GrpI->second.Xmean, common));                     
                s2Ri = variance(GrpI->second.Fitted, GrpI->second.Fitted, 1);
                switch(stat_choice)
                     {   case 0:  s2Ri -= ME.second; break;   //OLS
                         case 1:  s2Ri -= ME.first + common*common*ME.second; break;//SMA
                         case 2:  s2Ri -= ME.first*common*common + ME.second; break;}//MA 
                         
                gsl_vector_set (s2R, i, s2Ri *(GrpI->second.n-1)/(GrpI->second.n-2)/GrpI->second.n);    
                if(stat_choice == 1) gsl_vector_set (X, i, GrpI->second.Xmean);        //SMA 
                else if(stat_choice == 2) gsl_vector_set (X, i, GrpI->second.Ymean);        //MA
                break;
           case 3:  //Compare Resdiuals for testing shift in elevation
                gsl_vector_set (A, i,  GrpI->second.Ymean - common*GrpI->second.Xmean);     
                s2Ri = variance(GrpI->second.Resid, GrpI->second.Resid, 1) - ME.first- common*common*ME.second;
                gsl_vector_set (s2R, i, s2Ri*(GrpI->second.n-1)/(GrpI->second.n-2)/GrpI->second.n);    
                gsl_vector_set (X, i, GrpI->second.Xmean);break;}
       i++;}
     
      
     s2B_com=1/s2B_com;     
    //covariance matrix  - equation 51   
    gsl_matrix_set_identity(s2A);
    for(i=0; i<g; i++) for(j=0; j<g; j++) gsl_matrix_set (s2A, i, j,  gsl_matrix_get(s2A,i,j)*gsl_vector_get(s2R,i)+gsl_vector_get(X,i)*gsl_vector_get(X,j)*s2B_com);        
     //Make L matrix 
    for(i=0; i<g-1; i++) gsl_matrix_set (L, i, 0,1); 
    for(i=0; i<g-1; i++) for(j=0; j<g-1; j++) gsl_matrix_set (L, i, j+1,0-gsl_matrix_get(I,i,j)); 
    //L.A
    gsl_vector *LA = gsl_vector_alloc(g-1);
    gsl_blas_dgemv(CblasNoTrans, 1.0, L, A,0, LA);  //LA =L.A  
    //s2A.L'
    gsl_matrix *s2A_Ltr = gsl_matrix_alloc(g,g-1);    
    gsl_blas_dgemm(CblasNoTrans,CblasTrans, 1.0, s2A, L,0, s2A_Ltr);
    //L.s2A.L'
    gsl_matrix *L_s2A_Ltr = gsl_matrix_alloc(g-1,g-1);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1.0, L, s2A_Ltr,0,L_s2A_Ltr);
    //(L.s2A.L')-1
    gsl_matrix *Inv= gsl_matrix_alloc(g-1,g-1);   
    int s;   gsl_permutation *perm = gsl_permutation_alloc (g-1);
    gsl_linalg_LU_decomp (L_s2A_Ltr, perm, &s);
    gsl_linalg_LU_invert (L_s2A_Ltr, perm, Inv);
    //Right product  = Inv(L.s2A.L').LA
    gsl_vector *Prod = gsl_vector_alloc(g-1);
    gsl_blas_dgemv(CblasNoTrans, 1, Inv, LA, 0, Prod );
   //left product gives teststat
    gsl_blas_ddot(LA, Prod, &teststat);
    //pvalue
     #if using_GNU //Using GSL library
     p = gsl_cdf_chisq_Q(teststat, g-1); 
     #else        //Using R math library
     p= pchisq(teststat, g-1,0,0);
     #endif

    MSd =teststat;
    //caluclate common intercept - App D Eq 54 
    //inverse s2A
    if(p>p_crit) //calulcate commmon estaimate
        {gsl_matrix_free (Inv);  Inv= gsl_matrix_alloc(g,g);   
        gsl_permutation_free(perm); perm = gsl_permutation_alloc (g);
        gsl_linalg_LU_decomp (s2A, perm, &s); gsl_linalg_LU_invert (s2A, perm, Inv);
        //Vector 1(1,g)
        gsl_vector *one = gsl_vector_alloc(g); gsl_vector_set_all (one, 1);
        //Eqn54 top - product  = Inv(s2A).A
        gsl_vector_free (Prod); Prod = gsl_vector_alloc(g);
        gsl_blas_dgemv(CblasNoTrans, 1, Inv, A, 0, Prod);
        for(i=0; i< g; i++) acom+=gsl_vector_get(Prod, i);
        //Eqn 54 bototm 
        gsl_blas_dgemv(CblasNoTrans, 1, Inv, one, 0, Prod);
        for(i=0; i< g; i++) temp+=gsl_vector_get(Prod, i);
        acom=acom/temp;}
   //print results
   if(print) {
             file <<"df\tstat\tp\tA_com"<<endl;          cout <<"df\tstat\tp\tA_com"<<endl;
             file<<setiosflags(ios::left)<<setiosflags(ios::fixed)<<setprecision(0)<<g-1<<"\t"<<setprecision(3)<<setw(7)<<teststat<<"\t"<<setw(7)<<p<<"\t";
             cout<<setiosflags(ios::left)<<setiosflags(ios::fixed)<<setprecision(0)<<g-1<<"\t"<<setprecision(4)<<setw(7)<<teststat<<"\t"<<setw(7)<<p<<"\t";
             if(p>p_crit) {file <<acom; cout<<acom;} 
             else {file <<"na"; cout<<"na";} 
             cout<<endl; file<<endl;}         
    //free memeory
    gsl_vector_free(A); gsl_vector_free(s2R); gsl_vector_free(X); gsl_vector_free(LA); gsl_vector_free(Prod);
    gsl_matrix_free(Inv);   gsl_matrix_free (s2A);   gsl_matrix_free(L);   gsl_matrix_free(I); gsl_matrix_free(s2A_Ltr); gsl_matrix_free(L_s2A_Ltr); 
    gsl_permutation_free(perm);
    return p;}
     
//Run post-hoc pairwise comparisons on among groups. Outputs results in matrix,
//and in sorted list in acending order of group means
void SMAstats::ANCOVA_post(Group& new_groups, int WALD, int column, double MSd, double p_crit, ofstream &file)
      {
      int i,j,n, precis=4, letter='1', max_namelength=7;  
      char resp;
      double Si2, Sj2;
      
      i=n=0;
      PH_data temp;
      // data structure storing pairwise group comparion reults
      struct testdata {double diff; double df; double SE; double p;};testdata testmatrix[new_groups.size()][new_groups.size()];
      Grouplist grouplist;
      Grouplist::iterator Mstart, MSig, Mwork;
      Group::iterator It1, It2;
      Group::iterator GrpI;
      Group Pair; 
      //which method to use. If Fstat, run levenes test for equality of variances     
      resp = 'w';
      #if allow_ftest
      if(WALD!=1) resp = levenes(new_groups, column, p_crit, file);
      #endif
      //count total n    
      for(GrpI=new_groups.begin(); GrpI!=new_groups.end(); GrpI++) n+=GrpI->second.n;
         
      cout<<"processing........"<<endl;
      for(It1=new_groups.begin(); It1!=new_groups.end(); It1++)
          {j=0;
          if(It1->first.size() > max_namelength) max_namelength= It1->first.size();
          for(It2=new_groups.begin(); It2!=It1; It2++)
                {testmatrix[j][i].diff=testmatrix[i][j].diff=fabs(            //mean difference between two groups
                                      XYFR(column, It1->second.Xmean, It1->second.Ymean, It1->second.Fmean, It1->second.Rmean)
                                     -XYFR(column, It2->second.Xmean, It2->second.Ymean, It2->second.Fmean, It2->second.Rmean));
                switch(resp)
                          {
                          case 'w':  //wald statistic
                                         {Pair.insert(make_pair(It1-> first, It1-> second)); Pair.insert(make_pair(It2-> first, It2-> second));
                                         testmatrix[j][i].df=testmatrix[i][j].df=1;
                                         testmatrix[j][i].p=testmatrix[i][j].p=walds(Pair, column, testmatrix[i][j].SE, 0, file, p_crit);
                                         testmatrix[j][i].SE=testmatrix[i][j].SE;
                                         Pair.clear(); break;}
                          #if allow_ftest //Can only use ftest with Rmath library           
                          case 'g': //unequal variances, use Games-Howell Method for post-hoc tests
                                         Si2=Sj2=0.0;
                                         if(column==0) {Si2= SumSq(It1->second.rawX, It1->second.rawX, 1)/(It1->second.n-1); Sj2= SumSq(It2->second.rawX, It2->second.rawX, 1)/(It2->second.n-1);} //use x data
                                         else if (column==1){Si2= SumSq(It1->second.rawY, It2->second.rawY, 1)/(It1->second.n-1); Sj2= SumSq(It2->second.rawY, It2->second.rawY, 1)/(It2->second.n-1);} //use y data
                                         else if(column==2) {Si2= SumSq(It1->second.Fitted, It1->second.Fitted, 1)/(It1->second.n-1); Sj2= SumSq(It2->second.Fitted, It2->second.Fitted, 1)/(It2->second.n-1);} //use fitted data
                                         else if (column==3){Si2= SumSq(It1->second.Resid, It2->second.Resid, 1)/(It1->second.n-1); Sj2= SumSq(It2->second.Resid, It2->second.Resid, 1)/(It2->second.n-1);} //use Residual data
                                         testmatrix[j][i].SE=testmatrix[i][j].SE=sqrt(Si2 /It1->second.n + Sj2/It2->second.n);
                                         testmatrix[j][i].df=testmatrix[i][j].df=floor(   pow(Si2/It1->second.n+Sj2/It2->second.n,2)  /  (   pow(Si2,2)/(pow(It1->second.n,2)*(It1->second.n-1)) + pow(Sj2,2)/(pow(It2->second.n,2)*(It2->second.n-1))));
                                         testmatrix[j][i].p=testmatrix[i][j].p=ptukey(testmatrix[j][i].diff/(testmatrix[j][i].SE/sqrt(2.0)),1,(double)new_groups.size(),testmatrix[j][i].df,0,0);
                                         break;
                          case 't':  //equal variance - use Tukey-Kramer method for post-hoc tests
                                         testmatrix[j][i].SE=testmatrix[i][j].SE=sqrt(MSd*(1.0/(double)It1->second.n+1.0/(double)It2->second.n));
                                         testmatrix[j][i].df=testmatrix[i][j].df=floor(n-new_groups.size());
                                         testmatrix[j][i].p=testmatrix[i][j].p=ptukey(testmatrix[j][i].diff/(testmatrix[j][i].SE/sqrt(2.0)),1,(double)new_groups.size(),testmatrix[j][i].df,0,0);
                                         break;
                          #endif
                          }
                j++;}
          i++;}
     //print matix to file
     //full matrix
     if(resp=='w')   file<<"Difference, standard error, df, p-value\n"<<setw(max_namelength)<<setiosflags(ios::left)<<"Group"<<"\t";
     else file<<"Difference, standard error, df, p-value\n"<<setw(max_namelength)<<setiosflags(ios::left)<<"Group"<<"\t";
     for(It1=new_groups.begin(); It1!=new_groups.end(); It1++) {Print_String(It1->first, file, 21+precis); file<<"\t";} file<<endl;
     It1=new_groups.begin();
     for(i=0; i<new_groups.size(); i++)
              {Print_String(It1->first, file, max_namelength); file<<"\t"; It1++;
              for(j=0; j<new_groups.size(); j++)
                      {if(j==i) file<<"("<<setw(19+precis)<<setprecision(1)<<1.0<<")\t";
                       else file<<setiosflags(ios::fixed)<<setprecision(3)<<"("<<setw(6)<<testmatrix[i][j].diff<<","<<setw(6)<<testmatrix[i][j].SE<<","<<setprecision(0)<<setw(3)<<testmatrix[j][i].df<<","<<setprecision(3)<<testmatrix[i][j].p<<")\t";}
              file<<endl;}
     //p-values matrix
     file<<"\np-values only\n"<<setw(max_namelength)<<"Group"<<"\t";
     for(It1=new_groups.begin(); It1!=new_groups.end(); It1++) {Print_String(It1->first, file, max_namelength); file<<"\t";} file<<endl;
     It1=new_groups.begin();
     for(i=0; i<new_groups.size(); i++)
             {Print_String(It1->first, file,max_namelength); file<<"\t"; It1++;
             for(j=0; j<new_groups.size(); j++)
                     {if(j==i) file<<setiosflags(ios::fixed)<<setprecision(1)<<setw(max_namelength)<<1.0<<"\t";
                     else file<<setiosflags(ios::fixed)<<setprecision(precis)<<setw(max_namelength)<<testmatrix[i][j].p<<"\t";}
             file<<endl;}
     //make significance lists - prints matrix with columns and rows ordered by slope
    i=0;
    for(It1=new_groups.begin(); It1!=new_groups.end(); It1++)   //put new_groups into list sorted by means
          {temp.groupNo=i; temp.group=It1->first;
          temp.n=It1->second.n; 
          Mstart=grouplist.begin();           
          temp.slope=XYFR(column, It1->second.Xmean, It1->second.Ymean, It1->second.Fmean, It1->second.Rmean);
          while(Mstart!=grouplist.end()&& XYFR(column, It1->second.Xmean, It1->second.Ymean, It1->second.Fmean, It1->second.Rmean) > Mstart->slope) 
                {Mstart++;}
          Mstart=grouplist.insert(Mstart, temp);  i++;}
     //print matrix  - indicates not-signifacntly different groups
     file <<"\nSorted pair-wise significance matrix for critical p = "<<setprecision(3)<<p_crit<<"; x is indicated where p > "<<p_crit<<"\nCode\t"<<setw(max_namelength)<<setiosflags(ios::left)<<"Group"<<"\t"<<setw(8)<<"mean"<<"\t"<<"n"<<"\t|";
     for(Mwork=grouplist.begin(); Mwork!=grouplist.end(); Mwork++)
            {file<<letter; letter++; } file<<"|";
     letter='1';
     for(Mwork=grouplist.begin(); Mwork!=grouplist.end(); Mwork++)
            {file<<"\n"<<letter<<"\t"; Print_String(Mwork->group, file, max_namelength);  file<<"\t"<<setiosflags(ios::fixed)<<setprecision(3)<<setw(8)<<Mwork->slope<<"\t"<<Mwork->n<<"\t|";
            letter++;
            for(Mstart=grouplist.begin(); Mstart!=grouplist.end(); Mstart++)
                   if(Mstart==Mwork) file<<'0';
                   else {if(testmatrix[Mstart->groupNo][Mwork->groupNo].p >= p_crit) file <<'x';   else file <<' ';}
            file<<"|";}
     file<<endl;
     return;}

//calculate slope
double SMAstats::B_est(double s2y, double s2x, double sxy)
      {switch(stat_choice)
           {   case 0:  return (sxy/s2x);//OLS
               case 1:  return (signDF(sxy)*sqrt(s2y/s2x));//SMA
               case 2:  return (1/(2*sxy)*(s2y - s2x + sqrt( pow(s2y - s2x,2) + 4*pow(sxy,2))));   //MA 
               }}

double SMAstats::residual_axis(double Y, double X, double B){return Y-B*X;}
 
double SMAstats::fitted_axis(double Y, double X, double B)
     {switch(stat_choice) 
           {case 0: return X; //OLS
            case 1: return Y+B*X; //SMA
            case 2: return B*Y+X; }//MA              
     }   
     
//caluclates variance of residual values - formuals page Appednix D of Warton et al 2006
double SMAstats::var_resid(double s2y, double s2x, double sxy, double b, int n, bool from_data)
      {double temp = s2y - 2*b*sxy + b*b*s2x;   //OLS, SMA OR MA
       if(from_data) temp*= (n-1)/(n-2);
       return temp;}
            
//caluclates variance of FITTED values - formuals page Appednix D of Warton et al 2006
double SMAstats::var_fitted(double s2y, double s2x, double sxy, double b, int n, bool from_data)
      {double temp;
      switch(stat_choice)
           { case 0:  return (s2x);//OLS
             case 1:  temp = s2y + 2*b*sxy + b*b*s2x;   //SMA
                      if(from_data) temp*= (n-1)/(n-2);
                      return temp;
             case 2:  temp = b*b*s2y + 2*b*sxy + s2x;   //MA
                      if(from_data) temp*= (n-1)/(n-2);
                       return temp;}
      }   

//caluclates variance of slope  - formuals page Table 4 of Warton et al 2006
double SMAstats::var_slope(double s2y, double s2x, double sxy, double s2f, double s2r, double b, int N)
      {switch(stat_choice)
           { case 0: return s2y/s2x * (1.0-sxy*sxy/s2y/s2x)/(N-2.0); //OLS / SMA 
             case 1: return s2y/s2x * (1.0-sxy*sxy/s2y/s2x)/(N-2.0); //OLS / SMA 
             case 2: return pow(1.0 + b*b,2)/(N-2.0) /(s2f/s2r + s2r/s2f - 2.0);} //MA           
      }            
                 
//caluclates covarince of residual & fitted values - formuals page Appednix D of Warton et al 2006
double SMAstats::covar_rf(double s2y, double s2x, double sxy, double b, int n, bool from_data)
      {double temp;
      switch(stat_choice)
           {   case 0:  temp = sxy - b*s2x;   //OLS
                        if(from_data) temp*= (n-1)/(n-2);
                        return temp;
               case 1:  temp = s2y - b*b*s2x;   //SMA
                        if(from_data) temp*= (n-1)/(n-2);
                        return temp;
               case 2:  temp = sxy +b*(s2y-s2x) - sxy*b*b;   //MA
                        if(from_data) temp*= (n-1)/(n-2);
                        return temp;}
      }   

//Equation 46 from Warton et al
double SMAstats::kb(double b){
     if(stat_choice ==0) return 1;
     else if(stat_choice ==1) return 2*fabs(b);
     else return 1+b*b;
     }

//returns X, Y, Fitted or Resdiual score according to flag 
double SMAstats::XYFR(int flag, double X, double Y, double F, double R)
      {switch(flag)  {case 0: return X;   case 1: return Y;  case 2: return F;  case 3: return R;}}


double SMAstats::variance(vector<double> dataX, vector<double> dataY, bool int_flag)
      {int n =dataX.size();
      if(int_flag) n=n-1;
      return (SumSq(dataX, dataY, int_flag)/n);}

double SMAstats::SumSq(vector<double> dataX, vector<double> dataY, bool int_flag)
      {double xy=0, X, Y;
      if(int_flag) {X=sum(dataX)/dataX.size(); Y=sum(dataY)/dataX.size();}//standard 
      else         {X=0; Y=0;}//means set to zero - used when forcing through origin
      for(int i=0; i<dataX.size(); i++)  xy+=(dataX[i]-X)*(dataY[i]-Y);
      return (xy); }
                            
double SMAstats::sum(vector<double> data)
      {double X=0;  for(int i=0; i<data.size(); i++)   X+=data[i];
      return X;}    


#if allow_ftest
//Runs one-sample ANOVA on data stored in "column" (0=X, 1=Y, 2 = fitted, 3 = residual)
double SMAstats::ANOVA(Group& new_groups, int column, double& MSd, bool print, ofstream &file, double p_crit)
      {
      double x, SST, SSTr, SSE, x_dotdot, x_idot, MSE, MSTr,p; int n;
       SST=SSTr=SSE=x_dotdot=x_idot=MSE=MSTr=n=0;      
      
      if(print) {file<< " using ANOVA"<<endl; cout<< " using ANOVA"<<endl;}                            
      //CALCULATE Group SUMS OF SQUARES
      for (Group::iterator GrpI=new_groups.begin(); GrpI!=new_groups.end(); GrpI++)
          {n+=GrpI->second.n;
          x_idot=0;
          for(int i=0;i < GrpI->second.n; i++)
               {x= XYFR(column, GrpI->second.rawX[i], GrpI->second.rawY[i], GrpI->second.Fitted[i], GrpI->second.Resid[i]);
                x_dotdot+=x;
                SST+=pow(x,2);
                x_idot+=x;}
          SSTr+=pow(x_idot,2)/GrpI->second.n;
          }
      SSTr-=pow(x_dotdot,2)/n;
      SST-=pow(x_dotdot,2)/n;
      SSE=SST-SSTr;
      //degrees of freedom
      int df_Tr = new_groups.size()-1;
      int df_Er = n- new_groups.size();
          if(column>1) df_Er= df_Er-1;        //1 less degree of freedom for tests on fitted and residual axes
      
      MSTr=SSTr/df_Tr;   MSE=SSE/df_Er;
      //Get probability using inverse function for f distribution
      p=pf(MSTr/MSE, df_Tr, df_Er,0,0);
      //print to file
      if(print) {
            file <<"\nSource\tdf\tSumSq\tMeanSq\tf\tp"<<endl;
            file<<"Group\t"<<setprecision(1)<<df_Tr <<setprecision(3)<<"\t"<<SSTr<<"\t"<<MSTr<<"\t"<<MSTr/MSE<<"\t"<<p<<endl;
            file<<"Error\t"<<setprecision(1)<<df_Er<<setprecision(3)<<"\t"<<SSE<<"\t"<<MSE<<endl;
            file<<"Total\t"<<setprecision(1)<<n-1<<"\t"<<setprecision(3)<<SST<<endl;
            cout <<"\nSource\tdf\tSumSq\tMeanSq\tf\tp"<<endl;
            cout <<"Group\t"<<setprecision(1)<<df_Tr<<"\t"<<setprecision(4)<<SSTr<<"\t"<<MSTr<<"\t"<<MSTr/MSE<<"\t"<<p<<endl;
            cout <<"Error\t"<<setprecision(1)<<df_Er<<"\t"<<setprecision(4)<<SSE<<"\t"<<MSE<<endl;
            cout <<"Total\t"<<setprecision(1)<<n-1<<"\t"<<setprecision(4)<<SST<<endl;}
      MSd= MSE;
      return p;
      }

  
//LEVENES TEST FOR EQUALITY OF VARIANCES
char SMAstats::levenes(Group& new_groups, int column, double p_crit, ofstream &file)
     {double Z_dotdot,Z_idot,W,W2,p, diff;
     int i, n;  char resp;
     Group::iterator GrpI;
     W=W2=Z_dotdot=Z_idot=n=0.0;
     
     for(GrpI=new_groups.begin(); GrpI!=new_groups.end(); GrpI++) //FIND GRANDMEAN of ZIJ
          {n+=GrpI->second.n;
          for(i=0;i < GrpI->second.n; i++)
               Z_dotdot+=fabs(XYFR(column, GrpI->second.rawX[i], GrpI->second.rawY[i], GrpI->second.Fitted[i], GrpI->second.Resid[i])
                              - XYFR(column, GrpI->second.Xmean, GrpI->second.Ymean, GrpI->second.Fmean, GrpI->second.Rmean));
          }
      Z_dotdot/=n;
      for(GrpI=new_groups.begin(); GrpI!=new_groups.end(); GrpI++) //caluclate levenes statistic W
          {Z_idot=0.0;
          for(i=0;i < GrpI->second.n; i++)//find group mean
               {diff = fabs(XYFR(column, GrpI->second.rawX[i], GrpI->second.rawY[i], GrpI->second.Fitted[i], GrpI->second.Resid[i])
                              - XYFR(column, GrpI->second.Xmean, GrpI->second.Ymean, GrpI->second.Fmean, GrpI->second.Rmean));
               Z_idot+= diff;} 
          Z_idot=Z_idot/GrpI->second.n;
          W+=(GrpI->second.n*pow((Z_idot-Z_dotdot),2)); //numerator
          for(i=0;i<GrpI->second.n; i++)//calculate denomitor
               {diff = fabs(XYFR(column, GrpI->second.rawX[i], GrpI->second.rawY[i], GrpI->second.Fitted[i], GrpI->second.Resid[i])
                              - XYFR(column, GrpI->second.Xmean, GrpI->second.Ymean, GrpI->second.Fmean, GrpI->second.Rmean));
               W2+=pow(diff-Z_idot,2); }
          }
      //degrees of freedom
      int df_Tr = new_groups.size()-1;
      int df_Er = n- new_groups.size();
          if(column>1) df_Er= df_Er-1;        //1 less degree of freedom for tests on fitted and residual axes    
      W=W*(df_Er)/((df_Tr)*W2);
      //Get probability using inverse function for f distribution from R mathlibrary
      p=pf(W,df_Tr,df_Er,0,0);
      file <<"Levene's test for equal variance across groups"<<endl;
      file<<"W\tdf1\tdf2\tp\n"<<setprecision(3)<<W<<"\t"<<setprecision(0)<<df_Tr<<"\t"<<df_Er<<"\t"<<setprecision(3)<<p<<endl;
      cout <<"\nLevene's test for equal variance across groups"<<endl;
      cout<<"W\tdf1\tdf2\tp\n"<<setprecision(3)<<W<<"\t"<<setprecision(0)<<df_Tr<<"\t"<<df_Er<<"\t"<<setprecision(4)<<p<<endl;

      //run post hoc tests - method depends on whether variance was signifcantly heterogeneous or not
      if(p<p_crit)
         {file<<"Variance significantly heterogeneous"<<endl;
          cout<<"Variance significantly heterogeneous. Proceed with: \n\tGames-Howell method (recommended) - y \n\tor Tukey-Kramer method            - n"<<endl;
          cin >> resp;
          while(resp!='y' && resp!='n') {cout <<"\nInvalid choice, try again"<<endl; cin >> resp;}}
     else {cout<<"Homogeneity of variance across groups"<<endl; file<<"Homogeneity of variance across groups"<<endl; resp='n';}
     switch(resp)
         {case('y'): file<< "\nGames-Howell method"<<endl; break;
          case('n'):file<< "\nTukey-Kramer method"<<endl; break;}
     if(resp=='y') return 'g'; //use games howell
     else return 't';  //use tuker-Kramer
     }
#endif
