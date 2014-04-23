/*******************************************************************************
    * SMATR 2.O :  Standardised Makjor Axis Tests and Routines
    * http://www.bio.mq.edu.au/ecology/SMATR/
    * Definition of class SMAstats                                   *
    * Date   : 24/10/06                                 
    * Copyright (C) 2006 Daniel Falster
     *
    * This program is free software; you can redistribute it and/or
    * modify it under the terms of the GNU General Public License
    * as published by the Free Software Foundation, version 2, as 
    * stated at http://www.gnu.org/licenses/gpl.html
     
    * This program is distributed in the hope that it will be useful,
    * but WITHOUT ANY WARRANTY; without even the implied warranty of
    * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    * GNU General Public License for more details.
    * 
    * Contact details: Daniel Falster, dfalster [at] bio.mq.edu.au
    * Dept Biological Sciences, Macquarie University 2109, Australia
*****************************************************************************/
#include "menu.h"
#include <algorithm>
using std::random_shuffle;

//DATA STRUCTURES USED
//store data for inidividual groups, including raw data, fitted line etc.
// rawX/Y, Fitted, Resid = X,Y, Fitted and Residual scores,  n = number datapoints, X/Y/F/Rmean = average of X/Y/F/R, s2_X/Y/F/R/B/A = variance of s2_X/Y/F/R/slopeB/interceptA, covariance of XY, RF, R2 & p of correlation, slope & slopeCI, intercept &interceptCI
struct summary_stats {vector<double> rawX; vector<double> rawY; vector<double> Fitted; vector<double> Resid; 
                     int n; double Xmean; double Ymean;  double Rmean; double Fmean; double s2_X; double s2_Y; double s_XY; double s2_B; double s2_A; double s2_F; double s2_R; double s_RF;
                     double R2; double p; double slope; pair<double,double> slopeCI; double intercept; pair<double,double> interceptCI;};
typedef map<string, summary_stats> Group;      // makes a map of individual groups data, key is group name

class SMAstats {
  public:
     // constructor    
     SMAstats(); 
     /*retrieves data stored in the profile data structure from menu class (vector of
     strings), and divides it into groups. Individual groups are stored in a map with key = groupvariable. */          
     void input(profile& data, ofstream &file);
     //display input columns to screen
     void displayData();
     //fits lines to individual groups and record summary data (r2, p-value, Xmean, Ymean, CI's, etc) 
     void fit_individual(profile& data); 
     //print details of individual groups; run comparions to comparison to hypothesised slope / intercept
     void print_individual_details(profile& data, ofstream &file);
     //fit common slope & test significance
     bool test_common(profile& data, ofstream &file);
     //run post host comparisons of slope on pairwise basis
     void slope_postHoc(profile& data, ofstream &file);
     //Run comparisons for lines with common slope
     void ANCOVA(profile& data, ofstream &file);
     
   private:
//menu choices
     int stat_choice;         //0 = OLS, 1 = SMA, 2= MA 
     bool intercept_flag;     //include intercept in equation
     bool ME_flag;            //include measurement error
     
//VARIABLES
     Group groups;             //all data stored in a map, keyed by group
     double common;           //COMMON SLOP ESTIMATE
     pair<double, double> ME; //Measurement Errors
     
//For formatting
     int Dec_B, Dec_A, W_A, W_B;//for decimal point formatting
     int max_namelength;      //used to set column widths in formatting
    
//ROUTINES APPLYING TO ALL GROUPS
     //fits lines to individual groups and record summary data  - rotines consistent with what is required for resampling. 
     void fit_lines(Group& datagroup, bool resid_flag);  
     //calculates CIs, R2, variance of slope & intercept and other for individual lines
     void fit_CIs(Group& datagroup, double CI); 
     //fit common slope
     bool fit_common(Group& new_groups,int iterations, double &SLOPE); 
     //estimate CI of common slope
     pair<double, double> commonCI(Group& new_groups, double CI, double slope, int iterations); 
     //teststat for common slope test
     double common_TestStat(Group& new_groups, double slope, double H0_b);       
     //test for slope heterogenity
     double common_pval(Group& new_groups, int iterations, double common_slope, int resample);  
     //compare fitted or residual values using walds statistic
     double walds(Group& new_groups, int column, double& MSd, bool print, ofstream &file, double p_crit);  
     //Post goc comparions for lines of equal slope  
     void ANCOVA_post(Group& new_groups, int WALD, int column, double MSd, double p_crit, ofstream &file); 
     //ANOVA and levenes test for ANCOVA comparisons using ftest instead of walds (dsiabled).
     #if allow_ftest 
     double ANOVA(Group& new_groups, int column, double& MSd, bool print, ofstream &file, double p_crit);  //compare groups using ANOVA
     char levenes(Group& new_groups, int column, double p_crit, ofstream &file);                           //test for equality of variance
     #endif
     
//USEFUL SUBROUTINES 
     //return X,Y,F or R depedning on flag
     double XYFR(int flag, double X, double Y, double F, double R);           
     //estimate slope
     double B_est(double s2y, double s2x, double sxy);  
     double kb(double b);                                                
     //calculate fitted score
     double fitted_axis(double Y, double X, double B);                     
     //calculate residual score
     double residual_axis(double Y, double X, double B);                   
     //Calculate covaraince of X and Y (variance if X =Y)
     double variance(vector<double> dataX, vector<double> dataY, bool int_flag);
     //Calculate cross product of X and Y (sum of squares if X =Y)
     double SumSq(vector<double> dataX, vector<double> dataY, bool int_flag);
     //sum all elements in a vector
     double sum(vector<double> data);
     //calculate variance of fitted score
     double var_fitted(double s2y, double s2x, double sxy, double b, int n, bool from_data); 
     //calculate variance of residual score
     double var_resid(double s2y, double s2x, double sxy, double b, int n, bool from_data);  
     //calculate covariance of residual & fitted score
     double covar_rf(double s2y, double s2x, double sxy, double b, int n, bool from_data); 
     //Calculate variance of slope
     double var_slope(double s2y, double s2x, double sxy, double s2f, double s2r, double b, int N);
};
