# hmmn
Codes and data for the paper "A Proposal for Finite But Unbounded Human Lifetimes" coauthored by 
Fei Huang, Ross Maller, Brandon Milholland, Xu Ning

## Code Availability
Code used in analysing the data and constructing the figures and tables are provided. 

“Figure1_Lifetime densities.m” – Figure 1

“acceleration females.R” and “acceleration males.R” – Tables 1, 2, 3.

“combined cohorts CDF.R” – Figure 2, Figure 4

“alive deduction.R” – Table 6 (this script is to generate the numbers of alive population at threshold age for females and males “F_N_pop” and “M_N_pop” for “Table6_expected and CI.m”) 

“Table 6_expected and CI.m” – Table 6

“Figure3_CIplots.m”-- Figure 3

“Figure 5_Distribution of Gamma.m” – Figure 5

“Figure 6_Centenarians analysis_IDL.m” – Figure 6

Tables 4 and 5 are generated using the STLT package published in GitHub and available at https://github.com/u5838836/STLT based on the Australian data obtained from the Human Mortality Database (https://www.mortality.org/). 

## Data Availability 
The data used in this paper are all publicly available. The Netherlands data is obtained from the online data source provided at the Harvard Dataverse (Einmahl, 2018) and the Human Mortality Database (HMD, https://www.mortality.org/). We extracted the 1x1 cohort life tables and exposures-to-risk for males and females of Netherlands separately from the HMD. For the data combination and quality check of the augmented Netherlands data used in this paper, please refer to Section 3 in Huang et al. (2020).  The Australian data is also obtained from HMD, please refer to Section 2 in Fu et al. (2021) for more details of the data.  The England and Wales supercentenarians data is obtained from the International Database on Longevity (https://www.supercentenarians.org/). 

## References:
F. Huang, R. Maller, X. Ning, Modelling life tables with advanced ages: An extreme value theory approach, Insurance: Mathematics and Economics 93, 95–115 (2020).

B. Fu, F. Huang, R. Maller, Modelling Australian Life Tables with Advanced Ages – A Report Prepared for the Australian Government Actuary. Available at SSRN: http://dx.doi.org/10.2139/ssrn.4020309  (2021)

J. Einmahl, Replication Data for: Limits to human life span through extreme value theory, https://doi.org/10.7910/DVN/RNZA5D, Harvard Dataverse, V1, UNF:6:hY4ZMepNlEF7lzYY1QJ98g== [fileUNF] (2018)

