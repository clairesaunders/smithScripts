# Background
The auto.py script parses a melting csv file and creates proper xml files for each run. It then runs the dgsreductionmantid script for every xml file. The original excel file should be converted into a csv file. The csv file here is Gedata.csv. A template xml file is given. The template here is Ge_ramp_II.xml. 

auto.py finds every ramping instance in the csv file. It procures the temperature range that the ramping was done over. It then searches for the corresponding background (empty) sample that was ramped up between the same temperature limits. It uses the template xml to update the scan attributes to the correct values for both the metal and empty based on the previous search. It writes these updates to a new xml file. It does this for all ramping instances. 

It imports the dgsreductionmantid script and runs the script for each of the newly created xml files. 

# How to run

Download all the contents of this repository into one directory and run auto.py. 
If you do not wish to run the dgsreductionmantid script, as this takes a long time, and simply want to view the output xml files, comment out the last for-loop in auto.py.
There may be some further dependecies needed for the dgsreductionmantid script to work. 

# Robustness

The script has been optimized to be as robust as I could make it. There is very little hardcoding. The csv file can be replaced with any other melting csv file. The template xml can be replaced. If future requirements need to be met, like including more than just the ramping runs, these requirements can be easily met with changes to < 5 lines of code. Because auto.py is robust, it can be easily modified to meet further requirements. 

