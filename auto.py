# Shiva Mudide

# Making sure to import the reduction script
import xml.etree.ElementTree as etree
import csv
import dgsreductionmantid

mylist = []
with open('Gedata.csv', 'r') as csv_file:
    # Use csv.reader to read the csv line by line
    # Convert each line into its own list
    # Append each line to a master list
    csv_reader = csv.reader(csv_file)
    for line in csv_reader:
        mylist.append(list(line))
mylist = mylist[0:]

# Import the template xml file as a tree
mytree = etree.parse('Ge_ramp_II.xml')
# Set the root
root = mytree.getroot()

allruns = []

# Iterate through the lines and extract only those which are ramping runs
for i in range(len(mylist)):
    run = mylist[i]
    # Select for ramping runs
    if run[1] != 'Ramping':
        continue
    # A trial is a list with elements: run number, lowTemp, highTemp
    trial = []
    number = run[0]
    # Split about the -> arrow and append to a master list
    lowhigh = run[2].split('-')
    low = lowhigh[0]
    high = lowhigh[1][1:]
    trial.append(number)
    trial.append(low)
    trial.append(high)
    allruns.append(trial)

# print(allruns)
files = []

# Iterate through the ramping runs to find the corresponding background run
# Update the xml with the metal and corresponding empty run
# Write to a new xml and store the file name 
for i in range(len(allruns)):
    trial1 = allruns[i]
    trial1low = trial1[1]
    trial1high = trial1[2]
    # For each ramping run, search through the rest of the master list to find the
    # corresponding background
    for j in range(i + 1, len(allruns)):
        # If found a match
        if(allruns[j][1] == trial1low and allruns[j][2] == trial1high):
            count = 0
            # Need to set both scan attributes
            # Setting the scan attribute for the metal
            for scan in root.iter('scan'):
                if count == 0:
                    scan.set('logvaluemin', trial1low)
                    scan.set('logvaluemax', trial1high)
                    scan.set('runs', trial1[0])
                    # print(trial1[0], trial1low, trial1high)

                else:
                    # Setting the scan attribute for the background
                    scan.set('logvaluemin', allruns[j][1])
                    scan.set('logvaluemax', allruns[j][2])
                    scan.set('runs', allruns[j][0])
                    # print(allruns[j])
                # Counter ensures we don't double set
                count = count + 1
            # Write to a new xml
            mytree.write('output' + str(i) + '.xml')
            print('WROTE: ' + 'output' + str(i) + '.xml' + ' !!!')
            # Store the file name
            files.append('output' + str(i) + '.xml')
            # We must break to ensure we only create xml files for the desired runs
            break
    

# Run the dgsreductionmantid script for each of the ramping runs 
# by passing each of the newly created xml files as an argument
# The second counter variable can be uncommented to allow for a set number of reductions

# count2 = 0
for file in files:
    dgsreductionmantid.dgsreduction(XMLfile=file)
    print("FINISHED RUNNING: " + file)
    # count2 = count2 + 1
    # if count == 2:
    #     break
    

