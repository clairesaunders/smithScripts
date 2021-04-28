# Shiva Mudide

import xml.etree.ElementTree as etree
import csv
import dgsreductionmantid

mylist = []
with open('Gedata.csv', 'r') as csv_file:
    csv_reader = csv.reader(csv_file)
    for line in csv_reader:
        mylist.append(list(line))
mylist = mylist[0:]

mytree = etree.parse('Ge_ramp_II.xml')
root = mytree.getroot()

allruns = []

for i in range(len(mylist)):
    run = mylist[i]
    if run[1] != 'Ramping':
        continue
    trial = []
    number = run[0]
    lowhigh = run[2].split('-')
    low = lowhigh[0]
    high = lowhigh[1][1:]
    trial.append(number)
    trial.append(low)
    trial.append(high)
    allruns.append(trial)

# print(allruns)
files = []

for i in range(len(allruns)):
    trial1 = allruns[i]
    trial1low = trial1[1]
    trial1high = trial1[2]
    for j in range(i + 1, len(allruns)):
        if(allruns[j][1] == trial1low and allruns[j][2] == trial1high):
            count = 0
            for scan in root.iter('scan'):
                if count == 0:
                    scan.set('logvaluemin', trial1low)
                    scan.set('logvaluemax', trial1high)
                    scan.set('runs', trial1[0])
                    # print(trial1[0], trial1low, trial1high)

                else:
                    scan.set('logvaluemin', allruns[j][1])
                    scan.set('logvaluemax', allruns[j][2])
                    scan.set('runs', allruns[j][0])
                    # print(allruns[j])
                count = count + 1
            mytree.write('output' + str(i) + '.xml')
            print('WROTE: ' + 'output' + str(i) + '.xml' + ' !!!')
            files.append('output' + str(i) + '.xml')
            break
    


# count2 = 0
for file in files:
    dgsreductionmantid.dgsreduction(XMLfile=file)
    print("FINISHED RUNNING: " + file)
    # count2 = count2 + 1
    # if count == 2:
    #     break
    

