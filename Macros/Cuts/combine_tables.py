#!/usr/bin/env python2
import sys


#Setup_______________ first line (5)
periods =["W07", "W08", "W09", "W10", "W11", "W12", "W13", "W14", "W15", "WAll"]

mass_type ="HMDY" #"HMDY", "LowM_AMDY"
production ="t3" #"t3"=t3, "slot1"=t5
which_cuts ="FinalCuts" #"PhastCuts", "FinalCuts"




#Setup_______________ last line (15)

#Output settings
if len(sys.argv) < 2:
    print "\nScript combines tables made from tableOfCuts.C"
    print "\nCurrent periods to combine:  "
    print "    " + str(periods)
    print "Mass type considered:     " + mass_type
    print "Production considered:    " + production
    print "Which cuts considered:    " + which_cuts
    print "\n\nUsage:"
    print "\tpython combine_tables.py 1"
    print "\n\nView cuts file using jupyter notebook with script:"
    print "      ViewCutTable.ipynb"
    sys.exit()

#Basic Checks
if periods[0] != "W07":
    print "First period should be \"W07\""
    sys.exit()

#Determine cut_names and period_cuts from files
cut_names =[]
per_cuts = {}

for per in periods:
    per_cuts[per] = []
    
    #Open period file
    table_per ="Data/TableOfCuts/tableOfCuts_" + per +"_"+ mass_type
    table_per +="_"+ production +"_"+ which_cuts +".csv"
    
    try:
        f_per = open(table_per)
    except IOError:
        print "File for period  " + per + "  did not open"
        sys.exit()

    #Loop through period file
    line_num=0
    for line in f_per:
        line = line.rstrip("\n")
        separate = line.split()
        separate[0] = separate[0].rstrip(",")
        separate[0] = separate[0].replace(",", "_")
        
        if per == "W07":
            cut_names.append(separate[0])
        elif separate[0] != cut_names[line_num]:
            print "Cut names do not match"
            print "\t"+separate[0] +"\n\t"+ cut_names[line_num]
            sys.exit()

        per_cuts[per].append(separate[-1])
        line_num +=1

#Output file setup
output_name = "Data/CombineTables/combineTables_"
output_name += mass_type +"_"+ production +"_"+ which_cuts +".csv"
f_output = open(output_name, 'w')

#Write periods as header
for per in periods:
    f_output.write(",          ")
    f_output.write(per)
    
f_output.write("\n")

#Write cut names followed by cuts per period    
for icut in range(len(per_cuts["W07"])):
    f_output.write(cut_names[icut])
    
    for per in periods:
        f_output.write(",   ")
        f_output.write(per_cuts[per][icut])

    f_output.write("\n")
    
f_output.close()
#
#final settings
print "\n\nSettings________"
print "\nTables combined:  "
print "\nCurrent periods to combine:  "
print "    " + str(periods)
print "Mass type considered:     " + mass_type
print "Production considered:    " + production
print "Which cuts considered:    " + which_cuts
print "\nFile:"
print "\t"+output_name
print "was written"
