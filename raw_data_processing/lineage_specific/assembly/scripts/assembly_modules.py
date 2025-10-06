## This file contains all necessary python functions for running the ASSEMBLY snakemake
## Changes to any of these functions should only be done on the dev branch and thoroughly tested prior to merging

## functions
def read_samplesCSV(spls):
    # reads in samples.csv file, format: Batch, Sample,Alignments,ProviderName,Patient
    hdr_check = ['Batch', 'Sample', 'Alignments', 'ProviderName', 'Patient']
    switch = "on"
    file = open(spls, 'r')
    list_path = []
    list_splID = []
    list_providerNames = []
    list_refG = []
    list_patient = []
    for line in file:
        line = line.strip('\n').split(',')
        # Test Header. Note: Even when header wrong code continues (w/ warning), but first line not read.
        if switch == "on":
            if (line == hdr_check):
                print("Passed CSV header check")
            else:
                Warning("CSV did NOT pass header check! Code continues, but first line ignored")
            switch = "off"
            continue
        # build lists
        list_path.append(line[0])
        list_splID.append(line[1])
        list_refG.append(line[2])
        list_providerNames.append(line[3])
        list_patient.append(line[4])
    return [list_path,list_splID,list_refG,list_providerNames,list_patient] # 2019.02.08, Arolyn: removed set(list_patient) provides only unique subject IDs
