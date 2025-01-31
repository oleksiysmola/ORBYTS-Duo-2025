import numpy as np
import os
import glob
import shutil

cwd = os.getcwd()

## create folder in a given directory
def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print ('Error: Creating directory. ' +  directory)


## find all objects to grep within a file
def find_obj(file, point_flag):
    print("Finding MolPro Objects...")
    with open(file) as f:
        targets = [line for line in f if ("!MRCI trans" in line)or("!MRCI expec" in line)]
    #
    temp = []
    for obj in targets:
        temp.append(obj.split()[0]+" "+obj.split()[2])
    #
    if point_flag == "single":
        print()
        print(f"Found all objects in {file}")
        return np.unique(temp)
    #
    if point_flag == "many": 
        print()
        print(f"Found all objects in {file}")
        return np.unique(temp)

## Loading in files, extracting objects to grep
files = glob.glob("./*.out")                             # find all .out files in wd
for fname in files:
    createFolder(fname[0:-4])                            # create folders where we will save grep output to
    shutil.copy(fname,fname[0:-4]+"/"+fname[2:len(fname)])

errors = []
def grep(file,geom_Flag):
    os.chdir(cwd)
    # seperation vector
    print(file)
    with open(file[2:len(file)],'r') as f:
        temp = [line for line in f if "r=[" in line]
        if len(temp)==0:
            print(f"{file} is probably a single point calculation. Please refer to output file for geometry... setting seperation to -1 for identification")
            spFLAG = True
            seperation = [-1]
        else:
            seperation = temp[0][temp[0].index("[")+1:-2].split(',')

    np.savetxt(file[0:-4]+'/Nuclear_seperation', seperation, fmt='%s')

    ## grep objects
    obj = find_obj(file,geom_Flag)

    for i in range(len(obj)):
        os.chdir(cwd)
        #
        with open(file) as f:
            targets = [line for line in f if (obj[i].split()[0] in line)&(obj[i].split()[1] in line)]
        #
        os.chdir(os.getcwd()+"/"+file[0:-4])
        #
        coupling = obj[i].split()[1]
        #
        keyword =coupling #name of object.txt
        #
        #name of folder based on curve type (DMZ, DMY...)
        foldername = obj[i][obj[i].index("|")+1:-obj[i][::-1].index("|")-1]
        createFolder(obj[i][obj[i].index("|")+1:-obj[i][::-1].index("|")-1])
        #
        # name file
        filename = "./"+foldername+str('/')+keyword+"_"+foldername
        #
        #Save file
        np.savetxt(filename+str(".txt"), targets, fmt='%s')
        #
        if ("DMX" in coupling) or ("DMY" in coupling) or ("DMZ" in coupling):
            #
            unitkeyword =coupling+str("_Debye") 
            #
            unitfoldername = str("./")+foldername+str("/Debye")
            #
            createFolder(unitfoldername)
            #
            unitfilename = unitfoldername + str('/')+unitkeyword+"_"+foldername
            #
            sliced_targets = [w.split()[6] for w in targets]
            
            if len(targets)==len(seperation):
                np.savetxt(unitfilename+str(".txt"), np.stack([np.array(seperation), np.array(sliced_targets)],axis=1), fmt='%s')

            if len(targets)!=len(seperation):

                err_conditions = [((len(targets)<len(seperation)&(seperation[0]!=-1))),((len(targets)>len(seperation))&(seperation[0]!=-1)),(seperation[0]==-1)]
                err_values = ["less", "more","single calculated"]

                print(f"There are some missing points/repeats/single point calculation for {obj[i].split()[1]} in file {file}. Saving record in MolPro.e and moving on...")
                errors.append(f"{obj[i].split()[1]} from file {file} has {np.select(err_conditions,err_values)} points ({len(sliced_targets)}) than geometries ({len(seperation)})!\n")

                np.savetxt(unitfilename+str(".txt"), np.array(sliced_targets), fmt='%s')

        if ("LX" in coupling) or ("LY" in coupling) or ("LZ" in coupling):
            unitkeyword =coupling+str("_i") 
            unitfoldername = str("./")+foldername+str("/i")
            createFolder(unitfoldername)
            # name file
            unitfilename = unitfoldername + str('/')+unitkeyword+"_"+foldername
            #Save file
            sliced_targets = [w.split()[3] for w in targets]
            
            if len(targets)==len(seperation):
                np.savetxt(unitfilename+str(".txt"), np.stack([np.array(seperation), np.array(sliced_targets)],axis=1), fmt='%s')

            if len(targets)!=len(seperation):

                err_conditions = [((len(targets)<len(seperation)&(seperation[0]!=-1))),((len(targets)>len(seperation))&(seperation[0]!=-1)),(seperation[0]==-1)]
                err_values = ["less", "more","single calculated"]

                print(f"There are some missing points/repeats for {obj[i].split()[1]} in file {file}. Saving record in MolPro.e and moving on...")
                errors.append(f"{obj[i].split()[1]} from file {file} has {np.select(err_conditions,err_values)} points ({len(sliced_targets)}) than geometries ({len(seperation)})!\n")
                
                np.savetxt(unitfilename+str(".txt"), np.array(sliced_targets), fmt='%s')


        if ("LSX" in coupling) or ("LSY" in coupling) or ("LSZ" in coupling):
            ## create folder for different hamiltonians
            createFolder("Mean_Field")
            #
            unitkeyword =coupling+str("_cm-1") 
            #
            unitfoldername = str("./")+foldername+str("/Mean_Field/cm-1")
            #
            createFolder(unitfoldername)
            # name file
            unitfilename = unitfoldername + str('/')+unitkeyword+"_MF"+"_"+foldername
            #Save file
            # print(list(a[b]))
            #
            sliced_wavn = list(np.array([w.split()[6] for w in targets])[list(np.arange(0,len(targets),2))])
            #
            #
            ## create folder for different hamiltonians
            createFolder("Breit_Pauli")
            #
            unitfoldername2 = str("./")+foldername+str("/Breit_Pauli/cm-1")
            #
            createFolder(unitfoldername2)
            # name file
            unitfilename2 = unitfoldername2 + str('/')+unitkeyword+"_BP"+"_"+foldername
            #Save file
            sliced_wavn2 = list(np.array([w.split()[6] for w in targets])[list(np.arange(1,len(targets),2))])
            #
            
            #
            if float(len(targets))==float(len(seperation)/2):
                np.savetxt(unitfilename2+".txt", np.stack([np.array(seperation), np.array(sliced_wavn2)],axis=1), fmt='%s')
                np.savetxt(unitfilename+".txt",  np.stack([np.array(seperation), np.array(sliced_wavn)],axis=1), fmt='%s')

            if float(len(targets))!=float(len(seperation)/2):

                err_conditions = [((len(targets)<len(seperation)&(seperation[0]!=-1))),((len(targets)>len(seperation))&(seperation[0]!=-1)),(seperation[0]==-1)]
                err_values = ["less", "more","single calculated"]

                print(f"There are some missing points/repeats for {obj[i].split()[1]} in file {file}. Saving record in MolPro.e and moving on...")
                errors.append(f"{obj[i].split()[1]} from file {file} has {np.select(err_conditions,err_values)} points than geometries!\n")
                np.savetxt(unitfilename2+str(".txt"), np.array(sliced_wavn2), fmt='%s')
                np.savetxt(unitfilename+str(".txt"), np.array(sliced_wavn), fmt='%s')

for i in files:
    grep(i,"many")
    # path = i[:-4]
    # objects = ["DMZ", "DMX", "DMY", "LX", "LY", "LZ", "LSX", "LSY", "LSZ"]
    # for obj in objects:
    #     filenames = []
    #     for root, dirs, files in os.walk(path):
	#         for file in files:
	#     	    if(file.endswith(f"{obj}.out")):
	#     		    filenames.append(os.path.join(root,file))
    #                 print(filenames)
        # with open(f'{obj}_df.txt', 'w') as writer:
        #     readers = [open(filename) for filename in filenames]
        #     for lines in zip(*readers):
        #         print(' '.join([line.strip() for line in lines]), file=writer)

np.savetxt(cwd+"/MolPro.e",errors,'%s')