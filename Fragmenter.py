"""
Script to generate a set of profmark benchmark files where the inserted homologous domain in target sequences is 
fragmented either by a % coverage of the original or by only including some set quantity of residues in the inserted domain

The script will identify the .pos, .neg, .tbl, .test.fa, .train.msa files given the path to the "basename" of the files. 
It will then perform the "fragmenting" on the .test.fa sequences, and adjust values in the .pos file accordingly. 
.neg, .tbl, .train.msa files will not be edited. 

All new files will take the form [basename][Fragment_Size/Coverage] and will be written to the [Path_to_outdirectory] location. 


General Usage:

python3 [filename].py [-h] [Path_to_Basename] [Path_to_outdirectory] [Fragment or coverage size] [C/S for fragment size in residues or fragment size in coverage]

Example Usage:

python3 [filename].py C:/User/Home/Directory/Profmark_Benchmarks/Pmark C:/User/Home/Directory/Split_Benchmarks 10 C

Each input read as:

[Path_to_Basename] = C:/User/Home/Directory/Profmark_Benchmarks/Pmark
[Path_to_outdirectory] = C:/User/Home/Directory/Split_Benchmarks
[Fragment or coverage size] = 10
[C/S for fragment size in residues or fragment size in coverage] = C

Thus the script would look in C:/User/Home/Directory/Profmark_Benchmarks/ 
for benchmarks with a base of Pmark.neg, Pmark.pos, Pmark.tbl, Pmark.test.fa, Pmark.train.msa. 

Copies of these files would be written to the location C:/User/Home/Directory/Split_Benchmarks as files with names Pmark10.neg, Pmark10.pos, ect...
Additionally, these files would have fragmented inserted sequences with 10 percent of the original size.

"""

#Imports for usage
import sys
import Bio
from Bio import SeqIO
from Bio import AlignIO
from colorama import Fore
import pandas as pd
import numpy as np
import random
import os
import shutil
import os
import copy
import shutil
import argparse
import sys
from argparse import RawTextHelpFormatter
import os

"""
Created 9/21/2024 
@Author: Ezra Taub
"""


"""
Parsing arguments using argparse package:
"""
parser = argparse.ArgumentParser(description='Generate a set of benchmark files (.neg, .pos, .test.fa, .tbl, .train.msa) corresponding to fragmented inserted domains. These fragments can be either percent coverage or literal sizes. Chosen uniformly at random depending upon the Frag_Type parameter.\n\nExample Usage:\n\n python3 [filename] C:/User/Home/Directory/Profmark_Benchmarks/Pmark C:/User/Home/Directory/Split_Benchmarks 10 C\n\nThus the script would look in C:/User/Home/Directory/Profmark_Benchmarks/ for benchmarks with a base of Pmark.neg, Pmark.pos, Pmark.tbl, Pmark.test.fa, Pmark.train.msa. Copies of these files would be written to the location C:/User/Home/Directory/Split_Benchmarks as files with names Pmark10.neg, Pmark10.pos, ect... Additionally, these files would have fragmented inserted sequences with 10 percent of the original size.',formatter_class=RawTextHelpFormatter)
parser.add_argument('Path_to_basename',help="Search for an input set of benchmark files at this directory, with the final part referring to the basename name", type=str)
parser.add_argument('Path_to_outdirectory',help="Search for an output directory to write the new set of benchmark files to",type=str)
parser.add_argument('Frag_Size',help='Input fragment size as a percentage coverage or a number of residues', type=int)
parser.add_argument('Frag_type',help = 'Type of fragment, C = coverage, S = Size in amino acid residues. Default is C, coverage',type=str)
args = parser.parse_args()

#Assigning parsed arguments below. These are used in main() at the very bottom of this script.
Path_to_basename = args.Path_to_basename #Path to the basename
Path_to_outdirectory = args.Path_to_outdirectory #Path to write out to 
Frag_Size = args.Frag_Size #Size of desired fragments
Frag_type = args.Frag_type #Type of fragmenting to perform





#
# Find_Min_Length(path_to_pos):
# Objective: This function is called to return the minimum length of all inserted domains. It helps to check 
# whether this dataset is compatible with the input size of a fragment. If the input size of the fragment 
# exceeds that of the smallest inserted domain, then the script cannot run. 
#
def Find_Min_Length(path_to_pos):
    #find the minimum of the lengths of inserted domains found in the .pos file
    return np.min(pd.read_fwf(path_to_pos, header=None)[11])

#
# CheckFragType(Frag_type, Frag_Size, pos_path)
# Objective: Check if frag_type input is valid (C or S) and return default (C) if input is not S or C
# In addition, checks if the input length is valid given the type. If C is input and size exceeds 100, this is impossible and quits
# Alternatively, if S is the type of fragment and there is a inserted domain whose size is greater than the length, this is impossible and the program quits
#
def CheckFragType_Validity(Frag_type, Frag_Size, pos_path):

    if Frag_type == "S":
        min = Find_Min_Length(pos_path)
        if Frag_Size > min:
            print(f"The fragment size provided, {Frag_Size}, exceeds the length of a domain in the dataset (length {min}). Please specify a size smaller than this. Quitting")
            quit()
        else:
            return Frag_type
        
    else:
        if Frag_Size > 100:
            print(f"The Fragment Coverage provided, {Frag_Size} Exceeds 100 percent. This is not possible. Try again")
            quit()
        else:
            return "C"
#
# Get_Basename(Path_to_basename)
# Objective:
# Function takes in the whole path to the basename provided as input and grabs just the pmark basename string from the end
# Example: if input is C:\User\Directory\Folder\Basename then this just returns "Basename" 
#
def Get_Basename(Path_to_basename):
    #assume path is of the form dir/dir/dir/Basename
    return Path_to_basename.split('/')[-1]

#
# get_random_substring(s, num_chars, type):
#
# Objective: 
# Function to get the random substring from the input inserted domain 
#
# IN: s, num_chars, type
# s = string, inserted domain
# num_chars = int, desired amount of characters OR coverage
# type = string, "C" or "S" corresponding to coverage or size respectively
# OUT: substring, start_index, end_index
# substring = string, the selected substring
# start_index = first index relative to whole input (s) where the substring was "spliced"
# end_index = last index relative to whole input (s) where the subtring was "spliced"
#
# EXAMPLE:
#
# get_random_substring("EEMAAAAAAAAAMEMDEE", 9, "S")
# >> "AAAAAAAAA", 3, 12
#
def get_random_substring(s, num_chars, type):
    string_length = len(s)
    # CHECK TYPE, assign number of characters accordingly
    if type == "S":
        num_chars = max(1, min(len(s), num_chars))
    else:
        num_chars = max(1, min(string_length, int((num_chars / 100) * string_length)))

   
    if num_chars == 0 or string_length == 0:
        return "", 0, 0

    # Randomly select a start index for the substring
    start_index = random.randint(0, string_length - num_chars)
    end_index = start_index + num_chars
   
    # Get the substring
    substring = s[start_index:end_index]
   
    return substring, start_index, end_index


#
# format_row(row)
# Objective: This function is called to take the new values in the .pos file corresponding to fragmented sequences 
# and then write them out with correct spacing.
# IN: pandas dataframe row (series)
# OUT: Single-line string + \n formatted appropriately. 
# Example:
"""
IN:

0          14-3-3/10/40-151
1                       678
2                        39
3                       111
4                       527
5     sp|Q9NTX5|ECHD1_HUMAN
6                       307
7                       240
8                       278
9                         .
10      A2F1A4_TRIV3/12-123
11                      111
12                        1
13                      111
14                        .
15    sp|Q99J09|MEP50_MOUSE
16                      342
17                      196
18                       38
19                        c

OUT:
14-3-3/10/40-90                       617     39     50    527 sp|Q9NTX5|ECHD1_HUMAN                         307    240    278 .  A2F1A4_TRIV3/92-142             50      1     50 .  sp|Q99J09|MEP50_MOUSE        342    196     38 c 

"""
def format_row(row):
    return f"{row.iloc[0]:35s} {row.iloc[1]:5d} {row.iloc[2]:6d} {row.iloc[3]:6d} {row.iloc[4]:6d} {row.iloc[5]:40s} {row.iloc[6]:8d} {row.iloc[7]:6d} {row.iloc[8]:6d} {row.iloc[9]:2s} {row.iloc[10]:25s} {row.iloc[11]:8d} {row.iloc[12]:6d} {row.iloc[13]:6d} {row.iloc[14]:2s} {row.iloc[15]:25s} {row.iloc[16]:6d} {row.iloc[17]:6d} {row.iloc[18]:6d} {row.iloc[19]:2s}\n"

"""
Process_TESTFA

The function actually doing the splitting, the heavy lifting. 
This function operates by 
1) reading in the .pos and .test.fa files 
2) reading each entry of the .test.fa file as a FASTA entry
3) Getting indicies of the inserted domain
4) Splitting the inserted domain
5) Splicing the inserted domain back in 
6) Updating the information about the split locations
7) Updating corresponding information in the .pos file
8) return new text file (test.fa file) and new pos file

Additional note:
This function also tracks whether there are any duplicates detected, and will return the value of such duplicates / whether
they are present
"""
def Process_TESTFA(Path_to_TestFA, Path_to_pos,Frag_type, Frag_Size, path_to_outdirectory, basename):

    #Declare a variable to track the quantity of duplicated sequences
    duplicate_count = 0 
    seen = set() #to track duplicates

    #Declare a dictionary to store existing values 
    

    #This is the out index that helps track what the position is in the .pos file
    index= 0

    #We will write to this repeatedly as our new test.fa file it is a string, which once
    #Completely processed, will be written out to the new .test.fa file
    txt = ''

    #This is our new text file for the .pos file that we will write to in this function
    outpf = open(path_to_outdirectory+"/"+basename+str(Frag_Size)+".pos", 'w')
    #This is our new text file for the .test.fa file that we will write to in this function
    outfa = open(path_to_outdirectory+"/"+basename+str(Frag_Size)+".test.fa", 'w')

    #Will be written out as BasenameFrag_Size.pos or BasenameFrag_Size.test.fa



    try:
        Fasta = SeqIO.parse(Path_to_TestFA, "fasta") #try parsing the .test.fa file in as a FASTA file
    except:
        print("Error in parsing test.fa file") #quit if this goes wrong
        quit()
    
    try:
        pos_tbl = pd.read_fwf(Path_to_pos, header=None) #try parsing the .pos file in as a pandas dataframe
    except:
        print("Error in parsing .pos file") #quit if this goes wrong
        quit()

    #NOW, we do the actual processing and fragmenting

    for record in Fasta: #Looping through all sequences in the test file

        if 'decoy' in record.description:     #If this is a decoy
            


            txt+=(">"+record.description+"\n"+record.seq+"\n")  #Write to the text file the same, unchanged content

            if ">"+record.description+"\n"+record.seq+"\n" in seen:
                duplicate_count+=1
            else:
                seen.append(">"+record.description+"\n"+record.seq+"\n")  

        else:    #and if this is NOT a decoy

            #get split indicies from the fasta header
            indicies = record.description.split(' ')[0].split('/')[-1].split('-')

            #Now that we know where the embedded domain is, get that segment
            domain = record.seq[int(indicies[0]):int(indicies[1])]

            #Call helper function to randomly split homologous segment.
            domain, idxs, idxf = get_random_substring(domain, Frag_Size, Frag_type)

            #Note - idxs refers to start index, and idxf to final index

            #We splice the new fragment into the flanking regions as follows:
            newseq = str(record.seq[0:int(indicies[0])])+str(domain)+str(record.seq[int(indicies[1])::])

            #We need to relabel the header/description now with the new split locations:
            newsplits = indicies[0]+"-"+str(int(indicies[0])+len(domain))
            
            #Get the original description split to recalculate indicies
            newline = record.description.split(' ')

            #Recalculate start index
            idxs+=int(newline[2].split('/')[1].split('-')[0])
            #Recalculate finish index
            idxf+=int(newline[2].split('/')[1].split('-')[0])

            #reformat relevant parts of the description:
            newline[0] = newline[0].split('/')[0]+"/"+newline[0].split('/')[1]+"/"+newsplits

            #Finally, change the description:
            record.description = ">"+newline[0]+" "+newline[1]+" "+newline[2].split('/')[0]+f"/{idxs}-{idxf}"
            

            #Given that we have the record description as a unique identifier, we need to check if this has resulted in a duplicate
            if record.description in seen:
                duplicate_count+=1
            else:
                seen.append(record.description)

            #Now add that new description to the new test.fa file (txt string)
            txt+=str(record.description)+"\n"+str(newseq)+"\n"


            #NOW we have to update the positive table entry (the one with 16 columns)
            #With our relevant new information

            #We update the description (first column)
            pos_tbl.loc[index,0] = newline[0]

            
            pos_tbl.loc[index,1] = int((int(indicies[0])+len(domain))+int(len(record.seq[int(indicies[1]):])))

            #Update length of inserted domain
            pos_tbl.loc[index,3] = len(domain)
            #Update length of inserted domain (redundant but still needs changing)
            pos_tbl.loc[index,11] = len(domain)

            pos_tbl.loc[index,10] = pos_tbl.loc[index,10].split('/')[0] + f"/{idxs}-{idxf}"
            #Update length of inserted domain (redundant but still needs changing)
            pos_tbl.loc[index,13] = len(domain)  

    

            outpf.write(format_row(pos_tbl.loc[index]))

            index+=1



    #now we write out the new information, converting the descriptions and sequences to 
    #a string
    outfa.write(str(txt))
    
    #close both files
    outpf.close()
    outfa.close()
    
    return duplicate_count

"""
Check_files(Path_to_basename, Path_to_splitdir)

Objective:

Function to check all the filenames for the benchmark files, 
these will be returned as strings in main() to be used by other functions. This
will throw an error (prior to processing) if any of these files or directories cannot be located

In: Path_to_basename, Path_to_splitdir

Path_to_basename  = the input path as argument 1
Path_to_splitdir = the desired output directory

Out: Test_FA_path, pos_path, neg_path, tbl_path, train_msa_path

Test_FA_path = path to basename.test.fa
pos_path = path to basename.pos
neg_path = path to basename.neg
tbl_path = path to basename.tbl
train_msa_path = path to basename.train.msa

"""
def Check_files(Path_to_basename, Path_to_splitdir):

    #We append ".test.fa" to the basename as this is where we expect to find the .test.fa file
    Test_FA_path = Path_to_basename+".test.fa"
    #We try to parse the .test.fa file...
    if os.path.exists(Test_FA_path):
        None
    else:
        print(f"Could not find default test.fa file at {Test_FA_path}. Quitting")
        quit()
    #Construct the presumed path to the .pos file
    pos_path = Path_to_basename+".pos"
    #Now we do the same for the .pos, trying to read it in as a fwf file
    
    if os.path.exists(pos_path):
        None
    else:  
        print(f"Could not find default test.fa file at {Test_FA_path}. Quitting")
        quit()
    #Construct the presumed path to the .neg file
    neg_path = Path_to_basename+".neg"
    #Now we do the same for the .neg, trying to read it in as a fwf file
    if os.path.exists(neg_path):
        None
    else:
        print(f"Could not find default neg file at {neg_path}. Quitting")
        quit()
    tbl_path = Path_to_basename+".tbl"
    if os.path.exists(tbl_path):
        None
    else:
        print(f"Could not find default .tbl file at {tbl_path}. Quitting")
        quit()
    train_msa_path = Path_to_basename+".train.msa"
    if os.path.exists(train_msa_path):
        None
    else:
        print(f"Could not find .train.msa file at {train_msa_path}. Quitting")
        quit()
    if os.path.exists(Path_to_splitdir):
        None
    else:
        print(f"Could not find directory to write fragmented copy files to at {Path_to_splitdir}. Quitting")
        quit()


    return Test_FA_path, pos_path, neg_path, tbl_path, train_msa_path


"""
copy_files(path_to_outdirectory, fragment_size, neg_path, tbl_path, train_msa_path, basename)

Objective:
Function to make copies of all non-modified files associated with the profmark benchmark
looks for the files at specified paths, makes a copy written out to the path_to_outdirectory path specified

In: path_to_outdirectory, fragment_size, neg_path, tbl_path, train_msa_path, basename

path_to_outdirectory = path provided as argument to write new fragmented pmark benchmark files to
fragment_size = size of fragments as integer
neg_path = path to basename.neg 
tbl_path = path to basename.tbl
train_msa_path = path to basename.train.msa 
basename = the basename itself, not a path

Out:
None

"""
def copy_files(path_to_outdirectory, fragment_size, neg_path, tbl_path, train_msa_path, basename):
    

    #Describe the path to search for the default profmark files...   
   
    #Search for the default .tbl file in the search_path specified
    tbl = open(tbl_path,'r')
    #Make a copy of that table
    tblr = copy.deepcopy(tbl.read())
   
    #open a NEW file at the specified "split" location with the same name as the basename but "basename[cov].tbl"
    tblrout = open(path_to_outdirectory+"/"+basename+str(fragment_size)+".tbl",'w')
    tblrout.write(tbl.read())
    tblrout.write(tblr)
    tblrout.close()
   
    #Search for the default .neg file in the search_path specified
    neg = open(neg_path, 'r')
    negcopy = copy.deepcopy(neg.read())
    neg.close()
    #open a NEW file at the specified "split" location with the same name as the basename but "basename[cov].train.msa"
    negcopyout = open(path_to_outdirectory+"/"+basename+str(fragment_size)+".neg",'w')
    negcopyout.write(negcopy)
   

    #Search for the default .train.msa file in the search_path specified
    train = open(train_msa_path,'r')
    trainr = copy.deepcopy(train.read())
    train.close()
    #open a NEW file at the specified "split" location with the same name as the basename but "basename[cov].train.msa"
    trainrout = open(path_to_outdirectory+"/"+basename+str(fragment_size)+".train.msa",'w')
    trainrout.write(trainr)


"""
Main(Path_to_basename, path_to_outdirectory, Frag_Size, Frag_type)


MAIN FUNCTION, highest level part of the script. Follows these steps:

1) Check all files can be located
2) Check fragment type, set to default if something is off
3) Get basename
4) Actually do the fragmenting. Write files to specified outdirectory
5) Copy the files that are not altered to the outdirectory with the correct names. 

"""
def Main(Path_to_basename, path_to_outdirectory, Frag_Size, Frag_type):
    
    #First we check that all basename files can be located, and the splitdirectory folder exists
    Test_FA_path, pos_path, neg_path, tbl_path, train_msa_path = Check_files(Path_to_basename, path_to_outdirectory)


    #Check the validity of the fragment type and size. If looking for coverage > 100
    #or given a size in AA that exceeds the minimum size of a sequence, the program will end.
    Frag_type = CheckFragType_Validity(Frag_type, Frag_Size, pos_path)

    #Get basename (if input is C:\User\Directory\Folder\Basename then this just returns "Basename"):
    basename = Get_Basename(Path_to_basename)

    #Next we call a function to reconstruct our test.fa file
    duplicate_count = Process_TESTFA(Test_FA_path, pos_path, Frag_type, Frag_Size, path_to_outdirectory, basename)
    

    copy_files(path_to_outdirectory, Frag_Size, neg_path, tbl_path, train_msa_path, basename)

    
    if duplicate_count > 0:
        print(f"{duplicate_count} duplicates in test.fa file")
    else:
        print("No duplicates in the test.fa file")


#CALL MAIN
Main(Path_to_basename, Path_to_outdirectory, Frag_Size, Frag_type)
