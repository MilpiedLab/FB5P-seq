#!/usr/bin/env python
'''
Created on May 23, 2019

@author: cervera(modified by Dong)
'''
import pysam
import csv
import sys
from re import sub


if len(sys.argv) != 4:
	print ("############################################")
	print ("\ni.e." ,sys.argv[0], "file.bam barcode.txt output.bam")
	print ("############################################")
	sys.exit()


def modifyXC(read):
    ''' replace XC tag, tags is from a pysam.AlignmentFile.tags, returns list of tags '''
    tag_value = read.get_tag(tag="XC")
    if(tag_value not in BCsFromProtocol):
        for BC in BCsFromProtocol:
            mism_counter=0
            for j in range(0,8):
                if tag_value[j] != BC[j] : 
                    mism_counter+=1
                    if mism_counter > 1:
                        break
            if(mism_counter==1):
                # this BC has 1 mismatch so we merge it with a known BC
                read.set_tag(tag="XC", value=BC)
                break       
    

# take the barcodes used in the protocol: BCsFromProtocol
# #######################################################
BCsFromProtocol_filename = sys.argv[2]
#print ("###########", sys.argv[2])
#sys.exit()
BCsFromProtocol=[]
with open(BCsFromProtocol_filename, 'r') as f:
    reader = csv.reader(f, delimiter=';')
    for row in reader:
        BCsFromProtocol.append(row[0])


# ########################################
# ########################################
# ########################################

if len(sys.argv) == 4:
    assert sys.argv[1].endswith('.bam')
    inbamfn = sys.argv[1]
    outbamfn = sys.argv[3] #sub('.bam$', '_modifiedBC.bam', inbamfn)

    inbam = pysam.AlignmentFile(inbamfn, 'rb',check_sq=False)
    outbam = pysam.AlignmentFile(outbamfn, 'wb', template=inbam)
    
    n = 0
    for read in inbam:
        modifyXC(read)
        outbam.write(read)
        n+=1
        if n % 100000 == 0:
            print (n)
    
    outbam.close()
    inbam.close()
    print ("Finished. Reads were written to", outbamfn)
