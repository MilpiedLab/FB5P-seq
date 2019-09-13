#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on Jan 18, 2017

@author: cervera (motified by chuang)
'''
import pysam
import csv
import sys
import os
from re import sub

if len(sys.argv) != 3:
	print ("############################################")
	print ("\ni.e." ,sys.argv[0], "file.bam barcode.txt")
	print ("############################################")
	sys.exit()



def modifyXC(read):
    ''' replace XC tag, tags is from a pysam.AlignmentFile.tags, returns list of tags '''
    tag_value = read.get_tag(tag="XC")
    if(tag_value in BCsFromProtocol):
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
    




# ########################################
# ########################################
# ########################################

if len(sys.argv) == 3:
    assert sys.argv[1].endswith('.bam')
    inbamfn = sys.argv[1]
    inbamfn_modify=inbamfn[0:inbamfn.rindex('/')]+'/BAMbyBC'+inbamfn[inbamfn.rindex('/'):len(inbamfn)]
#    print (inbamfn_modify)
    directory = os.path.dirname(inbamfn_modify)
#    print (directory)
    if not os.path.exists(directory):
        os.makedirs(directory)
    inbam = pysam.AlignmentFile(inbamfn, 'rb', check_sq=False)
    
    # take the barcodes used in the protocol: BCsFromProtocol
    # #######################################################
    BCsFromProtocol_filename = sys.argv[2]
    #BCsFromProtocol_filename = '/Scripts/barcode_seq.txt'
    BCsFromProtocol=[]
    outbam_BCsFromProtocol_dict=dict()
    with open(BCsFromProtocol_filename, 'r') as f:
        reader = csv.reader(f, delimiter=';')
        for row in reader:
            BCsFromProtocol.append(row[0])
            outbamfn = sub('.bam$', '_'+row[0]+'.bam', inbamfn_modify)
#	    print(outbamfn)
#	    exit()
            outbam_BCsFromProtocol_dict[row[0]] =  pysam.AlignmentFile(outbamfn, 'wb', template=inbam) 
    
    n = 0
    for read in inbam:
        tag_value = read.get_tag(tag="XC")
        if(tag_value in BCsFromProtocol):
            outbam_BCsFromProtocol_dict[tag_value].write(read)
        
        n+=1
        if n % 100000 == 0:
            print (n)
    
    for BC in BCsFromProtocol:
        outbam_BCsFromProtocol_dict[BC].close()
    inbam.close()
    #print "Finished. Reads were written in multiple different files."
    
    #plate_name=inbamfn[(inbamfn.rindex('/')+1):len(inbamfn)]
    #plate_name=plate_name[0:plate_name.index('_')]
    #snakemake_output_file = open('/Data/'+plate_name+'/snakemake_output_split_bam_by_barcode.txt',”w”) 
    #snakemake_output_file.write("Finished. Reads were written in multiple different files.")
    #snakemake_output_file.close()
else:
    print ('you need to give a bam file')
