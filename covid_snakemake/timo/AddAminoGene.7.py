## Written by: Kate Johnson, kej310@nyu.edu
"""
Adding amino acid and gene information

Requires the UPDATED snplist files and csv file with features interested in

Snplist files are generated with the readreport.py script.
Updated snplist files are generated with the ConsensusFasta.Coverage.v4.py script.
Amino acid snplist files are generated with the AddAminoGene.5.py script.

The csv file with features that are of interest should have the following columns:
SEGMENT,START,END,NAME

python3 AddAminoGene.7.py --freqcut 0.001 \
                           --ref ../../reference/sars/sars-cov-2.fasta \
                           --var ../../marc_cryptic/varfiles/MCoV-116919_S289.SARS.NC_045512.2.0.001.snplist.csv \
                           --strain SARS \
                           --features ../../reference/sars/sars-cov-2-features.5.csv \
                           --save_dir ../../marc_cryptic/varfiles


Updated: 10.14.2022
Changed to work on one sample
Will need to provide segment information now (so will have to be ran 8 times for flu)

UPDATED: 03.05.2023: 
- Adjusted the printed outputs
- Fixed a copying error with pandas
- Adjusted the naming

Still need to: 
- INCLUDE SEGMENT INFORMATION

"""
import os
import glob
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--ref','-r', required=True,help='Indicate path and reference fasta file') #args.ref
parser.add_argument('--var','-v', required=True,help='Indicate variant file name')
#parser.add_argument('--path','-p', required=True,help='Indicate variant file name')
parser.add_argument('--strain','-s', required=True,help='Indicate strain')
#parser.add_argument('--segment','-c', required=True,help='Indicate segment/chrom information')
parser.add_argument('--features','-f', required=True,help='Indicate csv file with features of interest')
parser.add_argument('--save_dir','-sd', required=True,help='Indicate where new file should save')
parser.add_argument('--freqcut','-q', default=0.001,help='Indicate frequency cutoff')
args = parser.parse_args()

#amino acid dictionary
aminoacid = {'ttt': 'F', 'tct': 'S', 'tat': 'Y', 'tgt': 'C',
'ttc': 'F', 'tcc': 'S', 'tac': 'Y', 'tgc': 'C',
'tta': 'L', 'tca': 'S', 'taa': '*', 'tga': '*',
'ttg': 'L', 'tcg': 'S', 'tag': '*', 'tgg': 'W',
'ctt': 'L', 'cct': 'P', 'cat': 'H', 'cgt': 'R',
'ctc': 'L', 'ccc': 'P', 'cac': 'H', 'cgc': 'R',
'cta': 'L', 'cca': 'P', 'caa': 'Q', 'cga': 'R',
'ctg': 'L', 'ccg': 'P', 'cag': 'Q', 'cgg': 'R',
'att': 'I', 'act': 'T', 'aat': 'N', 'agt': 'S',
'atc': 'I', 'acc': 'T', 'aac': 'N', 'agc': 'S',
'ata': 'I', 'aca': 'T', 'aaa': 'K', 'aga': 'R',
'atg': 'M', 'acg': 'T', 'aag': 'K', 'agg': 'R',
'gtt': 'V', 'gct': 'A', 'gat': 'D', 'ggt': 'G',
'gtc': 'V', 'gcc': 'A', 'gac': 'D', 'ggc': 'G',
'gta': 'V', 'gca': 'A', 'gaa': 'E', 'gga': 'G',
'gtg': 'V', 'gcg': 'A', 'gag': 'E', 'ggg': 'G'}

def list_duplicates(ntpos_list):
    """
    INPUT: Tims code will sometimes duplicate positions that have only 1 read
    input list of ntpos to see if there are duplicates
    OUTPUT: return list of positions that have been duplicated
    """
    seen = set()
    seen_add = seen.add
    # adds all elements it doesn't know yet to seen and all other to seen_twice
    seen_twice = set( x for x in ntpos_list if x in seen or seen_add(x) )
    print("     Checking for dup positions. Positions observed multiple times: ",list(seen_twice))
    return list(seen_twice )

def read_fasta(fp):
    """
    INPUT: Fasta file reference
    """
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))

def open_fasta(filename):
    """
    INPUT: fasta file reference
    OUTPUT: dictionary with segment name as key and sequence of segment as value
    """
    segdict = {}
    with open(filename) as fp:
        for name, seq in read_fasta(fp):
            segdict[name[1:]] = seq
    return segdict

def creatList(r1,r2):
    """
    INPUT: the start and end of a gene position
    OUTPUT: list of positions within the range of the start and end of a gene
    """
    return list(range(r1,r2+1)) #python up to number you don't want included

def getindex(codon):
    """
    Identifying where the nucleotide falls in codon
    INPUT: major codon from snplist file
    OUTPUT: Index location of nucleotide
    """
    return [i for i, c in enumerate(codon) if c.isupper()]

def AdjustMinor(minordf):
    """
    INPUT: a codon df generated from translate function
    OUTPUT: a minor codon, minor amino acid, and nonysn info if present
    """
    print("")
    print("Gathering minor information")
    minordf = minordf.reset_index(drop=True)

    for index, row in minordf.iterrows():
        codon = row['majorcodon']
        aa = row['majoraa']
        #print(codon)

        if pd.isna(codon):  # will be 'na' if utr
            continue

        else:
            if '-' in codon: #if there is a deletion
                continue

            elif 'n' in codon or 'N' in codon: #if there is an N- will happen with primer trimming
                continue

            else:
                upperIndex = getindex(codon)[0] #identify where nucleotide is
                minorcodon = list(codon) #make codon a list of three nucleotides
                minorcodon[upperIndex] = row['minor'].upper() #put minor nt in index position and make uppercase
                minorcodon = ''.join(minorcodon) #join into string (no longer list)

                if '-' in minorcodon: #if there is a deletion for minor nt
                    continue

                else: #grab aa info for min
                    #print(minorcodon)
                    minoraa = aminoacid[minorcodon.lower()] #amino acid letter
                    if minoraa == aa:#if minor aa is the same as the major amino acid
                        nonsyn = 'syn'
                    elif minoraa != aa:
                        nonsyn = 'nonsyn'

                    #add all info to the main dataframe
                    minordf.loc[index,"minorcodon"] = minorcodon
                    minordf.loc[index,"minoraa"] = minoraa
                    minordf.loc[index,"nonsyn"] = nonsyn

    return minordf # return main dataframe of minor information



def translate(region_df, regionList, gene_name):
    """
    INPUT: subseted dataframe for region that we are interested
    OUTPUT: outputs a dataframe with updated information such as the major
    codon, amino acid, and gene id used with the region
    """
    print("Translating codons for: ", gene_name)
    region_df = region_df.copy() ## ADDED ON 03.06.2023
    region_df["major"] = region_df["major"].fillna("N") # change any missing data to N
    region_df = region_df.sort_values(by=["ntpos"])
    regionLength = len(regionList) #length of gene - must be divisible by three to translate
    aalist = list(range(0, regionLength, 3))  #generate a list going multiples of three
    dfntposLen = len(list(region_df['ntpos'])) #incase pos repeated due to poor alignments
    ntpos_list = list(region_df['ntpos'])

    if dfntposLen % 3 == 0:
        #print(dfntposLen%3)
        full_df = [] #empty df to append to
        for i in range(0, regionLength, 3): #iterate by 3's
            aaposition = aalist.index(i) + 1 #amino acid position - add 1 as python is 0based
            cod_df = region_df[i:i+3] #grab rows that correspond to amino acid
            cod_df = cod_df.reset_index(drop=True) #reset index (now 0-2)
            cList = list(region_df[i:i + 3]['major']) #grab the major nts for codon
            #print(cList)

            codonList = [x.lower() for x in cList] #make it lower to grab amino acid in dict
            codon = (''.join(codonList)).lower() #join the list of major nts to make string

            rlist = list(region_df[i:i + 3]['refnt'])
            refcodlist = [x.lower() for x in rlist]
            rcodon = (''.join(refcodlist)).lower()

            cod_df['gene_id'] = gene_name #add gene name for subset data
            cod_df['aapos'] = aaposition

            if 'n' in codon or '-' in codon:
                cod_df['gene_id'] = gene_name #add gene name for subset data
                cod_df['aapos'] = aaposition #add position
                full_df.append(cod_df) #but don't translate


            elif len(codon)==3:

                refaa = aminoacid[rcodon]
                r = refcodlist

                r[0] = r[0].upper()
                r1 = "".join(r) #join as string

                r[0] = r[0].lower() #change the first nt to lowercase
                r[1] = r[1].upper() #make second position in codon uppercase
                r2 = "".join(r)

                r[1] = r[1].lower() #make second position lower case
                r[2] = r[2].upper() #make third position uppercase
                r3 = "".join(r)

                cod_df['refaa'] = refaa  # add reference aa
                cod_df.loc[0,"refcodon"] = r1
                cod_df.loc[1,"refcodon"] = r2
                cod_df.loc[2,"refcodon"] = r3

                AA = aminoacid[codon] #grab the amino acid letter using the codon string
                c = codonList
                c[0] = c[0].upper() #for first nt position in codon make it upppercase
                first = "".join(c) #join as string

                c[0] = c[0].lower() #change the first nt to lowercase
                c[1] = c[1].upper() #make second position in codon uppercase
                second = "".join(c)

                c[1] = c[1].lower() #make second position lower case
                c[2] = c[2].upper() #make third position uppercase
                third = "".join(c)

                cod_df['majoraa'] = AA #add major amino acid information to dataframe
                cod_df.loc[0,"majorcodon"] = first
                cod_df.loc[1,"majorcodon"] = second
                cod_df.loc[2,"majorcodon"] = third

                full_df.append(cod_df) #append to empty list

        full_df=pd.concat(full_df,sort=False) #concatenate empty list into full dataframe
        return full_df

    elif dfntposLen % 3 != 0: #if not divisible by 3, still want to add info
        print("     {0} not divisible by 3. Length: {1}".format(gene_name, dfntposLen))
        #print(gene_name)
        x = list_duplicates(ntpos_list)  # check to make sure that there aren't multiple positions
        cod_df = region_df #just add gene id column name
        cod_df = cod_df.assign(gene_id = gene_name)
        return cod_df

def NoGen(subsetDF):
    """
    if there isn't a gene id/name in a region (5'utrs/3'utrs) just add empty gene id
    to the dataframe
    """
    print("")
    print('Attaching untranslated regions')
    subsetDF= subsetDF.assign(gene_id="")
    #print(subsetDF.head)
    return subsetDF

def markRegions(feature_df, df):
    """
    INPUT: Input a features dataframe and snplist variant file
    OUTPUT: Output info to add to new snplist file with aa information
    """
    print("Subsetting dataframe by specified region")
    print("")
    total_df = [] #empty df list to append and concatenate later
    tot_Regions = [] #looking at all ntpos covered through iterations
    
    for index, row in feature_df.iterrows(): #iterate through featers df
        region = creatList(row['START'],row['END']) #generate list of ntpos
        tot_Regions = tot_Regions + region #add ntpos to list
        df_region = df[df.ntpos.isin(region)] #subset dataframe for ntpos only in region of interest
        
        #check to make sure the df_region is the same length as the region we are interested in
        #print(len(region))
        #print(len(list(df_region['ntpos'])))

        gene_id = row['NAME'] #name of gene
        #print(gene_id)
        if len(region) == len(list(df_region['ntpos'])):
            trans_df = translate(df_region,region,gene_id) #pipe through translate function to add info
            total_df.append(trans_df) #append the output df to total df

    tot_Regions = list(set(tot_Regions)) #take the set of the ntpos list to remove duplicates
    noRegion = NoGen(df[~df.ntpos.isin(tot_Regions)]) #tilde indicates not in

    total_df.append(noRegion)
    total_df=pd.concat(total_df, sort=False)
    total_df.sort_values(by=['ntpos'])

    return total_df

def AddRefSeq(df,refseqdf):
    """
    INPUT: Reference sequence
    OUTPUT: Data frame with reference nt attached
    """
    print("Adding Reference Seq")

    m = pd.merge(df, refseqdf, on="ntpos", how="outer")

    return m


#### RUN CODE BELOW ####
infile = args.var
feat = pd.read_csv(args.features) #pandas will read the features csv file
STRAIN = args.strain
ref_dict = open_fasta(args.ref)
freqcut = args.freqcut
NT_LIST = ['A','T','G','C']


fname = str(infile)  # file name
df = pd.read_csv(infile)  # pandas will read the csv file create dataframe


SEGMENT = list(set(df['segment']))[0]  # pull segment info from file
seg_feats = feat[feat.SEGMENT == SEGMENT]  # filter features for segment
refseq = list(ref_dict[SEGMENT])  # make the ref seq a list
SegLength = list(range(1, len(refseq) + 1))  # generate list of ntpos for ref
refd = {"refnt":refseq,"ntpos":SegLength}  # build dict to build df
refdf = pd.DataFrame(refd)  # build ref df with ntpos and nt

# drop duplicate ntpos from snplist files (happens if low coverage)
df = df.drop_duplicates(subset=['ntpos'],keep='first')
df = AddRefSeq(df,refdf) #add ref info to dataframe
name = list(set(df['name']))[0]  # pull name from df - used to name csv
print("Input segment: ", SEGMENT)
print("Input name: ", name)
print("")

UpdatedDF = markRegions(seg_feats, df) # add features - aa info to majors
dfmin = UpdatedDF[(UpdatedDF.minor.isin(NT_LIST))]
nomin = UpdatedDF.merge(dfmin, how = 'outer' ,indicator=True).loc[lambda x : x['_merge']=='left_only'] # pull out positions that did not have a min

minchange = AdjustMinor(dfmin)
#minchange.to_csv("{0}/checkit3.csv".format(args.save_dir,name,STRAIN,SEGMENT,freqcut),index=False)


frames = [nomin, minchange]
final_df = pd.concat(frames,sort=False)  # concatenate frames
final_df = final_df.drop(columns=['_merge'])  # remove the 'merge' column from df

# checking major and minor information before printing
final_df['majRefSame'] = final_df.apply(lambda x: (x['major'] == x['refnt']), axis=1)
final_df['minRefSame'] = final_df.apply(lambda x: (x['minor'] == x['refnt']), axis=1)

final_len = len(list(set(final_df['ntpos'])))  # check to make sure all ntpos are accounted for



if "nonsyn" not in final_df: #the nonsyn column won't be added if no minor var
    final_df['nonsyn'] = ""
    final_df.to_csv("{0}/{1}.{2}.{3}.{4}.aa.snplist.csv".format(args.save_dir,name,STRAIN,SEGMENT,freqcut),index=False)

else:
    final_df.to_csv("{0}/{1}.{2}.{3}.{4}.aa.snplist.csv".format(args.save_dir,name,STRAIN,SEGMENT,freqcut),index=False)


if final_len == len(refseq):
    print("")
    print("FINISHED: all ntpos accounted for")
    print("Output named: ", "{0}/{1}.{2}.{3}.{4}.aa.snplist.csv".format(args.save_dir,name,STRAIN,SEGMENT,freqcut))

elif final_len != len(refseq): 
    print("")
    print("FINISHED: Error: ntpos positions missing in final output. Check input & output. Number output ntpos:", final_len)
    print("Output named: ", "{0}/{1}.{2}.{3}.{4}.aa.snplist.csv".format(args.save_dir,name,STRAIN,SEGMENT,freqcut))
