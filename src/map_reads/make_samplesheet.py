import glob

# this script reflects how our files were named from the sequencing core,
# where nf-core was run separately for the heart samples.
# This will have to be altered if using files from GEO
# but idea to make a nf-core compataible sampleshhet
FNs = glob.glob("../221017_A00405_0629_AHF3YJDSX5/*R1*.gz")

mylines = []
mylines.append("sample,fastq_1,fastq_2,strandedness")
for aFN in FNs:
    base_name = aFN.strip().split("/")[-1]
    samp_name = base_name.split("_DavJes")[0]
    R2FN = aFN.replace("R1","R2")
    mylines.append(f"{samp_name},{aFN},{R2FN},forward")

with open("map_reads/samplesheet.csv","w") as thefile:
    for item in mylines:
        thefile.write("%s\n"%item)


    
