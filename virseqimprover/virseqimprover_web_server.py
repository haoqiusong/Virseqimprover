import subprocess
import os
import math
import datetime
import sys

args = sys.argv
sid = args[1]
outputDir = "/home/lqzhang/virseqimprover_server/Samples/" + sid
spadesKmerlen = "default"
minOverlapCircular = 5000
minIdentityCircular = 95
salmonReadFraction = 0
minSuspiciousLen = 1000
global avgReadLen
avgReadLen = 0

n_fastq = 0
n_fq = 0
tem_fi = []
for images in os.listdir(outputDir):
    if (images.endswith(".fastq") and images[0] != '.'):
        n_fastq += 1
        tem_fi.append(images)
    elif (images.endswith(".fq") and images[0] != '.'):
        n_fq += 1
        tem_fi.append(images)
if len(tem_fi) == 1:
    if n_fastq == 0 and n_fq == 1:
        read1 = "/home/lqzhang/virseqimprover_server/Samples/" + sid + '/' + tem_fi[0]
        read2 = ""
    elif n_fastq == 1 and n_fq == 0:
        read1 = "/home/lqzhang/virseqimprover_server/Samples/" + sid + '/' + tem_fi[0]
        read2 = ""
elif len(tem_fi) == 2:
    read1 = "/home/lqzhang/virseqimprover_server/Samples/" + sid + '/' + tem_fi[0]
    read2 = "/home/lqzhang/virseqimprover_server/Samples/" + sid + '/' + tem_fi[1]
tem_fi2 = ""
for images in os.listdir(outputDir):
    if (images.endswith(".fasta") or images.endswith(".fa")):
        if images[0] != '.':
            tem_fi2 = images
scaffold = "/home/lqzhang/virseqimprover_server/Samples/" + sid + '/' + tem_fi2


def getReadLen():
    cmd = ""
    if read1.endswith(".gz"):
        cmd = "gzip -dc " + read1 + \
        " | awk 'NR%4 == 2 {lenSum+=length($0); readCount++;} END {print lenSum/readCount}'"
    else:
        cmd = "awk 'NR%4 == 2 {lenSum+=length($0); readCount++;} END {print lenSum/readCount}' " + read1

    shellFileWriter = open(outputDir + "/run.sh",'w')
    shellFileWriter.write('#'+"!/bin/bash\n")
    shellFileWriter.write(cmd)
    shellFileWriter.close()

    cmd = outputDir + "/run.sh"
    reader = subprocess.check_output(cmd, shell=True)
    reader = reader.decode()
    str_val = reader.replace("\n", "")
    avgReadLen = float(str_val)

    if avgReadLen == 0:
        print("Could not extract average read length from read file.")
        return
    print("Estimated average read length: " + str(avgReadLen))


def createBed():
    # move contig.fasta and contig.fasta.fai from spades-res folder to outputDir
    scaffoldFile = outputDir + "/scaffold-truncated/tmp/spades-res/scaffold.fasta"
    if (os.path.isfile(scaffoldFile) == True) and (os.path.isdir(scaffoldFile) == False):
        cmd = "cd " + outputDir + "/scaffold-truncated\n" + "rm scaffold.fasta*\n" + "cd tmp/spades-res\n" + "mv scaffold.fasta " + outputDir + "/scaffold-truncated\n" + "mv scaffold.fasta.fai " + outputDir + "/scaffold-truncated\n"

        shellFileWriter = open(outputDir + "/run.sh",'w')
        shellFileWriter.write('#'+"!/bin/bash\n")
        shellFileWriter.write(cmd)
        shellFileWriter.close()

        cmd = outputDir + "/run.sh"
        subprocess.check_output(cmd, shell=True)

    br = open(outputDir + "/scaffold-truncated/scaffold.fasta.fai")
    bw = open(outputDir + "/scaffold-truncated/scaffold-start-end.bed",'w')
    bwOutLog = open(outputDir + "/output-log.txt",'a')
    str_val = br.readline()
    results = []
    results = str_val.split("\t")
    scaffoldLength = int(results[1].strip())

    if scaffoldLength >= avgReadLen * 2:
        scaffoldId = results[0].strip()
        bw.write(scaffoldId + "\t" + '0' + "\t" + str(math.ceil(avgReadLen * 1.5)) + "\n")
        bw.write(scaffoldId + "\t" + str(math.ceil(scaffoldLength - avgReadLen * 1.5)) + "\t" + str(scaffoldLength) + "\n")
        bwOutLog.write("Trying to grow scaffold " + scaffoldId + " with length " + str(scaffoldLength) + "\n")
        if scaffoldLength > 300000:
            bwOutLog.write("Length of " + scaffoldId + \
                           " is already greater than 300kbp, so stop extending this one.\n")

    br.close()
    bw.close()
    bwOutLog.close()


def runAlignment():
    cmd = ""
    if len(read2) == 0:
        cmd = str("cd " + outputDir + "/scaffold-truncated\n" \
                  + "bedtools getfasta -fi scaffold.fasta -bed scaffold-start-end.bed -fo scaffold-start-end.fasta\n" \
                  + "rm -r salmon-index\n" \
                  + "rm -r salmon-res\n" \
                  + "rm salmon-mapped.sam\n" \
                  + "LD_PRELOAD=/lib/x86_64-linux-gnu/libpthread.so.0:/lib/x86_64-linux-gnu/librt.so.1 salmon index -t scaffold-start-end.fasta -i salmon-index\n" \
                  + "LD_PRELOAD=/lib/x86_64-linux-gnu/libpthread.so.0:/lib/x86_64-linux-gnu/librt.so.1 salmon quant -i salmon-index -l A " \
                  + "-r " + read1 + " -o salmon-res --writeMappings -p 16 --quasiCoverage " \
                  + str(salmonReadFraction) \
                  + " | samtools view -bS - | samtools view -h -F 0x04 - > salmon-mapped.sam\n")
    else:
        cmd = str("cd " + outputDir + "/scaffold-truncated\n" \
                  + "bedtools getfasta -fi scaffold.fasta -bed scaffold-start-end.bed -fo scaffold-start-end.fasta\n" \
                  + "rm -r salmon-index\n" \
                  + "rm -r salmon-res\n" \
                  + "rm salmon-mapped.sam\n" \
                  + "LD_PRELOAD=/lib/x86_64-linux-gnu/libpthread.so.0:/lib/x86_64-linux-gnu/librt.so.1 salmon index -t scaffold-start-end.fasta -i salmon-index\n" \
                  + "LD_PRELOAD=/lib/x86_64-linux-gnu/libpthread.so.0:/lib/x86_64-linux-gnu/librt.so.1 salmon quant -i salmon-index -l A " \
                  + "-1 " + read1 + " -2 " + read2 + " -o salmon-res --writeMappings -p 16 --quasiCoverage " \
                  + str(salmonReadFraction) \
                  + "| samtools view -bS - | samtools view -h -F 0x04 - > salmon-mapped.sam\n")

    shellFileWriter = open(outputDir + "/run.sh",'w')
    shellFileWriter.write('#'+"!/bin/bash\n")
    shellFileWriter.write(cmd)
    shellFileWriter.close()

    cmd = outputDir + "/run.sh"
    subprocess.check_output(cmd, shell=True)


def getMappedReads():
    cmd = ""
    if len(read2) == 0:
        cmd = str("cd " + outputDir + "/scaffold-truncated\n" \
                  + "rm -r tmp\n" \
                  + "mkdir tmp\n" \
                  + "bash filterbyname.sh in=" + read1 \
                  + " out=tmp/mapped_reads_1.fastq names=" \
                  + "salmon-mapped.sam include=t\n")
    else:
        cmd = str("cd " + outputDir + "/scaffold-truncated\n" \
                  + "rm -r tmp\n" \
                  + "mkdir tmp\n" \
                  + "bash filterbyname.sh in=" + read1 + " in2=" + read2 \
                  + " out=tmp/mapped_reads_1.fastq out2=tmp/mapped_reads_2.fastq names=" \
                  + "salmon-mapped.sam include=t\n")

    shellFileWriter = open(outputDir + "/run.sh",'w')
    shellFileWriter.write('#'+"!/bin/bash\n")
    shellFileWriter.write(cmd)
    shellFileWriter.close()

    cmd = outputDir + "/run.sh"
    subprocess.check_output(cmd, shell=True)


def runSpades():
    cmd = ""
    if len(read2) == 0:
        if spadesKmerlen == "default":
            cmd = str("cd " + outputDir + "/scaffold-truncated/tmp\n" \
                      + "spades.py -o " \
                      + "spades-res -s mapped_reads_1.fastq " \
                      + "--trusted-contigs ../scaffold.fasta" \
                      + " --only-assembler\n" \
                      + "cd spades-res\n" \
                      + "samtools faidx scaffolds.fasta\n")
        else:
            cmd = str("cd " + outputDir + "/scaffold-truncated/tmp\n" \
                      + "spades.py -o " \
                      + "spades-res -s mapped_reads_1.fastq " \
                      + "--trusted-contigs ../scaffold.fasta -k " + str(spadesKmerlen) \
                      + " --only-assembler\n" \
                      + "cd spades-res\n" \
                      + "samtools faidx scaffolds.fasta\n")
    else:
        if spadesKmerlen == "default":
            cmd = str("cd " + outputDir + "/scaffold-truncated/tmp\n" \
                      + "spades.py -o " \
                      + "spades-res -1 mapped_reads_1.fastq -2 mapped_reads_2.fastq " \
                      + "--trusted-contigs ../scaffold.fasta" \
                      + " --only-assembler\n" \
                      + "cd spades-res\n" \
                      + "samtools faidx scaffolds.fasta\n")
        else:
            cmd = str("cd " + outputDir + "/scaffold-truncated/tmp\n" \
                      + "spades.py -o " \
                      + "spades-res -1 mapped_reads_1.fastq -2 mapped_reads_2.fastq " \
                      + "--trusted-contigs ../scaffold.fasta -k " + str(spadesKmerlen) \
                      + " --only-assembler\n" \
                      + "cd spades-res\n" \
                      + "samtools faidx scaffolds.fasta\n")

    shellFileWriter = open(outputDir + "/run.sh",'w')
    shellFileWriter.write('#'+"!/bin/bash\n")
    shellFileWriter.write(cmd)
    shellFileWriter.close()

    cmd = outputDir + "/run.sh"
    subprocess.check_output(cmd, shell=True)


def getScaffoldFromScaffolds():
    maxScaffoldLength = 0
    maxScaffoldId = ""

    scaffoldFile = outputDir + "/scaffold-truncated/tmp/spades-res/scaffolds.fasta"
    if (os.path.isfile(scaffoldFile) == True) and (os.path.isdir(scaffoldFile) == False):
        br = open(outputDir + "/scaffold-truncated/tmp/spades-res/scaffolds.fasta.fai")
        str_val = br.readline()
        str_val.strip()
        results = []
        results = str_val.split("\t")
        maxScaffoldId = results[0].strip()
        maxScaffoldLength = int(results[1].strip())
        br.close()
        if (maxScaffoldId != "") and (maxScaffoldLength != 0):
            cmd = str("cd " + outputDir + "/scaffold-truncated/tmp/spades-res\n" \
                      + "bash filterbyname.sh " \
                      + "in=scaffolds.fasta " \
                      + "out=scaffold.fasta names=" + maxScaffoldId \
                      + " include=t\n" \
                      + "samtools faidx scaffold.fasta\n")

            shellFileWriter = open(outputDir + "/run.sh",'w')
            shellFileWriter.write('#'+"!/bin/bash\n")
            shellFileWriter.write(cmd)
            shellFileWriter.close()

            cmd = outputDir + "/run.sh"
            subprocess.check_output(cmd, shell=True)  
    else:
        br = open(outputDir + "/scaffold-truncated/scaffold.fasta.fai")
        str_val = br.readline()
        str_val.strip()
        results = []
        results = str_val.split("\t")
        maxScaffoldLength = int(results[1].strip())
        br.close()

    return maxScaffoldLength


def getScaffoldLenFromTruncatedExtend():
    scaffoldLength = 0
    scaffoldFile = str(outputDir + "/scaffold-truncated/scaffold.fasta")
    if (os.path.exists(scaffoldFile)) == True and (os.path.isdir(scaffoldFile) == False):
        br = open(outputDir + "/scaffold-truncated/scaffold.fasta.fai")
        str_val = br.readline()
        str_val.strip()
        results = []
        results = str_val.split("\t")
        scaffoldLength = int(results[1].strip())
        br.close()
    else:
        br = open(outputDir + "/scaffold.fasta.fai")
        str_val = br.readline()
        str_val.strip()
        results = []
        results = str_val.split("\t")
        scaffoldLength = int(results[1].strip())
        br.close()

    return scaffoldLength


def updateCurrentScaffold():
    scaffoldFile = str(outputDir + "/scaffold-truncated/scaffold.fasta")
    if (os.path.exists(scaffoldFile) == True) and (os.path.isdir(scaffoldFile) == False):
        cmd = str("cd " + outputDir + "\n" \
                  + "rm scaffold.fasta*\n" \
                  + "cd scaffold-truncated\n" \
                  + "cp scaffold.fasta " + outputDir + "\n" \
                  + "cp scaffold.fasta.fai " + outputDir + "\n")

        shellFileWriter = open(outputDir + "/run.sh",'w')
        shellFileWriter.write('#'+"!/bin/bash\n")
        shellFileWriter.write(cmd)
        shellFileWriter.close()

        cmd = outputDir + "/run.sh"
        subprocess.check_output(cmd, shell=True)


def createBedForTruncatedScaffold(truncatedLen):
    br = open(outputDir + "/scaffold.fasta.fai")
    bw = open(outputDir + "/scaffold-truncated.bed",'w')
    str_val = br.readline()
    results = str_val.split("\t")
    scaffoldLength = int(results[1].strip())
    if scaffoldLength >= truncatedLen:
        scaffoldId = results[0].strip()
        bw.write(scaffoldId + "\t" + str(truncatedLen) + "\t" + str(scaffoldLength - truncatedLen) + "\n")
    br.close()
    bw.close()


def growScaffoldWithAssembly():
    iteration = 1
    extendContig = True
    prevLength = 0
    while(extendContig):
        currentLength = getScaffoldFromScaffolds()
        if currentLength > 300000:
            extendContig = False
        elif currentLength > prevLength:
            prevLength = currentLength
            createBed()
            runAlignment()
            getMappedReads()
            runSpades()

            scaffoldFile = outputDir + "/scaffold-truncated/tmp/spades-res/scaffolds.fasta"
            if (os.path.isfile(scaffoldFile) == True) and (os.path.isdir(scaffoldFile) == False):
                iteration += 1
                extendContig = True
            else:
                extendContig = False
        else:
            extendContig = False

        if extendContig and iteration > 100:
            extendContig = False


def getTruncatedScaffoldAndExtend(currentLength):
    cmd = ""

    # truncate 300bp
    createBedForTruncatedScaffold(300)
    cmd = str("cd " + outputDir + "\n" \
              + "rm -r scaffold-truncated\n" \
              + "mkdir scaffold-truncated\n" \
              + "bedtools getfasta -fi scaffold.fasta -bed scaffold-truncated.bed -fo scaffold-truncated/scaffold.fasta\n" \
              + "cd scaffold-truncated\n" \
              + "samtools faidx scaffold.fasta\n")

    shellFileWriter = open(outputDir + "/run.sh",'w')
    shellFileWriter.write('#'+"!/bin/bash\n")
    shellFileWriter.write(cmd)
    shellFileWriter.close()

    cmd = outputDir + "/run.sh"
    subprocess.check_output(cmd, shell=True)

    growScaffoldWithAssembly()
    lengthFromGrowingTruncatedScaffold = 0
    lengthFromGrowingTruncatedScaffold = getScaffoldLenFromTruncatedExtend()
    if lengthFromGrowingTruncatedScaffold > currentLength:
        return
    else:
        # truncate 500bp
        createBedForTruncatedScaffold(500)
        cmd = str("cd " + outputDir + "\n" \
                  + "rm -r scaffold-truncated\n" \
                  + "mkdir scaffold-truncated\n" \
                  + "bedtools getfasta -fi scaffold.fasta -bed scaffold-truncated.bed -fo scaffold-truncated/scaffold.fasta\n" \
                  + "cd scaffold-truncated\n" \
                  + "samtools faidx scaffold.fasta\n")

        shellFileWriter = open(outputDir + "/run.sh",'w')
        shellFileWriter.write('#'+"!/bin/bash\n")
        shellFileWriter.write(cmd)
        shellFileWriter.close()

        cmd = outputDir + "/run.sh"
        subprocess.check_output(cmd, shell=True)

        growScaffoldWithAssembly()
        lengthFromGrowingTruncatedScaffold = getScaffoldLenFromTruncatedExtend()
        if lengthFromGrowingTruncatedScaffold > currentLength:
            return
        else:
            # truncate 700bp
            createBedForTruncatedScaffold(700)
            cmd = str("cd " + outputDir + "\n" \
                      + "rm -r scaffold-truncated\n" \
                      + "mkdir scaffold-truncated\n" \
                      + "bedtools getfasta -fi scaffold.fasta -bed scaffold-truncated.bed -fo scaffold-truncated/scaffold.fasta\n" \
                      + "cd scaffold-truncated\n" \
                      + "samtools faidx scaffold.fasta\n")

            shellFileWriter = open(outputDir + "/run.sh",'w')
            shellFileWriter.write('#'+"!/bin/bash\n")
            shellFileWriter.write(cmd)
            shellFileWriter.close()

            cmd = outputDir + "/run.sh"
            subprocess.check_output(cmd, shell=True)

            growScaffoldWithAssembly()
            lengthFromGrowingTruncatedScaffold = getScaffoldLenFromTruncatedExtend()
            if lengthFromGrowingTruncatedScaffold > currentLength:
                return
            else:
                # truncate 1000bp
                createBedForTruncatedScaffold(1000)
                cmd = str("cd " + outputDir + "\n" \
                          + "rm -r scaffold-truncated\n" \
                          + "mkdir scaffold-truncated\n" \
                          + "bedtools getfasta -fi scaffold.fasta -bed scaffold-truncated.bed -fo scaffold-truncated/scaffold.fasta\n" \
                          + "cd scaffold-truncated\n" \
                          + "samtools faidx scaffold.fasta\n")

                shellFileWriter = open(outputDir + "/run.sh",'w')
                shellFileWriter.write('#'+"!/bin/bash\n")
                shellFileWriter.write(cmd)
                shellFileWriter.close()

                cmd = outputDir + "/run.sh"
                subprocess.check_output(cmd, shell=True)

                growScaffoldWithAssembly()
                lengthFromGrowingTruncatedScaffold = getScaffoldLenFromTruncatedExtend()
                if lengthFromGrowingTruncatedScaffold > currentLength:
                    return
                else:
                    # truncate 1300bp
                    createBedForTruncatedScaffold(1300)
                    cmd = str("cd " + outputDir + "\n" \
                              + "rm -r scaffold-truncated\n" \
                              + "mkdir scaffold-truncated\n" \
                              + "bedtools getfasta -fi scaffold.fasta -bed scaffold-truncated.bed -fo scaffold-truncated/scaffold.fasta\n" \
                              + "cd scaffold-truncated\n" \
                              + "samtools faidx scaffold.fasta\n")

                    shellFileWriter = open(outputDir + "/run.sh",'w')
                    shellFileWriter.write('#'+"!/bin/bash\n")
                    shellFileWriter.write(cmd)
                    shellFileWriter.close()

                    cmd = outputDir + "/run.sh"
                    subprocess.check_output(cmd, shell=True)

                    growScaffoldWithAssembly()
                    lengthFromGrowingTruncatedScaffold = getScaffoldLenFromTruncatedExtend()
                    if lengthFromGrowingTruncatedScaffold > currentLength:
                        return
                    else:
                        # truncate 1500bp
                        createBedForTruncatedScaffold(1500)
                        cmd = str("cd " + outputDir + "\n" \
                                  + "rm -r scaffold-truncated\n" \
                                  + "mkdir scaffold-truncated\n" \
                                  + "bedtools getfasta -fi scaffold.fasta -bed scaffold-truncated.bed -fo scaffold-truncated/scaffold.fasta\n" \
                                  + "cd scaffold-truncated\n" \
                                  + "samtools faidx scaffold.fasta\n")

                        shellFileWriter = open(outputDir + "/run.sh",'w')
                        shellFileWriter.write('#'+"!/bin/bash\n")
                        shellFileWriter.write(cmd)
                        shellFileWriter.close()

                        cmd = outputDir + "/run.sh"
                        subprocess.check_output(cmd, shell=True)

                        growScaffoldWithAssembly()
                        lengthFromGrowingTruncatedScaffold = getScaffoldLenFromTruncatedExtend()
                        if lengthFromGrowingTruncatedScaffold > currentLength:
                            return
                        else:
                            # truncate 1700bp
                            createBedForTruncatedScaffold(1700)
                            cmd = str("cd " + outputDir + "\n" \
                                      + "rm -r scaffold-truncated\n" \
                                      + "mkdir scaffold-truncated\n" \
                                      + "bedtools getfasta -fi scaffold.fasta -bed scaffold-truncated.bed -fo scaffold-truncated/scaffold.fasta\n" \
                                      + "cd scaffold-truncated\n" \
                                      + "samtools faidx scaffold.fasta\n")

                            shellFileWriter = open(outputDir + "/run.sh",'w')
                            shellFileWriter.write('#'+"!/bin/bash\n")
                            shellFileWriter.write(cmd)
                            shellFileWriter.close()

                            cmd = outputDir + "/run.sh"
                            subprocess.check_output(cmd, shell=True)

                            growScaffoldWithAssembly()
                            lengthFromGrowingTruncatedScaffold = getScaffoldLenFromTruncatedExtend()
                            if lengthFromGrowingTruncatedScaffold > currentLength:
                                return
                            else:
                                # truncate 2000bp
                                createBedForTruncatedScaffold(2000)
                                cmd = str("cd " + outputDir + "\n" \
                                          + "rm -r scaffold-truncated\n" \
                                          + "mkdir scaffold-truncated\n" \
                                          + "bedtools getfasta -fi scaffold.fasta -bed scaffold-truncated.bed -fo scaffold-truncated/scaffold.fasta\n" \
                                          + "cd scaffold-truncated\n" \
                                          + "samtools faidx scaffold.fasta\n")

                                shellFileWriter = open(outputDir + "/run.sh",'w')
                                shellFileWriter.write('#'+"!/bin/bash\n")
                                shellFileWriter.write(cmd)
                                shellFileWriter.close()

                                cmd = outputDir + "/run.sh"
                                subprocess.check_output(cmd, shell=True)

                                growScaffoldWithAssembly()


def checkCircularity():
    file_list = os.listdir(outputDir)
    has_doc = False
    for file in file_list:
        path = os.path.join(outputDir, file)
        if "blastn-" in path:
            has_doc = True
            break
    if has_doc:
        cmd = "cd " + outputDir + "\n" + "rm -r blastn-*\n"

        shellFileWriter = open(outputDir + "/run.sh",'w')
        shellFileWriter.write('#'+"!/bin/bash\n")
        shellFileWriter.write(cmd)
        shellFileWriter.close()

        cmd = outputDir + "/run.sh"
        subprocess.check_output(cmd, shell=True)

    # create query and subject fasta
    br = open(outputDir + "/scaffold.fasta.fai")
    bw1 = open(outputDir + "/blastn-subject-1stround.bed",'w')
    bw2 = open(outputDir + "/blastn-query-1stround.bed",'w')
    str_val = br.readline()
    results = []
    results = str_val.split("\t")
    scaffoldLength = int(results[1].strip())
    intReadLen = int(avgReadLen)
    if scaffoldLength > (intReadLen*2):
        scaffoldId = results[0].strip()
        bw1.write(str(scaffoldId) + "\t" + str(0) + "\t" + str(scaffoldLength - (intReadLen*2)) + "\n")
        bw2.write(str(scaffoldId) + "\t" + str(scaffoldLength - (intReadLen*2)) + "\t" + str(scaffoldLength-intReadLen) + "\n")
    br.close()
    bw1.close()
    bw2.close()

    cmd = str("cd " + outputDir + "\n" \
              + "bedtools getfasta -fi scaffold.fasta -bed blastn-subject-1stround.bed -fo blastn-subject-1stround.fasta\n" \
              + "bedtools getfasta -fi scaffold.fasta -bed blastn-query-1stround.bed -fo blastn-query-1stround.fasta\n" \
              + "samtools faidx blastn-subject-1stround.fasta\n" \
              + "samtools faidx blastn-query-1stround.fasta\n" \
              + "makeblastdb -in blastn-subject-1stround.fasta -dbtype nucl\n" \
              + "blastn -query blastn-query-1stround.fasta -db blastn-subject-1stround.fasta -num_threads 16 -outfmt '7' -out blastn-res-1stround.txt\n")

    shellFileWriter = open(outputDir + "/run.sh",'w')
    shellFileWriter.write('#'+"!/bin/bash\n")
    shellFileWriter.write(cmd)
    shellFileWriter.close()

    cmd = outputDir + "/run.sh"
    subprocess.check_output(cmd, shell=True)

    # parse blastn result
    isCircular = False
    subjectStart = 0
    queryStart = 0
    br = open(outputDir + "/blastn-res-1stround.txt")
    bwCircularOutputLog = open(outputDir + "/circularity-output-log.txt",'a')
    bwOutputLog = open(outputDir + "/output-log.txt",'a')
    str_val = br.readline()
    results = []
    while(str_val != "" and isCircular == False):
        if str_val[0] != "#":
            str_val = str_val.strip()
            results = str_val.split("\t")
            percentIden = float(results[2].strip())
            alignmentLen = int(results[3].strip())
            subjectStart = int(results[8].strip())
            queryStart = int(results[6].strip())
            minAlignmentLength = int(avgReadLen*0.95)

            if (percentIden >= minIdentityCircular) and (alignmentLen >= minAlignmentLength):
                bwCircularOutputLog.write("Scaffold seems circular. Scaffold position " + str(results[6]) \
                                          + " to " + str(results[7]) + " mapped to position " + str(results[8]) + " to " + str(results[9]) \
                                          + " with " + str(percentIden) + "% identity and " + str(alignmentLen) + " alignment length.\n")
                bwOutputLog.write("Scaffold seems circular. Scaffold position " + str(results[6]) \
                                  + " to " + str(results[7]) + " mapped to position " + str(results[8]) + " to " + str(results[9]) \
                                  + " with " + str(percentIden) + "% identity and " + str(alignmentLen) + " alignment length.\n")
                isCircular = True
            else:
                bwCircularOutputLog.write("Scaffold does not seem to be circular. Scaffold position " + str(results[6]) \
                                          + " to " + str(results[7]) + " mapped to position " + str(results[8]) + " to " + str(results[9]) \
                                          + " with " + str(percentIden) + "% identity and " + str(alignmentLen) + " alignment length.\n")
                bwOutputLog.write("Scaffold does not seem to be circular. Scaffold position " + str(results[6]) \
                                  + " to " + str(results[7]) + " mapped to position " + str(results[8]) + " to " + str(results[9]) \
                                  + " with " + str(percentIden) + "% identity and " + str(alignmentLen) + " alignment length.\n")

        str_val = br.readline()
    br.close()
    bwOutputLog.close()
    bwCircularOutputLog.close()

    # blastn subject query 2nd round
    if isCircular:
        br = open(outputDir + "/scaffold.fasta.fai")
        bw1 = open(outputDir + "/blastn-subject-2ndround.bed",'w')
        bw2 = open(outputDir + "/blastn-query-2ndround.bed",'w')
        str_val = br.readline()
        results = []
        results = str_val.split("\t")
        scaffoldLength = int(results[1].strip())
        intReadLen = int(avgReadLen)
        if(scaffoldLength > (intReadLen*2)):
            scaffoldId = results[0].strip()
            bw1.write(str(scaffoldId) + "\t" + str(0) + "\t" + str(subjectStart) + "\n")
            bw2.write(str(scaffoldId) + "\t" + str(scaffoldLength - (intReadLen*2) - subjectStart) + "\t" + str(scaffoldLength-(intReadLen*2)) + "\n")
        br.close()
        bw1.close()
        bw2.close()

        cmd = str("cd " + outputDir + "\n" \
                  + "bedtools getfasta -fi scaffold.fasta -bed blastn-subject-2ndround.bed -fo blastn-subject-2ndround.fasta\n" \
                  + "bedtools getfasta -fi scaffold.fasta -bed blastn-query-2ndround.bed -fo blastn-query-2ndround.fasta\n" \
                  + "samtools faidx blastn-subject-2ndround.fasta\n" \
                  + "samtools faidx blastn-query-2ndround.fasta\n" \
                  + "makeblastdb -in blastn-subject-2ndround.fasta -dbtype nucl\n" \
                  + "blastn -query blastn-query-2ndround.fasta -db blastn-subject-2ndround.fasta -num_threads 16 -outfmt '7' -out blastn-res-2ndround.txt\n")

        shellFileWriter = open(outputDir + "/run.sh",'w')
        shellFileWriter.write('#'+"!/bin/bash\n")
        shellFileWriter.write(cmd)
        shellFileWriter.close()

        cmd = outputDir + "/run.sh"
        subprocess.check_output(cmd, shell=True)

    return isCircular


def runAlignmentGetCoverage():
    cmd = ""

    file_list = os.listdir(outputDir)
    has_bowtie2 = False
    has_samtools = False
    for file in file_list:
        path = os.path.join(outputDir, file)
        if "bowtie2" in path:
            has_bowtie2 = True
            break
    if has_bowtie2:
        cmd = "cd " + outputDir + "\n" + "rm -r bowtie2*\n"

        shellFileWriter = open(outputDir + "/run.sh",'w')
        shellFileWriter.write('#'+"!/bin/bash\n")
        shellFileWriter.write(cmd)
        shellFileWriter.close()

        cmd = outputDir + "/run.sh"
        subprocess.check_output(cmd, shell=True)

    for file in file_list:
        path = os.path.join(outputDir, file)
        if "samtools" in path:
            has_samtools = True
            break
    if has_samtools:
        cmd = "cd " + outputDir + "\n" + "rm -r samtools*\n"

        shellFileWriter = open(outputDir + "/run.sh",'w')
        shellFileWriter.write('#'+"!/bin/bash\n")
        shellFileWriter.write(cmd)
        shellFileWriter.close()

        cmd = outputDir + "/run.sh"
        subprocess.check_output(cmd, shell=True)

    if len(read2) == 0:
        cmd = str("cd " + outputDir + "\n" \
                  + "bowtie2-build scaffold.fasta bowtie2-index\n" \
                  + "bowtie2 -x bowtie2-index " \
                  + "-U " + read1 \
                  + " | samtools view -bS - | samtools view -h -F 0x04 -b - | " \
                  + "samtools sort - -o bowtie2-mapped.bam\n" \
                  + "samtools depth -a bowtie2-mapped.bam > samtools-coverage.txt\n")
    else:
        cmd = str("cd " + outputDir + "\n" \
                  + "bowtie2-build scaffold.fasta bowtie2-index\n" \
                  + "bowtie2 -x bowtie2-index " \
                  + "-1 " + read1 \
                  + " -2 " + read2 \
                  + " | samtools view -bS - | samtools view -h -F 0x04 -b - | " \
                  + "samtools sort - -o bowtie2-mapped.bam\n" \
                  + "samtools depth -a bowtie2-mapped.bam > samtools-coverage.txt\n")

    shellFileWriter = open(outputDir + "/run.sh",'w')
    shellFileWriter.write('#'+"!/bin/bash\n")
    shellFileWriter.write(cmd)
    shellFileWriter.close()

    cmd = outputDir + "/run.sh"
    subprocess.check_output(cmd, shell=True)


def percentile(values, percentiles):
    percentileResults = []
    temp = values
    values = values.sort()
    for percentile in percentiles:
        index = int(math.ceil((percentile / 100.00) * len(temp)))
        percentileResults.append(temp[index - 1])
    return percentileResults


def readCoverageGetPercentile():
    coverages = []
    br = open(outputDir + "/samtools-coverage.txt")
    str_val = br.readline()
    results = []
    while(str_val != ""):
        str_val = str_val.strip()
        results = str_val.split("\t")
        coverage = int(results[2])
        if coverage != 0:
            coverages.append(coverage)
        str_val = br.readline()
    br.close()

    percentiles = []
    percentiles.append(15.00)
    percentiles.append(85.00)
    percentileResults = percentile(coverages, percentiles)

    return percentileResults


def writeCoverageQuantile(quantile15Percent, quantile85Percent):
    hasSuspiciousRegion = False
    scaffoldLength = 0
    scaffoldId = ""
    suspiciousStarts = []
    suspiciousEnds = []

    br = open(outputDir + "/scaffold.fasta.fai")
    str_val = br.readline()
    str_val = str_val.strip()
    results = []
    results = str_val.split("\t")
    scaffoldLength = int(results[1].strip())
    scaffoldId = results[0].strip()
    br.close()

    str_val = ""
    results = []
    coverage = 0
    startBase = 0
    currentBase = 0
    lowHigh = False

    br = open(outputDir + "/samtools-coverage.txt")
    str_val = br.readline()
    while(str_val != ""):
        str_val = str_val.strip()
        results = str_val.split("\t")
        coverage = int(results[2])
        currentBase = int(results[1])
        if (coverage >= quantile15Percent) and (coverage <= quantile85Percent):
            if lowHigh:
                if (currentBase - startBase) > minSuspiciousLen:
                    suspiciousStarts.append(startBase)
                    suspiciousEnds.append(currentBase-1)
            lowHigh = False
        else:
            if (startBase == 0) or (lowHigh == False):
                startBase = currentBase
            lowHigh = True
        str_val = br.readline()
    br.close()

    if len(suspiciousStarts) == 0:
        hasSuspiciousRegion = False
        bwCoverageOutputLog = open(outputDir + "/suspicious-regions.log",'w')
        bwCoverageOutputLog.write(scaffoldId + " has no suspicious region\n")
        bwCoverageOutputLog.close()
    else:
        hasSuspiciousRegion = True
        maxStartBase = 1
        maxLength = 0
        start = 0
        bwCoverageOutputLog = open(outputDir + "/suspicious-regions.log",'w')
        bwCoverageOutputLog.write(scaffoldId + " has " + str(len(suspiciousStarts)) + " suspicious regions\n")

        # get start and end of longest non-suspicious region
        for i in range(len(suspiciousStarts)):
            start = suspiciousStarts[i]
            if maxLength == 0:
                maxLength = start - 1
            else:
                if (start - suspiciousEnds[i-1]) > maxLength:
                    maxLength = start - suspiciousEnds[i-1]
                    maxStartBase = suspiciousEnds[i-1] + 1

        if (scaffoldLength - suspiciousEnds[len(suspiciousEnds)-1]) > maxLength:
            maxLength = scaffoldLength - suspiciousEnds[len(suspiciousEnds)-1]
            maxStartBase = suspiciousEnds[len(suspiciousEnds)-1] + 1

        bwCoverageOutputLog.close()

        # write bed file
        bw = open(outputDir + "/longest-non-suspicious.bed",'w')
        bw.write(scaffoldId + "\t" + str(maxStartBase) + "\t" + str(maxStartBase + maxLength - 1) + "\n")
        bw.close()

        # get fasta from bed file
        cmd = str("cd " + outputDir + "\n" \
                  + "rm -r longest-non-suspicious.fasta*\n" \
                  + "bedtools getfasta -fi scaffold.fasta -bed longest-non-suspicious.bed -fo longest-non-suspicious.fasta\n" \
                  + "samtools faidx longest-non-suspicious.fasta\n")

        shellFileWriter = open(outputDir + "/run.sh",'w')
        shellFileWriter.write('#'+"!/bin/bash\n")
        shellFileWriter.write(cmd)
        shellFileWriter.close()

        cmd = outputDir + "/run.sh"
        subprocess.check_output(cmd, shell=True)

    return hasSuspiciousRegion


def updateScaffoldWithLongest():
    scaffoldFile = outputDir + "/scaffold.fasta"
    if (os.path.isfile(scaffoldFile) == True) and (os.path.isdir(scaffoldFile) == False):
        cmd = str("cd " + outputDir + "\n" \
                  + "rm scaffold.fasta*\n" \
                  + "cp longest-non-suspicious.fasta scaffold.fasta\n" \
                  + "samtools faidx scaffold.fasta\n")

        shellFileWriter = open(outputDir + "/run.sh",'w')
        shellFileWriter.write('#'+"!/bin/bash\n")
        shellFileWriter.write(cmd)
        shellFileWriter.close()

        cmd = outputDir + "/run.sh"
        subprocess.check_output(cmd, shell=True)


def copyScaffoldForAssembly():
    scaffoldFile = outputDir + "/scaffold.fasta"
    if (os.path.isfile(scaffoldFile) == True) and (os.path.isdir(scaffoldFile) == False):
        cmd = str("cd " + outputDir + "\n" \
                    + "rm -r scaffold-truncated\n" \
                    + "mkdir scaffold-truncated\n" \
                    + "cp scaffold.fasta scaffold-truncated/scaffold.fasta\n" \
                    + "cd scaffold-truncated\n" \
                    + "samtools faidx scaffold.fasta\n")

        shellFileWriter = open(outputDir + "/run.sh",'w')
        shellFileWriter.write('#'+"!/bin/bash\n")
        shellFileWriter.write(cmd)
        shellFileWriter.close()

        cmd = outputDir + "/run.sh"
        subprocess.check_output(cmd, shell=True)


def checkCoverage():
    runAlignmentGetCoverage()
    percentiles = readCoverageGetPercentile()
    hasSuspiciousRegion = writeCoverageQuantile(percentiles[0], percentiles[1])

    bwOutLog = open(outputDir + "/output-log.txt",'a')
    if hasSuspiciousRegion:
        updateScaffoldWithLongest()
        bwOutLog.write("Scaffold got updated for suspicious region.\n")
        bwOutLog.close()
        copyScaffoldForAssembly()
        growScaffoldWithAssembly()
    else:
        bwOutLog.write("Scaffold didn't get updated as there is no suspicious region.\n")
        bwOutLog.close()

    return hasSuspiciousRegion


def extendOneScaffold():
    iteration = 1
    cmd = str("cd " + outputDir + "\n" \
              + "cp " + scaffold + " scaffold.fasta\n" \
              + "samtools faidx scaffold.fasta\n")

    shellFileWriter = open(outputDir + "/run.sh",'w')
    shellFileWriter.write('#'+"!/bin/bash\n")
    shellFileWriter.write(cmd)
    shellFileWriter.close()

    cmd = outputDir + "/run.sh"
    subprocess.check_output(cmd, shell=True)

    extendContig = True
    prevLength = 0
    needUpdate = False

    while(extendContig):
        currentLength = getScaffoldLenFromTruncatedExtend()
        if currentLength > 300000:
            extendContig = False
        elif currentLength > prevLength:
            prevLength = currentLength
            updateCurrentScaffold()
            if checkCircularity():
                extendContig = False
            else:
                needUpdate = checkCoverage()
                if needUpdate:
                    # update current length
                    newLength = getScaffoldLenFromTruncatedExtend()
                    updateCurrentScaffold()
                    prevLength = newLength
                    getTruncatedScaffoldAndExtend(newLength)
                else:
                    getTruncatedScaffoldAndExtend(currentLength)
                iteration += 1
        else:
            extendContig = False

        if extendContig and iteration > 200:
            extendContig = False


#print("Finished parsing input arguments.")
#getReadLen()


print("Started growing scaffold.")
extendOneScaffold()
print("Finished growing scaffold.")

if len(read2) == 0:
    cmd = str("cd " + outputDir + "\n" \
            + "bowtie2-build scaffold.fasta bowtie2-index\n" \
            + "bowtie2 -x bowtie2-index -U " + read1 + " | samtools view -bS - | samtools view -h -F 0x04 -b - | samtools sort - -o bowtie2-mapped.bam\n" \
            + "samtools index bowtie2-mapped.bam\n" \
            + "pilon --genome scaffold.fasta --bam bowtie2-mapped.bam --output pilon_out --outdir ./")
else:
    cmd = str("cd " + outputDir + "\n" \
            + "bowtie2-build scaffold.fasta bowtie2-index\n" \
            + "bowtie2 -x bowtie2-index -1 " + read1 + "-2 " + read2 + " | samtools view -bS - | samtools view -h -F 0x04 -b - | samtools sort - -o bowtie2-mapped.bam\n" \
            + "samtools index bowtie2-mapped.bam\n" \
            + "pilon --genome scaffold.fasta --bam bowtie2-mapped.bam --output pilon_out --outdir ./")
shellFileWriter = open(outputDir + "/run.sh",'w')
shellFileWriter.write('#'+"!/bin/bash\n")
shellFileWriter.write(cmd)
shellFileWriter.close()

cmd = outputDir + "/run.sh"
subprocess.check_output(cmd, shell=True)
print("Program finished.")