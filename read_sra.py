from __future__ import print_function
from optparse import OptionParser
import subprocess
import os

def main():
    usage = 'usage: %prog [options] <SRA file id>'

    parser = OptionParser(usage)
    parser.add_option('-o', dest='out_dir', default='./')
    parser.add_option('-i', dest = 'prefix',  default=None, type='str', help='Add prefix [Default: %default]')

    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide SRA number')
    else:
        sra_file = args[0]

    if not os.path.isdir(options.out_dir):
        os.mkdir(options.out_dir)

    if options.prefix != None:
       sra_file = '%s_%s' %(options.prefix, sra_file)

    print (sra_file)
    run_hisat(sra_file, options.out_dir)
    bam_file = '%s%s.sort.bam' %( options.out_dir, sra_file)

    if not os.path.exists(bame_file):
        print ('hisat2 run was probably unsuccessfull')
        raise NameError('hisat2 run was probably unsuccessfull')
    run_GATK(bam_file, options.out_dir)

    ###to do: clean after each run



def run_hisat(sra_file='', out_path = './'):
        cmd = 'hisat2 -p 8 -x /opt/grch38/genome --sra-acc %s | samtools sort -@ 4 > %s%s.sort.bam' %(sra_file, out_path, sra_file)
        print (cmd)
        subprocess.call(cmd, shell=True)

def run_GATK(sra_file, out_path):
        cmd = '/opt/gatk-4.0.3.0/gatk MarkDuplicates --INPUT %s.sort.bam --OUTPUT %s.sort.markd.bam --METRICS_FILE %s.markd.metrics.bam' %(sra_file, sra_file, sra_file)
        print (cmd)
        subprocess.call(cmd, shell=True)
# __main__
################################################################################
if __name__ == '__main__':
    main()

