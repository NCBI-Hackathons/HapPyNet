from __future__ import print_function
from optparse import OptionParser
import subprocess
import os

def main():
    usage = 'usage: %prog [options] <SRA_id or vcf file>'

    parser = OptionParser(usage)
    parser.add_option('-o', dest='out_dir', default='./')
    parser.add_option('-i', dest = 'prefix',  default=None, type='str', help='Add prefix [Default: %default]')
    parser.add_option('-v', dest = 'vcf', default=None,  type='str', help='VCF input file: %default]')

    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide SRA id or input vcf file using -v option')
    else:
        sra_file = args[0] #or sra_file_path
        save_name = os.path.basename(sra_file)

    if not os.path.isdir(options.out_dir):
        os.mkdir(options.out_dir)

    if options.prefix != None:
        save_name = '%s_%s' %(options.prefix, save_name)

    #check log file:
    print (save_name)
    log_file_name = '%s%s.log' %(options.out_dir, save_name)

    log_file = open(log_file_name, "a+")
    log_file.write('~~~~~~ started new run for : %s ~~~~~~~~\n' % save_name)
    if options.vcf != None:
        # TODO
        #process vcf file
        pass
    else:
        log_file.write(save_name)
        log_file.write("\n")

        print (sra_file, save_name)
        log_file.write('started running HISAT2 \n')
        bam_file = '%s%s.sort.bam' % (options.out_dir, save_name)


        if not os.path.exists(bam_file):
            run_hisat(sra_file, options.out_dir)

        if not os.path.exists(bam_file):
            log_file.write('HISAT2 run was unsuccessfull\n')
            raise NameError('HISAT2 run was  unsuccessfull')
        else:
            log_file.write('HISAT2 finished succesfuly\n')

        log_file.write('started running GATK \n')
        run_GATK(bam_file, options.out_dir)
        log_file.write('finished running GATK \n')

    ###to do: clean after each run


def run_hisat(sra_path, save_name  ='', out_dir = './'):

        cmd = 'hisat2 -p 8 -x /opt/grch38/genome --sra-acc %s | samtools sort -@ 4 > %s%s.sort.bam' %(sra_path, out_dir, save_name)
        print (cmd)
     #   subprocess.call(cmd, shell=True)

def run_GATK(save_name = '', out_dir = './'):
        cmd = '/opt/gatk-4.0.3.0/gatk MarkDuplicates --INPUT %s%s.sort.bam --OUTPUT %s%s.sort.markd.bam --METRICS_FILE %s%s.markd.metrics.bam' %(out_dir, save_name, out_dir, save_name, out_dir, save_name)
        print (cmd)
     #   subprocess.call(cmd, shell=True)
# __main__
################################################################################
if __name__ == '__main__':
    main()


## for parallel running:
