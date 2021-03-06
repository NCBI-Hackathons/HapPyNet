'''
Written by Delasa Aghamirzaie for NCBI Hackaton
process SRA files and creates LD counts
'''

from __future__ import print_function
from optparse import OptionParser
import subprocess
import os
import time


#########################
global ALIGN_THREADS
global SORT_THREADS


ALIGN_THREADS = 1
SORT_THREADS = 1


#########################

def main():
    usage = 'usage: %prog [options] <SRA_id or vcf file>. call %prog with -h to get more information'

    parser = OptionParser(usage)
    #.add_option('-h',  usage = 'usage: %prog [options] <SRA_id or vcf file>')
    parser.add_option('-o', dest='out_dir', default='../results/')
    parser.add_option('-i', dest = 'prefix', default=None, type='str', help='Add prefix [Default: %default]')
    parser.add_option('-v', dest = 'vcf', default=False,  action='store_true', help='boolean flag setting input type to VCF: %default]')
    parser.add_option('-s', dest = 'sra', default=False,  action='store_true', help='boolean flag setting input type to SRA : %default]')
    parser.add_option('-c', dest='clean', default=False, action='store_true', help='clean intermediate files: %default]')
    parser.add_option('-b', dest = 'batch', default=False, action='store_true', help = 'generate counts for multiple samples' )
    parser.add_option('-n', dest = 'nationality', default='EUR', type='str', help='nationality can be one of EUR, ASN, AFR [Default: %default]')

    ## TODO: add genome version and species


    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide SRA id with -s or input vcf file using -v option, or -b for a list of samples to run')
    else:
        input_file = args[0] #or sra_file_path
        save_name = os.path.basename(input_file)

    if (not options.sra and not options.vcf):
        parser.error('Must provide input type using -v or -s')

    if not os.path.isdir(options.out_dir):
        os.mkdir(options.out_dir)


    if options.batch:
        samples = open(input_file, 'r')
        samples_file = samples.readlines()
        map(lambda x: runner(x.strip(), options), samples_file)
    else:
        runner(input_file, options)

    ###to do: clean after each run

def runner(input_file, options):

    if options.sra:
        input_type = 'sra'
        vcf_run_stat = False
    elif options.vcf:
        input_type = 'vcf'
        vcf_run_stat = True

    if options.prefix != None:
        save_name = '%s_%s' %(options.prefix, input_file)
    else:
        save_name = input_file

    out_dir = options.out_dir.strip('/') + '/' + save_name + '/'
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
        print('created %s' % out_dir)


    counts_dir = options.out_dir.strip('/') + '/' + 'counts/'
    if not os.path.isdir(counts_dir):
        os.mkdir(counts_dir)

    naming_dic = set_naming_convention(save_name, out_dir, counts_dir)
    naming_dic['input_type'] = input_type
    naming_dic['input_file'] = input_file

    if options.nationality == 'AFR':
        naming_dic['nationality'] = 'AFR'
        naming_dic['ld_blocks'] = '../ref_data/lddetect_GRCh38/AFR_ldetect.bed'
    elif options.nationality == 'ASN':
        naming_dic['nationality'] = 'ASN'
        naming_dic['ld_blocks'] = '../ref_data/lddetect_GRCh38/ASN_ldetect.bed'
    else:
        naming_dic['nationality'] = 'EUR'
        naming_dic['ld_blocks'] = '../ref_data/lddetect_GRCh38/EUR_ldetect.bed'

    if options.vcf:
        naming_dic['sorted_markd_recal_vcf'] = input_file
    #    print (naming_dic)
    check_file_existence(file_path=naming_dic['log_file'], create=True)
    log_file = open(naming_dic['log_file'], 'a+')
    log_file.write('~~~~~~ started new run for : %s ~~~~~~~~\n' % save_name)

    # TODO clean after each successful run
    if options.vcf:
        # TODO
        # process vcf file via bedtools intersect
        vcf_run_stat = True

    else:
        print ('processing SRA file using custom pipeline')

        ### do steps to create SRA file from SRA id
        log_file.write(save_name)
        log_file.write("\n")

        print('save_name: %s' %save_name)
        log_file.write('started running HISAT2 \n')

        start = time.time()
        #### ~~~~~~~~~~~~~~ run HISAT2 ~~~~~~~~~~~~~~~~~~~
        if (not check_file_existence(naming_dic['sorted_bam'])) :
            run_hisat(naming_dic)

        else:
            log_file.write('bam file exsists. skipping alignemnt \n')
            print('bam file exsists. skipping alignemnt')

        if not check_file_existence(naming_dic['sorted_bam'], create=False):
            log_file.write('HISAT2 run was unsuccessfull\n')
            raise NameError('HISAT2 run was  unsuccessfull')
        else:
            log_file.write('HISAT2 finished succesfuly\n')

        end = time.time()
        duration = end - start
        log_file.write('took %s seconds, %s minutes\n\n' % (duration, duration / 60))
        ### ~~~~~~~~~~~~~~ mark_duplicates ~~~~~~~~~~~~~~
        log_file.write('started running GATK mark_duplicates\n')
        start = time.time()
        mark_duplicates(naming_dic)

        if (check_file_existence(naming_dic['sorted_markd_bam']) and (check_file_existence(naming_dic['sorted_marked_metrics']))):
            log_file.write('Marking duplicates finished \n')
        else:
            log_file.write('Marking duplicates was unsuccessfull\n')
            raise NameError('Marking duplicates was unsucessful')

        end = time.time()
        duration = end - start
        log_file.write('took %s seconds, %s minutes\n' % (duration, duration / 60))
        ## ~~~~~~~~~~~~~~ BaseRecalibrator ~~~~~~~~~~~~~~

        start = time.time()
        base_recalibrator(naming_dic)

        if not check_file_existence(file_path=naming_dic['sorted_markd_recal_table']):
            log_file.write('base recalibration was unsuccessfull\n')
            raise NameError('base recalibration was unsucessful')
        end = time.time()
        duration = end - start
        log_file.write('took %s seconds, %s minutes\n' % (duration, duration / 60))
        ## ~~~~~~~~~~~~~~ PrintReads ~~~~~~~~~~~~~~

        start = time.time()
        print_reads(naming_dic)

        if not check_file_existence(file_path=naming_dic['sorted_markd_recal_bam']):
            log_file.write('printing reads was unsuccessfull\n')
            raise NameError('printing reads was unsucessful')
        end = time.time()
        duration = end - start
        log_file.write('took %s seconds, %s minutes' % (duration, duration / 60))
        ## ~~~~~~~~~~~~~~ HaplotypeCaller ~~~~~~~~~~~~~~
        start = time.time()

        haplotype_caller(naming_dic)
        if not check_file_existence(file_path=naming_dic['sorted_markd_recal_vcf']):
            log_file.write('Haplotype calling was unsuccessfull\n')
            raise NameError('Haplotype calling was unsucessful')
        end = time.time()
        duration = end - start
        log_file.write('took %s seconds, %s minutes\n' % (duration, duration / 60))

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~``

        vcf_run_stat = True

    ##~~~~~~~~~~~~~~ Bedtools intersect ~~~~~~~~~~~~
    if vcf_run_stat:
        log_file.write('started running bedtools intersect\n')
        start = time.time()
        ### steps after having a VCF file
        bedtools_intersect(naming_dic)

        if not check_file_existence(file_path=naming_dic['ld_counts']):
            log_file.write('Generating LD counts was unsuccessfull\n')
            raise NameError('Generating LD counts was unsucessful')

        end =  time.time()
        duration = end - start
        log_file.write('took %s seconds, %s minutes' % (duration, duration / 60))
    ###
    log_file.write('~~~~~~~~~~~~~~~~  Preprocessing data finished successfully  ~~~~~~~~~~~~~~~~~~~~~~~~`\n')
    print ('Preprocessing data finished successfully')
    ### ~~~~~~~~~~~~~~~~~

def check_file_existence(file_path, create=False):
    exists = False
    if (os.path.exists(file_path)):
        if (os.stat(file_path).st_size > 0):
            return True
    elif create:
        new_file = open(file_path, "w")
        print ('created a new empty file just for test')
        new_file.close()
        return True
    else:
        return False

def set_naming_convention(save_name, out_dir, counts_dir):
    '''
    :param save_name:
    :return: a dictionary containing the expected file names in this pipeline if SRA file is inputted
    '''

    naming_dic = {
        'save_name': save_name,
        'out_dir': out_dir,
        'log_file': '%s%s.log' %(out_dir,save_name),
        'sorted_bam':'%s%s.sort.bam' %(out_dir,save_name),
        'sorted_markd_bam': '%s%s.sort.markd.bam' %(out_dir,save_name),
        'sorted_marked_metrics': '%s%s.sort.markd.metrics.txt' %(out_dir, save_name),
        'sorted_markd_recal_table':  '%s%s.sort.markd.recal.table' %(out_dir,save_name),
        'sorted_markd_recal_bam': '%s%s.sort.markd.recal.bam' %(out_dir,save_name),
        'sorted_markd_recal_vcf': '%s%s.sort.markd.recal.vcf.gz' %(out_dir,save_name),
        'ld_counts': '%s%s.ld_counts' %(counts_dir, save_name)

    }

    return (naming_dic)

def run_hisat(naming_dic):
        tmp_dir = naming_dic['out_dir'] + naming_dic['save_name'] +  '_tmp'
        os.mkdir(tmp_dir)
        #with tempfile.NamedTemporaryFile(dir=naming_dic['out_dir'], prefix=naming_dic['save_name'], suffix='tmp') as t:
        hisat_cmd = 'hisat2 -p %s -x /opt/grch38/Homo_sapiens_assembly38 --sra-acc %s --rg-id %s --rg SM:%s --rg PL:ILLUMINA  | samtools sort -@ %s > %s -T %s' %(ALIGN_THREADS, naming_dic['input_file'], naming_dic['input_file'], SORT_THREADS,naming_dic['input_file'], naming_dic['sorted_bam'], tmp_dir)
        print (hisat_cmd)
        subprocess.call(hisat_cmd, shell=True)
        os.rmdir(tmp_dir)


def mark_duplicates(naming_dic):
        print ('\n#### mark_duplicates ####')
        #/opt/gatk-4.0.3.0/gatk MarkDuplicates --INPUT ${SRR}.sort.bam --OUTPUT ${SRR}.sort.markd.bam --METRICS_FILE ${SRR}.sort.markd.metrics.bam
        cmd = '/opt/gatk-4.0.3.0/gatk MarkDuplicates --INPUT %s --OUTPUT %s --METRICS_FILE %s' %(naming_dic['sorted_bam'], naming_dic['sorted_markd_bam'], naming_dic['sorted_marked_metrics'])
        print (cmd)
        subprocess.call(cmd, shell=True)

def base_recalibrator(naming_dic):
    ## BaseRecalibrator
    #/opt/gatk-4.0.3.0/gatk BaseRecalibrator -R /opt/grch38/Homo_sapiens_assembly38.fasta --known-sites /opt/grch38/dbsnp_138.hg38.vcf.gz --known-sites /opt/grch38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz -I ${SRR}.sort.markd.bam -O ${SRR}.sort.markd.recal.table
    print('\n#### base_recalibrator ####')
    cmd = '/opt/gatk-4.0.3.0/gatk BaseRecalibrator -R /opt/grch38/Homo_sapiens_assembly38.fasta --known-sites /opt/grch38/dbsnp_138.hg38.vcf.gz --known-sites /opt/grch38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz -I %s -O %s' %(naming_dic['sorted_markd_bam'], naming_dic['sorted_markd_recal_table'])
    print (cmd)
    subprocess.call(cmd, shell=True)

def print_reads(naming_dic):
    print ('\n#### print_reads ####')
    #/opt/gatk-4.0.3.0/gatk ApplyBQSR -R /opt/grch38/Homo_sapiens_assembly38.fasta -I ${SRR}.sort.markd.bam --bqsr-recal-file ${SRR}.sort.markd.recal.table -O ${SRR}.sort.markd.recal.bam
    cmd = '/opt/gatk-4.0.3.0/gatk ApplyBQSR -R /opt/grch38/Homo_sapiens_assembly38.fasta -I %s --bqsr-recal-file %s -O %s'%(naming_dic['sorted_markd_bam'],naming_dic['sorted_markd_recal_table'],naming_dic['sorted_markd_recal_bam'])
    print (cmd)
    subprocess.call(cmd, shell=True)

def haplotype_caller(naming_dic):
    print('\n#### haplotype_caller ####')
    #/opt/gatk-4.0.3.0/gatk HaplotypeCaller -R /opt/grch38/Homo_sapiens_assembly38.fasta --dbsnp /opt/grch38/dbsnp_138.hg38.vcf.gz -I ${SRR}.sort.markd.recal.bam -O ${SRR}.sort.markd.recal.vcf.gz
    cmd = '/opt/gatk-4.0.3.0/gatk HaplotypeCaller -R /opt/grch38/Homo_sapiens_assembly38.fasta --dbsnp /opt/grch38/dbsnp_138.hg38.vcf.gz -I %s -O %s' %(naming_dic['sorted_markd_recal_bam'],naming_dic['sorted_markd_recal_vcf'])
    print (cmd)
    subprocess.call(cmd, shell=True)

def bedtools_intersect(naming_dic):
    print('\n#### bedtools_intersect ####')
    #'/home/ubuntu/bin/bedtools intersect -a /home/ubuntu/ldetect_GRCh38/EUR_ldetect.bed -b ${SRR}.sort.markd.recal.vcf.gz -c | sort -k1,1V -k2,2n > ${SRR}.count'
    cmd = '/home/ubuntu/bin/bedtools intersect -a %s -b %s -c | sort -k1,1V -k2,2n > %s' %(naming_dic['ld_blocks'],naming_dic['sorted_markd_recal_vcf'], naming_dic['ld_counts'])
    print (cmd)
    subprocess.call(cmd, shell=True)

# __main__
################################################################################

if __name__ == '__main__':
    main()

## for parallel running:
#parallel -j8 'echo {}' ::: $(seq 1 10)


#####


####