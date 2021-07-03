import os
import subprocess
import time


# go through hmm results, and for every pfam, make fasta for all genes that hit for that pfam
def hmm_fastas(result, dir1, proteome_address):
    # dictionary to store all gene results from pfam hmmsearch
    result_pfams = {}

    with open(result, 'r') as file1:
        prev = ''  # pfam of previous row
        for line1 in file1:
            first = line1[:2]
            # if this is true, line contains data
            if first == 'sp' or first == 'tr':
                result = line1.split()
                gene = result[0]
                pfam1 = result[4]
                filename1 = dir1 + pfam1 + '.fasta'
                if pfam1 == prev:
                    if not (result_pfams[pfam1].__contains__(gene)):
                        result_pfams.setdefault(pfam1, []).append(gene)
                        # write a to file
                        get_fasta(gene, proteome_address, filename1, 'a')

                else:
                    if not (pfam1 in result_pfams):
                        result_pfams.setdefault(pfam1, []).append(gene)
                        prev = pfam1
                        # current is a new pfam so write w to file
                        get_fasta(gene, proteome_address, filename1, 'w')

    return result_pfams


# uses shell commands to get fasta for gene
def get_fasta(gene_name, proteome_address, file_name, ch):
    arg1 = '/projects/MOLBIO/local/src/hmmer-3.1b2/easel/miniapps/esl-sfetch --index ' + proteome_address
    subprocess.call(arg1, shell=True)

    arg1 = '/projects/MOLBIO/local/src/hmmer-3.1b2/easel/miniapps/esl-sfetch ' + proteome_address + ' "' + gene_name + '" >'

    # if ch == a, then append to existing file
    if ch == 'a':
        arg1 += '> ' + file_name
        subprocess.call(arg1, shell=True)
    # ch == w, so rewrite existing file
    else:
        arg1 += ' ' + file_name
        subprocess.call(arg1, shell=True)


# Go through organism's pfam results to tabulate and make included list and excluded list --> determine which pfams
# to go through by going through accession list (to make sure we dont miss any pfams)
# Want to go through blast result, count up number of hits in blast and figure out which hmm results were not
# included in blast result
def tabulate(ac_list1, hmm_pfams, proteome_address):
    row = organism.ljust(15, ' ')

    with open(ac_list1, 'r') as doc:
        # each line = pfam name
        for pfam1 in doc:
            pfam1 = pfam1.strip('\n')
            pfam_result = '/tigress/takegami/pfam_hmmblast_biogen/blast_results/' + organism + '/' + pfam1 + '.txt'
            print(pfam_result)
            blast_list = []  # list of genes that hit with blast for pfam

            subdir1 = '/tigress/takegami/pfam_hmmblast_biogen/blast_results/' + organism + '/'
            file_name = pfam1 + '.txt'
            print(file_name)
            cwd = os.getcwd()
            print('location: ' + cwd)

            # read files only if pfam files exist in blast_results
            if os.listdir(subdir1).__contains__(file_name):
                # reading in blast result line by line
                with open(pfam_result, 'r') as file1:
                    cwd = os.getcwd()
                    print('location: should be blast result file    ' + cwd)
                    for line1 in file1:
                        first = line1[:2]
                        # if this is true, line contains data
                        if first == 'sp' or first == 'tr':
                            words = line1.split()
                            name = words[1]
                            if not blast_list.__contains__(name):
                                blast_list.append(name)

            # go through pfam genes, if not contained in blast_list get fasta and add to continous text file (fasta)
            cwd = os.getcwd()
            print('location:    ' + cwd)
            if pfam1 in hmm_pfams:
                hmm_genes = hmm_pfams[pfam1]
                for i in hmm_genes:
                    if not blast_list.__contains__(i):
                        filename1 = '/tigress/takegami/pfam_hmmblast_biogen/excluded.fasta'
                        get_fasta(i, proteome_address, filename1, 'a')

                data = str(len(blast_list)) + '/' + str(len(hmm_pfams[pfam1]))
                row += data.ljust(15, ' ')
            else:
                data = '0/0'
                row += data.ljust(15, ' ')

    return row


if __name__ == '__main__':
    start_time = time.time()
    # go through each organism and do HMM and BLAST
    directory = '/tigress/MOLBIO/databases/uniprot/Bacteria/'
    organism_counter = 0

    filename = '/tigress/takegami/pfam_hmmblast_biogen/excluded.fasta'
    f = open(filename, 'w')
    f.write('')
    f.close()

    for organism in os.listdir(directory):

        if not(organism == '.DS_Store') and not ('ssi' in organism):
            organism_counter += 1
            print('-------------organism: ' + organism + '-------------')

            # file name in directory for the proteome
            subdir = directory + organism + '/'
            organism_name = ''
            for filename in os.listdir(subdir):
                suffix = filename[-5:]
                if not ('ssi' in organism) and suffix == 'fasta' and not ('DNA' in filename) and not('additional' in organism):
                    organism_name = filename    # organism_name contains .fasta
                    break

            # address of where proteome fasta is stored
            proteome = directory + organism + '/' + organism_name

            # HMMer to get all gene hits and fastas
            hmm_result = "/tigress/takegami/pfam_hmmblast_biogen/hmm_results/" + organism + '_hmm.txt'
            arg = "hmmsearch --noali --domtblout " + hmm_result + " /tigress/takegami/pfam_hmmblast_biogen/fes_biogen.hmm " + proteome
            subprocess.call(arg, shell=True)

            # make directory to store hmm gene fastas
            arg = "mkdir /tigress/takegami/pfam_hmmblast_biogen/hmm_gene_fastas/" + organism + '/'
            subprocess.call(arg, shell=True)
            subdir = "/tigress/takegami/pfam_hmmblast_biogen/hmm_gene_fastas/" + organism + '/'
            hmm_pfam = hmm_fastas(hmm_result, subdir, proteome)

            # put results into blast_results folder under organism subfolder
            arg = "mkdir /tigress/takegami/pfam_hmmblast_biogen/blast_results/" + organism + "/"
            subprocess.call(arg, shell=True)

            # blast each pfam fasta against query gene fasta by going through all pfam fastas
            pfamdir = '/tigress/takegami/pfam_hmmblast_biogen/hmm_gene_fastas/' + organism + '/'
            for pfam in os.listdir(pfamdir):
                pfam_name = pfam[:-6]
                gene_file = "/tigress/takegami/pfam_hmmblast_biogen/query_genes.fasta"
                subject_file = pfamdir + pfam
                output_file = "/tigress/takegami/pfam_hmmblast_biogen/blast_results/" + organism + "/" + pfam_name + ".txt"
                print('subject file is ' + subject_file)
                print('output file is ' + output_file)



                arg = "blastp -query " + gene_file + " -subject " + subject_file + ' -evalue 1.0e-15 -outfmt "7 qacc sacc evalue" > ' + output_file
                subprocess.call(arg, shell=True)

            # tabulate results into master matrix
            ac_list = '/tigress/takegami/pfam_hmmblast_biogen/ac_list_biogen.txt'
            master_file = '/tigress/takegami/pfam_hmmblast_biogen/master_matrix.txt'

            # first row of master matrix (headers of the columns) if necessary
            if organism_counter == 1:
                header = "organism".ljust(15, ' ')
                f = open(master_file, 'w')
                with open(ac_list, 'r') as file2:
                    for line in file2:
                        line = line.strip('\n')
                        header += str(line.ljust(15, ' '))
                    header += '\n'
                f.write(header)
                f.close()

            # add line for organism summary
            master_row = tabulate(ac_list, hmm_pfam, proteome)
            print('row = ' + master_row)
            f = open(master_file, 'a')
            f.write(master_row)
            f.write('\n')
            f.close()

            print('number of organisms run: ' + str(organism_counter))
            print('\n')

    print("%s seconds total" % (time.time() - start_time))