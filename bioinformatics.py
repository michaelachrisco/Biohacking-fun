from Bio import SeqIO


# if __name__ == '__main__':
#     # idx = 0
#     # for sequence in SeqIO.parse('ls_orchid.fasta','fasta'):
#     #     print(sequence.id)
#     #     print(sequence.seq)
#     #     print('No of nucleotides: {}'.format(len(sequence)))
#     #     idx+=1

#     # print('Total Sequence found {}'.format(idx))

if __name__ == '__main__':
    idx = 0
    for sequence in SeqIO.parse('covid_sample.fasta','fasta'):
        print(sequence.id)
        print(sequence.seq)
        print('No of nucleotides: {}'.format(len(sequence)))
        idx+=1

    print('Total Sequence found {}'.format(idx))