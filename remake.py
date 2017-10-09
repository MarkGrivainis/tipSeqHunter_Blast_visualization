# TODO

import os
import re
import xml.etree.ElementTree as ET
from collections import defaultdict

import pandas as pd
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.SeqIO import parse

TAGS = ['Hsp_align-len', 'Hsp_gaps', 'Hsp_query-from',
        'Hsp_query-to', 'Hsp_query-from', 'Hsp_query-to', 'Hsp_evalue', 'Hsp_bit-score', 'Hsp_qseq', 'Hsp_midline', 'Hsp_hseq', 'Hsp_hit-to', 'Hsp_hit-from']
rev_nucleotides = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N', '-': '-'}

FLANKING = 10
sample_type = 'primary'


def reverse_comp(seq):
    new_seq = ''.join([rev_nucleotides[x] for x in seq])[::-1]
    return new_seq


# df_repred = pd.read_table('simu.frag-300.len-100.cov-100.error-10.snps-5.rpt-1.wsize100.regwsize1.minreads1.clip1.clipflk5.mindis150.FP.uniqgs.bed.csinfo.lm.l1hs.pred.txt.repred')
if sample_type == 'normal':
    df_repred = pd.read_table('A152-Normal.s5_9_SL31275.wsize100.regwsize1.minreads1.clip1.clipflk5.mindis150.FP.uniqgs.bed.csinfo.lm.l1hs.pred.txt.repred')
elif sample_type == 'met':
    df_repred = pd.read_table('A152-Met.s5_11_SL31277.wsize100.regwsize1.minreads1.clip1.clipflk5.mindis150.FP.uniqgs.bed.csinfo.lm.l1hs.pred.txt.repred')
elif sample_type == 'primary':
    df_repred = pd.read_table('A152-Primary.s5_10_SL31276.wsize100.regwsize1.minreads1.clip1.clipflk5.mindis150.FP.uniqgs.bed.csinfo.lm.l1hs.pred.txt.repred')
else:
    df_repred = pd.read_table(
        'simu.frag-300.len-100.cov-100.error-10.snps-5.rpt-1.wsize100.regwsize1.minreads1.clip1.clipflk5.mindis150.FP.uniqgs.bed.csinfo.lm.l1hs.pred.txt.repred')
df_200FP = pd.read_table('ucsc.rm327.rmsk.L1Hs-Ta.with-strand.matchFP.bed.l1prm', dtype=str)

output_json = {}
output_json['hit'] = defaultdict(list)
output_json['miss'] = defaultdict(list)

# current_chrom = ''
results = {}
# with open('simu.frag-300.len-100.cov-100.error-10.snps-5.rpt-1.fa', 'rU') as handle:
xml_files = os.listdir('./xml_' + sample_type + '/')
if sample_type == 'normal':
    filename = 'k55-A152-Normal.s5_9_SL31275.fa'
elif sample_type == 'met':
    filename = 'k55-A152-Met.s5_11_SL31277.fa'
elif sample_type == 'primary':
    filename = 'k55-A152-Primary.s5_10_SL31276.fa'
else:
    filename = 'simu.frag-300.len-100.cov-100.error-10.snps-5.rpt-1.fa'
with open(filename, 'rU') as handle:
    hg_19 = parse(handle, 'fasta')
    for each in hg_19:
        if True:
        # if each.name in ['chr8-48240865-48244565_NODE_1_length_2657_cov_195.635']:
            match = re.match('^(.+?)-(\d+?)-(\d+?)_', each.name)
            match_df = df_repred.loc[(df_repred['H5_TargChr'] == match.group(1)) & (df_repred['H6_TargS'] == int(match.group(2))) & (df_repred['H7_TargE'] == int(match.group(3)))]
            if not str(match_df.iloc[0]['H17_IfGS']) == 'nan':
                # print(each.name)
                results[each.name] = []
                gs = match_df.iloc[0]['H17_IfGS']
                GS_from_repred = dict(zip(['chrm', 'L1S', 'L1E'], match_df.iloc[0]['H17_IfGS'].split('-')))
                TargS = match_df.iloc[0]['H6_TargS']
                TargE = match_df.iloc[0]['H7_TargE']
                # print('TargS', TargS, 'TargE', TargE)
                GS_row = df_200FP.loc[
                    (df_200FP['chrNum'] == GS_from_repred['chrm']) & \
                    (df_200FP['L1Start'] == GS_from_repred['L1S']) & \
                    (df_200FP['L1End'] == GS_from_repred['L1E'])
                    ]
                if GS_row.iloc[0]['L1Strand'] == '+':
                    start = int(GS_row.iloc[0]['L1PrimerStart'])
                    end = TargE
                    taat_position = int(GS_row.iloc[0]['EndBeforePolyA']) - int(GS_row.iloc[0]['L1PrimerStart'])
                    junc_position = int(GS_row.iloc[0]['L1End']) - int(GS_row.iloc[0]['L1PrimerStart'])
                    junc_offset = TargS - start
                elif GS_row.iloc[0]['L1Strand'] == '-':
                    start = TargS
                    end = int(GS_row.iloc[0]['EndBeforePolyA'])
                    taat_position = int(GS_row.iloc[0]['L1PrimerStart']) - TargS
                    junc_position = int(GS_row.iloc[0]['L1Start']) - TargS
                    junc_offset = end - TargE
                else:
                    raise RuntimeError('Strand is neither + or -')
                if each.name + '.xml' not in xml_files:
                    with open('temp_sim.fa', 'w') as sim_fa:
                        sim_fa.write('>{}\n{}'.format(each.name, each.seq))
                    with open(os.path.expanduser('~/PycharmProjects/transposcope_backend/reference/chromFa/{}.fa'.format(GS_from_repred['chrm'])), 'rU') as handle:
                        ref = parse(handle, 'fasta')
                        with open('temp_ref.fa', 'w') as ref_fa:
                            # print(GS_row)
                            # print(GS_row.iloc[0]['L1Strand'])
                            theSeq = next(ref)
                            seq = theSeq.seq[start: end]
                            ref_fa.write('>{}\n{}'.format(match_df.iloc[0]['H17_IfGS'], seq))

                    blastn_cline = NcbiblastnCommandline(task="megablast",
                                                         query="temp_ref.fa",
                                                         subject="temp_sim.fa", max_target_seqs=100, word_size=28,
                                                         penalty=-2,
                                                         reward=1,
                                                         outfmt=5)
                    stdout, stderr = blastn_cline()
                    with open('xml_' + sample_type + '/' + each.name + '.xml', 'w') as fh:
                        fh.write(stdout)
                    root = ET.fromstring(stdout)
                else:
                    with open('xml_' + sample_type + '/' + each.name + '.xml', 'r') as fh:
                        root = ET.fromstring(fh.read())
                for iteration in root.iter('BlastOutput_iterations'):

                    for i in iteration:
                        a = i.find('Iteration_hits')
                        if len(a) > 0:
                            b = a.find('Hit')
                            # print(b.find('Iteration_query-def').text)
                            c = b.find('Hit_hsps')
                            for d in c.findall('Hsp'):
                                hit = {
                                    'taat_position': taat_position,
                                    'junc_position': junc_position,
                                    'strand': GS_row.iloc[0]['L1Strand'],
                                    'length': end - start,
                                    'GS': gs
                                }
                                for tag in TAGS:
                                    e = d.find(tag)
                                    hit[tag] = e.text
                                results[each.name].append(hit)

for k, v in results.items():
    print(k, len(v))
    for each in v:
        blast_start = int(each['Hsp_query-from']) if int(each['Hsp_query-from']) <= int(each['Hsp_query-to']) else each['length'] - int(each['Hsp_query-from'])
        blast_end = int(each['Hsp_query-to']) if int(each['Hsp_query-from']) <= int(each['Hsp_query-to']) else each['length'] - int(each['Hsp_query-to'])
        print(each['length'])
        hit_or_miss = 'miss'
        if each['strand'] == '+' and int(each['Hsp_query-from']) <= int(each['Hsp_query-to']):
            print('+ | +')
            print(blast_start, each['taat_position'], each['junc_position'], blast_end)
            start = -1
            end = -1
            if (blast_start <= each['taat_position'] and each['junc_position'] <= blast_end) and (
                            blast_start <= each['taat_position'] - FLANKING or each['junc_position'] + FLANKING <=
                        blast_end):
                hit_or_miss = 'hit'
                start = each['taat_position'] - blast_start + 1
                end = each['junc_position'] - blast_start + 1
                start_offset_correction = each['Hsp_qseq'][0:start - 1].count('-')
                # poly_a_offset_correction = start_offset_correction + each['Hsp_qseq'][start:end + start_offset_correction].count('-')
                poly_a_offset_correction = 0
                poly_a_offset = 0
                current_char = start + start_offset_correction
                print(each['Hsp_qseq'])
                print(each['Hsp_midline'])
                print(each['Hsp_hseq'])
                # add one so that it counts the first character at the taat or atta
                while poly_a_offset < (end + start_offset_correction) - start + 1:
                    if current_char == len(each['Hsp_qseq']):
                        break
                    if each['Hsp_qseq'][current_char] != '-':
                        poly_a_offset += 1
                    else:
                        poly_a_offset_correction += 1
                    current_char += 1
                poly_a_offset_correction += start_offset_correction

        if each['strand'] == '+' and int(each['Hsp_query-from']) >= int(each['Hsp_query-to']):
            print('+ | -')
            # print(blast_start, each['taat_position'], each['junc_position'], blast_end)
            # if blast_start <= each['taat_position'] and each['junc_position'] <= blast_end:
            #     hit.add(k)
            #
            # else:
            #     miss.add(k)
            # start = 0
            # end = 0
            raise AssertionError('strand has been reverse complemented')
        if each['strand'] == '-' and int(each['Hsp_query-from']) <= int(each['Hsp_query-to']):
            print('- | +')
            print(blast_start, each['junc_position'], each['taat_position'], blast_end)
            start = -1
            end = -1
            if (blast_start <= each['junc_position'] and each['taat_position'] <= blast_end) and (
                            blast_start <= each['junc_position'] - FLANKING or each['taat_position'] + FLANKING <=
                        blast_end):
                hit_or_miss = 'hit'
                start = each['junc_position'] - blast_start + 1
                end = each['taat_position'] - blast_start + 1
                start_offset_correction = each['Hsp_qseq'][0:start - 1].count('-')
                # poly_a_offset_correction = start_offset_correction + each['Hsp_qseq'][start:end + start_offset_correction].count('-')
                poly_a_offset_correction = 0
                poly_a_offset = 0
                current_char = start + start_offset_correction
                print(each['Hsp_qseq'])
                print(each['Hsp_midline'])
                print(each['Hsp_hseq'])
                # add one so that it counts the first character at the taat or atta
                while poly_a_offset < (end + start_offset_correction) - start + 1:
                    if current_char == len(each['Hsp_qseq']):
                        break
                    if each['Hsp_qseq'][current_char] != '-':
                        poly_a_offset += 1
                    else:
                        poly_a_offset_correction += 1
                    current_char += 1
                poly_a_offset_correction += start_offset_correction
        if each['strand'] == '-' and int(each['Hsp_query-from']) >= int(each['Hsp_query-to']):
            print('+ | -')
            # print(blast_start, each['taat_position'], each['junc_position'], blast_end)
            # if blast_start <= each['taat_position'] and each['junc_position'] <= blast_end:
            #     hit.add(k)
            #
            # else:
            #     miss.add(k)
            # start = 0
            # end = 0
            raise AssertionError('strand has been reverse complemented')
        if start != -1 and end != -1:
            # -1 so that the last t is not included as part of the poly a tail


            end_flanking_offset = 0
            end_flanking_offset_correction = 0
            while end_flanking_offset < FLANKING:
                if current_char == len(each['Hsp_qseq']):
                    break
                if each['Hsp_qseq'][current_char] != '-':
                    end_flanking_offset += 1
                else:
                    end_flanking_offset_correction += 1
                current_char += 1

            pad_sequence = ''
            extend_poly_a = 0
            hseq = each['Hsp_hseq'][start + start_offset_correction: end + poly_a_offset_correction]
            qseq = each['Hsp_qseq'][start + start_offset_correction: end + poly_a_offset_correction]
            if -10 <= start - FLANKING + start_offset_correction <= 0:
                pad_sequence = '*' * (-1 * (start - FLANKING + start_offset_correction))
                extend_poly_a = -1 * (start - FLANKING + start_offset_correction)
                start = 10
                print('needs correction', start - FLANKING + start_offset_correction, pad_sequence)
            print(each['Hsp_qseq'])
            print(each['Hsp_midline'])
            print(each['Hsp_hseq'])

            print('.' * (start - FLANKING + start_offset_correction) + '*' * FLANKING + 'A' * (extend_poly_a + end + poly_a_offset_correction - (start + start_offset_correction))
                  + '*' *
                  FLANKING)
            print('.' * (start - FLANKING + start_offset_correction) + each['Hsp_qseq'][start - FLANKING + start_offset_correction: end + FLANKING + poly_a_offset_correction])
            print('.' * (start - FLANKING + start_offset_correction) + each['Hsp_midline'][start - FLANKING + start_offset_correction: end + FLANKING + poly_a_offset_correction])
            print('.' * (start - FLANKING + start_offset_correction) + each['Hsp_hseq'][start - FLANKING + start_offset_correction: end + FLANKING + poly_a_offset_correction])
            # print('.' * (start-FLANKING + offset_correction) + '*' * FLANKING + 'A' * (end - start) +'*' * 10)

            print('*' * FLANKING + 'A' * (extend_poly_a + end + poly_a_offset_correction - (start + start_offset_correction)) + '*' * FLANKING)
            print('*' * FLANKING + 'A' * (end - (start)) + '*' * FLANKING)
            print(pad_sequence + each['Hsp_qseq'][start - FLANKING + start_offset_correction: end + FLANKING + poly_a_offset_correction])
            print(pad_sequence + each['Hsp_midline'][start - FLANKING + start_offset_correction: end + FLANKING + poly_a_offset_correction])
            print(pad_sequence + each['Hsp_hseq'][start - FLANKING + start_offset_correction: end + FLANKING + poly_a_offset_correction])
            # if output_json[hit_or_miss][k].get
            output_json[hit_or_miss][k].append({
                'ref_strand': each['strand'],
                'query-f': int(each['Hsp_query-from']),
                'query-t': int(each['Hsp_query-to']),
                'ref_len': each['length'],
                'query-len': len(each['Hsp_qseq']),
                'hit-t': int(each['Hsp_hit-to']),
                'hit-f': int(each['Hsp_hit-from']),
                'poly_a_template': '*' * FLANKING + 'A' * (extend_poly_a + end + poly_a_offset_correction - (start + start_offset_correction)) + '*' * FLANKING,
                'poly_a_ref': pad_sequence + each['Hsp_qseq'][start - FLANKING + start_offset_correction: end + FLANKING + poly_a_offset_correction],
                'poly_a_midline': pad_sequence + each['Hsp_midline'][start - FLANKING + start_offset_correction: end + FLANKING + poly_a_offset_correction],
                'poly_a_hit': pad_sequence + each['Hsp_hseq'][start - FLANKING + start_offset_correction: end + FLANKING + poly_a_offset_correction],
                'hit_poly_a_length': len(hseq) - hseq.count('-'),
                'query_poly_a_length': len(qseq) - qseq.count('-'),
                'seq_ref': each['Hsp_qseq'],
                'seq_midline': each['Hsp_midline'],
                'seq_query': each['Hsp_hseq'],
                # 'expected': expected,
                # 'expected_template': each['template'],
                'GS': each['GS']
            })
        else:
            output_json[hit_or_miss][k].append('a')

fh = open('output_' + sample_type + '.csv', 'w')
fh.write('{},{},{},{},{}\n'.format('id', 'patient_length', 'reference_length', 'type', 'GS'))
print('hits', len(output_json['hit']))
for k, value in output_json['hit'].items():
    print(k, len(value))
    hit = max(value, key=lambda x: x['query-len'])
    # for hit in value:
    print()
    # print('\texpected-\t:', hit['expected'])
    # print('\texp_temp-\t:', hit['expected_template'])
    print('\tmatchTemp\t:', hit['poly_a_template'])
    print('\treference\t:', hit['poly_a_ref'])
    print('\tmidline--\t:', hit['poly_a_midline'])
    print('\thit------\t:', hit['poly_a_hit'])
    print('\tmatchTemp\t:', hit['poly_a_template'])
    print('\thit len--\t:', hit['hit_poly_a_length'])
    print('\tquery len\t:', hit['query_poly_a_length'])
    fh.write('{},{},{},{},{}\n'.format(k, hit['hit_poly_a_length'], hit['query_poly_a_length'], sample_type, hit['GS']))
# for k, values in output_json['miss'].items():
#     print(k)
#     for hit in values:
#         print()
#         print('\texpected-\t:', hit['expected'])
#         print('\texp_temp-\t:', hit['expected_template'])
#         print('\treference\t:', hit['poly_a_ref'])
#         print('\tmidline--\t:', hit['poly_a_midline'])
#         print('\tquery----\t:', hit['poly_a_query'])
#         print('\tmatchTemp\t:', hit['poly_a_template'])
print('hits', len(output_json['hit']))
print('misses', len(output_json['miss']))
print(len(set(output_json['miss']).difference(set(output_json['hit']))))
miss_fh = open('misses_' + sample_type + '.txt', 'w')
for each in set(output_json['miss']).difference(set(output_json['hit'])):
    miss_fh.write('{}\n'.format(each))
miss_fh.close()
fh.close()
