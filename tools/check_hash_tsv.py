#!/usr/bin/env python3

import sys
import argparse


def pad_well_col(well_col, zero_pad, id_length):
  if zero_pad:
      template = '%%0%sd' % id_length
  else:
      template = '%s'
  col_id = template % (well_col)
  return col_id


def index_to_well(well_index, across_row_first):
  if(well_index < 0):
    return((0, 'none'))
  nrow = 8
  ncol = 12
  ipl = int( well_index / 96 )
  i96 = well_index - ipl * 96
  if across_row_first:
      well_row = chr(65 + int(i96 / ncol))
      well_col = (i96 % ncol) + 1
  else:
      well_row = chr(65 + (i96 % nrow))
      well_col = int(i96 / nrow) + 1

  well_id = '%s%s' % ( well_row, pad_well_col( well_col, True, 2 ) )

  return( (ipl, well_id ) )


def base_convert(i, b):
  if(i == 0):
    return([0])
  result = []
  while i > 0:
    result.insert(0, i % b)
    i = i // b
  return result


def index_base_4_to_encoded_index(index_base_4):
  encoder_base = ['A', 'C', 'G', 'T']
  encoded_index_list = list()
  for i in index_base_4:
    encoded_index_list.append(encoder_base[i])
  encoded_index = ''.join(encoded_index_list).rjust(7, 'A')
  return(encoded_index)


imax = 768

if __name__ == '__main__':

  parser = argparse.ArgumentParser(description='A program to write a well-index map.')
  parser.add_argument('-i', '--input', required=True, help='Input hash read TSV filename.')
  args = parser.parse_args()

  input_tsv = args.input

  across_row_first = 1
  well_to_encoded_index_true_dict = dict()
  well_to_encoded_index_true_dict['none'] = 'AAAAAAA'
  for i in range(imax):
    (ipl, well_id) = index_to_well(i, across_row_first)
    well_string = 'P%02d-%s' % (ipl+1, well_id)
    index_base_4 = base_convert(i+1, 4)
    encoded_index = index_base_4_to_encoded_index(index_base_4)
    well_to_encoded_index_true_dict[well_string] = encoded_index

  across_row_first = 0
  well_to_encoded_index_false_dict = dict()
  well_to_encoded_index_false_dict['none'] = 'AAAAAAA'
  for i in range(imax):
    (ipl, well_id) = index_to_well(i, across_row_first)
    well_string = 'P%02d-%s' % (ipl+1, well_id)
    index_base_4 = base_convert(i+1, 4)
    encoded_index = index_base_4_to_encoded_index(index_base_4)
    well_to_encoded_index_false_dict[well_string] = encoded_index

  lig_to_encoded_index_dict = dict()
  for i in range(imax):
    lig_string = 'LIG%d' % (i+1)
    index_base_4 = base_convert(i+1, 4)
    encoded_index = index_base_4_to_encoded_index(index_base_4)
    lig_to_encoded_index_dict[lig_string] = encoded_index


with open(input_tsv) as ifh:
  next(ifh)
  # 25.0340-P5P01-E01-P7none_1732|25.0340|P01-E01|none|P08-H12_LIG53|CGGGGTAC
  i = 0
  for line in ifh:
    i = i + 1
    parts = line.split('\t')
    qname = parts[0]
    encoded_index = parts[1]

    qname_parts = qname.split('|')
    p5_well = qname_parts[2]
    p7_well = qname_parts[3]
    rt_lig  = qname_parts[4]
    rt_well = rt_lig.split('_')[0]
    lig_nam = rt_lig.split('_')[1]
    
#     print('line: %s' % (line))
#     print('rt well: %s' % (rt_well))
#     print('lig nam: %s' % (lig_nam))
#     print('p7 well: %s' % (p7_well))
#     print('p5 well: %s' % (p5_well))

    encoded_index_translated = '%s%s%s%s' % (well_to_encoded_index_true_dict[rt_well],
                        lig_to_encoded_index_dict[lig_nam],
                        well_to_encoded_index_true_dict[p7_well],
                        well_to_encoded_index_false_dict[p5_well])

#     print(encoded_index)
#     print(encoded_index_translated)
    if(encoded_index != encoded_index_translated):
      print('uh-oh')

  print('%d rows checked' % (i))
