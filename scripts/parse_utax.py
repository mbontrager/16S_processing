#!/usr/bin/python

import sys, getopt, csv, re

############################################################
#
# Parse the utax output
#
# Author: Martin Bontrager
############################################################

def main():
    path = ''
    
    try:
        opts, args = getopt.getopt(sys.argv[1:],'hp:',['path='])
    except getopt.GetoptError:
        print('UPARSE_pipeline.py -p path/to/UTAX/file.tax')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('UPARSE_pipeline.py -p path/to/UTAX/file.tax')
            sys.exit()
        elif opt in ('-p', '--path'):
            path = arg

    d = {}
    fix_d = {}
    
    with open(path, 'r') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        for row in csv_reader:
            d[row[0]] = row[2]

    for k, v in d.items():
        new_key = re.sub(r'(OTU_.*);.*;', r'\1', k)
        tax_list = ["Unclassified"] * 7
        fix_d[new_key] = tax_list

        v = re.sub('"Bacteroidetes"_incertae_sedis', 'Bacteroidetes_incertae_sedis', v)

        dom = re.match('d:([^,]*).*$', v)
        if dom:
            fix_d[new_key][0] = dom.group(1).strip('"')
        
        phy = re.match('.*p:([^,]*).*$', v)
        if phy:
            fix_d[new_key][1] = phy.group(1).strip('"')

        cla = re.match('.*c:([^,]*).*$', v)
        if cla:
            fix_d[new_key][2] = cla.group(1).strip('"')

        order = re.match('.*o:([^,]*).*$', v)
        if order:
            fix_d[new_key][3] = order.group(1).strip('"')

        fam = re.match('.*f:([^,]*).*$', v)
        if fam:
            fix_d[new_key][4] = fam.group(1).strip('"')
            
        gen = re.match('.*g:([^,]*).*$', v)
        if gen:
            fix_d[new_key][5] = gen.group(1).strip('"')

        spe = re.match('.*s:([^,]*).*$', v)
        if spe:
            fix_d[new_key][6] = spe.group(1).strip('"')

    with open('utaxtable.tsv', 'w', newline='') as tsv_file:
        twrite = csv.writer(tsv_file, delimiter='\t')
        for k, v in fix_d.items():
            row = [k] + v
            twrite.writerow(row)    
            
if __name__ == "__main__":
    main()
