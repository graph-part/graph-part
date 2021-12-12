import collections
import os
import re

def parse_pdb(filepath):
    '''Adapted from biopython source. They have protein hardcoded.
    https://github.com/biopython/biopython/blob/master/Bio/SeqIO/PdbIO.py
    '''
    chains = collections.defaultdict(list)
    metadata = collections.defaultdict(list)
    rec_name = None
    with open(filepath) as f:
        for line in f:
            rec_name = line[0:6].strip()
            if rec_name == "SEQRES":
                chn_id = line[11]
                residues = [res for res in line[19:].split()]
                chains[chn_id].extend(residues)
            elif rec_name == "DBREF":
                pdb_id = line[7:11]
                chn_id = line[12]
                database = line[26:32].strip()
                db_acc = line[33:41].strip()
                db_id_code = line[42:54].strip()
                metadata[chn_id].append(
                    {
                        "pdb_id": pdb_id,
                        "database": database,
                        "db_acc": db_acc,
                        "db_id_code": db_id_code,
                    }
                )

        if rec_name is None:
            raise ValueError("Empty file.")
            

        out_list = []
        for chn_id, residues in sorted(chains.items()):
            rec = {'sequence': "".join(residues),
                   'chain': chn_id,
                  }
            if chn_id in metadata:
                m = metadata[chn_id][0]
                rec['id'] = f"{m['pdb_id']}:{chn_id}"

            # TODO this will cause duplicates
            else:
                pdb_id = filepath[:-4].split('/')[1]
                rec['id'] = f"{pdb_id}:{chn_id}"
            out_list.append(rec)
        return out_list


directory = 'raw'
recs = []
for file in os.listdir(directory):
    full_path = os.path.join(directory, file)
    if file[-3:] != 'pdb':
        continue

    recs = recs + parse_pdb(full_path)

filtered = []
for rec in recs:
    seq = rec['sequence']
    if  len(re.findall(r'[^AUGC]', seq)) != 0:
        continue
    else:
        filtered.append(rec)



with open('spot_rna_1d.fasta', 'w') as f:
    for rec in filtered:
        f.write('>'+rec['id'].replace(':','_')+'\n')
        f.write(rec['sequence']+'\n')