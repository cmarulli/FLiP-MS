import os
import json

from ost import io
from ost import seq

from var3d import uniprot
from var3d import pipeline
from var3d import var_importer
from var3d import struct_importer

if not os.path.exists("/scicore/home/schwede/pantol0000/lip_ms/FLiP/data/jsons_new"):
    os.makedirs("/scicore/home/schwede/pantol0000/lip_ms/FLiP/data/jsons_new")

s_imp = struct_importer.SMRStructImporter()
import_pipeline = pipeline.DataImportPipeline(struct_importer = [s_imp])

# Process remaining entries in uniprot proteome UP000464024
def Process(AC):
    print("processing", AC)
    print(AC)
    fn = os.path.join("/scicore/home/schwede/pantol0000/lip_ms/FLiP/data/jsons_new", AC + ".json")
    if os.path.exists(fn):
        print("already there...")
        return
    entry = uniprot.FetchUniprotEntry(AC)
    sequence = seq.CreateSequence(AC, uniprot.ParseUniprotSequence(entry))
    data = pipeline.DataContainer(sequence)
    data = import_pipeline.Import(data)
    with open(fn, 'w') as fh:
        json.dump(data.ToJSON(), fh)

id_list_path = '/scicore/home/schwede/pantol0000/lip_ms/FLiP/data/uniprot_id.txt'
acs = list()
with open(id_list_path, 'r') as f:
    for line in f:
        line = line.strip()
        acs.append(line)

#acs = ["A0A663DJA2", "P0DTC3", "P0DTC4"]
n=0
with open('unprocessed_ids.txt', 'w') as file:
    for ac in acs:
        try:
            print(n)
            n+=1
            Process(ac)
        except:
            line = ac + '\n'
            file.write(line)


