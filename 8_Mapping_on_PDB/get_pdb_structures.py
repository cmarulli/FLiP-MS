import os
import json

from ost import io
from ost import seq

from var3d import uniprot
from var3d import pipeline
from var3d import var_importer
from var3d import struct_importer


jsons_path = "data/jsons"
if not os.path.exists(jsons_path):
    os.makedirs(jsons_path)

s_imp = struct_importer.SMRStructImporter()
import_pipeline = pipeline.DataImportPipeline(struct_importer = [s_imp])

# Process remaining entries in uniprot proteome UP000464024
def Process(AC):
    print("processing", AC)
    print(AC)
    fn = os.path.join(jsons_path, AC + ".json")
    if os.path.exists(fn):
        print("already there...")
        return
    entry = uniprot.FetchUniprotEntry(AC)
    sequence = seq.CreateSequence(AC, uniprot.ParseUniprotSequence(entry))
    data = pipeline.DataContainer(sequence)
    data = import_pipeline.Import(data)
    with open(fn, 'w') as fh:
        json.dump(data.ToJSON(), fh)

id_list_path = "data/uniprot_id.txt"
acs = list()
with open(id_list_path, 'r') as f:
    for line in f:
        line = line.strip()
        acs.append(line)


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


