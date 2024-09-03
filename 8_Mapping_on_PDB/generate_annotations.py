import json
import os
from var3d import struct_anno
from var3d import pipeline

working_path = '/scicore/home/schwede/pantol0000/lip_ms/FLiP'
data_path = os.path.join(working_path, 'data')
jsons_path = os.path.join(data_path, 'jsons_new')
ann_error_path = os.path.join(data_path, 'annotation_erros.txt')

#acc_anno = struct_anno.AccessibilityAnno()
interface_anno = struct_anno.InterfaceAnno()
anno_pipeline = pipeline.AnnotationPipeline(struct_anno = [interface_anno])#, acc_anno])
n=0

data_jsons = os.listdir(jsons_path)
data_jsons = [item for item in data_jsons if "_ann" not in item]

for json_f in data_jsons:
    path_to_json = os.path.join(jsons_path, json_f)
    n+=1
    print('Annotating file: ', n)
    ann_file = json_f.split('.')[0] + '_ann.' + json_f.split('.')[1]
    ann_path = os.path.join(jsons_path, ann_file)

    if os.path.exists(ann_path):
        print("already there...")
        continue

    with open(path_to_json, 'r') as fh:
        #print(path_to_json)
        data = pipeline.DataContainer.FromJSON(json.load(fh))


    try:
        anno = anno_pipeline.Run(data, structure_hashes = data.structure_hashes)
        with open(ann_path, 'w') as fh:
            json.dump(anno, fh)
    
    except Exception as e:
        with open(ann_error_path, 'a') as error_file:
            line = json_f + '\n'
            error_file.write(line)
        print('Failed in: ', json_f)
        print(e)
 
