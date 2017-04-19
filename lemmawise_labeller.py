import csv, os, pickle
import bz2

def open_dsbz2(filename):
    with bz2.BZ2File(filename, 'r') as f:
        loader = pickle.load(f)
    
    conflicts_Dict_correct = loader['conflicts_Dict_correct']
    nodelist_to_correct_mapping = loader['nodelist_to_correct_mapping']
    nodelist_correct = loader['nodelist_correct']
    featVMat_correct = loader['featVMat_correct']
    featVMat = loader['featVMat']
    conflicts_Dict = loader['conflicts_Dict']
    nodelist = loader['nodelist']
    
    return (nodelist_correct, conflicts_Dict_correct, featVMat_correct, nodelist_to_correct_mapping,\
            nodelist, conflicts_Dict, featVMat)

files = os.listdir('outputs/dump_predictions/')

with open('outputs/dump_predictions/lemmawise_labelled.csv', 'w') as out_fh:
    out_fh_csv = csv.writer(out_fh)
    fi = 0
    for root_file in files:
        with open(os.path.join('outputs/dump_predictions/', root_file)) as fh:
            print('Processing File: ', root_file)
            fh_csv = csv.reader(fh)
            for lr in fh_csv:
                if fi % 100 == 0:
                    print('Files done: ', fi)
                fi += 1
                sent_id = lr[0]
                dcs_name = sent_id + '.ds.bz2'
                (nodelist_correct, _, _, nodelist_to_correct_mapping,\
                    _, _, _) = open_dsbz2(os.path.join('../NewData/skt_dcs_DS.bz2_1L_bigram_heldout/', dcs_name))
                for rx in range(5):
                    lr = next(fh_csv)[1:]
                    if rx == 3:
                        iam = [int(x) for x in lr]
                for i in range(len(nodelist_correct)):
                    out_fh_csv.writerow([sent_id, nodelist_correct[i].lemma, 1*(nodelist_to_correct_mapping[i] in iam)])