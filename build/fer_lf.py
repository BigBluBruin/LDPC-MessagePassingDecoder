import argparse
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser('FER Calculation')
    parser.add_argument('--snr_str', type=str, default='000', nargs='*', help='snr strings ')
    parser.add_argument('--chop', type=str, default="0", nargs='*', help='right precision')
    arglist = parser.parse_args()
    return arglist


def cal_fer(filename):
    pre = np.loadtxt(filename)
    result = np.zeros(5)

    
    if pre.ndim>1:
        ss = np.sum(pre, axis=0)
        dimm = np.shape(pre)[0]
        result[0] = pre[0][0]
    else :
        ss = pre
        dimm = 1
        result[0] = pre[0]
    result [1] = ss[1]
    result [2] = ss[2]
    result [3] = ss[2]/ss[1]
    result [4] = ss[3]/ dimm
    return result




def main(arglist):
    result = np.empty(0)
    for rr in arglist.chop:
        result = np.empty(0)
        for ii in arglist.snr_str :
            filename = "result_n_2193_k_1161_SNR_" + ii + "_Rate_0.50_check_" + rr + ".txt"
            result = np.append(result, cal_fer(filename))
            print(result)
        result = np.reshape(result,(-1,5))
        outputfile = "result_"+rr+".txt"
        np.savetxt(outputfile, result)
        




if __name__ == "__main__":
    arglist = parse_args()
    main(arglist)
