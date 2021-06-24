import re
from collections import Counter
import os
import math
import numpy as np
from importlib.resources import open_text
from . import data


FEAT_TYPES = ('aac', 'cksaap', 'tpc', 'dpc', 'dde', 'gaac', 'cksaagp', 'gtpc', 'gdpc', 'moran', 'geary', 'nmbroto',
              'ctdc', 'ctdt', 'ctdd', 'ctriad', 'ksctriad', 'socnumber', 'qsorder', 'paac', 'apaac')


def read_fasta(file):
    """
    Reads a fasta file and returns a list of lists where the inner list contains an enzyme name and sequence.
    :param file: filepath to the fasta file
    :return: nested list where each inner list contains an enzyme name and its sequence
    """
    if not os.path.exists(file):
        raise ValueError('Error: "' + file + '" does not exist.')

    with open(file) as f:
        records = f.read()

    if re.search('>', records) is None:
        raise TypeError('The input file does not seem to be in fasta format.')

    records = records.split('>')[1:]
    myFasta = []
    for fasta in records:
        array = fasta.split('\n')
        name, sequence = array[0].split()[0], re.sub('[^ARNDCQEGHILKMFPSTWYV-]', '-', ''.join(array[1:]).upper())
        myFasta.append([name, sequence])
    return myFasta


def min_seq_len(fastas):
    """
    A function to determine the length of the shortest sequence in the given fasta file
    :param fastas:
    :return:
    """
    minLen = 10000
    for i in fastas:
        if minLen > len(i[1]):
            minLen = len(i[1])
    return minLen


def min_seq_len_norm_aa(fastas):
    """
    A function to determine the length of the shortest sequence in the given fasta file where sequences can only be
    represented by the set of 20 "Normal" amino acids
    :param fastas:
    :return:
    """
    minLen = 10000
    for i in fastas:
        if minLen > len(re.sub('-', '', i[1])):
            minLen = len(re.sub('-', '', i[1]))
    return minLen


def save_csv(encodings, file='encoding.csv'):
    """
    Save the protein feature encodings in a csv file
    :param encodings:
    :param file:
    :return:
    """
    with open(file, 'w') as f:
        for lines in encodings:
            f.write(','.join(map(str, lines)))
            f.write('\n')
    return None


def aac(fastas, **kw):
    """
    Function to create Amino Acid Composition encoding of protein sequences
    :param fastas: protein sequences in fasta format
    :param kw: optional keyword arguments that can be provided by the user
    :return: Amino Acid Composition encoding
    """
    AA = kw['order'] if kw['order'] is not None else 'ACDEFGHIKLMNPQRSTVWY'
    encodings = []

    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])
        count = Counter(sequence)
        for key in count:
            count[key] = count[key] / len(sequence)
        code = [name]
        for aa in AA:
            code.append(count[aa])
        encodings.append(code)
    return encodings


def apaac(fastas, lambda_value=30, w=0.05, **kw):
    """
    Function to generate the pseudo amino acid encoding of protein sequences
    :param fastas:
    :param lambda_value:
    :param w:
    :param kw: Not required in this case. It is given to increase make function calling generalizable
    :return:
    """
    if min_seq_len_norm_aa(fastas) < lambda_value + 1:
        raise AssertionError('All sequences should have length larger than parameter lambda_value:' +
                             str(lambda_value + 1))

    # data_dir = os.path.dirname(os.path.realpath(__file__))
    # data_file = os.path.join(data_dir, 'data', 'PAAC.txt')

    # with open(data_file) as f:
    #     records = f.readlines()
    records = open_text(data, "PAAC.txt").readlines()

    AA = ''.join(records[0].rstrip().split()[1:])
    AADict = {}
    for i in range(len(AA)):
        AADict[AA[i]] = i
    AAProperty = []
    AAPropertyNames = []
    for i in range(1, len(records) - 1):
        array = records[i].rstrip().split() if records[i].rstrip() != '' else None
        AAProperty.append([float(j) for j in array[1:]])
        AAPropertyNames.append(array[0])

    AAProperty1 = []
    for i in AAProperty:
        meanI = sum(i) / 20
        fenmu = math.sqrt(sum([(j - meanI) ** 2 for j in i]) / 20)
        AAProperty1.append([(j - meanI) / fenmu for j in i])

    encodings = []

    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])
        code = [name]
        theta = []
        for n in range(1, lambda_value + 1):
            for j in range(len(AAProperty1)):
                theta.append(sum([AAProperty1[j][AADict[sequence[k]]] * AAProperty1[j][AADict[sequence[k + n]]] for k in
                                  range(len(sequence) - n)]) / (len(sequence) - n))
        myDict = {}
        for aa in AA:
            myDict[aa] = sequence.count(aa)

        code = code + [myDict[aa] / (1 + w * sum(theta)) for aa in AA]
        code = code + [w * value / (1 + w * sum(theta)) for value in theta]
        encodings.append(code)
    return encodings


def gen_group_pairs(group_key):
    """
    Helper function for cksaagp
    :param group_key:
    :return:
    """
    gPair = {}
    for key1 in group_key:
        for key2 in group_key:
            gPair[key1 + '.' + key2] = 0
    return gPair


def cksaagp(fastas, gap=5, **kw):
    """
    Function to generate cksaagp encoding of protein sequences
    :param fastas:
    :param gap:
    :param kw:
    :return:
    """
    if gap < 0:
        raise AssertionError('Gap should be equal or greater than zero')

    if min_seq_len(fastas) < gap + 2:
        raise AssertionError(f"All sequences should have length greater than the (gap value) + 2 = {str(gap + 2)}")

    group = {
        'alphaticr': 'GAVLMI',
        'aromatic': 'FYW',
        'postivecharger': 'KRH',
        'negativecharger': 'DE',
        'uncharger': 'STCPNQ'
    }

    AA = 'ARNDCQEGHILKMFPSTWYV'

    groupKey = group.keys()

    index = {}
    for key in groupKey:
        for aa in group[key]:
            index[aa] = key

    gPairIndex = []
    for key1 in groupKey:
        for key2 in groupKey:
            gPairIndex.append(key1 + '.' + key2)

    encodings = []


    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])
        code = [name]
        for g in range(gap + 1):
            gPair = gen_group_pairs(groupKey)
            sum_var = 0
            for p1 in range(len(sequence)):
                p2 = p1 + g + 1
                if p2 < len(sequence) and sequence[p1] in AA and sequence[p2] in AA:
                    gPair[index[sequence[p1]] + '.' + index[sequence[p2]]] = gPair[index[sequence[p1]] + '.' + index[
                        sequence[p2]]] + 1
                    sum_var = sum_var + 1

            if sum_var == 0:
                for _ in gPairIndex:
                    code.append(0)
            else:
                for gp in gPairIndex:
                    code.append(gPair[gp] / sum_var)

        encodings.append(code)

    return encodings


def cksaap(fastas, gap=5, **kw):
    """
    Function to generate cksaap encoding for protein sequences
    :param fastas:
    :param gap:
    :param kw:
    :return:
    """
    if gap < 0:
        raise AssertionError("Gap should be equal or greater than zero")

    if min_seq_len(fastas) < gap + 2:
        raise AssertionError(f"All sequences should have length greater than the (gap value) + 2 = {str(gap + 2)}")

    AA = kw['order'] if kw['order'] is not None else 'ACDEFGHIKLMNPQRSTVWY'
    encodings = []

    for i in fastas:
        name, sequence = i[0], i[1]
        code = [name]
        for g in range(gap + 1):
            myDict = {}
            for pair in aaPairs:
                myDict[pair] = 0
            sum_var = 0
            for index1 in range(len(sequence)):
                index2 = index1 + g + 1
                if index1 < len(sequence) and index2 < len(sequence) and sequence[index1] in AA \
                        and sequence[index2] in AA:
                    myDict[sequence[index1] + sequence[index2]] = myDict[sequence[index1] + sequence[index2]] + 1
                    sum_var = sum_var + 1
            for pair in aaPairs:
                code.append(myDict[pair] / sum_var)
        encodings.append(code)
    return encodings


def count_ctdc(seq1, seq2):
    """
    helper function for ctdc encoding
    :param seq1:
    :param seq2:
    :return:
    """
    sum_var = 0
    for aa in seq1:
        sum_var = sum_var + seq2.count(aa)
    return sum_var


def ctdc(fastas, **kw):
    """
    Function to generate ctdc encoding for protein sequences
    :param fastas:
    :param kw:
    :return:
    """
    group1 = {
        'hydrophobicity_PRAM900101': 'RKEDQN',
        'hydrophobicity_ARGP820101': 'QSTNGDE',
        'hydrophobicity_ZIMJ680101': 'QNGSWTDERA',
        'hydrophobicity_PONP930101': 'KPDESNQT',
        'hydrophobicity_CASG920101': 'KDEQPSRNTG',
        'hydrophobicity_ENGD860101': 'RDKENQHYP',
        'hydrophobicity_FASG890101': 'KERSQD',
        'normwaalsvolume': 'GASTPDC',
        'polarity': 'LIFWCMVY',
        'polarizability': 'GASDT',
        'charge': 'KR',
        'secondarystruct': 'EALMQKRH',
        'solventaccess': 'ALFCGIVW'
    }
    group2 = {
        'hydrophobicity_PRAM900101': 'GASTPHY',
        'hydrophobicity_ARGP820101': 'RAHCKMV',
        'hydrophobicity_ZIMJ680101': 'HMCKV',
        'hydrophobicity_PONP930101': 'GRHA',
        'hydrophobicity_CASG920101': 'AHYMLV',
        'hydrophobicity_ENGD860101': 'SGTAW',
        'hydrophobicity_FASG890101': 'NTPG',
        'normwaalsvolume': 'NVEQIL',
        'polarity': 'PATGS',
        'polarizability': 'CPNVEQIL',
        'charge': 'ANCQGHILMFPSTWYV',
        'secondarystruct': 'VIYCWFT',
        'solventaccess': 'RKQEND'
    }
    group3 = {
        'hydrophobicity_PRAM900101': 'CLVIMFW',
        'hydrophobicity_ARGP820101': 'LYPFIW',
        'hydrophobicity_ZIMJ680101': 'LPFYI',
        'hydrophobicity_PONP930101': 'YMFWLCVI',
        'hydrophobicity_CASG920101': 'FIWC',
        'hydrophobicity_ENGD860101': 'CVLIMF',
        'hydrophobicity_FASG890101': 'AYHWVMFLIC',
        'normwaalsvolume': 'MHKFRYW',
        'polarity': 'HQRKNED',
        'polarizability': 'KMHFRYW',
        'charge': 'DE',
        'secondarystruct': 'GNPSD',
        'solventaccess': 'MSPTHY'
    }

    groups = [group1, group2, group3]
    property_var = (
        'hydrophobicity_PRAM900101', 'hydrophobicity_ARGP820101', 'hydrophobicity_ZIMJ680101',
        'hydrophobicity_PONP930101',
        'hydrophobicity_CASG920101', 'hydrophobicity_ENGD860101', 'hydrophobicity_FASG890101', 'normwaalsvolume',
        'polarity', 'polarizability', 'charge', 'secondarystruct', 'solventaccess')

    encodings = []

    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])
        code = [name]
        for p in property_var:
            c1 = count_ctdc(group1[p], sequence) / len(sequence)
            c2 = count_ctdc(group2[p], sequence) / len(sequence)
            c3 = 1 - c1 - c2
            code = code + [c1, c2, c3]
        encodings.append(code)
    return encodings


def count_ctdd(aa_set, sequence):
    """
    helper function for ctdd encoding
    :param aa_set:
    :param sequence:
    :return:
    """
    number = 0
    for aa in sequence:
        if aa in aa_set:
            number = number + 1
    cutoffNums = [1, math.floor(0.25 * number), math.floor(0.50 * number), math.floor(0.75 * number), number]
    cutoffNums = [i if i >= 1 else 1 for i in cutoffNums]

    code = []
    for cutoff in cutoffNums:
        myCount = 0
        for i in range(len(sequence)):
            if sequence[i] in aa_set:
                myCount += 1
                if myCount == cutoff:
                    code.append((i + 1) / len(sequence) * 100)
                    break
        if myCount == 0:
            code.append(0)
    return code


def ctdd(fastas, **kw):
    """
    Function to generate ctdd encoding for protein sequences
    :param fastas:
    :param kw:
    :return:
    """
    group1 = {
        'hydrophobicity_PRAM900101': 'RKEDQN',
        'hydrophobicity_ARGP820101': 'QSTNGDE',
        'hydrophobicity_ZIMJ680101': 'QNGSWTDERA',
        'hydrophobicity_PONP930101': 'KPDESNQT',
        'hydrophobicity_CASG920101': 'KDEQPSRNTG',
        'hydrophobicity_ENGD860101': 'RDKENQHYP',
        'hydrophobicity_FASG890101': 'KERSQD',
        'normwaalsvolume': 'GASTPDC',
        'polarity': 'LIFWCMVY',
        'polarizability': 'GASDT',
        'charge': 'KR',
        'secondarystruct': 'EALMQKRH',
        'solventaccess': 'ALFCGIVW'
    }
    group2 = {
        'hydrophobicity_PRAM900101': 'GASTPHY',
        'hydrophobicity_ARGP820101': 'RAHCKMV',
        'hydrophobicity_ZIMJ680101': 'HMCKV',
        'hydrophobicity_PONP930101': 'GRHA',
        'hydrophobicity_CASG920101': 'AHYMLV',
        'hydrophobicity_ENGD860101': 'SGTAW',
        'hydrophobicity_FASG890101': 'NTPG',
        'normwaalsvolume': 'NVEQIL',
        'polarity': 'PATGS',
        'polarizability': 'CPNVEQIL',
        'charge': 'ANCQGHILMFPSTWYV',
        'secondarystruct': 'VIYCWFT',
        'solventaccess': 'RKQEND'
    }
    group3 = {
        'hydrophobicity_PRAM900101': 'CLVIMFW',
        'hydrophobicity_ARGP820101': 'LYPFIW',
        'hydrophobicity_ZIMJ680101': 'LPFYI',
        'hydrophobicity_PONP930101': 'YMFWLCVI',
        'hydrophobicity_CASG920101': 'FIWC',
        'hydrophobicity_ENGD860101': 'CVLIMF',
        'hydrophobicity_FASG890101': 'AYHWVMFLIC',
        'normwaalsvolume': 'MHKFRYW',
        'polarity': 'HQRKNED',
        'polarizability': 'KMHFRYW',
        'charge': 'DE',
        'secondarystruct': 'GNPSD',
        'solventaccess': 'MSPTHY'
    }

    property_var = (
        'hydrophobicity_PRAM900101', 'hydrophobicity_ARGP820101', 'hydrophobicity_ZIMJ680101',
        'hydrophobicity_PONP930101',
        'hydrophobicity_CASG920101', 'hydrophobicity_ENGD860101', 'hydrophobicity_FASG890101', 'normwaalsvolume',
        'polarity', 'polarizability', 'charge', 'secondarystruct', 'solventaccess')

    encodings = []

    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])
        code = [name]
        for p in property_var:
            code = code + count_ctdd(group1[p], sequence) + \
                   count_ctdd(group2[p], sequence) + count_ctdd(group3[p], sequence)
        encodings.append(code)
    return encodings


def ctdt(fastas, **kw):
    """
    Function to generate ctdt encoding for protein sequences
    :param fastas:
    :param kw:
    :return:
    """
    group1 = {
        'hydrophobicity_PRAM900101': 'RKEDQN',
        'hydrophobicity_ARGP820101': 'QSTNGDE',
        'hydrophobicity_ZIMJ680101': 'QNGSWTDERA',
        'hydrophobicity_PONP930101': 'KPDESNQT',
        'hydrophobicity_CASG920101': 'KDEQPSRNTG',
        'hydrophobicity_ENGD860101': 'RDKENQHYP',
        'hydrophobicity_FASG890101': 'KERSQD',
        'normwaalsvolume': 'GASTPDC',
        'polarity': 'LIFWCMVY',
        'polarizability': 'GASDT',
        'charge': 'KR',
        'secondarystruct': 'EALMQKRH',
        'solventaccess': 'ALFCGIVW'
    }
    group2 = {
        'hydrophobicity_PRAM900101': 'GASTPHY',
        'hydrophobicity_ARGP820101': 'RAHCKMV',
        'hydrophobicity_ZIMJ680101': 'HMCKV',
        'hydrophobicity_PONP930101': 'GRHA',
        'hydrophobicity_CASG920101': 'AHYMLV',
        'hydrophobicity_ENGD860101': 'SGTAW',
        'hydrophobicity_FASG890101': 'NTPG',
        'normwaalsvolume': 'NVEQIL',
        'polarity': 'PATGS',
        'polarizability': 'CPNVEQIL',
        'charge': 'ANCQGHILMFPSTWYV',
        'secondarystruct': 'VIYCWFT',
        'solventaccess': 'RKQEND'
    }
    group3 = {
        'hydrophobicity_PRAM900101': 'CLVIMFW',
        'hydrophobicity_ARGP820101': 'LYPFIW',
        'hydrophobicity_ZIMJ680101': 'LPFYI',
        'hydrophobicity_PONP930101': 'YMFWLCVI',
        'hydrophobicity_CASG920101': 'FIWC',
        'hydrophobicity_ENGD860101': 'CVLIMF',
        'hydrophobicity_FASG890101': 'AYHWVMFLIC',
        'normwaalsvolume': 'MHKFRYW',
        'polarity': 'HQRKNED',
        'polarizability': 'KMHFRYW',
        'charge': 'DE',
        'secondarystruct': 'GNPSD',
        'solventaccess': 'MSPTHY'
    }

    property_var = (
        'hydrophobicity_PRAM900101', 'hydrophobicity_ARGP820101', 'hydrophobicity_ZIMJ680101',
        'hydrophobicity_PONP930101',
        'hydrophobicity_CASG920101', 'hydrophobicity_ENGD860101', 'hydrophobicity_FASG890101', 'normwaalsvolume',
        'polarity', 'polarizability', 'charge', 'secondarystruct', 'solventaccess')

    encodings = []

    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])
        code = [name]
        aaPair = [sequence[j:j + 2] for j in range(len(sequence) - 1)]
        for p in property_var:
            c1221, c1331, c2332 = 0, 0, 0
            for pair in aaPair:
                if (pair[0] in group1[p] and pair[1] in group2[p]) or (pair[0] in group2[p] and pair[1] in group1[p]):
                    c1221 = c1221 + 1
                    continue
                if (pair[0] in group1[p] and pair[1] in group3[p]) or (pair[0] in group3[p] and pair[1] in group1[p]):
                    c1331 = c1331 + 1
                    continue
                if (pair[0] in group2[p] and pair[1] in group3[p]) or (pair[0] in group3[p] and pair[1] in group2[p]):
                    c2332 = c2332 + 1
            code = code + [c1221 / len(aaPair), c1331 / len(aaPair), c2332 / len(aaPair)]
        encodings.append(code)
    return encodings


def calc_ksctriad(sequence, gap, features, aa_dict):
    """
    Helper function for ctraid and ksctraid
    :param sequence:
    :param gap:
    :param features:
    :param aa_dict:
    :return:
    """
    res = []
    for g in range(gap + 1):
        myDict = {}
        for f in features:
            myDict[f] = 0

        for i in range(len(sequence)):
            if i + gap + 1 < len(sequence) and i + 2 * gap + 2 < len(sequence):
                fea = aa_dict[sequence[i]] + '.' + aa_dict[sequence[i + gap + 1]] + '.' + aa_dict[
                    sequence[i + 2 * gap + 2]]
                myDict[fea] = myDict[fea] + 1

        maxValue, minValue = max(myDict.values()), min(myDict.values())
        for f in features:
            res.append((myDict[f] - minValue) / maxValue)

    return res


def ctriad(fastas, **kw):
    """
    Function to generate ctriad encoding for protein sequences
    :param fastas: 
    :param kw: 
    :return: 
    """
    # TODO: Need to check if the gap parameter should be kept or is it only needed for ksctriad
    AAGroup = {
        'g1': 'AGV',
        'g2': 'ILFP',
        'g3': 'YMTS',
        'g4': 'HNQW',
        'g5': 'RK',
        'g6': 'DE',
        'g7': 'C'
    }

    myGroups = sorted(AAGroup.keys())

    AADict = {}
    for g in myGroups:
        for aa in AAGroup[g]:
            AADict[aa] = g

    features = [f1 + '.' + f2 + '.' + f3 for f1 in myGroups for f2 in myGroups for f3 in myGroups]

    encodings = []

    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])
        code = [name]
        if len(sequence) < 3:
            print('Error: for "CTriad" encoding, the input fasta sequences should be greater than 3. \n\n')
            return 0
        code = code + calc_ksctriad(sequence, 0, features, AADict)
        encodings.append(code)

    return encodings


def dde(fastas, **kw):
    """
    Function to generate dde encoding for protein sequences
    :param fastas:
    :param kw:
    :return:
    """
    AA = kw['order'] if kw['order'] is not None else 'ACDEFGHIKLMNPQRSTVWY'

    myCodons = {
        'A': 4,
        'C': 2,
        'D': 2,
        'E': 2,
        'F': 2,
        'G': 4,
        'H': 2,
        'I': 3,
        'K': 2,
        'L': 6,
        'M': 1,
        'N': 2,
        'P': 4,
        'Q': 2,
        'R': 6,
        'S': 6,
        'T': 4,
        'V': 4,
        'W': 1,
        'Y': 2
    }

    encodings = []

    myTM = []
    for pair in diPeptides:
        myTM.append((myCodons[pair[0]] / 61) * (myCodons[pair[1]] / 61))

    AADict = {}
    for i in range(len(AA)):
        AADict[AA[i]] = i

    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])
        code = [name]
        tmpCode = [0] * 400
        for j in range(len(sequence) - 2 + 1):
            tmpCode[AADict[sequence[j]] * 20 + AADict[sequence[j + 1]]] = tmpCode[AADict[sequence[j]] * 20 + AADict[
                sequence[j + 1]]] + 1
        if sum(tmpCode) != 0:
            tmpCode = [i / sum(tmpCode) for i in tmpCode]

        myTV = []
        for j in range(len(myTM)):
            myTV.append(myTM[j] * (1 - myTM[j]) / (len(sequence) - 1))

        for j in range(len(tmpCode)):
            tmpCode[j] = (tmpCode[j] - myTM[j]) / math.sqrt(myTV[j])

        code = code + tmpCode
        encodings.append(code)
    return encodings


def dpc(fastas, **kw):
    """
    Function to generate dpc encoding for protein sequences
    :param fastas:
    :param kw:
    :return:
    """
    AA = kw['order'] if kw['order'] is not None else 'ACDEFGHIKLMNPQRSTVWY'
    encodings = []

    AADict = {}
    for i in range(len(AA)):
        AADict[AA[i]] = i

    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])
        code = [name]
        tmpCode = [0] * 400
        for j in range(len(sequence) - 2 + 1):
            tmpCode[AADict[sequence[j]] * 20 + AADict[sequence[j + 1]]] = tmpCode[AADict[sequence[j]] * 20 + AADict[
                sequence[j + 1]]] + 1
        if sum(tmpCode) != 0:
            tmpCode = [i / sum(tmpCode) for i in tmpCode]
        code = code + tmpCode
        encodings.append(code)
    return encodings


def gaac(fastas, **kw):
    """
    Function to generate gaac encoding for protein sequences
    :param fastas:
    :param kw:
    :return:
    """
    group = {
        'alphatic': 'GAVLMI',
        'aromatic': 'FYW',
        'postivecharge': 'KRH',
        'negativecharge': 'DE',
        'uncharge': 'STCPNQ'
    }

    groupKey = group.keys()

    encodings = []

    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])
        code = [name]
        count = Counter(sequence)
        myDict = {}
        for key in groupKey:
            for aa in group[key]:
                myDict[key] = myDict.get(key, 0) + count[aa]

        for key in groupKey:
            code.append(myDict[key] / len(sequence))
        encodings.append(code)

    return encodings


def gdpc(fastas, **kw):
    """
    Function to generate gdpc encoding for protein sequences
    :param fastas:
    :param kw:
    :return:
    """
    group = {
        'alphaticr': 'GAVLMI',
        'aromatic': 'FYW',
        'postivecharger': 'KRH',
        'negativecharger': 'DE',
        'uncharger': 'STCPNQ'
    }

    groupKey = group.keys()
    dipeptide = [g1 + '.' + g2 for g1 in groupKey for g2 in groupKey]

    index = {}
    for key in groupKey:
        for aa in group[key]:
            index[aa] = key

    encodings = []

    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])

        code = [name]
        myDict = {}
        for t in dipeptide:
            myDict[t] = 0

        sum_var = 0
        for j in range(len(sequence) - 2 + 1):
            myDict[index[sequence[j]] + '.' + index[sequence[j + 1]]] = myDict[index[sequence[j]] + '.' + index[
                sequence[j + 1]]] + 1
            sum_var = sum_var + 1

        if sum_var == 0:
            for _ in dipeptide:
                code.append(0)
        else:
            for t in dipeptide:
                code.append(myDict[t] / sum_var)
        encodings.append(code)

    return encodings


def geary(fastas, props=('CIDH920105', 'BHAR880101', 'CHAM820101', 'CHAM820102',
                         'CHOC760101', 'BIGC670101', 'CHAM810101', 'DAYM780201'),
          nlag=30, **kw):
    """
    Function to generate geary encoding for protein sequences
    :param fastas:
    :param props:
    :param nlag:
    :param kw:
    :return:
    """
    if min_seq_len_norm_aa(fastas) < nlag + 1:
        raise AssertionError(f"All sequences should have length greater than nlag+1: {str(nlag + 1)}")

    AA = 'ARNDCQEGHILKMFPSTWYV'

    data_dir = os.path.dirname(os.path.realpath(__file__))
    data_file = os.path.join(data_dir, 'data', 'AAidx.txt')

    with open(data_file) as f:
        records = f.readlines()[1:]
    myDict = {}
    for i in records:
        array = i.rstrip().split('\t')
        myDict[array[0]] = array[1:]

    aa_idx = []
    aa_idx_name = []
    for i in props:
        if i in myDict:
            aa_idx.append(myDict[i])
            aa_idx_name.append(i)
        else:
            print('"' + i + '" properties not exist.')
            return None

    aa_idx1 = np.array([float(j) for i in aa_idx for j in i])
    aa_idx = aa_idx1.reshape((len(aa_idx), 20))

    propMean = np.mean(aa_idx, axis=1)
    propStd = np.std(aa_idx, axis=1)

    for i in range(len(aa_idx)):
        for j in range(len(aa_idx[i])):
            aa_idx[i][j] = (aa_idx[i][j] - propMean[i]) / propStd[i]

    index = {}
    for i in range(len(AA)):
        index[AA[i]] = i

    encodings = []


    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])
        code = [name]
        N = len(sequence)
        for prop in range(len(props)):
            xmean = sum([aa_idx[prop][index[aa]] for aa in sequence]) / N
            for n in range(1, nlag + 1):
                if len(sequence) > nlag:
                    # if key is '-', then the value is 0
                    rn = (N - 1) / (2 * (N - n)) * ((sum(
                        [(aa_idx[prop][index.get(sequence[j], 0)] - aa_idx[prop][index.get(sequence[j + n], 0)]) ** 2 for
                         j in range(len(sequence) - n)])) / (sum([(aa_idx[prop][index.get(sequence[j], 0)] - xmean) ** 2
                                                                  for j in range(len(sequence))])))
                else:
                    rn = 'NA'
                code.append(rn)
        encodings.append(code)
    return encodings


def gtpc(fastas, **kw):
    """
    Function to generate gtpc encoding for protein sequences
    :param fastas:
    :param kw:
    :return:
    """
    group = {
        'alphaticr': 'GAVLMI',
        'aromatic': 'FYW',
        'postivecharger': 'KRH',
        'negativecharger': 'DE',
        'uncharger': 'STCPNQ'
    }

    groupKey = group.keys()
    triple = [g1 + '.' + g2 + '.' + g3 for g1 in groupKey for g2 in groupKey for g3 in groupKey]

    index = {}
    for key in groupKey:
        for aa in group[key]:
            index[aa] = key

    encodings = []

    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])

        code = [name]
        myDict = {}
        for t in triple:
            myDict[t] = 0

        sum_var = 0
        for j in range(len(sequence) - 3 + 1):
            myDict[index[sequence[j]] + '.' + index[sequence[j + 1]] + '.' + index[sequence[j + 2]]] = \
                myDict[index[sequence[j]] + '.' + index[sequence[j + 1]] + '.' + index[sequence[j + 2]]] + 1
            sum_var = sum_var + 1

        if sum_var == 0:
            for _ in triple:
                code.append(0)
        else:
            for t in triple:
                code.append(myDict[t] / sum_var)
        encodings.append(code)

    return encodings


def ksctriad(fastas, gap=0, **kw):
    """
    Function to generate ksctriad encoding for protein sequences
    :param fastas:
    :param gap:
    :param kw:
    :return:
    """
    AAGroup = {
        'g1': 'AGV',
        'g2': 'ILFP',
        'g3': 'YMTS',
        'g4': 'HNQW',
        'g5': 'RK',
        'g6': 'DE',
        'g7': 'C'
    }

    myGroups = sorted(AAGroup.keys())

    AADict = {}
    for g in myGroups:
        for aa in AAGroup[g]:
            AADict[aa] = g

    features = [f1 + '.' + f2 + '.' + f3 for f1 in myGroups for f2 in myGroups for f3 in myGroups]

    encodings = []

    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])
        code = [name]
        if len(sequence) < 2 * gap + 3:
            raise AssertionError(f"All sequences should have length greater than (2*gap+3): {str(2*gap+3)}")
        code = code + calc_ksctriad(sequence, gap, features, AADict)
        encodings.append(code)

    return encodings


def moran(fastas, props=('CIDH920105', 'BHAR880101', 'CHAM820101', 'CHAM820102',
                         'CHOC760101', 'BIGC670101', 'CHAM810101', 'DAYM780201'),
          nlag=30, **kw):
    if min_seq_len_norm_aa(fastas) < nlag + 1:
        raise AssertionError(f"All sequences should have length greater than nlag+1: {str(nlag + 1)}")


    AA = 'ARNDCQEGHILKMFPSTWYV'

    data_dir = os.path.dirname(os.path.realpath(__file__))
    data_file = os.path.join(data_dir, 'data', 'AAidx.txt')

    with open(data_file) as f:
        records = f.readlines()[1:]
    myDict = {}
    for i in records:
        array = i.rstrip().split('\t')
        myDict[array[0]] = array[1:]

    aa_idx = []
    aa_idx_name = []
    for i in props:
        if i in myDict:
            aa_idx.append(myDict[i])
            aa_idx_name.append(i)
        else:
            print('"' + i + '" properties not exist.')
            return None

    aa_idx1 = np.array([float(j) for i in aa_idx for j in i])
    aa_idx = aa_idx1.reshape((len(aa_idx), 20))

    propMean = np.mean(aa_idx, axis=1)
    propStd = np.std(aa_idx, axis=1)

    for i in range(len(aa_idx)):
        for j in range(len(aa_idx[i])):
            aa_idx[i][j] = (aa_idx[i][j] - propMean[i]) / propStd[i]

    index = {}
    for i in range(len(AA)):
        index[AA[i]] = i

    encodings = []

    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])
        code = [name]
        N = len(sequence)
        for prop in range(len(props)):
            xmean = sum([aa_idx[prop][index[aa]] for aa in sequence]) / N
            for n in range(1, nlag + 1):
                if len(sequence) > nlag:
                    # if key is '-', then the value is 0
                    fenzi = sum([(aa_idx[prop][index.get(sequence[j], 0)] - xmean) * (
                            aa_idx[prop][index.get(sequence[j + n], 0)] - xmean) for j in
                                 range(len(sequence) - n)]) / (N - n)
                    fenmu = sum(
                        [(aa_idx[prop][index.get(sequence[j], 0)] - xmean) ** 2 for j in range(len(sequence))]) / N
                    rn = fenzi / fenmu
                else:
                    rn = 'NA'
                code.append(rn)
        encodings.append(code)
    return encodings


def nmbroto(fastas, props=('CIDH920105', 'BHAR880101', 'CHAM820101', 'CHAM820102',
                           'CHOC760101', 'BIGC670101', 'CHAM810101', 'DAYM780201'),
            nlag=30, **kw):
    """
    Function to generate nmbroto encoding for protein sequences
    :param fastas:
    :param props:
    :param nlag:
    :param kw:
    :return:
    """
    if min_seq_len_norm_aa(fastas) < nlag + 1:
        raise AssertionError(f"All sequences should have length greater than nlag+1: {str(nlag + 1)}")

    AA = 'ARNDCQEGHILKMFPSTWYV'
    data_dir = os.path.dirname(os.path.realpath(__file__))
    data_file = os.path.join(data_dir, 'data', 'AAidx.txt')

    with open(data_file) as f:
        records = f.readlines()[1:]
    myDict = {}
    for i in records:
        array = i.rstrip().split('\t')
        myDict[array[0]] = array[1:]

    aa_idx = []
    aa_idx_name = []
    for i in props:
        if i in myDict:
            aa_idx.append(myDict[i])
            aa_idx_name.append(i)
        else:
            print('"' + i + '" properties not exist.')
            return None

    aa_idx1 = np.array([float(j) for i in aa_idx for j in i])
    aa_idx = aa_idx1.reshape((len(aa_idx), 20))
    pstd = np.std(aa_idx, axis=1)
    pmean = np.average(aa_idx, axis=1)

    for i in range(len(aa_idx)):
        for j in range(len(aa_idx[i])):
            aa_idx[i][j] = (aa_idx[i][j] - pmean[i]) / pstd[i]

    index = {}
    for i in range(len(AA)):
        index[AA[i]] = i

    encodings = []

    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])
        code = [name]
        N = len(sequence)
        for prop in range(len(props)):
            for n in range(1, nlag + 1):
                if len(sequence) > nlag:
                    # if key is '-', then the value is 0
                    rn = sum(
                        [aa_idx[prop][index.get(sequence[j], 0)] * aa_idx[prop][index.get(sequence[j + n], 0)] for j in
                         range(len(sequence) - n)]) / (N - n)
                else:
                    rn = 'NA'
                code.append(rn)
        encodings.append(code)
    return encodings


def r_value(aa1, aa2, aa_dict, matrix):
    """

    :param aa1:
    :param aa2:
    :param aa_dict:
    :param matrix:
    :return:
    """
    return sum([(matrix[i][aa_dict[aa1]] - matrix[i][aa_dict[aa2]]) ** 2 for i in range(len(matrix))]) / len(matrix)


def paac(fastas, lambda_value=30, w=0.05, **kw):
    """
    Function to generate paac encoding for protein sequences
    :param fastas:
    :param lambda_value:
    :param w:
    :param kw:
    :return:
    """
    if min_seq_len_norm_aa(fastas) < lambda_value + 1:
        raise AssertionError(f'All sequences should have length greater than lambda_value+1: {str(lambda_value + 1)}')

    data_dir = os.path.dirname(os.path.realpath(__file__))
    data_file = os.path.join(data_dir, 'data', 'PAAC.txt')

    with open(data_file) as f:
        records = f.readlines()
    AA = ''.join(records[0].rstrip().split()[1:])
    AADict = {}
    for i in range(len(AA)):
        AADict[AA[i]] = i
    AAProperty = []
    AAPropertyNames = []
    for i in range(1, len(records)):
        array = records[i].rstrip().split() if records[i].rstrip() != '' else None
        AAProperty.append([float(j) for j in array[1:]])
        AAPropertyNames.append(array[0])

    AAProperty1 = []
    for i in AAProperty:
        meanI = sum(i) / 20
        fenmu = math.sqrt(sum([(j - meanI) ** 2 for j in i]) / 20)
        AAProperty1.append([(j - meanI) / fenmu for j in i])

    encodings = []

    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])
        code = [name]
        theta = []
        for n in range(1, lambda_value + 1):
            theta.append(
                sum([r_value(sequence[j], sequence[j + n], AADict, AAProperty1) for j in range(len(sequence) - n)]) / (
                        len(sequence) - n))
        myDict = {}
        for aa in AA:
            myDict[aa] = sequence.count(aa)
        code = code + [myDict[aa] / (1 + w * sum(theta)) for aa in AA]
        code = code + [(w * j) / (1 + w * sum(theta)) for j in theta]
        encodings.append(code)
    return encodings


def qsorder(fastas, nlag=30, w=0.1, **kw):
    """
    Function to generate qsorder encoding for protein sequences
    :param fastas:
    :param nlag:
    :param w:
    :param kw:
    :return:
    """
    if min_seq_len_norm_aa(fastas) < nlag + 1:
        raise AssertionError(f'All sequences should have length greater than the nlag+1: {str(nlag + 1)}')

    data_dir = os.path.dirname(os.path.realpath(__file__))
    data_file = os.path.join(data_dir, 'data', 'Schneider-Wrede.txt')
    data_file1 = os.path.join(data_dir, 'data', 'Grantham.txt')

    AA = 'ACDEFGHIKLMNPQRSTVWY'
    AA1 = 'ARNDCQEGHILKMFPSTWYV'

    DictAA = {}
    for i in range(len(AA)):
        DictAA[AA[i]] = i

    DictAA1 = {}
    for i in range(len(AA1)):
        DictAA1[AA1[i]] = i

    with open(data_file) as f:
        records = f.readlines()[1:]
    AADistance = []
    for i in records:
        array = i.rstrip().split()[1:] if i.rstrip() != '' else None
        AADistance.append(array)
    AADistance = np.array(
        [float(AADistance[i][j]) for i in range(len(AADistance)) for j in range(len(AADistance[i]))]).reshape((20, 20))

    with open(data_file1) as f:
        records = f.readlines()[1:]
    AADistance1 = []
    for i in records:
        array = i.rstrip().split()[1:] if i.rstrip() != '' else None
        AADistance1.append(array)
    AADistance1 = np.array(
        [float(AADistance1[i][j]) for i in range(len(AADistance1)) for j in range(len(AADistance1[i]))]).reshape(
        (20, 20))

    encodings = []

    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])
        code = [name]
        arraySW = []
        arrayGM = []
        for n in range(1, nlag + 1):
            arraySW.append(
                sum([AADistance[DictAA[sequence[j]]][DictAA[sequence[j + n]]] ** 2 for j in range(len(sequence) - n)]))
            arrayGM.append(sum(
                [AADistance1[DictAA1[sequence[j]]][DictAA1[sequence[j + n]]] ** 2 for j in range(len(sequence) - n)]))
        myDict = {}
        for aa in AA1:
            myDict[aa] = sequence.count(aa)
        for aa in AA1:
            code.append(myDict[aa] / (1 + w * sum(arraySW)))
        for aa in AA1:
            code.append(myDict[aa] / (1 + w * sum(arrayGM)))
        for num in arraySW:
            code.append((w * num) / (1 + w * sum(arraySW)))
        for num in arrayGM:
            code.append((w * num) / (1 + w * sum(arrayGM)))
        encodings.append(code)
    return encodings


def socnumber(fastas, nlag=30, **kw):
    """
    Function to generate socnumber encoding for protein sequences
    :param fastas:
    :param nlag:
    :param kw:
    :return:
    """
    if min_seq_len_norm_aa(fastas) < nlag + 1:
        print('Error: all the sequence length should be larger than the nlag+1: ' + str(nlag + 1) + '\n\n')
        return 0

    data_dir = os.path.dirname(os.path.realpath(__file__))
    data_file = os.path.join(data_dir, 'data', 'Schneider-Wrede.txt')
    data_file1 = os.path.join(data_dir, 'data', 'Grantham.txt')

    AA = 'ACDEFGHIKLMNPQRSTVWY'
    AA1 = 'ARNDCQEGHILKMFPSTWYV'

    DictAA = {}
    for i in range(len(AA)):
        DictAA[AA[i]] = i

    DictAA1 = {}
    for i in range(len(AA1)):
        DictAA1[AA1[i]] = i

    with open(data_file) as f:
        records = f.readlines()[1:]
    AADistance = []
    for i in records:
        array = i.rstrip().split()[1:] if i.rstrip() != '' else None
        AADistance.append(array)
    AADistance = np.array(
        [float(AADistance[i][j]) for i in range(len(AADistance)) for j in range(len(AADistance[i]))]).reshape((20, 20))

    with open(data_file1) as f:
        records = f.readlines()[1:]
    AADistance1 = []
    for i in records:
        array = i.rstrip().split()[1:] if i.rstrip() != '' else None
        AADistance1.append(array)
    AADistance1 = np.array(
        [float(AADistance1[i][j]) for i in range(len(AADistance1)) for j in range(len(AADistance1[i]))]).reshape(
        (20, 20))

    encodings = []

    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])
        code = [name]
        for n in range(1, nlag + 1):
            code.append(sum(
                [AADistance[DictAA[sequence[j]]][DictAA[sequence[j + n]]] ** 2 for j in range(len(sequence) - n)]) / (
                                len(sequence) - n))

        for n in range(1, nlag + 1):
            code.append(sum([AADistance1[DictAA1[sequence[j]]][DictAA1[sequence[j + n]]] ** 2 for j in
                             range(len(sequence) - n)]) / (len(sequence) - n))
        encodings.append(code)
    return encodings


def tpc(fastas, **kw):
    """
    Function to generate tpc encoding for protein sequences
    :param fastas:
    :param kw:
    :return:
    """
    AA = kw['order'] if kw['order'] is not None else 'ACDEFGHIKLMNPQRSTVWY'
    encodings = []

    AADict = {}
    for i in range(len(AA)):
        AADict[AA[i]] = i

    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])
        code = [name]
        tmpCode = [0] * 8000
        for j in range(len(sequence) - 3 + 1):
            tmpCode[AADict[sequence[j]] * 400 + AADict[sequence[j + 1]] * 20 + AADict[sequence[j + 2]]] = \
                tmpCode[AADict[sequence[j]] * 400 + AADict[sequence[j + 1]] * 20 + AADict[sequence[j + 2]]] + 1
        if sum(tmpCode) != 0:
            tmpCode = [i / sum(tmpCode) for i in tmpCode]
        code = code + tmpCode
        encodings.append(code)
    return encodings


def get_feature(protein_fasta_file, feature_type, output_dir):
    """
    A function to create numerically encoded feature of a specific type for protein sequences
    :param protein_fasta_file: The path to a file that contains all the protein sequences in fasta format
    :param feature_type: The feature encoding type. Must be one of the 21 types mentioned in README
    :param output_dir: The path to a directory where the feature encoded files will be stored
    :return: Name of the output file where the features are stored
    """
    if feature_type not in FEAT_TYPES:
        raise ValueError(f"feature_type must be one of {FEAT_TYPES}")

    fastas = read_fasta(protein_fasta_file)
    myOrder = 'ACDEFGHIKLMNPQRSTVWY'
    kw = {'order': myOrder}

    myFun = f"{feature_type}(fastas, **kw)"
    print('Descriptor type: ' + feature_type)
    encodings = eval(myFun)

    if os.path.exists(output_dir):
        if os.path.isdir(output_dir):
            output_file = os.path.join(output_dir, feature_type+'.csv')
            save_csv(encodings, output_file)
            return output_file
        else:
            raise OSError(f"output_dir: {output_dir} is not a directory!")

    else:
        raise OSError(f"output_dir: {output_dir} does not exist!")

    pass


def get_all_features(fasta_file, output_dir):
    """
    A function to create 21 numerically encoded features for protein sequences
    :param fasta_file: The path to a file that contains all the protein sequences in fasta format
    :param output_dir: The path to a directory where the feature encoded files will be stored
    :return: None
    """
    for ft in FEAT_TYPES:
        get_feature(fasta_file, ft, output_dir)
    return
