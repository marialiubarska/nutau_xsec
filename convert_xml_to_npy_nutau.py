#!/usr/bin/env python 

import numpy as np
import xml.etree.ElementTree as ET
import pickle
import os
import argparse

def get_args():
    parser = argparse.ArgumentParser(description='get path to xml cross section file')
    parser.add_argument('--file', type=str, help='')
    args = parser.parse_args()
    return args

def conv_to_dict_v2(root):
    xsec = {}
    E = {}
    for elem in root:
        nm = elem.get('name')
        print(nm)
        xsec[nm] = []
        E[nm] = []
        for subelem in elem:
            xsec[nm].append( float(subelem[1].text) )
            E[nm].append( float(subelem[0].text) )
    return E, xsec

def conv_to_dict_v3(root):
    xsec = {}
    E =	{}
    for elem in root:
        for subelem in elem:
            nm = subelem.get('name')
            print(nm)
            xsec[nm] = []
            E[nm] = []
            for subsubelem in subelem:
                xsec[nm].append( float(subsubelem[1].text) )
                E[nm].append( float(subsubelem[0].text) )
    return E, xsec

def make_fits(E, xsec):

    E_fit = np.logspace(0.,4.,num=200)
    
    xsec_fits = {}
    for k in E.keys():
        x = np.array(E[k])
        y = np.array(xsec[k])/x
        yinterp = np.interp(E_fit, x, y)
        xsec_fits[k] = yinterp
#        fitpar = np.polyfit(np.log10(x), y, 8)
#        xsec_fits[k] = np.polyval(fitpar, logE_fit)

    return E_fit, xsec_fits

def sum_fits_all_regular_genie(E_fit, xsec_fits, E, xsec):

    genie_xml = {'E_fit':E_fit, 'xsec_fits':[], 'full_name':[], 'E_init':[], 'xsec_init':[],
                 'pdg':[], 'tgt':[], 'hitN':[], 'hitqrk':[], 'sea':[], 'CC':[], 'qe':[], 'charm_incl':[]}

    for PDG in [-16,-14,-12,12,14,16]:

        for TGT in [1000000010, 1000010010]:
            if TGT == 1000000010:
                HITN = 2112
            elif TGT == 1000010010:
                HITN = 2212
            else:
                continue
            for CURR in ['CC', 'NC']:

                passed = False
                xsec_this_step = np.zeros(len(E_fit))
                
                for k in xsec_fits.keys():
                
                    st_pdg = k.find('/nu:') + len('/nu:')
                    end_pdg = st_pdg + k[st_pdg:].find(';')
                    pdg = int(k[st_pdg:end_pdg])
                    if pdg != PDG: continue
                    print pdg

                    st_tgt = k.find(';tgt:') +len(';tgt:')
                    end_tgt = st_tgt + k[st_tgt:].find(';')
                    tgt = int(k[st_tgt:end_tgt])
                    if tgt != TGT: continue
                    print tgt

                    if ';N:' not in k:
                        print 'no hitN', k
                        continue
                    st_N = k.find(';N:') +len(';N:')
                    end_N = st_N + k[st_N:].find(';')
                    hitN = int(k[st_N:end_N])
                    if hitN != HITN: continue
                    print hitN

                    st_curr = k.find('proc:Weak[') + len('proc:Weak[')
                    end_curr =st_curr +  k[st_curr:].find(']')
                    if  k[st_curr:end_curr] != CURR: continue
                    print k[st_curr:end_curr]
                    
                    xsec_this_step += np.array(xsec_fits[k])
                    passed = True

                if passed:
                    if CURR == 'CC':
                        cc = 1
                    elif CURR == 'NC':
                        cc = 0
                    else:
                        cc = -1
                
                    genie_xml['pdg'].append(PDG)
                    genie_xml['tgt'].append(TGT)
                    genie_xml['hitN'].append(HITN)
                    genie_xml['CC'].append(cc)

                    genie_xml['full_name'].append(k)
                    genie_xml['xsec_fits'].append(xsec_this_step)

#                    genie_xml['E_init'].append(E[k])
#                    genie_xml['xsec_init'].append(xsec[k])        
    
    return genie_xml

def conv_xml_to_dict(filename):

    tree = ET.parse(filename)
    root = tree.getroot()

    genie_version = float(root.get('version'))
    is_hedis = False
    
    if genie_version == 2:
        print('detected GENIE version = 2')
        E, xsec = conv_to_dict_v2(root)
       
    elif genie_version == 3:
        print('detected GENIE version = 3')
        E, xsec = conv_to_dict_v3(root)
        
        genie_tune = root.getchildren()[0].get('name')
        if genie_tune[:3] == 'GHE':
            is_hedis = True
            print('detected GENIE-HEDIS tune!')

    else:
        print('Failed to detect GENIE version. Returning 1.')
        return 1

    E_fit, xsec_fits = make_fits(E, xsec)
            
#    if is_hedis:
#        dict_xml = sort_DIS_genie_hedis(E, xsec)
#        dict_fits_xml = sort_fits_DIS_genie_hedis(E_fit, xsec_fits, E, xsec)
#    else:
#        dict_xml = sort_DIS_regular_genie(E, xsec)
#        dict_fits_xml = sort_fits_DIS_regular_genie(E_fit, xsec_fits, E, xsec)
#            
    dict_fits_xml = sum_fits_all_regular_genie(E_fit, xsec_fits, E, xsec)
        
    return dict_fits_xml
        
def main():
    args = get_args()
    xml_file = args.file
    print('opening file:', xml_file)
    dict_fits_xml = conv_xml_to_dict(xml_file)

    outname_fits = xml_file[:xml_file.find('.xml')] + '_fits_SUM.pckl'
    with open(outname_fits, 'wb') as pfile:
        pickle.dump(dict_fits_xml, pfile)
        
if __name__=='__main__':
    main()
