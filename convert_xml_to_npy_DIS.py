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

def sort_fits_DIS_regular_genie(E_fit, xsec_fits, E, xsec):

    genie_xml = {'E_fit':E_fit, 'xsec_fits':[], 'full_name':[], 'E_init':[], 'xsec_init':[],
                 'pdg':[], 'tgt':[], 'hitN':[], 'hitqrk':[], 'sea':[], 'CC':[], 'dis':[], 'charm_incl':[]}

    for k in xsec_fits.keys():
        
        if 'QPMDISPXSec/Default/nu:' in k:
            model = 'QPMDISPXSec/Default/nu:'
        elif 'AivazisCharmPXSecLO/CC-Default/nu:' in k:
            model = 'AivazisCharmPXSecLO/CC-Default/nu:'
        else:
            continue

        genie_xml['full_name'].append(k)
        genie_xml['xsec_fits'].append(xsec_fits[k])

        genie_xml['E_init'].append(E[k])
        genie_xml['xsec_init'].append(xsec[k])
        
        st_pdg = k.find(model) + len(model)
        end_pdg = k.find(';tgt:')
        pdg = int(k[st_pdg:end_pdg])
        genie_xml['pdg'].append(pdg)
        print pdg

        st_tgt = end_pdg +len(';tgt:')
        end_tgt = k.find(';N:')
        tgt = int(k[st_tgt:end_tgt])
        genie_xml['tgt'].append(tgt)
        print tgt

        st_N = end_tgt +len(';N:')
        end_N = k.find(';q:')
        hitN = int(k[st_N:end_N])
        genie_xml['hitN'].append(hitN)
        print hitN

        st_q = end_N +len(';q:')
        end_q = k.find('(')
        hitqrk = int(k[st_q:end_q])
        genie_xml['hitqrk'].append(hitqrk)
        print hitqrk

        st_sea = end_q + 1
        end_sea = k.find(')')
        if k[st_sea:end_sea] == 's':
            sea = 1
        elif k[st_sea:end_sea] == 'v':
            sea = 0
        else:
            sea = -1
        genie_xml['sea'].append(sea)
        print k[st_sea:end_sea], sea

        st_curr = k.find('proc:Weak[') + len('proc:Weak[')
        end_curr = k.find(']')
        if k[st_curr:end_curr] == 'CC':
            cc = 1
        elif k[st_curr:end_curr] == 'NC':
            cc = 0
        else:
            cc = -1
        genie_xml['CC'].append(cc)        
        print k[st_curr:end_curr], cc

        st_int = k.find('],') + len('],')
        end_int = st_int + 3 #len(k) - 1
        if k[st_int:end_int] == 'DIS':
            dis = 1
        else:
            dis = 0
        genie_xml['dis'].append(dis)
        print k[st_int:end_int], dis
        
        if end_int == len(k) - 1:
            charm_incl = 0
        else:
            charm_incl = 1
            print '~~~~~~~~~~~~~~~> incl charm:', k[st_int:len(k) - 1]
        genie_xml['charm_incl'].append(charm_incl)
        print 'ch', charm_incl  
    
    return genie_xml

def sort_fits_DIS_genie_hedis(E_fit, xsec_fits, E, xsec):

    hedis_xml = {'E_fit':E_fit, 'xsec_fits':[], 'full_name':[], 'E_init':[], 'xsec_init':[],
                 'pdg':[], 'tgt':[], 'hitN':[], 'hitqrk':[], 'sea':[], 'CC':[], 'int_hedis':[], 'finalqrk':[]}

    for k in xsec_fits.keys():

        if not 'HEDISPXSec/Default/nu:' in k: continue
        
        hedis_xml['full_name'].append(k)
        hedis_xml['xsec_fits'].append(xsec_fits[k])

        hedis_xml['E_init'].append(E[k])
        hedis_xml['xsec_init'].append(xsec[k])

        st_pdg = k.find('HEDISPXSec/Default/nu:') + len('HEDISPXSec/Default/nu:')
        end_pdg = k.find(';tgt:')
        pdg = int(k[st_pdg:end_pdg])
        hedis_xml['pdg'].append(pdg)
        print pdg

        st_tgt = end_pdg +len(';tgt:')
        end_tgt = k.find(';N:')
        tgt = int(k[st_tgt:end_tgt])
        hedis_xml['tgt'].append(tgt)
        print tgt

        st_N = end_tgt +len(';N:')
        end_N = k.find(';q:')
        hitN = int(k[st_N:end_N])
        hedis_xml['hitN'].append(hitN)
        print hitN

        st_q = end_N +len(';q:')
        end_q = k.find('(')
        hitqrk = int(k[st_q:end_q])
        hedis_xml['hitqrk'].append(hitqrk)
        print hitqrk

        st_sea = end_q + 1
        end_sea = k.find(')')
        if k[st_sea:end_sea] == 's':
            sea = 1
        elif k[st_sea:end_sea] == 'v':
            sea = 0
        else:
            sea = -1
        hedis_xml['sea'].append(sea)
        print k[st_sea:end_sea], sea

        st_curr = k.find('proc:Weak[') + len('proc:Weak[')
        end_curr = k.find(']')
        if k[st_curr:end_curr] == 'CC':
            cc = 1
        elif k[st_curr:end_curr] == 'NC':
            cc = 0
        else:
            cc = -1
        hedis_xml['CC'].append(cc)        
        print k[st_curr:end_curr], cc

        st_int = k.find('],') + len('],')
        end_int = k.find(';finalquark:')
        if k[st_int:end_int] == 'HEDIS':
            hedis = 1
        else:
            hedis = 0
        hedis_xml['int_hedis'].append(hedis)
        print k[st_int:end_int], hedis

        st_fq = end_int +len(';finalquark:')
        end_fq = len(k) - 1
        finalqrk = int(k[st_fq:end_fq])
        hedis_xml['finalqrk'].append(finalqrk)
        print finalqrk

    return hedis_xml

def sort_DIS_regular_genie(E, xsec):

    genie_xml = {'E':[], 'xsec':[], 
                 'pdg':[], 'tgt':[], 'hitN':[], 'hitqrk':[], 'sea':[], 'CC':[], 'dis':[], 'charm_incl':[]}

    for k in E.keys():
        
        if 'QPMDISPXSec/Default/nu:' in k:
            model = 'QPMDISPXSec/Default/nu:'
        elif 'AivazisCharmPXSecLO/CC-Default/nu:' in k:
            model = 'AivazisCharmPXSecLO/CC-Default/nu:'
        else:
            continue

        l = len(E[k])
        arr_ones = np.ones(l)        
        
        genie_xml['E'] = np.concatenate((genie_xml['E'],E[k]))
        genie_xml['xsec'] = np.concatenate((genie_xml['xsec'],xsec[k]))

        st_pdg = k.find(model) + len(model)
        end_pdg = k.find(';tgt:')
        pdg = int(k[st_pdg:end_pdg])
        genie_xml['pdg'] = np.concatenate((genie_xml['pdg'], arr_ones*pdg))
        print pdg

        st_tgt = end_pdg +len(';tgt:')
        end_tgt = k.find(';N:')
        tgt = int(k[st_tgt:end_tgt])
        genie_xml['tgt'] = np.concatenate((genie_xml['tgt'], arr_ones*tgt))
        print tgt

        st_N = end_tgt +len(';N:')
        end_N = k.find(';q:')
        hitN = int(k[st_N:end_N])
        genie_xml['hitN'] = np.concatenate((genie_xml['hitN'], arr_ones*hitN))
        print hitN

        st_q = end_N +len(';q:')
        end_q = k.find('(')
        hitqrk = int(k[st_q:end_q])
        genie_xml['hitqrk'] = np.concatenate((genie_xml['hitqrk'], arr_ones*hitqrk))
        print hitqrk

        st_sea = end_q + 1
        end_sea = k.find(')')
        if k[st_sea:end_sea] == 's':
            sea = 1
        elif k[st_sea:end_sea] == 'v':
            sea = 0
        else:
            sea = -1
        genie_xml['sea'] = np.concatenate((genie_xml['sea'], arr_ones*sea))
        print k[st_sea:end_sea], sea

        st_curr = k.find('proc:Weak[') + len('proc:Weak[')
        end_curr = k.find(']')
        if k[st_curr:end_curr] == 'CC':
            cc = 1
        elif k[st_curr:end_curr] == 'NC':
            cc = 0
        else:
            cc = -1
        genie_xml['CC'] = np.concatenate((genie_xml['CC'], arr_ones*cc))        
        print k[st_curr:end_curr], cc

        st_int = k.find('],') + len('],')
        end_int = st_int + 3 #len(k) - 1
        if k[st_int:end_int] == 'DIS':
            dis = 1
        else:
            dis = 0
        genie_xml['dis'] = np.concatenate((genie_xml['dis'], arr_ones*dis))
        print k[st_int:end_int], dis
        
        if end_int == len(k) - 1:
            charm_incl = 0
        else:
            charm_incl = 1
            print '~~~~~~~~~~~~~~~> incl charm:', k[st_int:len(k) - 1]
        genie_xml['charm_incl'] = np.concatenate((genie_xml['charm_incl'], arr_ones*charm_incl))
        print 'ch', charm_incl  
    
    return genie_xml

def sort_DIS_genie_hedis(E, xsec):

    hedis_xml = {'E':[], 'xsec':[], 
                 'pdg':[], 'tgt':[], 'hitN':[], 'hitqrk':[], 'sea':[], 'CC':[], 'int_hedis':[], 'finalqrk':[]}

    for k in E.keys():

        if not 'genie::HEDISPXSec/Default/nu:' in k: continue
        
        l = len(E[k])
        arr_ones = np.ones(l)

        hedis_xml['E'] = np.concatenate((hedis_xml['E'],E[k]))
        hedis_xml['xsec'] = np.concatenate((hedis_xml['xsec'],xsec[k]))

        st_pdg = k.find('genie::HEDISPXSec/Default/nu:') + len('genie::HEDISPXSec/Default/nu:')
        end_pdg = k.find(';tgt:')
        pdg = int(k[st_pdg:end_pdg])
        hedis_xml['pdg'] = np.concatenate((hedis_xml['pdg'], arr_ones*pdg))
        print pdg

        st_tgt = end_pdg +len(';tgt:')
        end_tgt = k.find(';N:')
        tgt = int(k[st_tgt:end_tgt])
        hedis_xml['tgt'] = np.concatenate((hedis_xml['tgt'], arr_ones*tgt))
        print tgt

        st_N = end_tgt +len(';N:')
        end_N = k.find(';q:')
        hitN = int(k[st_N:end_N])
        hedis_xml['hitN'] = np.concatenate((hedis_xml['hitN'], arr_ones*hitN))
        print hitN

        st_q = end_N +len(';q:')
        end_q = k.find('(')
        hitqrk = int(k[st_q:end_q])
        hedis_xml['hitqrk'] = np.concatenate((hedis_xml['hitqrk'], arr_ones*hitqrk))
        print hitqrk

        st_sea = end_q + 1
        end_sea = k.find(')')
        if k[st_sea:end_sea] == 's':
            sea = 1
        elif k[st_sea:end_sea] == 'v':
            sea = 0
        else:
            sea = -1
        hedis_xml['sea'] = np.concatenate((hedis_xml['sea'], arr_ones*sea))
        print k[st_sea:end_sea], sea

        st_curr = k.find('proc:Weak[') + len('proc:Weak[')
        end_curr = k.find(']')
        if k[st_curr:end_curr] == 'CC':
            cc = 1
        elif k[st_curr:end_curr] == 'NC':
            cc = 0
        else:
            cc = -1
        hedis_xml['CC'] = np.concatenate((hedis_xml['CC'], arr_ones*cc))        
        print k[st_curr:end_curr], cc

        st_int = k.find('],') + len('],')
        end_int = k.find(';finalquark:')
        if k[st_int:end_int] == 'HEDIS':
            hedis = 1
        else:
            hedis = 0
        hedis_xml['int_hedis'] = np.concatenate((hedis_xml['int_hedis'], arr_ones*hedis))
        print k[st_int:end_int], hedis

        st_fq = end_int +len(';finalquark:')
        end_fq = len(k) - 1
        finalqrk = int(k[st_fq:end_fq])
        hedis_xml['finalqrk'] = np.concatenate((hedis_xml['finalqrk'], arr_ones*finalqrk))
        print finalqrk

    return hedis_xml

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
            
    if is_hedis:
        dict_xml = sort_DIS_genie_hedis(E, xsec)
        dict_fits_xml = sort_fits_DIS_genie_hedis(E_fit, xsec_fits, E, xsec)
    else:
        dict_xml = sort_DIS_regular_genie(E, xsec)
        dict_fits_xml = sort_fits_DIS_regular_genie(E_fit, xsec_fits, E, xsec)
            

    return dict_xml, dict_fits_xml
        
def main():
    args = get_args()
    xml_file = args.file
    print('opening file:', xml_file)
    dict_xml, dict_fits_xml = conv_xml_to_dict(xml_file)

    outname = xml_file[:xml_file.find('.xml')] + '.pckl'
    with open(outname, 'wb') as pfile:
        pickle.dump(dict_xml, pfile)

    outname_fits = xml_file[:xml_file.find('.xml')] + '_fits.pckl'
    with open(outname_fits, 'wb') as pfile:
        pickle.dump(dict_fits_xml, pfile)
        
if __name__=='__main__':
    main()
