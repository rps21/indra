from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import os
import re
import csv
import xml.etree.ElementTree as ET
import requests
# Python3
try:
    from functools import lru_cache
# Python2
except ImportError:
    from functools32 import lru_cache
from indra.util import read_unicode_csv, UnicodeXMLTreeBuilder as UTB

hgnc_url = 'http://rest.genenames.org/fetch/'


def get_uniprot_id(hgnc_id):
    """Return the UniProt ID corresponding to the given HGNC ID.

    Parameters
    ----------
    hgnc_id : str
        The HGNC ID to be converted. Note that the HGNC ID is a number that is
        passed as a string. It is not the same as the HGNC gene symbol.

    Returns
    -------
    uniprot_id : str
        The UniProt ID corresponding to the given HGNC ID.
    """
    uniprot_id = uniprot_ids.get(hgnc_id)
    # The lookup can yield an empty string. Instead return None.
    if not uniprot_id:
        return None
    return uniprot_id

def get_entrez_id(hgnc_id):
    """Return the Entrez ID corresponding to the given HGNC ID.

    Parameters
    ----------
    hgnc_id : str
        The HGNC ID to be converted. Note that the HGNC ID is a number that is
        passed as a string. It is not the same as the HGNC gene symbol.

    Returns
    -------
    entrez_id : str
        The Entrez ID corresponding to the given HGNC ID.
    """
    entrez_id = entrez_ids.get(hgnc_id)
    # The lookup can yield an empty string. Instead return None.
    if not entrez_id:
        return None
    return entrez_id

def get_hgnc_from_entrez(entrez_id):
    """Return the HGNC ID corresponding to the given Entrez ID.

    Parameters
    ----------
    entrez_id : str
        The EntrezC ID to be converted, a number passed as a strig.

    Returns
    -------
    hgnc_id : str
        The HGNC ID corresponding to the given Entrez ID.
    """
    hgnc_id = entrez_ids_reverse.get(entrez_id)
    return hgnc_id

def get_hgnc_name(hgnc_id):
    """Return the HGNC symbol corresponding to the given HGNC ID.

    Parameters
    ----------
    hgnc_id : str
        The HGNC ID to be converted.

    Returns
    -------
    hgnc_name : str
        The HGNC symbol corresponding to the given HGNC ID.
    """
    try:
        hgnc_name = hgnc_names[hgnc_id]
    except KeyError:
        xml_tree = get_hgnc_entry(hgnc_id)
        if xml_tree is None:
            return None
        hgnc_name_tag =\
            xml_tree.find("result/doc/str[@name='symbol']")
        if hgnc_name_tag is None:
            return None
        hgnc_name = hgnc_name_tag.text.strip()
    return hgnc_name

def get_hgnc_id(hgnc_name):
    """Return the HGNC ID corresponding to the given HGNC symbol.

    Parameters
    ----------
    hgnc_name : str
        The HGNC symbol to be converted. Example: BRAF

    Returns
    -------
    hgnc_id : str
        The HGNC ID corresponding to the given HGNC symbol.
    """
    return hgnc_ids.get(hgnc_name)

def get_hgnc_from_mouse(mgi_id):
    """Return the HGNC ID corresponding to the given MGI mouse gene ID.

    Parameters
    ----------
    mgi_id : str
        The MGI ID to be converted. Example: "2444934"

    Returns
    -------
    hgnc_id : str
        The HGNC ID corresponding to the given MGI ID.
    """
    if mgi_id.startswith('MGI:'):
        mgi_id = mgi_id[4:]
    return mouse_map.get(mgi_id)

def get_hgnc_from_rat(rgd_id):
    """Return the HGNC ID corresponding to the given RGD rat gene ID.

    Parameters
    ----------
    rgd_id : str
        The RGD ID to be converted. Example: "1564928"

    Returns
    -------
    hgnc_id : str
        The HGNC ID corresponding to the given RGD ID.
    """
    if rgd_id.startswith('RGD:'):
        rgd_id = rgd_id[4:]
    return rat_map.get(rgd_id)


def get_rat_id(hgnc_id):
    """Return the RGD rat ID corresponding to the given HGNC ID.

    Parameters
    ----------
    hgnc_id : str
        The HGNC ID to be converted. Example: ""

    Returns
    -------
    rgd_id : str
        The RGD ID corresponding to the given HGNC ID.
    """
    for k, v in rat_map.items():
        if v == hgnc_id:
            return k

def get_mouse_id(hgnc_id):
    """Return the MGI mouse ID corresponding to the given HGNC ID.

    Parameters
    ----------
    hgnc_id : str
        The HGNC ID to be converted. Example: ""

    Returns
    -------
    mgi_id : str
        The MGI ID corresponding to the given HGNC ID.
    """
    for k, v in mouse_map.items():
        if v == hgnc_id:
            return k

@lru_cache(maxsize=1000)
def get_hgnc_entry(hgnc_id):
    """Return the HGNC entry for the given HGNC ID from the web service.

    Parameters
    ----------
    hgnc_id : str
        The HGNC ID to be converted.

    Returns
    -------
    xml_tree : ElementTree
        The XML ElementTree corresponding to the entry for the
        given HGNC ID.
    """
    url = hgnc_url + 'hgnc_id/%s' % hgnc_id
    headers = {'Accept': '*/*'}
    res = requests.get(url, headers=headers)
    if not res.status_code == 200:
        return None
    xml_tree = ET.XML(res.content, parser=UTB())
    return xml_tree


def _read_hgnc_maps():
    hgnc_file = os.path.dirname(os.path.abspath(__file__)) + \
                '/../resources/hgnc_entries.tsv'
    csv_rows = read_unicode_csv(hgnc_file, delimiter='\t', encoding='utf-8')
    hgnc_names = {}
    hgnc_ids = {}
    hgnc_withdrawn = []
    uniprot_ids = {}
    entrez_ids = {}
    entrez_ids_reverse = {}
    mouse_map = {}
    rat_map = {}
    for row in csv_rows:
        hgnc_id = row[0][5:]
        hgnc_status = row[3]
        if hgnc_status == 'Approved':
            hgnc_name = row[1]
            hgnc_names[hgnc_id] = hgnc_name
            hgnc_ids[hgnc_name] = hgnc_id
        elif hgnc_status == 'Symbol Withdrawn':
            descr = row[2]
            m = re.match(r'symbol withdrawn, see ([^ ]*)', descr)
            new_name = m.groups()[0]
            hgnc_withdrawn.append(hgnc_id)
            hgnc_names[hgnc_id] = new_name
        # Uniprot
        uniprot_id = row[6]
        uniprot_ids[hgnc_id] = uniprot_id
        # Entrez
        entrez_id = row[5]
        entrez_ids[hgnc_id] = entrez_id
        entrez_ids_reverse[entrez_id] = hgnc_id
        # Mouse
        mgi_id = row[7]
        if mgi_id:
            mgi_ids = mgi_id.split(', ')
            for mgi_id in mgi_ids:
                if mgi_id.startswith('MGI:'):
                    mgi_id = mgi_id[4:]
                mouse_map[mgi_id] = hgnc_id
        # Rat
        rgd_id = row[8]
        if rgd_id:
            rgd_ids = rgd_id.split(', ')
            for rgd_id in rgd_ids:
                if rgd_id.startswith('RGD:'):
                    rgd_id = rgd_id[4:]
                rat_map[rgd_id] = hgnc_id
    return (hgnc_names, hgnc_ids, hgnc_withdrawn,
            uniprot_ids, entrez_ids, entrez_ids_reverse, mouse_map, rat_map)

(hgnc_names, hgnc_ids, hgnc_withdrawn, uniprot_ids, entrez_ids,
 entrez_ids_reverse, mouse_map, rat_map) = \
    _read_hgnc_maps()
