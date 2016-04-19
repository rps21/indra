from indra.literature import literature_client

def test_id_lookup_pmid():
    ids = literature_client.id_lookup('25361007')
    assert(ids['doi'] == '10.18632/oncotarget.2555')
    assert(ids['pmid'] == '25361007')
    assert(ids['pmcid'] == 'PMC4322985')

def test_id_lookup_pmid2():
    ids = literature_client.id_lookup('PMID25361007')
    assert(ids['doi'] == '10.18632/oncotarget.2555')
    assert(ids['pmid'] == '25361007')
    assert(ids['pmcid'] == 'PMC4322985')

def test_id_lookup_pmcid():
    ids = literature_client.id_lookup('PMC4322985')
    assert(ids['doi'] == '10.18632/oncotarget.2555')
    assert(ids['pmid'] == '25361007')
    assert(ids['pmcid'] == 'PMC4322985')

def test_id_lookup_doi():
    ids = literature_client.id_lookup('10.18632/oncotarget.2555')
    assert(ids['doi'] == '10.18632/oncotarget.2555')
    assert(ids['pmid'] == '25361007')
    assert(ids['pmcid'] == 'PMC4322985')

def test_id_lookup_doi2():
    ids = literature_client.id_lookup('DOI10.18632/oncotarget.2555')
    assert(ids['doi'] == '10.18632/oncotarget.2555')
    assert(ids['pmid'] == '25361007')
    assert(ids['pmcid'] == 'PMC4322985')

def test_get_full_text_pmc():
    txt, txt_format = literature_client.get_full_text('PMC4322985')
    assert(txt_format == 'nxml')
    assert(len(txt) > 300000)

def test_get_full_text_doi():
    txt, txt_format = literature_client.get_full_text('10.18632/oncotarget.2555')
    assert(txt_format == 'nxml')
    assert(len(txt) > 300000)

def test_get_full_text_pubmed_abstract():
    txt, txt_format = literature_client.get_full_text('27075779')
    assert(txt_format == 'abstract')
    assert(len(txt) > 800)