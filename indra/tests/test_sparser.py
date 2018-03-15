import json
from indra.sources import sparser
from indra.sources.sparser.processor import _fix_agent
from indra.statements import Agent, Phosphorylation

# ############################
# XML processing tests
# ############################

def test_invalid_xml():
    sp = sparser.process_xml('xyz')
    assert(sp is None)


def test_phosphorylation():
    sp = sparser.process_xml(xml_str1)
    assert(len(sp.statements) == 1)
    assert(sp.statements[0].enz.name == 'MAP2K1')
    assert(sp.statements[0].sub.name == 'MAPK3')
    assert(len(sp.statements[0].evidence) == 1)
    ev = sp.statements[0].evidence[0]
    assert(ev.pmid == '54321')
    assert(ev.text)
    assert(ev.source_api == 'sparser')


def test_phosphorylation2():
    sp = sparser.process_xml(xml_str2)
    assert(len(sp.statements) == 1)
    assert(sp.statements[0].enz.name == 'MAPK1')
    assert(sp.statements[0].sub.name == 'TP53BP2')
    assert(sp.statements[0].residue == 'S')
    assert(sp.statements[0].position == '827')
    assert (len(sp.statements[0].evidence) == 1)
    ev = sp.statements[0].evidence[0]
    assert (ev.pmid == '12345')
    assert (ev.text)
    assert (ev.source_api == 'sparser')


def test_fix_agent_be_name():
    a = Agent('XXX', db_refs={'FPLX': 'CDK'})
    _fix_agent(a)
    assert(a.name == 'CDK')


def test_fix_agent_hgnc_only():
    a = Agent('XXX', db_refs={'HGNC': '7199'})
    _fix_agent(a)
    assert(a.name == 'MOS')
    assert(a.db_refs.get('UP') == 'P00540')


def test_fix_agent_fa_only():
    a = Agent('XXX', db_refs={'FA': '00815'})
    _fix_agent(a)
    assert(a.name == 'Cyclin')
    assert(a.db_refs.get('FPLX') == 'Cyclin')
    assert(a.db_refs.get('NXPFA') == '00815')
    assert('FA' not in a.db_refs)


def test_fix_agent_ncit_only():
    a = Agent('XXX', db_refs={'NCIT': 'C25785'})
    _fix_agent(a)
    assert(a.name == 'KRAS')
    assert(a.db_refs.get('HGNC') == '6407')
    assert(a.db_refs.get('UP') == 'P01116')


def test_fix_agent_ncit_only():
    a = Agent('XXX', db_refs={'NCIT': 'C104166'})
    _fix_agent(a)
    assert(a.name == 'TUBB')
    assert(a.db_refs.get('FPLX') == 'TUBB')


# ############################
# JSON processing tests
# ############################

def test_process_json_str():
    sp = sparser.process_json_dict(json.loads(json_str1))
    assert sp is not None
    assert len(sp.statements) == 1
    assert isinstance(sp.statements[0], Phosphorylation)
    sp.set_statements_pmid('1234567')
    assert sp.statements[0].evidence[0].pmid == '1234567'
    assert sp.json_stmts[0]['evidence'][0]['pmid'] == '1234567'


xml_str1 = '''
<article pmid="54321">
 <interpretation>
 <sentence-text>MEK1 phosphorylates ERK1</sentence-text>
 <sem>
     <ref category="phosphorylate">
         <var name="agent">
         <ref category="protein">
             <var name="name">MP2K1_HUMAN</var>
             <var name="uid">UP:MP2K1_HUMAN</var>
         </ref>
         </var>
         <var name="substrate">
            <ref category="protein">
                <var name="name">MK03_HUMAN</var>
                <var name="uid">UP:MK03_HUMAN</var>
            </ref>
         </var>
     <var name="present"><ref category="present"></ref></var>
     </ref>
 </sem>
</interpretation>
</article>
'''

xml_str2 = '''
<article pmid="12345">
<interpretation>
  <sentence-text>Hence ASPP2 can be phosphorylated at serine 827 by MAPK1 in vitro</sentence-text>
  <sem>
    <ref category="phosphorylate">
      <var name="subordinate-conjunction">
          <ref category="subordinate-conjunction"><var name="word">hence</var></ref></var>
      <var name="substrate">
          <ref category="protein">
              <var name="name">ASPP2_HUMAN</var>
              <var name="uid">UP:ASPP2_HUMAN</var>
          </ref>
      </var>
      <var name="agent">
        <ref category="protein">
          <var name="context">
            <ref category="in-vitro"></ref>
          </var>
          <var name="uid">UP:MK01_HUMAN</var>
          <var name="name">MK01_HUMAN</var>
        </ref>
      </var>
      <var name="site">
        <ref category="residue-on-protein">
          <var name="amino-acid">
            <ref category="amino-acid"><var name="name">serine</var></ref>
          </var>
          <var name="position"> 827</var>
        </ref>
      </var>
      <var name="modal"><ref category="can"></ref></var>
    </ref>
  </sem>
</interpretation>
</article>
'''

json_str1 = '''
[
 {
  "type": "Phosphorylation",
  "evidence": [
  {
    "source_api": "sparser",
    "text": "MEK phosphorylates ERK",
    "pmid": "PMC_3500"}],
  "sub": {
    "name": "ERK",
    "db_refs": {
      "NCIT": "C26360",
      "TEXT": "ERK"},
    "TEXT": "ERK"},
  "enz": {
    "name": "MEK",
    "db_refs": {
      "FPLX": "MEK",
      "TEXT": "MEK"},
    "TEXT": "MEK"}
 }
]'''
