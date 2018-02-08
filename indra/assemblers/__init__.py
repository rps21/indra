try:
    from indra.assemblers.pysb_assembler import PysbAssembler
except ImportError:
    pass
try:
    from indra.assemblers.graph_assembler import GraphAssembler
except ImportError:
    pass
try:
    from indra.assemblers.sif_assembler import SifAssembler
except ImportError:
    pass
try:
    from indra.assemblers.cx_assembler import CxAssembler
except ImportError:
    pass
try:
    from indra.assemblers.english_assembler import EnglishAssembler
except ImportError:
    pass
try:
    from indra.assemblers.sbgn_assembler import SBGNAssembler
except ImportError:
    pass
try:
    from indra.assemblers.index_card_assembler import IndexCardAssembler
except ImportError:
    pass
try:
    from indra.assemblers.cyjs_assembler import CyJSAssembler
except ImportError:
    pass
try:
    from indra.assemblers.kami_assembler import KamiAssembler
except ImportError:
    pass
try:
    from indra.assemblers.pybel_assembler import PybelAssembler
except ImportError:
    pass
try:
    from indra.assemblers.figaro_assembler import FigaroAssembler
except ImportError:
    pass

try:
    from indra.assemblers.cag_assembler import CAGAssembler
except ImportError:
    pass
