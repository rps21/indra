from indra.preassembler import Preassembler
from indra.preassembler.hierarchy_manager import hierarchies

def get_is_direct(stmt):
    '''Returns true if there is evidence that the statement is a direct
    interaction. If any of the evidences associated with the statement
    indicates a direct interatcion then we assume the interaction
    is direct. If there is no evidence for the interaction being indirect
    then we default to direct.'''
    any_indirect = False
    for ev in stmt.evidence:
        if ev.epistemics.get('direct') is True:
            return True
        elif ev.epistemics.get('direct') is False:
            # This guarantees that we have seen at least
            # some evidence that the statement is indirect
            any_indirect = True
    if any_indirect:
        return False
    return True


def filter_direct(stmts):
    direct_stmts = []
    for stmt in stmts:
        if get_is_direct(stmt):
            direct_stmts.append(stmt)
    return direct_stmts

    pa = Preassembler(hierarchies, stmts_valid)
    stmts_unique = pa.combine_duplicates()
    logger.info('%d unique phosphorylations in RAS Machine' % len(stmts_unique))

    stmts_unique = pa.combine_related()
    logger.info('%d top-level phosphorylations in RAS Machine' % len(stmts_unique))
