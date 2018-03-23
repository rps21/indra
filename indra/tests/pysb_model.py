# exported from PySB model 'None'

from pysb import Model, Monomer, Parameter, Expression, Compartment, Rule, Observable, Initial, MatchOnce, Annotation, ANY, WILD

Model()

Monomer('MDM2')
Monomer('TP53')

Parameter('kf_tm_synth_1', 4.0)
Parameter('Ka_tm_synth_1', 10000.0)
Parameter('n_tm_synth_1', 1.0)
Parameter('MDM2_0', 10000.0)
Parameter('TP53_0', 10000.0)

Observable('TP53_synthesizes_MDM2_subj_obs', TP53())

Expression('TP53_synthesizes_MDM2_rate', kf_tm_synth_1*TP53_synthesizes_MDM2_subj_obs**(-1 + n_tm_synth_1)*(TP53_synthesizes_MDM2_subj_obs**n_tm_synth_1 + Ka_tm_synth_1**n_tm_synth_1)**(-1))

Rule('TP53_synthesizes_MDM2', TP53() >> TP53() + MDM2(), TP53_synthesizes_MDM2_rate)

Initial(MDM2(), MDM2_0)
Initial(TP53(), TP53_0)

Annotation(TP53_synthesizes_MDM2, '7f46e853-c87d-4e6d-aa2f-3da7be4292cd', 'from_indra_statement')
Annotation(TP53_synthesizes_MDM2, 'MDM2', 'rule_has_object')
Annotation(TP53_synthesizes_MDM2, 'TP53', 'rule_has_subject')

