----> 1 ml.gather_implicit_activities()

/home/bobby/Dropbox/Sorger_Lab/indra/indra/mechlinker/__init__.py in gather_implicit_activities(self)
     86                     enz_base = self._get_base(stmt.enz)
     87                     enz_base.add_activity('kinase')
---> 88                     enz_base.add_active_state('kinase', stmt.enz.mods)
     89             elif isinstance(stmt, Dephosphorylation):
     90                 if stmt.enz is not None:

/home/bobby/Dropbox/Sorger_Lab/indra/indra/mechlinker/__init__.py in add_active_state(self, activity_type, agent)
    632 
    633     def add_active_state(self, activity_type, agent):
--> 634         agent_state = AgentState(agent)
    635         if activity_type in self.active_states:
    636             self.active_states[activity_type].append(agent_state)

/home/bobby/Dropbox/Sorger_Lab/indra/indra/mechlinker/__init__.py in __init__(self, agent)
    697     
    698     def __init__(self, agent):
--> 699         self.bound_conditions = agent.bound_conditions
    700         self.mods = agent.mods
    701         self.mutations = agent.mutations


enz_base.add_active_state('kinase', stmt.enz.mods) #in mech_linger, adding active state with activity type ('kinase') and list of modifications on the enzyme (stmts.enz.mods)

/home/bobby/Dropbox/Sorger_Lab/indra/indra/mechlinker/__init__.py in add_active_state(self, activity_type, agent)  #add_active_state method expects activity type and agent. This looks like a problem
 agent_state = AgentState(agent) #within add_active_state call AgentState with agent. Same problem, seem to be calling with list of modifications instead of agent 


 self.bound_conditions = agent.bound_conditions #setting bound conditions based on bound conditions of agent. Seems like if we went back to implicit activities and made call with stmts.enz instead of stmts.enz.mods this may work 
#Still have to see exactly what the output is if it's called with an agent instead of a list and theres not error 


#############
#EDITED IMPLICIT ACTIVITIES FUNCTION DUE TO ERRORS, WILL NEED TO SEE IF THIS IS CORRECT OR I MISSED SOMETHING
#aLSO SHOULD LOOK MORE INTO HOW THIS WORKS AND IF ITS USEFUL, DOESN'T SEEM TO DO MUCH NOW
