#!/usr/bin/python
#! python3.7
#try to specify that will use python version 3.7
__author__      = "Wheaton Schroeder"

#Written by; Wheaton L Schroeder
#Latest Version: 07/13/2020

#Written to test the OptFill library which I have built

import cobra
from cobra import Model, Reaction, Metabolite
from ConnectingProblems import ConnectingProblems
from TICFindingProblem import TICFindingProblem
import re

#read the model with cobrapy
model = cobra.io.read_sbml_model("test_model_1.xml")
database = cobra.io.read_sbml_model("test_db_1.xml")

#find exchange reactions, binary value
ex_rxns = list()

for rxn in model.reactions:
    
    #in this model exchange reactions are noted by "ex" and the end of the identifer
    if bool(re.search('ex$',rxn.id)):
    
        #add exchange reaction to the list
        ex_rxns.extend([rxn.id])
    
#by this point we should have the full list of exchange reactions model, and database
#now hopefully we can test the TIC-Finding Problem
TFP = TICFindingProblem(model,database,ex_rxns)

#solve TIC-finding
TFP_solution = TFP.find(epsilon=0.001)

#retrieve the TFP solutions so that can pass to connecting problems
prev_alphas = TFP_solution.prev_alphas
prev_betas = TFP_solution.prev_betas
solution_list = TFP_solution.solution_numbers

#now need to solve connecting problems
#initialize connecting problems
#note that must prod should match a string name of one or more metabolites
CP = ConnectingProblems(model,database,["biomass_wt"],prev_alphas,prev_betas,solution_list)

#solve connecting problems
CP_solution = CP.fill(iterations=25)