#!/usr/bin/python

#Written by; Wheaton L Schroeder
#Latest Version: 07/13/2020

#don't know if I need all these, gapfill used them or they are required for things
#that I am trying to do here, so importing them for now
from __future__ import absolute_import

from optlang.interface import OPTIMAL, FEASIBLE, INFEASIBLE, ITERATION_LIMIT, NUMERIC, OPTIMAL, SUBOPTIMAL, TIME_LIMIT
from optlang.symbolics import Zero, add

import cobra
import sys
import warnings
import time
import re

from cobra.core import Model, Reaction, Metabolite, Solution
from cobra.util import fix_objective_as_constraint
from copy import copy, deepcopy


#so I think how this needs to go is that I define a class to perform optfill
#since OptFill is four seperate problems (1 tic-finding, 3 connecting), probably need four seperate classes
#in programming, classes are defined as something which creates an object

#class for solving the Connecting Problems
class ConnectingProblems(object):

    def __init__(self,model,database,must_prod,TFP_alphas,TFP_betas,TFPsol_list):
    
        #model is the model with the holes
        #database is reactions with which the holes will be filled
        #must_prod is the list of metabolites which must be produced in the CPs solutions
        
        self.original_model = model.copy()              #stores the input, unfilled model
        self.database = database.copy()                 #stores the input database
        self.must_prod = must_prod                      #list of metabolites which must be produced
        self.solution_numbers = list()                  #list of TFP solution numbers
        self.num_mets = list()                          #list of TFP objective values
        self.num_rxns = list()                          #list of TFP objective values
        self.num_rev = list()                           #list of TFP objective values
        self.fluxes = dict()                            #list of flux values for each solution, each solutions' flux values are in a nested list indexed by reaction
        self.deltas_1 = dict()                          #list of each solutions' delta values are in a nested list indexed by reaction
        self.rhos_1 = dict()                            #list of each solutions' rho values are in a nested list indexed by reaction
        self.omegas_1 = dict()                          #list of each solutions' omega values are in a nested list indexed by reaction
        self.sigmas_1 = dict()                          #list of each solutions' sigma values are in a nested list indexed by reaction
        self.taus_1 = dict()                            #list of each solutions' tau values are in a nested list indexed by reaction
        self.xs_1 = dict()                              #list of each solutions' x values are in a nested list indexed by reaction
        self.zetas_1 = dict()                           #list of each solutions' zeta values are in a nested list indexed by reaction
        self.thetas_1 = dict()                          #list of each solutions' theta values are in a nested list indexed by reaction
        self.lambdas_1 = dict()                         #list of each solutions' lambda values are in a nested list indexed by reaction
        self.deltas_2 = dict()                          #list of each solutions' delta values are in a nested list indexed by reaction
        self.rhos_2 = dict()                            #list of each solutions' rho values are in a nested list indexed by reaction
        self.omegas_2 = dict()                          #list of each solutions' omega values are in a nested list indexed by reaction
        self.sigmas_2 = dict()                          #list of each solutions' sigma values are in a nested list indexed by reaction
        self.taus_2 = dict()                            #list of each solutions' tau values are in a nested list indexed by reaction
        self.xs_2 = dict()                              #list of each solutions' x values are in a nested list indexed by reaction
        self.zetas_2 = dict()                           #list of each solutions' zeta values are in a nested list indexed by reaction
        self.thetas_2 = dict()                          #list of each solutions' theta values are in a nested list indexed by reaction
        self.lambdas_2 = dict()                         #list of each solutions' lambda values are in a nested list indexed by reaction
        self.deltas_3 = dict()                          #list of each solutions' delta values are in a nested list indexed by reaction
        self.rhos_3 = dict()                            #list of each solutions' rho values are in a nested list indexed by reaction
        self.omegas_3 = dict()                          #list of each solutions' omega values are in a nested list indexed by reaction
        self.sigmas_3 = dict()                          #list of each solutions' sigma values are in a nested list indexed by reaction
        self.taus_3 = dict()                            #list of each solutions' tau values are in a nested list indexed by reaction
        self.xs_3 = dict()                              #list of each solutions' x values are in a nested list indexed by reaction
        self.zetas_3 = dict()                           #list of each solutions' zeta values are in a nested list indexed by reaction
        self.thetas_3 = dict()                          #list of each solutions' theta values are in a nested list indexed by reaction
        self.lambdas_3 = dict()                         #list of each solutions' lambda values are in a nested list indexed by reaction
        self.xis = dict()                               #stores 1 if metabolite i has positive stoichiometry, 0 otherwise
        self.psis = dict()                              #stores 1 if metabolite i has negative stoichiometry, 0 otherwise
        self.TFP_alphas = TFP_alphas                    #stores reactions participating in TFP solutions in the forward direction indexed first by solution number then by reaction
        self.TFP_betas = TFP_betas                      #stores reactions participating in TFP solutions in the backward direction indexed first by solution number then by reaction
        self.TFPsol_list = TFPsol_list                  #list of previous solutions

        #instead of using the native merge function, will add an "origin" property to the combined model by iterating through them
        #first we can start with adding all model reactions, defining their origin as "model"
        #according to documentation setattr() should create an attribute if it does not already exist
        self.combined_model_1 = self.original_model.copy()            #combination of the model and the database used to solve CP1
        self.combined_model_2 = self.original_model.copy()            #combination of the model and the database used to solve CP2
        self.combined_model_3 = self.original_model.copy()            #combination of the model and the database used to solve CP3
        
        for rxn in self.combined_model_1.reactions:
        
            #add origin attribute
            setattr(rxn,'origin','model')
            
        for rxn in self.combined_model_2.reactions:
        
            #add origin attribute
            setattr(rxn,'origin','model')
            
        for rxn in self.combined_model_3.reactions:
        
            #add origin attribute
            setattr(rxn,'origin','model')
            
        #make origins for metabolites
        for mets in self.combined_model_1.metabolites:
        
            #add origin attribute
            setattr(mets,'origin','model')
            
        for mets in self.combined_model_2.metabolites:
        
            #add origin attribute
            setattr(mets,'origin','model')
            
        for mets in self.combined_model_3.metabolites:
        
            #add origin attribute
            setattr(mets,'origin','model')
            
        #add database reactions
        self.combined_model_1.merge(self.database.copy())           #add the database reactions to the combination model 1
        self.combined_model_2.merge(self.database.copy())           #add the database reactions to the combination model 2
        self.combined_model_3.merge(self.database.copy())           #add the database reactions to the combination model 3
            
        #all model reactions already have the 'origin' attribute, so reactions that  do not yet have this attribute are database reactions
        for rxn in self.combined_model_1.reactions:
        
            if not hasattr(rxn,'origin'):
            
                #we found a database reaction
                setattr(rxn,'origin','database')
                
            #otherwise do nothing, we already have this attribute
        
        #repeat for other combined models
        for rxn in self.combined_model_2.reactions:
        
            if not hasattr(rxn,'origin'):
            
                #we found a database reaction
                setattr(rxn,'origin','database')
                
            #otherwise do nothing, we already have this attribute
            
        for rxn in self.combined_model_3.reactions:
        
            if not hasattr(rxn,'origin'):
            
                #we found a database reaction
                setattr(rxn,'origin','database')
                
            #otherwise do nothing, we already have this attribute
        
        #assign origins for metaboltes
        for met in self.combined_model_1.metabolites:
        
            if not hasattr(met,'origin'):
            
                #we found a database reaction
                setattr(met,'origin','database')
                
            #otherwise do nothing, we already have this attribute
        
        #repeat for other combined models
        for met in self.combined_model_2.metabolites:
        
            if not hasattr(met,'origin'):
            
                #we found a database reaction
                setattr(met,'origin','database')
                
            #otherwise do nothing, we already have this attribute
            
        for met in self.combined_model_3.metabolites:
        
            if not hasattr(met,'origin'):
            
                #we found a database reaction
                setattr(met,'origin','database')
        
        #In the initialization we can define xi and psi
        #xi and psi are indexed by reaction and metabolite
        for rxn in self.combined_model_1.reactions:

            #get S_ij values
            for met in self.combined_model_1.metabolites:
                
                #check if metabolite a part of the reaction
                if met in rxn.metabolites.keys():
                
                    #check the stoichiometry
                    if rxn.get_coefficient(met.id) > 0:
                        
                        #if stoichiometric coefficient is greater than zero then xi = 1, psi = 0
                        self.xis.update({str(rxn.id+met.id): 1})
                        self.psis.update({str(rxn.id+met.id): 0})
                            
                    elif rxn.get_coefficient(met.id) < 0:
                        
                        #if stoichiometric coeffici ent is greater than zero then xi = 0, psi = 1
                        self.xis.update({str(rxn.id+met.id): 0})
                        self.psis.update({str(rxn.id+met.id): 1})
                            
                    else:
                        
                        #otherwise xi and phi are zero
                        self.xis.update({str(rxn.id+met.id): 0})
                        self.psis.update({str(rxn.id+met.id): 0})
                        
                #if metabolite is not a part of the reaction, xi and phi are zero
                else:
                
                    #otherwise xi and phi are zero
                    self.xis.update({str(rxn.id+met.id): 0})
                    self.psis.update({str(rxn.id+met.id): 0})

    def fill(self,iterations=1,epsilon=0.001,tolerance=1E-8):
    
        output=open('CP_results.txt','w',buffering=1)
        output.write("Connecting Problems Solutions\n")
        output.write("Model: "+self.original_model.name+"\n")
        output.write("Database: "+self.database.name+"\n\n")
        
        #initialize the objective questions
        self.combined_model_1.objective = self.combined_model_1.problem.Objective(Zero, direction='max')
        self.combined_model_2.objective = self.combined_model_2.problem.Objective(Zero, direction='min')
        self.combined_model_3.objective = self.combined_model_3.problem.Objective(Zero, direction='max')
        
        #keeps track of the total algorithm time
        alg_start = time.time()
            
        #need to define variables and constraints related to all reaction
        for rxn in self.combined_model_1.reactions:
        
            #define variables related to all reactions
            theta_1 = self.combined_model_1.problem.Variable(name='theta_1_{}'.format(rxn.id), lb=0, ub=1, type='binary')
            lambda_var_1 = self.combined_model_1.problem.Variable(name='lambda_var_1_{}'.format(rxn.id), lb=0, ub=1, type='binary')
            
            #add variables to model
            self.combined_model_1.add_cons_vars([theta_1,lambda_var_1], sloppy=False)
            
            #have handle to reference these variables
            self.thetas_1.update({theta_1.name: theta_1})
            self.lambdas_1.update({lambda_var_1.name: lambda_var_1})
            
            #define constraints related to all reactions
            const_9_1 = self.combined_model_1.problem.Constraint(rxn.flux_expression + (rxn.upper_bound + epsilon) * theta_1, ub=rxn.upper_bound, name='const_9_1_{}'.format(rxn.id))
            const_10_1 = self.combined_model_1.problem.Constraint(rxn.flux_expression - (theta_1 * rxn.lower_bound), lb=0, name='const_10_1_{}'.format(rxn.id))
            const_11_1 = self.combined_model_1.problem.Constraint(rxn.flux_expression - (lambda_var_1 * rxn.upper_bound), ub=0, name='const_11_1_{}'.format(rxn.id))
            const_12_1 = self.combined_model_1.problem.Constraint(rxn.flux_expression + (rxn.lower_bound - epsilon) * lambda_var_1, lb=rxn.lower_bound, name='const_12_1_{}'.format(rxn.id))
            
            #add constraints to list of constraints
            self.combined_model_1.add_cons_vars([const_9_1,const_10_1,const_11_1,const_12_1],sloppy=False)
            
        #need to define variables and constraints related to all reaction
        for rxn in self.combined_model_2.reactions:
        
            #define variables related to all reactions
            theta_2 = self.combined_model_2.problem.Variable(name='theta_2_{}'.format(rxn.id), lb=0, ub=1, type='binary')
            lambda_var_2 = self.combined_model_2.problem.Variable(name='lambda_var_2_{}'.format(rxn.id), lb=0, ub=1, type='binary')
            
            #add variables to model
            self.combined_model_2.add_cons_vars([theta_2,lambda_var_2], sloppy=False)
            
            #have handle to reference these variables
            self.thetas_2.update({theta_2.name: theta_2})
            self.lambdas_2.update({lambda_var_2.name: lambda_var_2})
            
            #define constraints related to all reactions
            const_9_2 = self.combined_model_2.problem.Constraint(rxn.flux_expression + (rxn.upper_bound + epsilon) * theta_2, ub=rxn.upper_bound, name='const_9_2_{}'.format(rxn.id))
            const_10_2 = self.combined_model_2.problem.Constraint(rxn.flux_expression - (theta_2 * rxn.lower_bound), lb=0, name='const_10_2_{}'.format(rxn.id))
            const_11_2 = self.combined_model_2.problem.Constraint(rxn.flux_expression - (lambda_var_2 * rxn.upper_bound), ub=0, name='const_11_2_{}'.format(rxn.id))
            const_12_2 = self.combined_model_2.problem.Constraint(rxn.flux_expression + (rxn.lower_bound - epsilon) * lambda_var_2, lb=rxn.lower_bound, name='const_12_2_{}'.format(rxn.id))
            
            #add constraints to list of constraints
            self.combined_model_2.add_cons_vars([const_9_2,const_10_2,const_11_2,const_12_2],sloppy=False)
            
        #need to define variables and constraints related to all reaction
        for rxn in self.combined_model_3.reactions:
        
            #define variables related to all reactions
            theta_3 = self.combined_model_3.problem.Variable(name='theta_3_{}'.format(rxn.id), lb=0, ub=1, type='binary')
            lambda_var_3 = self.combined_model_3.problem.Variable(name='lambda_var_3_{}'.format(rxn.id), lb=0, ub=1, type='binary')
            
            #add variables to model
            self.combined_model_3.add_cons_vars([theta_3,lambda_var_3], sloppy=False)
            
            #have handle to reference these variables
            self.thetas_3.update({theta_3.name: theta_3})
            self.lambdas_3.update({lambda_var_3.name: lambda_var_3})
            
            #define constraints related to all reactions
            const_9_3 = self.combined_model_3.problem.Constraint(rxn.flux_expression + (rxn.upper_bound + epsilon) * theta_3, ub=rxn.upper_bound, name='const_9_3_{}'.format(rxn.id))
            const_10_3 = self.combined_model_3.problem.Constraint(rxn.flux_expression - (theta_3 * rxn.lower_bound), lb=0, name='const_10_3_{}'.format(rxn.id))
            const_11_3 = self.combined_model_3.problem.Constraint(rxn.flux_expression - (lambda_var_3 * rxn.upper_bound), ub=0, name='const_11_3_{}'.format(rxn.id))
            const_12_3 = self.combined_model_3.problem.Constraint(rxn.flux_expression + (rxn.lower_bound - epsilon) * lambda_var_3, lb=rxn.lower_bound, name='const_12_3_{}'.format(rxn.id))
            
            #add constraints to list of constraints
            self.combined_model_3.add_cons_vars([const_9_3,const_10_3,const_11_3,const_12_3],sloppy=False)
        
        #make metabolite variables for model 1
        for met in self.combined_model_1.metabolites:
        
            #define variables for all metabolites
            x_1 = self.combined_model_1.problem.Variable(name='x_1_{}'.format(met.id),lb=0, ub=1, type='binary')
            
            #add these variables to the model
            self.combined_model_1.add_cons_vars([x_1], sloppy=False)
            
            #have handle to reference these variables
            self.xs_1.update({x_1.name: x_1})
            
            #build constraint 13
            const_13_1 = self.combined_model_1.problem.Constraint(x_1,ub=0,name="const_13_1_{}".format(met.id),sloppy=False)
            self.combined_model_1.add_cons_vars([const_13_1], sloppy=False)
            self.combined_model_1.solver.update()
            
            #loop for each raection so can sum over each reaction
            for rxn in self.combined_model_1.reactions:
            
                const_13_1.set_linear_coefficients({self.lambdas_1['lambda_var_1_{}'.format(rxn.id)]: -self.xis[str(rxn.id+met.id)]})
                const_13_1.set_linear_coefficients({self.thetas_1['theta_1_{}'.format(rxn.id)]: -self.psis[str(rxn.id+met.id)]})
        
        #make metabolite variables for model 2
        for met in self.combined_model_2.metabolites:
        
            #define variables for all metabolites
            x_2 = self.combined_model_2.problem.Variable(name='x_2_{}'.format(met.id),lb=0, ub=1, type='binary')
            
            #add these variables to the model
            self.combined_model_2.add_cons_vars([x_2], sloppy=False)
            
            #have handle to reference these variables
            self.xs_2.update({x_2.name: x_2})
            
            #build constraint 13
            const_13_2 = self.combined_model_2.problem.Constraint(x_2,ub=0,name="const_13_2_{}".format(met.id),sloppy=False)
            self.combined_model_2.add_cons_vars([const_13_2], sloppy=False)
            self.combined_model_2.solver.update()
            
            #loop for each raection so can sum over each reaction
            for rxn in self.combined_model_2.reactions:
            
                const_13_2.set_linear_coefficients({self.lambdas_2['lambda_var_2_{}'.format(rxn.id)]: -self.xis[str(rxn.id+met.id)]})
                const_13_2.set_linear_coefficients({self.thetas_2['theta_2_{}'.format(rxn.id)]: -self.psis[str(rxn.id+met.id)]})
                    
        #make metabolite variables for model 3
        for met in self.combined_model_3.metabolites:
        
            #define variables for all metabolites
            x_3 = self.combined_model_3.problem.Variable(name='x_3_{}'.format(met.id),lb=0, ub=1, type='binary')
            
            x_3.met_id = met.id
            
            #add these variables to the model
            self.combined_model_3.add_cons_vars([x_3], sloppy=False)
            
            #have handle to reference these variables
            self.xs_3.update({x_3.name: x_3})
            
            #build constraint 13
            const_13_3 = self.combined_model_3.problem.Constraint(x_3,ub=0,name="const_13_3_{}".format(met.id),sloppy=False)
            self.combined_model_3.add_cons_vars([const_13_3], sloppy=False)
            self.combined_model_3.solver.update()
            
            #loop for each raection so can sum over each reaction
            for rxn in self.combined_model_3.reactions:
            
                const_13_3.set_linear_coefficients({self.lambdas_3['lambda_var_3_{}'.format(rxn.id)]: -self.xis[str(rxn.id+met.id)]})
                const_13_3.set_linear_coefficients({self.thetas_3['theta_3_{}'.format(rxn.id)]: -self.psis[str(rxn.id+met.id)]})
                
        #also need to define constraints to ensure that the number of metabolites is equal
        const_num_mets_2 = self.combined_model_2.problem.Constraint(Zero, lb=0, ub=0, name='const_num_mets_2')
        const_num_mets_3 = self.combined_model_3.problem.Constraint(Zero, lb=0, ub=0, name='const_num_mets_3')
        self.combined_model_2.add_cons_vars([const_num_mets_2], sloppy=False)
        self.combined_model_3.add_cons_vars([const_num_mets_3], sloppy=False)
        self.combined_model_2.solver.update()
        self.combined_model_3.solver.update()
        
        for met in self.combined_model_1.metabolites:
        
            #now that x is defined, can define first objective sum
            #summing over model metabolites
            if bool(re.search('^model$',met.origin)):
            
                self.combined_model_1.objective.set_linear_coefficients({self.xs_1['x_1_{}'.format(met.id)]: 1})
        
        for met in self.combined_model_2.metabolites:
        
            #summing over model metabolites
            if bool(re.search('^model$',met.origin)):
            
                const_num_mets_2.set_linear_coefficients({self.xs_2['x_2_{}'.format(met.id)]: 1})
        
        for met in self.combined_model_3.metabolites:
        
            #summing over model metabolites
            if bool(re.search('^model$',met.origin)):
            
                const_num_mets_3.set_linear_coefficients({self.xs_3['x_3_{}'.format(met.id)]: 1})
        
        #there is one constraints which sum over all database reactions
        #this constraint ensures taht at least one 
        const_1_1 = self.combined_model_1.problem.Constraint(Zero,lb=1,name="const_1_1",sloppy=False)
        const_1_2 = self.combined_model_2.problem.Constraint(Zero,lb=1,name="const_1_2",sloppy=False)
        const_1_3 = self.combined_model_3.problem.Constraint(Zero,lb=1,name="const_1_3",sloppy=False)
        self.combined_model_1.add_cons_vars([const_1_1], sloppy=False)
        self.combined_model_2.add_cons_vars([const_1_2], sloppy=False)
        self.combined_model_3.add_cons_vars([const_1_3], sloppy=False)
        self.combined_model_1.solver.update()
        self.combined_model_2.solver.update()
        self.combined_model_3.solver.update()
        
        #the results of the second objective form a constraint in the third connecting problem
        #note that the upper and lower bounds will be updated based on the solution to the second connecting problem
        const_num_rxns = self.combined_model_3.problem.Constraint(Zero, lb=0, ub=0, name='const_num_rxns')
        self.combined_model_3.add_cons_vars([const_num_rxns], sloppy=False)
        self.combined_model_3.solver.update()
        
        #need to define variables and constraints related to database reactions
        for rxn in self.combined_model_1.reactions:
        
            #check if has origin of database
            if bool(re.search('^database$',rxn.origin)):
        
                #if so then do all the defining
                #define variables which are dependent on just database reactions
                delta_1 = self.combined_model_1.problem.Variable(name='delta_1_{}'.format(rxn.id), lb=0, ub=1, type='binary')
                rho_1 = self.combined_model_1.problem.Variable(name='rho_1_{}'.format(rxn.id), lb=0, ub=1, type='binary')
                omega_1 = self.combined_model_1.problem.Variable(name='omega_1_{}'.format(rxn.id), lb=0, ub=1, type='binary')
                zeta_1 = self.combined_model_1.problem.Variable(name='zeta_1_{}'.format(rxn.id), lb=0, ub=1, type='binary')
                
                #have handle to reference these variables later
                self.deltas_1.update({delta_1.name: delta_1})
                self.rhos_1.update({rho_1.name: rho_1})
                self.omegas_1.update({omega_1.name: omega_1})
                self.zetas_1.update({zeta_1.name: zeta_1})
                
                #add these to the model
                self.combined_model_1.add_cons_vars([delta_1,rho_1,omega_1,zeta_1], sloppy=False)
                
                #make linear coefficients for constraint 1 so that can sum zeta over all database reactions
                const_1_1.set_linear_coefficients({self.zetas_1['zeta_1_{}'.format(rxn.id)]: 1})
                
                #write the constraints dependent on just database reactions
                #force database reactions to either have rho or delta = 1 if it has flu
                const_2_1 = self.combined_model_1.problem.Constraint(rxn.flux_expression - delta_1 * rxn.upper_bound, ub=0, name='const_2_1_{}'.format(rxn.id))
                const_3_1 = self.combined_model_1.problem.Constraint(rxn.flux_expression - rho_1 * rxn.lower_bound, lb=0, name='const_3_1_{}'.format(rxn.id))
                
                #Ensure that if a rection is added that it has a value of zeta = 1
                const_6_1 = self.combined_model_1.problem.Constraint(delta_1 + rho_1 - omega_1 - zeta_1, lb=0, ub=0, name='const_6_1_{}'.format(rxn.id))
                
                #ensure that omega can equal 1 iff delta and rho equal 2
                const_7_1 = self.combined_model_1.problem.Constraint(omega_1 - delta_1, ub=0, name='const_7_1_{}'.format(rxn.id))
                const_8_1 = self.combined_model_1.problem.Constraint(omega_1 - rho_1, ub=0, name='const_8_1_{}'.format(rxn.id))
                
                #ensure that omega rho and delta only occur if the database reaction has flux in the indicated direction
                #note that these limit flux to [-1,1], generally not an issue because then can run FBA later to deterime what growth rate is
                const_85_1 = self.combined_model_1.problem.Constraint(omega_1 - self.thetas_1['theta_1_{}'.format(rxn.id)] - self.lambdas_1['lambda_var_1_{}'.format(rxn.id)], ub=0, name='const_85_1_{}'.format(rxn.id))
                const_86_1 = self.combined_model_1.problem.Constraint(rho_1 - self.thetas_1['theta_1_{}'.format(rxn.id)] - self.lambdas_1['lambda_var_1_{}'.format(rxn.id)], ub=0, name='const_86_1_{}'.format(rxn.id))
                const_87_1 = self.combined_model_1.problem.Constraint(delta_1 - self.thetas_1['theta_1_{}'.format(rxn.id)] - self.lambdas_1['lambda_var_1_{}'.format(rxn.id)], ub=0, name='const_87_1_{}'.format(rxn.id))
                const_88_1 = self.combined_model_1.problem.Constraint(delta_1 - epsilon * rxn.flux_expression, lb=0, name='const_88_1_{}'.format(rxn.id))
                const_89_1 = self.combined_model_1.problem.Constraint(rho_1 + epsilon * rxn.flux_expression, lb=0, name='const_89_1_{}'.format(rxn.id))
                
                #constraint 14 has summation terms, used to determine if database reaction partcipating in the solution
                const_14_1 = self.combined_model_1.problem.Constraint(zeta_1,ub=0,name="const_14_1_{}".format(rxn.id),sloppy=False)
                
                self.combined_model_1.add_cons_vars([const_2_1,const_3_1,const_6_1,const_7_1,const_8_1,const_85_1,const_86_1,const_87_1,const_88_1,const_89_1,const_14_1], sloppy=False)
                self.combined_model_1.solver.update()
                
                for met in self.combined_model_1.metabolites:
                
                    #if xi for the current metabolite and ractions is 1, then need to include lambda in the sum
                    if self.xis[str(rxn.id+met.id)] == 1:
                        
                        const_14_1.set_linear_coefficients({self.lambdas_1['lambda_var_1_{}'.format(rxn.id)]: -1})
                        
                    elif self.psis[str(rxn.id+met.id)] == 1:
                    
                        const_14_1.set_linear_coefficients({self.thetas_1['theta_1_{}'.format(rxn.id)]: -1})
        
        #need to define variables and constraints related to database reactions
        for rxn in self.combined_model_2.reactions:
        
            #check if has origin of database
            if bool(re.search('^database$',rxn.origin)):
        
                #if so then do all the defining
                #define variables which are dependent on just database reactions
                delta_2 = self.combined_model_2.problem.Variable(name='delta_2_{}'.format(rxn.id), lb=0, ub=1, type='binary')
                rho_2 = self.combined_model_2.problem.Variable(name='rho_2_{}'.format(rxn.id), lb=0, ub=1, type='binary')
                omega_2 = self.combined_model_2.problem.Variable(name='omega_2_{}'.format(rxn.id), lb=0, ub=1, type='binary')
                zeta_2 = self.combined_model_2.problem.Variable(name='zeta_2_{}'.format(rxn.id), lb=0, ub=1, type='binary')
                
                #have handle to reference these variables later
                self.deltas_2.update({delta_2.name: delta_2})
                self.rhos_2.update({rho_2.name: rho_2})
                self.omegas_2.update({omega_2.name: omega_2})
                self.zetas_2.update({zeta_2.name: zeta_2})
                
                #add these to the model
                self.combined_model_2.add_cons_vars([delta_2,rho_2,omega_2,zeta_2], sloppy=False)
                
                #make linear coefficients for constraint 1 so that can sum zeta over all database reactions
                const_1_2.set_linear_coefficients({self.zetas_2['zeta_2_{}'.format(rxn.id)]: 1})
                
                #zeta values summed over database reactions is the second objective
                self.combined_model_2.objective.set_linear_coefficients({self.zetas_2['zeta_2_{}'.format(rxn.id)]: 1})
                
                #write the constraints dependent on just database reactions
                #force database reactions to either have rho or delta = 1 if it has flu
                const_2_2 = self.combined_model_2.problem.Constraint(rxn.flux_expression - delta_2 * rxn.upper_bound, ub=0, name='const_2_2_{}'.format(rxn.id))
                const_3_2 = self.combined_model_2.problem.Constraint(rxn.flux_expression - rho_2 * rxn.lower_bound, lb=0, name='const_3_2_{}'.format(rxn.id))
                
                #Ensure that if a rection is added that it has a value of zeta = 1
                const_6_2 = self.combined_model_2.problem.Constraint(delta_2 + rho_2 - omega_2 - zeta_2, lb=0, ub=0, name='const_6_2_{}'.format(rxn.id))
                
                #ensure that omega can equal 1 iff delta and rho equal 2
                const_7_2 = self.combined_model_2.problem.Constraint(omega_2 - delta_2, ub=0, name='const_7_2_{}'.format(rxn.id))
                const_8_2 = self.combined_model_2.problem.Constraint(omega_2 - rho_2, ub=0, name='const_8_2_{}'.format(rxn.id))
                
                #ensure that omega rho and delta only occur if the database reaction has flux in the indicated direction
                #note that these limit flux to [-1,1], generally not an issue because then can run FBA later to deterime what growth rate is
                const_85_2 = self.combined_model_2.problem.Constraint(omega_2 - self.thetas_2['theta_2_{}'.format(rxn.id)] - self.lambdas_2['lambda_var_2_{}'.format(rxn.id)], ub=0, name='const_85_2_{}'.format(rxn.id))
                const_86_2 = self.combined_model_2.problem.Constraint(rho_2 - self.thetas_2['theta_2_{}'.format(rxn.id)] - self.lambdas_2['lambda_var_2_{}'.format(rxn.id)], ub=0, name='const_86_2_{}'.format(rxn.id))
                const_87_2 = self.combined_model_2.problem.Constraint(delta_2 - self.thetas_2['theta_2_{}'.format(rxn.id)] - self.lambdas_2['lambda_var_2_{}'.format(rxn.id)], ub=0, name='const_87_2_{}'.format(rxn.id))
                const_88_2 = self.combined_model_2.problem.Constraint(delta_2 - epsilon * rxn.flux_expression, lb=0, name='const_88_2_{}'.format(rxn.id))
                const_89_2 = self.combined_model_2.problem.Constraint(rho_2 + epsilon * rxn.flux_expression, lb=0, name='const_89_2_{}'.format(rxn.id))
                
                #constraint 14 has summation terms, used to determine if database reaction partcipating in the solution
                const_14_2 = self.combined_model_2.problem.Constraint(zeta_2,ub=0,name="const_14_2_{}".format(rxn.id),sloppy=False)
                
                self.combined_model_2.add_cons_vars([const_2_2,const_3_2,const_6_2,const_7_2,const_8_2,const_85_2,const_86_2,const_87_2,const_88_2,const_89_2,const_14_2], sloppy=False)
                self.combined_model_2.solver.update()
                
                for met in self.combined_model_2.metabolites:
                
                    #if xi for the current metabolite and ractions is 1, then need to include lambda in the sum
                    if self.xis[str(rxn.id+met.id)] == 1:
                        
                        const_14_2.set_linear_coefficients({self.lambdas_2['lambda_var_2_{}'.format(rxn.id)]: -1})
                        
                    elif self.psis[str(rxn.id+met.id)] == 1:
                    
                        const_14_2.set_linear_coefficients({self.thetas_2['theta_2_{}'.format(rxn.id)]: -1})
                        
        #need to define variables and constraints related to database reactions
        for rxn in self.combined_model_3.reactions:
        
            #check if has origin of database
            if bool(re.search('^database$',rxn.origin)):
        
                #if so then do all the defining
                #define variables which are dependent on just database reactions
                delta_3 = self.combined_model_3.problem.Variable(name='delta_3_{}'.format(rxn.id), lb=0, ub=1, type='binary')
                rho_3 = self.combined_model_3.problem.Variable(name='rho_3_{}'.format(rxn.id), lb=0, ub=1, type='binary')
                omega_3 = self.combined_model_3.problem.Variable(name='omega_3_{}'.format(rxn.id), lb=0, ub=1, type='binary')
                zeta_3 = self.combined_model_3.problem.Variable(name='zeta_3_{}'.format(rxn.id), lb=0, ub=1, type='binary')
                
                #have handle to reference these variables later
                self.deltas_3.update({delta_3.name: delta_3})
                self.rhos_3.update({rho_3.name: rho_3})
                self.omegas_3.update({omega_3.name: omega_3})
                self.zetas_3.update({zeta_3.name: zeta_3})
                
                #add these to the model
                self.combined_model_3.add_cons_vars([delta_3,rho_3,omega_3,zeta_3], sloppy=False)
                
                #make linear coefficients for constraint 1 so that can sum zeta over all database reactions
                const_1_3.set_linear_coefficients({self.zetas_3['zeta_3_{}'.format(rxn.id)]: 1})
                
                #enforce the optimal value of CP2 in CP3
                const_num_rxns.set_linear_coefficients({self.zetas_3['zeta_3_{}'.format(rxn.id)]: 1})
                
                #omega values summed over databse reactions is the third objective 
                self.combined_model_3.objective.set_linear_coefficients({self.omegas_3['omega_3_{}'.format(rxn.id)]: 1})
                
                #write the constraints dependent on just database reactions
                #force database reactions to either have rho or delta = 1 if it has flu
                const_2_3 = self.combined_model_3.problem.Constraint(rxn.flux_expression - delta_3 * rxn.upper_bound, ub=0, name='const_2_3_{}'.format(rxn.id))
                const_3_3 = self.combined_model_3.problem.Constraint(rxn.flux_expression - rho_3 * rxn.lower_bound, lb=0, name='const_3_3_{}'.format(rxn.id))
                
                #Ensure that if a rection is added that it has a value of zeta = 1
                const_6_3 = self.combined_model_3.problem.Constraint(delta_3 + rho_3 - omega_3 - zeta_3, lb=0, ub=0, name='const_6_3_{}'.format(rxn.id))
                
                #ensure that omega can equal 1 iff delta and rho equal 2
                const_7_3 = self.combined_model_3.problem.Constraint(omega_3 - delta_3, ub=0, name='const_7_3_{}'.format(rxn.id))
                const_8_3 = self.combined_model_3.problem.Constraint(omega_3 - rho_3, ub=0, name='const_8_3_{}'.format(rxn.id))
                
                #ensure that omega rho and delta only occur if the database reaction has flux in the indicated direction
                #note that these limit flux to [-1,1], generally not an issue because then can run FBA later to deterime what growth rate is
                const_85_3 = self.combined_model_3.problem.Constraint(omega_3 - self.thetas_3['theta_3_{}'.format(rxn.id)] - self.lambdas_3['lambda_var_3_{}'.format(rxn.id)], ub=0, name='const_85_3_{}'.format(rxn.id))
                const_86_3 = self.combined_model_3.problem.Constraint(rho_3 - self.thetas_3['theta_3_{}'.format(rxn.id)] - self.lambdas_3['lambda_var_3_{}'.format(rxn.id)], ub=0, name='const_86_3_{}'.format(rxn.id))
                const_87_3 = self.combined_model_3.problem.Constraint(delta_3 - self.thetas_3['theta_3_{}'.format(rxn.id)] - self.lambdas_3['lambda_var_3_{}'.format(rxn.id)], ub=0, name='const_87_3_{}'.format(rxn.id))
                const_88_3 = self.combined_model_3.problem.Constraint(delta_3 - epsilon * rxn.flux_expression, lb=0, name='const_88_3_{}'.format(rxn.id))
                const_89_3 = self.combined_model_3.problem.Constraint(rho_3 + epsilon * rxn.flux_expression, lb=0, name='const_89_3_{}'.format(rxn.id))
                
                #constraint 14 has summation terms, used to determine if database reaction partcipating in the solution
                const_14_3 = self.combined_model_3.problem.Constraint(zeta_3,ub=0,name="const_14_3_{}".format(rxn.id),sloppy=False)
                
                self.combined_model_3.add_cons_vars([const_2_3,const_3_3,const_6_3,const_7_3,const_8_3,const_85_3,const_86_3,const_87_3,const_88_3,const_89_3,const_14_3], sloppy=False)
                self.combined_model_3.solver.update()
                
                for met in self.combined_model_3.metabolites:
                
                    #if xi for the current metabolite and ractions is 1, then need to include lambda in the sum
                    if self.xis[str(rxn.id+met.id)] == 1:
                        
                        const_14_3.set_linear_coefficients({self.lambdas_3['lambda_var_3_{}'.format(rxn.id)]: -1})
                        
                    elif self.psis[str(rxn.id+met.id)] == 1:
                    
                        const_14_3.set_linear_coefficients({self.thetas_3['theta_3_{}'.format(rxn.id)]: -1})
    
        #ensure that biomass or other target metabolite is produced
        const_4_1 = self.combined_model_1.problem.Constraint(Zero,lb=1,name="const_4_1",sloppy=False)
        const_4_2 = self.combined_model_2.problem.Constraint(Zero,lb=1,name="const_4_2",sloppy=False)
        const_4_3 = self.combined_model_3.problem.Constraint(Zero,lb=1,name="const_4_3",sloppy=False)
        self.combined_model_1.add_cons_vars([const_4_1], sloppy=False)
        self.combined_model_2.add_cons_vars([const_4_2], sloppy=False)
        self.combined_model_3.add_cons_vars([const_4_3], sloppy=False)
        self.combined_model_1.solver.update()
        self.combined_model_2.solver.update()
        self.combined_model_3.solver.update()
        
        for met in self.must_prod:
        
            const_4_1.set_linear_coefficients({self.xs_1["x_1_{}".format(met)]: 1})
            const_4_2.set_linear_coefficients({self.xs_2["x_2_{}".format(met)]: 1})
            const_4_3.set_linear_coefficients({self.xs_3["x_3_{}".format(met)]: 1})
        
        #note that in GAMS constraint 5 is mass balance, which is automatically defined in COBRApy
        
        #finally, need to add integer cuts based on the TIC-Finding Problem        
        #need to build in a loop for each TFp solution
        for solution in self.TFPsol_list:
        
            #initialize the integer cut variable
            tau_1 = self.combined_model_1.problem.Variable(name='tau_1_{}'.format(solution), lb=0, ub=1, type='binary')
            tau_2 = self.combined_model_2.problem.Variable(name='tau_2_{}'.format(solution), lb=0, ub=1, type='binary')
            tau_3 = self.combined_model_3.problem.Variable(name='tau_3_{}'.format(solution), lb=0, ub=1, type='binary')
            
            #first need to get the sum of alphas and betas over all database reactions for the given TFP solution
            sum_alphas = 0
            sum_betas = 0
            
            for rxn in self.database.reactions:
            
                sum_alphas = sum_alphas + self.TFP_alphas["{}_{}".format(solution,rxn.id)]
                sum_betas = sum_betas + self.TFP_betas["{}_{}".format(solution,rxn.id)]
        
            #initialize the integer cuts for the TFP solutions
            #integer cut 4 ensures that that forward TIC elements not added in entirety
            int_cut_4_1 = self.combined_model_1.problem.Constraint(tau_1,ub=sum_alphas,name="int_cut_4_1_{}".format(solution),sloppy=False)
            int_cut_4_2 = self.combined_model_2.problem.Constraint(tau_2,ub=sum_alphas,name="int_cut_4_2_{}".format(solution),sloppy=False)
            int_cut_4_3 = self.combined_model_3.problem.Constraint(tau_3,ub=sum_alphas,name="int_cut_4_3_{}".format(solution),sloppy=False)
            self.combined_model_1.add_cons_vars([int_cut_4_1], sloppy=False)
            self.combined_model_2.add_cons_vars([int_cut_4_2], sloppy=False)
            self.combined_model_3.add_cons_vars([int_cut_4_3], sloppy=False)
            self.combined_model_1.repair()
            self.combined_model_2.repair()
            self.combined_model_3.repair()
            self.combined_model_1.solver.update()
            self.combined_model_2.solver.update()
            self.combined_model_3.solver.update()
            
            #integer cut 5 ensures that that backward TIC elements not added in entirety
            int_cut_5_1 = self.combined_model_1.problem.Constraint(-tau_1,ub=(sum_betas-1),name="int_cut_5_1_{}".format(solution),sloppy=False)
            int_cut_5_2 = self.combined_model_2.problem.Constraint(-tau_2,ub=(sum_betas-1),name="int_cut_5_2_{}".format(solution),sloppy=False)
            int_cut_5_3 = self.combined_model_3.problem.Constraint(-tau_3,ub=(sum_betas-1),name="int_cut_5_3_{}".format(solution),sloppy=False)
            self.combined_model_1.add_cons_vars([int_cut_5_1], sloppy=False)
            self.combined_model_2.add_cons_vars([int_cut_5_2], sloppy=False)
            self.combined_model_3.add_cons_vars([int_cut_5_3], sloppy=False)
            self.combined_model_1.repair()
            self.combined_model_2.repair()
            self.combined_model_3.repair()
            self.combined_model_1.solver.update()
            self.combined_model_2.solver.update()
            self.combined_model_3.solver.update()
            
            #need to add coefficients for rho and delta to the integer cuts
            #just for database reactions 
            for rxn in self.combined_model_1.reactions:
            
                #check if has origin of database
                if bool(re.search('^database$',rxn.origin)):
            
                    int_cut_4_1.set_linear_coefficients({self.deltas_1['delta_1_{}'.format(rxn.id)]: self.TFP_alphas["{}_{}".format(solution,rxn.id)]})
                    int_cut_5_1.set_linear_coefficients({self.rhos_1['rho_1_{}'.format(rxn.id)]: self.TFP_betas["{}_{}".format(solution,rxn.id)]})
        
            #need to add coefficients for rho and delta to the integer cuts
            for rxn in self.combined_model_2.reactions:
            
                #check if has origin of database
                if bool(re.search('^database$',rxn.origin)):
            
                    int_cut_4_2.set_linear_coefficients({self.deltas_2['delta_2_{}'.format(rxn.id)]: self.TFP_alphas["{}_{}".format(solution,rxn.id)]})
                    int_cut_5_2.set_linear_coefficients({self.rhos_2['rho_2_{}'.format(rxn.id)]: self.TFP_betas["{}_{}".format(solution,rxn.id)]})
                
            #need to add coefficients for rho and delta to the integer cuts
            for rxn in self.combined_model_3.reactions:
            
                #check if has origin of database
                if bool(re.search('^database$',rxn.origin)):
            
                    int_cut_4_3.set_linear_coefficients({self.deltas_3['delta_3_{}'.format(rxn.id)]: self.TFP_alphas["{}_{}".format(solution,rxn.id)]})
                    int_cut_5_3.set_linear_coefficients({self.rhos_3['rho_3_{}'.format(rxn.id)]: self.TFP_betas["{}_{}".format(solution,rxn.id)]})
        
        #adjust the tolerance on the solvers
        self.combined_model_1.tolerance = tolerance
        self.combined_model_2.tolerance = tolerance
        self.combined_model_3.tolerance = 5E-8
        
        #make sure solvers are up-to-date
        self.combined_model_1.repair()
        self.combined_model_2.repair()
        self.combined_model_3.repair()
        self.combined_model_1.solver.update() 
        self.combined_model_2.solver.update()
        self.combined_model_3.solver.update()
        
        #iterate to get multiple solutions
        for iteration in range(iterations):
        
            #solve the first connecting problem
            #try to solve 
            try:
            
                #catch warnings related to infeasible
                with warnings.catch_warnings(): 
                
                    warnings.simplefilter('ignore',category=UserWarning)
                    
                    #attempt to solve
                    self.combined_model_1.repair()
                    self.combined_model_1.solver.update()
                    
                    print("attempting to solve CP1...")
                    
                    #time how long things take
                    start_time = time.time()
                    
                    CP1_solution = self.combined_model_1.optimize()
                    
                    end_time = time.time()
                    
                    CP1_time = end_time - start_time
                    
                    print("complete, objective value: "+str(CP1_solution.objective_value))
            
            except:
            
                end_time = time.time()
                
                CP1_time = end_time - start_time
                
                print("error occurred: "+str(sys.exec_info()[0]))
                print("no more solutions possible")
                
                output.write("error occurred: "+str(sys.exec_info()[0])+"\n")
                output.write("no more solutions possible\n")
                output.write("solution time: "+str(CP1_time)+"\n")
                
                #define an empty solution so errors don't get thrown later
                CP1_solution = cobra.core.Solution(objective_value=0,status=INFEASIBLE,fluxes=None)
            
            if CP1_solution.status == INFEASIBLE:
            
                #no solution found so all solutions must have been found by now
                print("no more solutions possible\n")
                output.write("no more solutions possible\n")
                output.write("solution time: "+str(CP1_time))
                
                #break the iteration loops
                break
            
            else:
            
                #has found a solution, now solve CP2 and CP3
                #if there is a solution then there likely won't be errors when solving the second and third problems
                
                #enforce CP1 solution level
                #errors will get thrown if ub < lb and visa versa, so will quickly reset the bounds to make sure
                #this is not the case
                const_num_mets_2.lb = 0
                const_num_mets_2.ub = 0
                const_num_mets_3.lb = 0
                const_num_mets_3.ub = 0
                
                const_num_mets_2.ub = CP1_solution.objective_value
                const_num_mets_2.lb = CP1_solution.objective_value
                const_num_mets_3.ub = CP1_solution.objective_value
                const_num_mets_3.lb = CP1_solution.objective_value
                
                output.write("\n\nCONNECTING PROBLEMS SOLUTION #"+str(iteration)+"\n")
                output.write("CP1 Objective Value: "+str(CP1_solution.objective_value)+" (solve time: "+str(CP1_time)+")\n")
                
                print("attempting to solve CP2...")
                
                #time how long things take
                start_time = time.time()
                
                self.combined_model_2.repair()
                self.combined_model_2.solver.update()
                
                CP2_solution = self.combined_model_2.optimize()
                
                end_time = time.time()
                
                CP2_time = end_time - start_time
                
                print("complete, objective value: "+str(CP2_solution.objective_value))
                
                #enforce CP2 ojbective value
                #reset to small numbers first so that we don't get a 
                const_num_rxns.lb = 0
                const_num_rxns.ub = 0
                const_num_rxns.ub = CP2_solution.objective_value
                const_num_rxns.lb = CP2_solution.objective_value
                
                output.write("CP2 Objective Value: "+str(CP2_solution.objective_value)+" (solve time: "+str(CP2_time)+")\n")
                
                print("attempting to solve CP3...")
                
                #time how long things take
                start_time = time.time()
                
                self.combined_model_3.repair()
                self.combined_model_3.solver.update()
                
                CP3_solution = self.combined_model_3.optimize()
                
                end_time = time.time()
                
                CP3_time = end_time - start_time
                
                print("complete, objective value: "+str(CP3_solution.objective_value))
                print("successfully solved all three connecting problems")
                print("solution #"+str(iteration)+" identifed")
                print("reporting results...")
                
                output.write("CP3 Objective Value: "+str(CP3_solution.objective_value)+" (solve time: "+str(CP3_time)+")\n\n")
                output.write("DB Rxn to add\t\tDirection\t\tDelta\t\tRho\t\tOmega\t\tZeta\t\ttheta\t\tlambda\t\tFlux\n")
                output.write("------------------------------------------------------------------------------------------------------------\n")
                
                #now write table for results
                #let's have it a bit abbreviated, just show those reactions where zeta = 1
                for rxn in self.database.reactions:
                
                    #if the database reaction was used in the solution report on it
                    #otherwise do nothing
                    
                    #get a string version of direction for human-readability
                    direction = "XXX"
                        
                    if (self.combined_model_3.solver.primal_values.get('delta_3_{}'.format(rxn.id)) == 1.0) and (self.combined_model_3.solver.primal_values.get('rho_3_{}'.format(rxn.id)) == 1.0): 
                        
                        #reaction added reversibly
                        direction = "<->"
                            
                    elif (self.combined_model_3.solver.primal_values.get('delta_3_{}'.format(rxn.id)) == 1.0) and (self.combined_model_3.solver.primal_values.get('rho_3_{}'.format(rxn.id)) == 0.0):
                        
                        #reaction added irreversibly forwards
                        direction = " ->"
                        
                    elif (self.combined_model_3.solver.primal_values.get('rho_3_{}'.format(rxn.id)) == 1.0) and (self.combined_model_3.solver.primal_values.get('delta_3_{}'.format(rxn.id)) == 0.0):
                        
                        #reaction added irreversibly forwards
                        direction = "<- "
                    
                    #get the net reaction flux
                    net_flux = self.combined_model_3.solver.primal_values.get('{}'.format(rxn.id)) - self.combined_model_3.solver.primal_values.get('{}'.format(rxn.reverse_id))
                    
                    #write the information about the 
                    output.write(str(rxn.id)+"\t\t\t"+str(direction)+"\t\t\t\t"+str(self.combined_model_3.solver.primal_values.get('delta_3_{}'.format(rxn.id)))+"\t\t\t"+str(self.combined_model_3.solver.primal_values.get('rho_3_{}'.format(rxn.id)))+"\t\t"+str(self.combined_model_3.solver.primal_values.get('omega_3_{}'.format(rxn.id)))+"\t\t\t"+str(self.combined_model_3.solver.primal_values.get('zeta_3_{}'.format(rxn.id)))+"\t\t\t"+str(self.combined_model_3.solver.primal_values.get('theta_3_{}'.format(rxn.id)))+"\t\t\t"+str(self.combined_model_3.solver.primal_values.get('lambda_var_3_{}'.format(rxn.id)))+"\t\t\t"+str(net_flux)+"\n")
                
                #write the reactions to add
                output.write("\nReactions to add for OptFill Solution\n---------------------------------------\n")
                
                for rxn in self.database.reactions:
                
                    if self.combined_model_3.solver.primal_values.get('omega_3_{}'.format(rxn.id)) == 1.0:
                    
                        output.write(str(rxn)+" reversible\n")
                        
                    elif self.combined_model_3.solver.primal_values.get('delta_3_{}'.format(rxn.id)) == 1.0:
                    
                        output.write(str(rxn)+" forward\n")
                    
                    elif self.combined_model_3.solver.primal_values.get('rho_3_{}'.format(rxn.id)) == 1.0:
                    
                        output.write(str(rxn)+" backward\n")
                
                #check status of model reactions
                output.write("\nMODEL REACTIONS\t\tTHETA\t\tLAMBDA\t\tFLUX\n")
                output.write("--------------------------------------------------------\n")                
                
                for rxn in self.original_model.reactions:
                
                    net_flux = self.combined_model_3.solver.primal_values.get('{}'.format(rxn.id)) - self.combined_model_3.solver.primal_values.get('{}'.format(rxn.reverse_id))
                    output.write(str(rxn.id)+"\t\t\t"+str(self.combined_model_3.solver.primal_values.get('theta_3_{}'.format(rxn.id)))+"\t\t\t"+str(self.combined_model_3.solver.primal_values.get('lambda_var_3_{}'.format(rxn.id)))+"\t\t"+str(net_flux)+"\n")
                
                output.write("\n\nMET\t\t\t\t\tX\n-----------------------\n");
                
                #now we have output all reactions which the connecting problems add, let us now report on the metabolites that can be produced
                for met in self.combined_model_3.metabolites:
                
                    output.write(str(met.id)+"\t\t\t"+str(self.combined_model_3.solver.primal_values.get('x_3_{}'.format(met.id)))+"\n")
                
                #now add some newlines to ensure adequate spacing in results
                output.write("\n\n")
                
                #next, build integer cuts for the three integer cuts for the CPs solutions
                #should create a cut for each combined model
                
                #for the integer cuts, need sum of: omegas, rhos, and deltas_1
                sum_omega = 0
                sum_delta = 0
                sum_rho = 0
                
                for rxn in self.database.reactions:
                
                    #rounding to ignore slight error/wiggle in binary values
                    sum_omega = sum_omega + (self.combined_model_3.solver.primal_values.get('omega_3_{}'.format(rxn.id)))
                    sum_rho = sum_rho + (self.combined_model_3.solver.primal_values.get('rho_3_{}'.format(rxn.id)))
                    sum_delta = sum_delta + (self.combined_model_3.solver.primal_values.get('delta_3_{}'.format(rxn.id)))
                    
                #define binary variable for integer cuts
                sigma_1 = self.combined_model_1.problem.Variable(name='sigma_1_{}'.format(iteration), lb=0, ub=1, type='binary')
                sigma_2 = self.combined_model_2.problem.Variable(name='sigma_2_{}'.format(iteration), lb=0, ub=1, type='binary')
                sigma_3 = self.combined_model_3.problem.Variable(name='sigma_3_{}'.format(iteration), lb=0, ub=1, type='binary')
                
                self.combined_model_1.add_cons_vars([sigma_1], sloppy=False)
                self.combined_model_2.add_cons_vars([sigma_2], sloppy=False)
                self.combined_model_3.add_cons_vars([sigma_3], sloppy=False)
                
                #now define integer cuts
                #first integer cut: make sure that deltas are different
                int_cut_1_1 = self.combined_model_1.problem.Constraint(Zero,ub=sum_delta,name="int_cut_1_1_{}".format(iteration),sloppy=False)
                int_cut_1_2 = self.combined_model_2.problem.Constraint(Zero,ub=sum_delta,name="int_cut_1_2_{}".format(iteration),sloppy=False)
                int_cut_1_3 = self.combined_model_3.problem.Constraint(Zero,ub=sum_delta,name="int_cut_1_3_{}".format(iteration),sloppy=False)
                self.combined_model_1.add_cons_vars([int_cut_1_1], sloppy=False)
                self.combined_model_2.add_cons_vars([int_cut_1_2], sloppy=False)
                self.combined_model_3.add_cons_vars([int_cut_1_3], sloppy=False)
                self.combined_model_1.repair()
                self.combined_model_2.repair()
                self.combined_model_3.repair()
                self.combined_model_1.solver.update()
                self.combined_model_2.solver.update()
                self.combined_model_3.solver.update()
                int_cut_1_1.set_linear_coefficients({sigma_1: 1})
                int_cut_1_2.set_linear_coefficients({sigma_2: 1})
                int_cut_1_3.set_linear_coefficients({sigma_3: 1})
                
                #second integer cut: make sure that rhos are different
                int_cut_2_1 = self.combined_model_1.problem.Constraint(Zero,ub=(sum_rho-1),name="int_cut_2_1_{}".format(iteration),sloppy=False)
                int_cut_2_2 = self.combined_model_2.problem.Constraint(Zero,ub=(sum_rho-1),name="int_cut_2_2_{}".format(iteration),sloppy=False)
                int_cut_2_3 = self.combined_model_3.problem.Constraint(Zero,ub=(sum_rho-1),name="int_cut_2_3_{}".format(iteration),sloppy=False)
                self.combined_model_1.add_cons_vars([int_cut_2_1], sloppy=False)
                self.combined_model_2.add_cons_vars([int_cut_2_2], sloppy=False)
                self.combined_model_3.add_cons_vars([int_cut_2_3], sloppy=False)
                self.combined_model_1.repair()
                self.combined_model_2.repair()
                self.combined_model_3.repair()
                self.combined_model_1.solver.update()
                self.combined_model_2.solver.update()
                self.combined_model_3.solver.update()
                int_cut_2_1.set_linear_coefficients({sigma_1: -1})
                int_cut_2_2.set_linear_coefficients({sigma_2: -1})
                int_cut_2_3.set_linear_coefficients({sigma_3: -1})
                
                #third integer cut: make sure that omegas are different
                int_cut_3_1 = self.combined_model_1.problem.Constraint(Zero,lb=(sum_omega-sum_delta-sum_rho+1),name="int_cut_3_1_{}".format(iteration),sloppy=False)
                int_cut_3_2 = self.combined_model_2.problem.Constraint(Zero,lb=(sum_omega-sum_delta-sum_rho+1),name="int_cut_3_2_{}".format(iteration),sloppy=False)
                int_cut_3_3 = self.combined_model_3.problem.Constraint(Zero,lb=(sum_omega-sum_delta-sum_rho+1),name="int_cut_3_3_{}".format(iteration),sloppy=False)
                self.combined_model_1.add_cons_vars([int_cut_3_1], sloppy=False)
                self.combined_model_2.add_cons_vars([int_cut_3_2], sloppy=False)
                self.combined_model_3.add_cons_vars([int_cut_3_3], sloppy=False)
                self.combined_model_1.repair()
                self.combined_model_2.repair()
                self.combined_model_3.repair()
                self.combined_model_1.solver.update()
                self.combined_model_2.solver.update()
                self.combined_model_3.solver.update()
                
                #now add rho and delta variables to the integer cuts
                for rxn in self.database.reactions:
                
                    #build backwards flux integer cuts
                    if self.combined_model_3.solver.primal_values.get('rho_3_{}'.format(rxn.id)) == 1:
                    
                        int_cut_2_1.set_linear_coefficients({self.rhos_1['rho_1_{}'.format(rxn.id)]: 1})
                        int_cut_2_2.set_linear_coefficients({self.rhos_2['rho_2_{}'.format(rxn.id)]: 1})
                        int_cut_2_3.set_linear_coefficients({self.rhos_3['rho_3_{}'.format(rxn.id)]: 1})
                        int_cut_3_1.set_linear_coefficients({self.rhos_1['rho_1_{}'.format(rxn.id)]: -1})
                        int_cut_3_2.set_linear_coefficients({self.rhos_2['rho_2_{}'.format(rxn.id)]: -1})
                        int_cut_3_3.set_linear_coefficients({self.rhos_3['rho_3_{}'.format(rxn.id)]: -1})
                    
                    #build forward flux integer cuts
                    if self.combined_model_3.solver.primal_values.get('delta_3_{}'.format(rxn.id)) == 1:
                    
                        int_cut_1_1.set_linear_coefficients({self.deltas_1['delta_1_{}'.format(rxn.id)]: 1})
                        int_cut_1_2.set_linear_coefficients({self.deltas_2['delta_2_{}'.format(rxn.id)]: 1})
                        int_cut_1_3.set_linear_coefficients({self.deltas_3['delta_3_{}'.format(rxn.id)]: 1})
                        int_cut_3_1.set_linear_coefficients({self.deltas_1['delta_1_{}'.format(rxn.id)]: -1})
                        int_cut_3_2.set_linear_coefficients({self.deltas_2['delta_2_{}'.format(rxn.id)]: -1})
                        int_cut_3_3.set_linear_coefficients({self.deltas_3['delta_3_{}'.format(rxn.id)]: -1})
                        
                #by this point should have the integer cuts done
                output.write("Integer cut 1: "+str(int_cut_1_1)+"\n")
                output.write("Integer cut 2: "+str(int_cut_2_1)+"\n")
                output.write("Integer cut 3: "+str(int_cut_3_1)+"\n")
                
            #then go to the next iteration
            #make sure solvers are up-to-date
            #update the model solvers
            self.combined_model_1.repair()
            self.combined_model_2.repair()
            self.combined_model_3.repair()
            self.combined_model_1.solver.update()
            self.combined_model_2.solver.update()
            self.combined_model_3.solver.update()
            
        #report on total algorithm time
        alg_end = time.time()
        total_alg_time = alg_end - alg_start
        output.write("\n\ntotal CPs runtime: "+str(total_alg_time)+" s")
       
        #return just the objective value of the 
        return self