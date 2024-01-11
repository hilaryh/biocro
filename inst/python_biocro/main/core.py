import pickle as pkl
import cobra
from cobra import Reaction, Metabolite, Model,DictList
import pandas as pd
import math
import signal
class WholePlantModel:
    ## Create whole plant model with leaves, stems, fruits, and root
    ## Based on Sanu Shameer's original BioCro-FBA code, based on Tarit Konuntakiet's code
    def __init__(self,id='soybeanPlant'):
        self.id                 = id
        self.growthCycle        = 1
        self.leaves             = DictList()
        self.stems              = DictList()
        self.fruits             = DictList()
        self.roots              = DictList()
        self.seeds              = DictList()
        self.shells              = DictList()
        self.symbionts          = DictList()
        self.parameters         = {}
        self.age                = 0
        self.mode               = 'daily'
        self.MetaboliteStores   = {}
        self.UnmetNGAM          = {}
        self.maxFlux            = 1000000
        self.maxDebt            = 10

    def addRoot(self,newOrgan):
        self.roots +=newOrgan
    def addLeaf(self, newOrgan):
        self.leaves += newOrgan
    # def addShell(self, newOrgan):
    #     self.fruits += newOrgan
    def addStem(self, newOrgan):
        self.stems += newOrgan
    def addSeed(self, newOrgan):
        self.seeds += newOrgan
    def addShell(self, newOrgan):
        self.shells += newOrgan
    def addSymbiont(self, newOrgan):
        self.symbionts += newOrgan
        
    def convertUnits(self, value, orig="umol/m2/s",final="mmol/gDW/day",organAreatoMass=0.036919917):
        if orig=="umol/m2/s" and final=="mmol/gDW/day":
            value2 = (value*60*60*self.parameters["photoperiod"]*organAreatoMass)/1000
        elif orig=="umol/m2/s" and final=="mmol/gDW/hr":
            value2 = (value*60*60*organAreatoMass)/(1000)
        elif orig=="mmol/gDW/hr" and final=="umol/m2/s":
            value2 = (value*1000)/(60*60*organAreatoMass)
        elif orig=="Mg/Ha/hr" and final=="mmol/plant/hr":
            value2 = value*1e6*1e-4*(0.38/20)
        return value2

    def save(self,filename):
        with open(filename,'wb') as file:
            pkl.dump(self,file)

    def load(self,filename):
        with open(filename,'rb') as file:
            temp_self = pkl.load(file)
        return temp_self
    
    def generateStoichiometricModels(self, sbmlFile = "PlantCoreMetabolismAA.xml",biomassEquations = "./Parameters/soyBiomassEquations.csv", biomassEquationsFruit = "./Parameters/shellBiomassComposition.csv",sink_tissues = ["R", "ST","SE","SH"]):
        model = cobra.io.read_sbml_model(sbmlFile)
        model.solver = 'glpk'
        for rxn in model.reactions:
            if rxn.objective_coefficient !=0:
                rxn.objective_coefficient = 0

        # Prevent exchange reactions
        exch_rxns = ["GLC_tx", "Sucrose_tx", "NH4_tx", "NADPH_Dehydrogenase_p", "Plastoquinol_Oxidase_p", "ATP_pc","Photon_tx", "Pi_tx", "SO4_tx", "Nitrate_tx", "Ca_tx", "Mg_tx", "K_tx"]
        for rxn in exch_rxns:
            model.reactions.get_by_id(rxn).lower_bound=0
            model.reactions.get_by_id(rxn).upper_bound=0
        unnec_rxns = [rxn.id for rxn in model.reactions if any([y in rxn.id for y in ['Biomass','unlProtHYPO']])]
        model.remove_reactions(unnec_rxns)
        
        # Increase flux bounds for 'unbound' reactions
        for rxn in model.reactions:
            if rxn.lower_bound == -1000:
                rxn.lower_bound = -self.maxFlux
            if rxn.upper_bound == 1000:
                rxn.upper_bound = self.maxFlux

        # Maintenance ratios
        ATP_NADPH_f = cobra.Metabolite('ATP_NADPH_f', name='ATP NADPH pseudo', compartment='f')
        model.reactions.ATPase_tx.add_metabolites({ATP_NADPH_f: 1.0})
        model.reactions.NADPHoxc_tx.add_metabolites({ATP_NADPH_f: -3.0})
        model.reactions.NADPHoxm_tx.add_metabolites({ATP_NADPH_f: -3.0})
        model.reactions.NADPHoxp_tx.add_metabolites({ATP_NADPH_f: -3.0})

        # # Total CO2 exchange
        # total_co2 = model.problem.Constraint(0,lb=-self.maxFlux,ub=0,name='Total_CO2')
        # model.add_cons_vars(total_co2)
        # cons_coeffs={}
        # for rxn in [x for x in model.reactions if 'CO2_tx' in x.id]:
        #     cons_coeffs[rxn.forward_variable]=1
        #     cons_coeffs[rxn.reverse_variable]=-1
        # model.solver.update()
        # total_co2.set_linear_coefficients(coefficients=cons_coeffs)

        # NOTE Chance this causes problems after model duplication...
        NAD_balcp_constr = model.problem.Constraint(model.reactions.NADPHoxc_tx.forward_variable-model.reactions.NADPHoxm_tx.forward_variable,lb=0,ub=0,name='NAD_balcp')
        NAD_balcm_constr = model.problem.Constraint(model.reactions.NADPHoxc_tx.forward_variable-model.reactions.NADPHoxp_tx.forward_variable,lb=0,ub=0,name='NAD_balcm')
        model.add_cons_vars([NAD_balcp_constr,NAD_balcm_constr])
        model.solver.update()

        #Adding biomass reactions for new biomass components 
        ## Double check this
        for met in ["FRU", "FRUCTOSE_6P", "GLC", "MALTOSE", "Pi", "PYRUVATE"]:
            met_b = model.metabolites.get_by_id(met + "_c").copy()
            met_b.id = "s" + met + "_b"
            met_b.compartment = "b"
            met_b.name = "s" + met + "_b"
            rxn = Reaction("s" + met + "_biomass", upper_bound = 0, lower_bound = -self.maxFlux, name = met + "_biomass")
            model.add_reactions([rxn])
            model.reactions.get_by_id("s" + met + "_biomass").add_metabolites({met_b: -1, model.metabolites.get_by_id(met + "_c"): 1})

        # Phloem
        ## Adjust phloem to include cytosolic sucrose
        sucrose_stoic = abs(model.reactions.Phloem_output_tx.metabolites[model.metabolites.sSUCROSE_b])
        model.reactions.Phloem_output_tx.add_metabolites({model.metabolites.sSUCROSE_b: sucrose_stoic,model.metabolites.SUCROSE_c: -sucrose_stoic})

        ## Removing amino acids without degradation pathways from phloem composition
        for met in ["MET", "HIS", "CYS", "PHE", "TRP", "TYR"]:
            metStoic = abs(model.reactions.Phloem_output_tx.metabolites[model.metabolites.get_by_id(met+"_c")])
            model.reactions.Phloem_output_tx.add_metabolites({model.metabolites.get_by_id(met+"_c"): metStoic, model.metabolites.PROTON_c: -metStoic, model.metabolites.PROTON_e: metStoic})
                        
        ## Make separate reaction for each phloem component and disable defined phloem reaction or
        if self.phloemComposition == "free":
            for met in model.reactions.Phloem_output_tx.metabolites:
                if not((met.id.startswith("PROTON"))):
                    met_name = met.id[:-2]
                    rxn = Reaction(met_name + "_ph", name = met_name + "_ph", subsystem = "Phloem", lower_bound = 0, upper_bound = self.maxFlux)
                    new_met = model.metabolites.get_by_id(met.id).copy()
                    new_met.name = met_name + "_ph"
                    new_met.id = met_name + "_ph"
                    new_met.compartment = "ph"
                    model.add_reactions([rxn])
                    model.reactions.get_by_id(met_name + "_ph").add_metabolites({model.metabolites.PROTON_c: 1,
                                                                                 model.metabolites.PROTON_e: -1,
                                                                                 met: -1,
                                                                                 new_met: 1})
            ###Remove phloem reaction
            model.reactions.Phloem_output_tx.remove_from_model()

        ## Establish fixed phloem composition
        elif self.phloemComposition == "Shameer2018":
            for met in model.reactions.Phloem_output_tx.metabolites:
                if "PROTON" not in met.id:
                    newMet = met.copy()
                    newMet.compartment = "ph"
                    newMet.name = met.id[:-2] + "_ph"
                    newMet.id = met.id[:-2] + "_ph"
                    model.reactions.Phloem_output_tx.add_metabolites({newMet: -model.reactions.Phloem_output_tx.metabolites[met]})

        # Xylem
        xylem_metabolites = ["KI", "CAII", "MGII", "NITRATE", "Pi", "SULFATE","AMMONIUM"]
        for met in xylem_metabolites:
            rxn = Reaction(met + "_xy")
            rxn.name = met + "_xy"
            rxn.subsystem = "Xylem"
            rxn.lower_bound = 0
            rxn.upper_bound = self.maxFlux
            new_met = model.metabolites.get_by_id(met + "_c").copy()
            new_met.name = met + "_xy"
            new_met.id = met + "_xy"
            new_met.compartment = "xy"
            if met == "Pi":
                rxn.add_metabolites({model.metabolites.get_by_id(met + "_c"): 0.7, model.metabolites.get_by_id("a" + met + "_c"): 0.3,
                                     new_met: -1, model.metabolites.PROTON_e: -3, model.metabolites.PROTON_c: 4.7})
            else:
                rxn.add_metabolites({model.metabolites.get_by_id(met + "_c"): 1, new_met: -1})
            if met == "NITRATE":
                rxn.add_metabolites({model.metabolites.PROTON_e: -2, model.metabolites.PROTON_c: 2})
            elif met == "SULFATE":
                rxn.add_metabolites({model.metabolites.PROTON_e: -3, model.metabolites.PROTON_c: 3})
            elif met == "KI":
                rxn.add_metabolites({model.metabolites.PROTON_e: -1, model.metabolites.PROTON_c: 1})
            model.add_reactions([rxn])

        # Save as attribute
        self.PlantCoreMetabolism_v1_2 = model

        # Set refresh rate
        if self.mode=="daily":
            self.leafModel = self.generateDielLeafModel()
        elif self.mode=="hourly":
            self.leafModel = self.generateHourlyLeafModel()

        # Create sink tissues
        df = pd.read_csv(biomassEquations,index_col="Unnamed: 0")
        for tissue in sink_tissues:
            sink_model = self.PlantCoreMetabolism_v1_2.copy()
            # Make unique identifiers for each reaction adn metabolite. 
            # Reverse direction of phloem and add uptake cost
            for rxn in sink_model.reactions:
                if rxn.id.endswith('_ph'):
                    rxn.upper_bound = 0
                    rxn.lower_bound = -self.maxFlux
                    rxn.add_metabolites({'PROTON_e':2,'PROTON_c':-2})
                if 'Phloem_output_tx' in rxn.id:
                    rxn.upper_bound = 0
                    rxn.lower_bound = -self.maxFlux
                    rxn.add_metabolites({'PROTON_e':1.9207920792,'PROTON_c':-1.9207920792})
                rxn.id = rxn.id + '_' +tissue
            for met in sink_model.metabolites:
                if self.parameters["singlePhloem"] and (met.id.endswith("_xy") or met.id.endswith("_ph")):
                        continue
                met.id = met.id + '_' + tissue
                met.compartment = met.compartment + '_' + tissue

            # Add starch storage reaction
            ## NOTE: added
            starch_storage_rxn = Reaction('Starch_storage_p_'+tissue,name='Starch storage '+tissue,lower_bound = -self.maxDebt,upper_bound=self.maxFlux)
            sink_model.add_reactions([starch_storage_rxn])
            starch_storage_rxn.add_metabolites({'STARCH_p_'+tissue:-1})
            # if tissue == 'FR':
            #     biomass_ratios = df['fruit']
            if tissue == 'ST':
                biomass_ratios = df['stem']
            elif tissue == 'SE':
                biomass_ratios = df['seed']
            elif tissue == 'SH':
                biomass_ratios = df['shell']
            elif tissue == 'R':
                biomass_ratios = df['root']
            biomass_rxn = Reaction("Biomass_tx_" + tissue, name = "Biomass_tx_" + tissue, lower_bound = 0, upper_bound = self.maxFlux)
            sink_model.add_reactions([biomass_rxn])
            for met in df.index:
                if float(biomass_ratios[met])>1e-5:
                    biomass_rxn.add_metabolites({sink_model.metabolites.get_by_id(met+'_'+tissue):-float(biomass_ratios[met])})
            # if tissue == 'FR':
            #     self.shellModel = sink_model
            if tissue == 'ST':
                self.stemModel = sink_model
            elif tissue == 'SE':
                self.seedModel = sink_model
            elif tissue == 'SH':
                self.shellModel = sink_model
            elif tissue == 'R':
                # Allow for soil exchange in root model
                ## Allow for nutrient uptake
                for rxn in ["Pi", "Nitrate", "SO4", "Ca", "Mg", "K"]:
                    sink_model.reactions.get_by_id(rxn + "_tx_R").upper_bound = self.maxFlux
                    sink_model.reactions.get_by_id(rxn + "_tx_R").lower_bound = -self.maxFlux

                # Reverse xylem reaction direction for output into the xylem
                for rxn in sink_model.reactions:
                    if rxn.id.endswith("_xy_R"):
                        rxn.upper_bound = 0
                        rxn.lower_bound = -self.maxFlux
                        #Remove cost of xylem loading
                        for met in rxn.metabolites:
                            if "PROTON" in met.id:
                                rxn.add_metabolites({met: -rxn.metabolites[met]})
                
                # Add soil symbiont
                for type in self.symbiosis:
                    if type == "N":
                        rhizo = cobra.io.read_sbml_model("./../Data/iCC541.xml")
                        tissue = "SY"
                        for met in rhizo.metabolites:
                            met.id = met.id+"_"+tissue
                            met.compartment = met.compartment+"_"+tissue
                        for rxn in rhizo.reactions:
                            rxn.id = rxn.id+"_"+tissue
                            rxn.objective_coefficient=0
                        sink_model = sink_model.merge(rhizo)
                        rxn = Reaction("Alanine_exchange_symbiont")
                        rxn.name = rxn.id.replace("_"," ")
                        rxn.add_metabolites({sink_model.metabolites.get_by_id("L_ALPHA_ALANINE_c_R"):-1,sink_model.metabolites.get_by_id("ala__L[e]_"+tissue):1})
                        rxn.lower_bound = -self.maxFlux
                        rxn.upper_bound = self.maxFlux
                        sink_model.add_reactions([rxn])

                        rxn = Reaction("Aspartate_exchange_symbiont")
                        rxn.name = rxn.id.replace("_"," ")
                        rxn.add_metabolites({sink_model.metabolites.get_by_id("L_ASPARTATE_c_R"):-1,sink_model.metabolites.get_by_id("asp__L[e]_"+tissue):1})
                        rxn.lower_bound = -self.maxFlux
                        rxn.upper_bound = self.maxFlux
                        sink_model.add_reactions([rxn])

                        rxn = Reaction("Glutamate_exchange_symbiont")
                        rxn.name = rxn.id.replace("_"," ")
                        rxn.add_metabolites({sink_model.metabolites.get_by_id("GLT_c_R"):-1,sink_model.metabolites.get_by_id("glu__L[e]_"+tissue):1})
                        rxn.lower_bound = -self.maxFlux
                        rxn.upper_bound = self.maxFlux
                        sink_model.add_reactions([rxn])

                        rxn = Reaction("Malate_exchange_symbiont")
                        rxn.name = rxn.id.replace("_"," ")
                        rxn.add_metabolites({sink_model.metabolites.get_by_id("MAL_c_R"):-1,sink_model.metabolites.get_by_id("mal__L[e]_"+tissue):1})
                        rxn.lower_bound = -self.maxFlux
                        rxn.upper_bound = self.maxFlux
                        sink_model.add_reactions([rxn])

                        rxn = Reaction("Succinate_exchange_symbiont")
                        rxn.name = rxn.id.replace("_"," ")
                        rxn.add_metabolites({sink_model.metabolites.get_by_id("SUC_c_R"):-1,sink_model.metabolites.get_by_id("succ[e]_"+tissue):1})
                        rxn.lower_bound = -self.maxFlux
                        rxn.upper_bound = self.maxFlux
                        sink_model.add_reactions([rxn])

                        rxn = Reaction("Ammonium_exchange_symbiont")
                        rxn.name = rxn.id.replace("_"," ")
                        rxn.add_metabolites({sink_model.metabolites.get_by_id("AMMONIUM_c_R"):-1,sink_model.metabolites.get_by_id("fixedNH3[e]_"+tissue):1})
                        rxn.lower_bound = -self.maxFlux
                        rxn.upper_bound = self.maxFlux
                        sink_model.add_reactions([rxn])

                        sink_model.reactions.get_by_id("EX_oxygen_"+tissue).lower_bound = -self.maxFlux
                        for rxn in sink_model.reactions:
                            if "EX_" in rxn.id and rxn.id not in ["EX_fe2_SY","EX_nitrogen_SY","EX_co2_SY","EX_orthophosphate_SY","EX_sulfate_SY","EX_oxygen_SY","EX_h2o_SY","EX_H_SY","EX_mg_SY","EX_cbl1_SY","EX_molybdate_SY","EX_thiamine_SY","EX_zn_SY","EX_h2_SY"]:
                                rxn.lower_bound = 0
                                rxn.upper_bound = 0
                self.rootModel = sink_model

    def generateDielLeafModel(self,biomassEquations = "./Parameters/tomatoBiomassEquations.csv", biomassEquationsFruit = "./Parameters/shellBiomassComposition.csv"):
        """
        This function generates a diel Leaf model
        """

        #Create diel source leaf model
        diel_model = Model('Leaf_model')
        model = self.PlantCoreMetabolism_v1_2.copy()
        for met in model.metabolites:
            if self.parameters["singlePhloem"] and (met.id.endswith("_xy") or met.id.endswith("_ph")):
                    continue
            met.id = met.id + "_LD"
            met.compartment = met.compartment + "_LD"
        for rxn in model.reactions:
            rxn_copy = rxn.copy()
            rxn_copy.id = rxn_copy.id + "_LD"
            diel_model.add_reactions([rxn_copy])
        model = self.PlantCoreMetabolism_v1_2.copy()
        for met in model.metabolites:
            met.id = met.id + "_LN"
            met.compartment = met.compartment + "_LN"
        for rxn in model.reactions:
            rxn_copy = rxn.copy()
            rxn_copy.id = rxn_copy.id + "_LN"
            diel_model.add_reactions([rxn_copy])

        #Adding linker reactions between day and night metabolites that are allowed to accumulate.
        linker_reactions = ["STARCH", "NITRATE", "GLC", "FRU"]
        for linker in linker_reactions:
            if linker == "STARCH":
                met_day = diel_model.metabolites.get_by_id(linker + "_p_LD")
                met_night = diel_model.metabolites.get_by_id(linker + "_p_LN")
                rxn = Reaction(linker + "_p_linker_LD")
            else:
                met_day = diel_model.metabolites.get_by_id(linker + "_v_LD")
                met_night = diel_model.metabolites.get_by_id(linker + "_v_LN")
                rxn = Reaction(linker + "_v_linker_LD")
            rxn.name = linker + " day to night linker"
            rxn.subsystem = "Day-Night linker reaction"
            rxn.lower_bound = -self.maxFlux
            rxn.upper_bound = self.maxFlux
            rxn.add_metabolites({met_day: -1.0, met_night: 1.0})
            diel_model.add_reactions([rxn])
        rxn = Reaction("MAL_v_linker_LD", name = "MAL day to night linker", lower_bound = -self.maxFlux, upper_bound = self.maxFlux)
        rxn.add_metabolites({diel_model.metabolites.MAL_v_LD: -0.7, diel_model.metabolites.aMAL_v_LD: -0.3, diel_model.metabolites.MAL_v_LN: 0.7, diel_model.metabolites.aMAL_v_LN: 0.3})
        diel_model.add_reactions([rxn])
        rxn = Reaction("CIT_v_linker_LD", name = "CIT day to night linker", lower_bound = -self.maxFlux, upper_bound = self.maxFlux)
        rxn.add_metabolites({diel_model.metabolites.CIT_v_LD: -0.5, diel_model.metabolites.aCIT_v_LD: -0.5, diel_model.metabolites.CIT_v_LN: 0.5, diel_model.metabolites.aCIT_v_LN: 0.5})
        diel_model.add_reactions([rxn])
        Amino_acids = ["L_ALPHA_ALANINE", "L_ASPARTATE", "ARG", "ASN", "CYS", "GLN", "GLT", "GLY", "ILE",
                       "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "bHIS"]
        for aa in Amino_acids:
            rxn = Reaction(aa + "_v_linker_LD", name = aa + " day to night linker", subsystem = "Day-Night linker reaction", lower_bound = 0, upper_bound = self.maxFlux)
            rxn.add_metabolites({diel_model.metabolites.get_by_id(aa + "_v_LD"): -1.0, diel_model.metabolites.get_by_id(aa + "_v_LN"): 1.0})
            diel_model.add_reactions([rxn])

        #Fix the nitrate uptake ratio to 3:2 (day:night).
        Nitrate_f = Metabolite('Nitrate_f_L', name='Nitrate pseudo-metabolite', compartment='f_L')
        diel_model.reactions.NITRATE_xy_LD.add_metabolites({Nitrate_f: 2.0})
        diel_model.reactions.NITRATE_xy_LN.add_metabolites({Nitrate_f: -3.0})

        #Fix the Vc/Vo ratio to 3:1 using pseudo-metabolites for both day and night.
        Rubisco_f = Metabolite('Rubisco_f_LD', name='Rubisco Vc/Vo day pseudo-metabolite', compartment='f_LD')
        diel_model.reactions.RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p_LD.add_metabolites({Rubisco_f: 1.0})
        diel_model.reactions.RXN_961_p_LD.add_metabolites({Rubisco_f: -3.0})

        #Add biomass reaction
        biomass_file = pd.read_csv(biomassEquations,header=None)
        diel_model.add_reactions([Reaction("Biomass_tx_LD", name = "Biomass_tx_LD", lower_bound = 0, upper_bound = self.maxFlux)])
        for it in range(len(biomass_file)):
            met_id = biomass_file[0][it]+"_LD"
            stoich = float(biomass_file[1][it])
            diel_model.reactions.get_by_id("Biomass_tx_LD").add_metabolites({diel_model.metabolites.get_by_id(met_id): -stoich})
        diel_model.add_reactions([Reaction("Biomass_tx_LN", name = "Biomass_tx_LN", lower_bound = 0, upper_bound = self.maxFlux)])
        for it in range(len(biomass_file)):
            met_id = biomass_file[0][it]+"_LN"
            stoich = float(biomass_file[1][it])
            diel_model.reactions.get_by_id("Biomass_tx_LN").add_metabolites({diel_model.metabolites.get_by_id(met_id): -stoich})

        #Constrain starch phosphorylase reaction in leaf
        diel_model.reactions.G6P_Pi_pc_LD.upper_bound = diel_model.reactions.G6P_Pi_pc_LD.lower_bound = 0
        diel_model.reactions.G6P_Pi_pc_LN.upper_bound = diel_model.reactions.G6P_Pi_pc_LN.lower_bound = 0

        #Allow photon uptake during the day
        diel_model.reactions.Photon_tx_LD.upper_bound = self.maxFlux

        #Create reverse phloem reactions for uptake from the phloem
        if self.phloemComposition == "free":
            revRxnDict = dict()
            for rxn in diel_model.reactions:
                if "_ph" in rxn.id:
                    revRxn = Reaction(rxn.id.split("_ph")[0] + "_rev" + rxn.id.split(rxn.id.split("_ph")[0])[1], name = rxn.id.split("_ph")[0] + "_rev" + rxn.id.split(rxn.id.split("_ph")[0])[1], upper_bound = self.maxFlux, lower_bound = 0, subsystem = "Phloem")
                    for met in rxn.metabolites:
                        if "PROTON" not in met.id:
                            revRxn.add_metabolites({met: -rxn.metabolites[met]})
                        else:
                            revRxn.add_metabolites({met: rxn.metabolites[met]})
                    revRxnDict[revRxn] = revRxn
            for rxn in revRxnDict:
                diel_model.add_reactions([revRxnDict[rxn]])

            #Phloem day-night 3:1 constraint
            Phloem_output_day = Metabolite("Phloem_LD_f", name = "Phloem output LD", compartment = "f")
            Phloem_output_night = Metabolite("Phloem_LN_f", name = "Phloem output LN", compartment = "f")
            Phloem_output_rxn = Reaction("Phloem_LD_LN_constraint_LD", name = "Phloem_LD_LN_constraint", upper_bound = self.maxFlux, lower_bound = -self.maxFlux)
            diel_model.add_reactions([Phloem_output_rxn])
            diel_model.reactions.get_by_id("Phloem_LD_LN_constraint_LD").add_metabolites({Phloem_output_day: -3, Phloem_output_night: -1})
            diel_model.reactions.Phloem_output_tx_LD.add_metabolites({Phloem_output_day: 1})
            diel_model.reactions.Phloem_output_rev_tx_LD.add_metabolites({Phloem_output_day: -1})
            diel_model.reactions.Phloem_output_tx_LN.add_metabolites({Phloem_output_night: 1})
            diel_model.reactions.Phloem_output_rev_tx_LN.add_metabolites({Phloem_output_night: -1})

        elif self.phloemComposition == "Shameer2018":
            for rxn in [diel_model.reactions.Phloem_output_tx_LD, diel_model.reactions.Phloem_output_tx_LN]:
                revRxn = Reaction(rxn.id.split("_tx")[0] + "_rev" + rxn.id.split(rxn.id.split("_tx")[0])[1], name = rxn.id.split("_tx")[0] + "_rev" + rxn.id.split(rxn.id.split("_tx")[0])[1], upper_bound = self.maxFlux, lower_bound = 0, subsystem = "Phloem")
                for met in rxn.metabolites:
                    if "PROTON" not in met.id:
                        revRxn.add_metabolites({met: -rxn.metabolites[met]})
                    else:
                        revRxn.add_metabolites({met: rxn.metabolites[met]})
                diel_model.add_reactions([revRxn])

            #Phloem day-night 3:1 constraint
            Phloem_output_day = Metabolite("Phloem_LD_f", name = "Phloem output LD", compartment = "f")
            Phloem_output_night = Metabolite("Phloem_LN_f", name = "Phloem output LN", compartment = "f")
            diel_model.add_reactions([Reaction("Phloem_LD_LN_constraint_LD", name = "Phloem_LD_LN_constraint", upper_bound = self.maxFlux, lower_bound = -self.maxFlux)])
            diel_model.reactions.Phloem_LD_LN_constraint_LD.add_metabolites({Phloem_output_day: -3, Phloem_output_night: -1})
            diel_model.reactions.Phloem_output_tx_LD.add_metabolites({Phloem_output_day: 1})
            diel_model.reactions.Phloem_output_rev_tx_LD.add_metabolites({Phloem_output_day: -1})
            diel_model.reactions.Phloem_output_tx_LN.add_metabolites({Phloem_output_night: 1})
            diel_model.reactions.Phloem_output_rev_tx_LN.add_metabolites({Phloem_output_night: -1})

        return diel_model

    def generateHourlyLeafModel(self,biomassEquations = "./Parameters/tomatoBiomassEquations.csv", biomassEquationsFruit = "./Parameters/shellBiomassComposition.csv"):
        """
        This function generate an hourly leaf model
        """
        hour_model = Model("Leaf_model")
        model = self.PlantCoreMetabolism_v1_2.copy()
        for met in model.metabolites:
            if self.parameters["singlePhloem"] and (met.id.endswith("_xy") or met.id.endswith("_ph")):
                    continue
            met.id = met.id + "_L"
            met.compartment = met.compartment + "_L"
        for rxn in model.reactions:
            rxn_copy = rxn.copy()
            rxn_copy.id = rxn_copy.id + "_L"
            hour_model.add_reactions([rxn_copy])

        linker_reactions = ["STARCH", "NITRATE", "GLC", "FRU"]
        for linker in linker_reactions:
            if linker == "STARCH":
                met = hour_model.metabolites.get_by_id(linker + "_p_L")
                rxn = Reaction(linker + "_p_linker")
            else:
                met = hour_model.metabolites.get_by_id(linker + "_v_L")
                rxn = Reaction(linker + "_v_linker")
            rxn.name = linker + " accumulation"
            rxn.subsystem = "Day-Night linker reaction"
            rxn.lower_bound = 0
            rxn.upper_bound = 0
            rxn.add_metabolites({met: -1.0})
            hour_model.add_reactions([rxn])
        rxn = Reaction("MAL_v_linker", name = "MAL day to night linker", lower_bound = 0, upper_bound = 0)
        rxn.add_metabolites({hour_model.metabolites.MAL_v_L: -0.7, hour_model.metabolites.aMAL_v_L: -0.3})
        hour_model.add_reactions([rxn])
        rxn = Reaction("CIT_v_linker", name = "CIT day to night linker", lower_bound = 0, upper_bound = 0)
        rxn.add_metabolites({hour_model.metabolites.CIT_v_L: -0.5, hour_model.metabolites.aCIT_v_L: -0.5})
        hour_model.add_reactions([rxn])
        Amino_acids = ["L_ALPHA_ALANINE", "L_ASPARTATE", "ARG", "ASN", "CYS", "GLN", "GLT", "GLY", "ILE",
                       "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "bHIS"]
        for aa in Amino_acids:
            rxn = Reaction(aa + "_v_linker", name = aa + " day to night linker", subsystem = "Day-Night linker reaction", lower_bound = 0, upper_bound = 0)
            rxn.add_metabolites({hour_model.metabolites.get_by_id(aa + "_v_L"): -1.0})
            hour_model.add_reactions([rxn])

        #Fix the Vc/Vo ratio to 3:1 using pseudo-metabolites for both day and night.
        Rubisco_f = Metabolite('Rubisco_f_L', name='Rubisco Vc/Vo day pseudo-metabolite', compartment='f_L')
        hour_model.reactions.RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p_L.add_metabolites({Rubisco_f: 1.0})
        hour_model.reactions.RXN_961_p_L.add_metabolites({Rubisco_f: -3.0})

        #Add biomass reaction
        biomass_file = pd.read_csv(biomassEquations,header=None)
        hour_model.add_reactions([Reaction("Biomass_tx_L", name = "Biomass_tx_L", lower_bound = 0, upper_bound = self.maxFlux)])
        for it in range(len(biomass_file)):
            met_id = biomass_file[0][it]+"_L"
            stoich = float(biomass_file[1][it])
            hour_model.reactions.get_by_id("Biomass_tx_L").add_metabolites({hour_model.metabolites.get_by_id(met_id): -stoich})

        #Constrain starch phosphorylase reaction in leaf
        hour_model.reactions.G6P_Pi_pc_L.upper_bound = hour_model.reactions.G6P_Pi_pc_L.lower_bound = 0

        #Allow photon uptake during the day
        hour_model.reactions.Photon_tx_L.upper_bound = self.maxFlux

        #Create reverse phloem reactions for uptake from the phloem
        ## NOTE: this might need to be rethought...
        if self.phloemComposition == "free":
            revRxnDict = dict()
            for rxn in hour_model.reactions:
                if "_ph" in rxn.id:
                    revRxn = Reaction(rxn.id.split("_ph")[0] + "_rev" + rxn.id.split(rxn.id.split("_ph")[0])[1], name = rxn.id.split("_ph")[0] + "_rev" + rxn.id.split(rxn.id.split("_ph")[0])[1], upper_bound = self.maxFlux, lower_bound = 0, subsystem = "Phloem")
                    for met in rxn.metabolites:
                        if "PROTON" not in met.id:
                            revRxn.add_metabolites({met: -rxn.metabolites[met]})
                        else:
                            revRxn.add_metabolites({met: rxn.metabolites[met]})
                    revRxnDict[revRxn] = revRxn
            for rxn in revRxnDict:
                hour_model.add_reactions([revRxnDict[rxn]])

        elif self.phloemComposition == "Shameer2018":
            for rxn in [hour_model.reactions.Phloem_output_tx_L]:
                revRxn = Reaction(rxn.id.split("_tx")[0] + "_rev" + rxn.id.split(rxn.id.split("_tx")[0])[1], name = rxn.id.split("_tx")[0] + "_rev" + rxn.id.split(rxn.id.split("_tx")[0])[1], upper_bound = self.maxFlux, lower_bound = 0, subsystem = "Phloem")
                for met in rxn.metabolites:
                    if "PROTON" not in met.id:
                        revRxn.add_metabolites({met: -rxn.metabolites[met]})
                    else:
                        revRxn.add_metabolites({met: rxn.metabolites[met]})
                hour_model.add_reactions([revRxn])
        return hour_model

    def calculateBiomassWeight(self, rxnID = "Biomass_tx_"):

        """Calculate the molecular weight of biomass to convert flux from mols to g biomass

        Parameters
        ----------
        - Organ specific stoichiometric model
        - Reaction ID of biomass reaction

        Returns
        -------
        - Molar mass of total biomass components for the organ of interest (g/mmol)

        """
        if self.mode == "daily":
            self.leafBiomassWeight = sum([abs(met.formula_weight*self.leafModel.reactions.get_by_id(rxnID + "LD").metabolites[met]) for met in self.leafModel.reactions.get_by_id(rxnID + "LD").metabolites])/1000
        else:
            self.leafBiomassWeight = sum([abs(met.formula_weight*self.leafModel.reactions.get_by_id(rxnID + "L").metabolites[met]) for met in self.leafModel.reactions.get_by_id(rxnID + "L").metabolites])/1000

        self.stemBiomassWeight = sum([abs(met.formula_weight*self.stemModel.reactions.get_by_id(rxnID + "ST").metabolites[met]) for met in self.stemModel.reactions.get_by_id(rxnID + "ST").metabolites])/1000

        self.rootBiomassWeight = sum([abs(met.formula_weight*self.rootModel.reactions.get_by_id(rxnID + "R").metabolites[met]) for met in self.rootModel.reactions.get_by_id(rxnID + "R").metabolites])/1000

        self.shellBiomassWeight = sum([abs(met.formula_weight*self.shellModel.reactions.get_by_id(rxnID + "SH").metabolites[met]) for met in self.shellModel.reactions.get_by_id(rxnID + "SH").metabolites])/1000

    def readStrucutrefromMask(self, datafile = "./Parameters/soyMask_singleLeafStem.csv"):

            """
            Function to import and read mask file for plant structure, specifying organ
            position along the developmental axis.

            Parameters
            ----------
            - datafile : GREENLAB mask file
                A GREENLAB mask file containing information on organ emergence in thermal time

            Returns
            -------
            - A list of time points of leaf emergence
            - A list of time points of stem emergence
            - A list of time points of fruit emergence
            - A list of time points of seed emergence
            - A list of expansion time scaling of leaves at each phytomer unit
            - A list of expansion time scaling of stems at each phytomer unit
            - A list of expansion time scaling of fruits at each phytomer unit
            
            """

            #Import positions
            df = pd.read_csv(datafile, delimiter=",")

            self.parameters["leafPosition"] = [int(df["leaf"][n]) for n in range(0, len(df))]
            self.parameters["stemPosition"] = [int(df["stem"][n]) for n in range(0, len(df))]
            self.parameters["shellPosition"] = [int(df["shell"][n]) for n in range(0, len(df))]
            self.parameters["seedPosition"] = [int(df["seed"][n]) for n in range(0, len(df))]

            self.parameters["leafExpScale"] = [float(df["leafExpScale"][n]) for n in range(0, len(df))]
            self.parameters["stemExpScale"] = [float(df["stemExpScale"][n]) for n in range(0, len(df))]

    def readConnections(self, datafile = "./Parameters/Organ_connections.csv"):
        #Import connections
        df = pd.read_csv(datafile, delimiter=",",header=None)
        tempDict = dict()
        for i in range(0,len(df)):
            tempDict[df.iloc[i][0]]=df.iloc[i][1]
        self.connections=tempDict

    def generateModelAtEmergence(self):

        """
        Generate whole plant model for soybean at emergence (1 leaf, 1 stem, 1 root)

        Parameters
        ----------
        - Organ-specific stoichiometric models (leaf, root, stem)
        - Initial organ biomasses (leaf, root, stem)

        Returns
        -------
        - Whole-plant model at emergence

        """
        plantModel = Model()
        plantModel.solver="glpk"
        for rxn in plantModel.reactions:
            if rxn.objective_coefficient!=0:
                rxn.objective_coefficient = 0
        t = 1
        plantModel.add_reactions([Reaction("Total_biomass_tx", name = "Total_biomass_tx", upper_bound = self.maxFlux, lower_bound= 0)])

        for model in [self.leafModel, self.rootModel, self.stemModel]:
            # Rename models to make it easier to add organs
            if model == self.leafModel:
                if self.mode == "daily":
                    organTag = "LD"
                else:
                    organTag = "L"
                model = self.leafModel.copy()
            elif model == self.rootModel:
                organTag = "R"
                model = self.rootModel.copy()
            elif model == self.stemModel:
                organTag = "ST"
                model = self.stemModel.copy()
            for rxn in model.reactions:
                rxn.id = rxn.id + str(t)
            for met in model.metabolites:
                if self.parameters["singlePhloem"]:
                    if met.id.endswith("_xy") or met.id.endswith("_ph"):
                        continue
                met.id = met.id + str(t)
                met.compartment = met.compartment + str(t)
            BiomassMet = Metabolite("Biomass_" + organTag + str(t), name = "Biomass_" + organTag + str(t), compartment="b_" + organTag + str(t))
            model.reactions.get_by_id("Biomass_tx_" + organTag + str(t)).add_metabolites({BiomassMet: 1})
            if organTag == "LD":
                BiomassMetLN = Metabolite("Biomass_" + "LN" + str(t), name = "Biomass_" + "LN" + str(t), compartment="b_" + "LN" + str(t))
                model.reactions.get_by_id("Biomass_tx_" + "LN" + str(t)).add_metabolites({BiomassMetLN: 1})
            plantModel = plantModel.merge(model)
            plantModel.reactions.Total_biomass_tx.add_metabolites({BiomassMet: -1})
            if organTag == "LD":
                plantModel.reactions.Total_biomass_tx.add_metabolites({BiomassMetLN: -1})

        # plantModel.reactions.Total_biomass_tx.objective_coefficient = 1

        self.stoichiometricModel = plantModel
        self.addRoot([Organ("root" + str(t), position=t, mass = self.parameters["RbInit"], age = 12, tag = "R"+str(t))])
        self.addStem([Organ("stem" + str(t), position=t, mass = self.parameters["SbInit"], age = 12, connection=self.connections["ST"+str(t)], tag = "ST"+str(t))])


        if not self.parameters["singlePhloem"]:
            for met in plantModel.metabolites:
                if met.compartment == "ph_ST1":
                    met2 = plantModel.metabolites.get_by_id(met.id.replace("ST"+str(t),self.connections["ST"+str(t)]))
                    rxn = Reaction(met.id+"_"+met2.id+"_exchange")
                    rxn.add_metabolites({met:-1,met2:1})
                    rxn.upper_bound = self.maxFlux
                    rxn.lower_bound = -self.maxFlux
                    plantModel.add_reactions([rxn])
                if met.compartment == "xy_ST1":
                    met2 = plantModel.metabolites.get_by_id(met.id.replace("ST"+str(t),self.connections["ST"+str(t)]))
                    rxn = Reaction(met.id+"_"+met2.id+"_exchange")
                    rxn.add_metabolites({met:-1,met2:1})
                    rxn.upper_bound = self.maxFlux
                    rxn.lower_bound = -self.maxFlux
                    plantModel.add_reactions([rxn])
        self.addLeaf([Organ("leaf" + str(t), position=t, mass = self.parameters["LbInit"], age = 12, connection=self.connections["L"+str(t)], tag = "L"+str(t))])
        print('Leaf added: ',list(self.leaves))
        if not self.parameters["singlePhloem"]:
            for met in plantModel.metabolites:
                if met.compartment == "ph_L1":
                    met2 = plantModel.metabolites.get_by_id(met.id.replace("L"+str(t),self.connections["L"+str(t)]))
                    rxn = Reaction(met.id+"_"+met2.id+"_exchange")
                    rxn.add_metabolites({met:-1,met2:1})
                    rxn.upper_bound = self.maxFlux
                    rxn.lower_bound = -self.maxFlux
                    plantModel.add_reactions([rxn])
                if met.compartment == "xy_L1":
                    met2 = plantModel.metabolites.get_by_id(met.id.replace("L"+str(t),self.connections["L"+str(t)]))
                    rxn = Reaction(met.id+"_"+met2.id+"_exchange")
                    rxn.add_metabolites({met:-1,met2:1})
                    rxn.upper_bound = self.maxFlux
                    rxn.lower_bound = -self.maxFlux
                    plantModel.add_reactions([rxn])

        #Initiate total biomass for organ type attribute
        self.leafBiomass = self.leaves.leaf1.mass
        self.rootBiomass = self.roots.root1.mass
        self.stemBiomass = self.stems.stem1.mass
        self.reproductiveSwitch = 0
        # self.shellBiomass = 0

    def addPhloemCtracker(self,CSinks=["R1","ST1",'SH1']):
        for tag in CSinks:
            Cmet = Metabolite("Cmet_"+tag)
            print("------------ADD Cmet SINK REACTION-----------------")
            self.stoichiometricModel.add_reactions([Reaction("Cmet_sink_"+tag,name=Cmet.id+" sink",lower_bound=-self.maxFlux,upper_bound=self.maxFlux)])
            self.stoichiometricModel.reactions.get_by_id("Cmet_sink_"+tag).add_metabolites({Cmet:1})

            for rxn in self.stoichiometricModel.reactions.query("_ph_"+tag):
                CmetN = 0
                for met in rxn.reactants:
                    if "C" in met.formula:
                        CmetN = CmetN+abs(float((met.formula.split("H")[0][1:]))*rxn.metabolites[met])
                rxn.add_metabolites({Cmet:CmetN})

        Cmet = Metabolite("Cmet_SE1")
        self.stoichiometricModel.add_reactions([Reaction("Cmet_sink_SE1",name=Cmet.id+" sink",lower_bound=-self.maxFlux,upper_bound=self.maxFlux)])
        self.stoichiometricModel.reactions.get_by_id("Cmet_sink_SE1").add_metabolites({Cmet:1})

        CmetN = 0
        for rxn in self.stoichiometricModel.reactions.query("_ph_"+tag):
            for met in rxn.reactants:
                if "C" in met.formula:
                    CmetN = CmetN+abs(float((met.formula.split("H")[0][1:]))*rxn.metabolites[met])
            rxn.add_metabolites({Cmet:CmetN})


    ## I'm not sure what's happening with this reaciton...
    def replaceBiomassMetabolitewithGDW(self):
        if len(self.stoichiometricModel.metabolites.query("Biomass_"))==0:
            return
        for leaf in self.leaves:
            self.stoichiometricModel.metabolites.get_by_id("Biomass_L"+str(leaf.position)).remove_from_model()
            Cmet = Metabolite("Cmet_L"+str(leaf.position),compartment="c_L"+str(leaf.position))
            leafBiomassCN = abs(self.stoichiometricModel.reactions.Biomass_tx_L1.check_mass_balance()["C"])
            self.stoichiometricModel.reactions.Biomass_tx_L1.add_metabolites({Cmet:leafBiomassCN})

            rxn = Reaction(Cmet.id+"_sink")
            rxn.add_metabolites({Cmet:-1})
            rxn.lower_bound = 0
            rxn.upper_bound = self.maxFlux #10000
            self.stoichiometricModel.add_reactions([rxn])

        self.stoichiometricModel.metabolites.get_by_id("Biomass_R1").remove_from_model()
        #self.stoichiometricModel.reactions.get_by_id("Biomass_tx_R1").objective_coefficient=1
        self.stoichiometricModel.metabolites.get_by_id("Biomass_ST1").remove_from_model()
        #self.stoichiometricModel.reactions.get_by_id("Biomass_tx_ST1").objective_coefficient=1
    
    def permitBiomassRemobilization(self,tags=["L1","R1","ST1","SE1",'SH1'],protein_degradation_ATPcost=1):
        self.replaceBiomassMetabolitewithGDW()
        for tag in tags:
            Rxn = self.stoichiometricModel.reactions.get_by_id("Biomass_tx_"+tag)
            print(Rxn) ## NOTE
            #add a new metabolite to represent biomass gDW
            MW = 0
            for met in Rxn.metabolites:
                MW = MW+(met.formula_weight*Rxn.metabolites[met]*-1)
            GDW = Metabolite("gDW_"+tag,name="1 gWD",compartment="c_"+tag)
            Rxn.add_metabolites({GDW:MW/1000})

            GDWRXN = Reaction("GDWsink_"+tag,name="GDWsink")
            GDWRXN.add_metabolites({GDW:-1})
            GDWRXN.lower_bound = 0
            GDWRXN.upper_bound = 1000
            self.stoichiometricModel.add_reactions([GDWRXN])

            # Add sink for every biomass metabolite
            for met in Rxn.metabolites.keys():
                if met == GDW or "Cmet" in met.id:
                    continue
                metID = met.id.replace("_"+tag,"")
                rxn = Reaction(metID+"_sink_"+tag,name=met.id+" sink")
                rxn.add_metabolites({self.stoichiometricModel.metabolites.get_by_id(met.id):-1, GDW:self.stoichiometricModel.metabolites.get_by_id(met.id).formula_weight/1000})
                rxn.lower_bound = 0
                rxn.upper_bound = 0
                self.stoichiometricModel.add_reactions([rxn])
                # NOTE: this is weird. Why is GLC_c treated differently to sGLC_b
                if any([y == metID[:len(y)] for y in ["GLC_c","FRU_c","L_PHOSPHATIDATE_p","LYS_c","THR_c","CYS_c","GLY_c","PRO_c","VAL_c","ILE_c","LEU_c","MET_c","PHE_c","TYR_c","SUC_c","MYO_INOSITOL_c",
                            "HIS_c","TRP_c","Biomass","Cmet","gDW","PALMITATE_p","CPD_9245_p","STEARIC_ACID_p","OLEATE_CPD_p","Octadecadienoate_p","ARACHIDIC_ACID_p","CPD_16709_p","DOCOSANOATE_p"]]):
                    continue
                elif any([y == metID[:len(y)] for y in ["sSUCROSE_b","Starch_b","Cellulose_b","Xylan_b","sMAL_b","sCIT_b","sFUM_b","sASP_b","sGLU_b","sGLN_b","sSER_b","sALA_b","sGABA_b","sGLC_b","sFRU_b"]+['sTHR_b', 'sCYS_b', 'sGLY_b', 'sLEU_b', 'sTYR_b', 'sSUC_b','sASN_b','sFRUCTOSE_6P_b','sGLN_b','sHIS_b','sILE_b','sLYS_b','sPRO_b','sVAL_b']]):
                    self.stoichiometricModel.reactions.get_by_id(met.id.replace("_b_","_biomass_")).lower_bound = -1000
                    self.stoichiometricModel.reactions.get_by_id(met.id.replace("_b_","_biomass_")).upper_bound = 0
                elif any([y in metID for y in ["FattyAcid_b",]]):
                    product = self.stoichiometricModel.metabolites.get_by_id("PALMITATE_c_"+tag)
                    rxn = Reaction("FattyAcid_degradation_"+tag)
                    rxn.add_metabolites({met:-1,product:1})
                    rxn.lower_bound = 0
                    self.stoichiometricModel.add_reactions([rxn])
                else:
                    METID = "s"+metID[1:]+"_"+tag
                    Met = self.stoichiometricModel.metabolites.get_by_id(METID)
                    rxn = Reaction(metID+"_degradation_"+tag)
                    self.stoichiometricModel.add_reactions([rxn])
                    rxn.add_metabolites({met.id:-1,Met.id:1})
                    Arxn = self.stoichiometricModel.reactions.get_by_id("ATPase_tx_"+tag)
                    for Amet in Arxn.metabolites.keys():
                        if "ATP_NADPH_f" in Amet.id:
                            continue
                        cost = protein_degradation_ATPcost
                        rxn.add_metabolites({Amet:cost*Arxn.metabolites[Amet]})
                    rxn.lower_bound = 0
                    rxn.upper_bound = 1000
                    # self.stoichiometricModel.add_reactions([rxn])
                    self.stoichiometricModel.reactions.get_by_id(METID.replace("_b_","_biomass_")).lower_bound = -1000
                    self.stoichiometricModel.reactions.get_by_id(METID.replace("_b_","_biomass_")).upper_bound = 1000

        self.stoichiometricModel.reactions._generate_index()
        self.stoichiometricModel.metabolites._generate_index()

    def inspectModel(self):
        ## Prints out a list of all reactions with non-default bounds or missing compartments
        for rxn in self.stoichiometricModel.reactions:
            if not(rxn.lower_bound == 0 or rxn.lower_bound == -self.maxFlux or rxn.lower_bound == -1000):
                print(rxn.id)
                print("LB="+str(rxn.lower_bound))
            if not(rxn.upper_bound == 0 or rxn.upper_bound == self.maxFlux or rxn.upper_bound == 1000):
                print(rxn.id)
                print("UB="+str(rxn.upper_bound))
        for met in self.stoichiometricModel.metabolites:
            if met.compartment == "" or met.compartment is None:
                metIDparts = met.id.split("_")
                met.compartment = "f_"+metIDparts[len(metIDparts)-1]
                print(met.id)
                print(met.compartment)
                #exit()
        self.stoichiometricModel.reactions.Total_biomass_tx.upper_bound = 1000
        print("GDW sink reactions !!!")
        print(self.stoichiometricModel.reactions.query("GDWsink"))

    def createNewOrgansAlongAxis(self):

        """
        Function to add organs to self.stoichiometricModel based on modeling time

        Parameters
        ----------
        - A list of time points of emergence for organ of interest
        - The whole plant model
        - The cobra.core.Model object of organ of interest
        - Thermal time point

        Returns
        -------
        - Updated whole plant model

        Notes: changed PlantModel to self.stoichiometricModel
        
        """
        t = int(self.growthCycle)

        # If the model mask has a new organ, add a new organ model
        print(t)
        print('Leaf: ',self.parameters["leafPosition"][t-1],' Stem: ',self.parameters["stemPosition"][t-1],' Shell: ',self.parameters["shellPosition"][t-1],' Seed: ',self.parameters["seedPosition"][t-1])
        if self.parameters["leafPosition"][t-1]:
            modelNew = self.leafModel.copy()
            organTag = "LD"
            for met in modelNew.metabolites:
                if not(met.id.endswith("_xy")) and not(met.id.endswith("_ph")):
                    met.id = met.id + str(t)
                    if met.compartment != "f":
                        met.compartment = met.compartment + str(t)
            BiomassMetLD = Metabolite("Biomass_" + organTag + str(t), name = "Biomass_" + organTag + str(t), compartment="b_" + organTag + str(t))
            BiomassMetLN = Metabolite("Biomass_" + "LN" + str(t), name = "Biomass_" + "LN" + str(t), compartment="b_" + "LN" + str(t))
            for rxn in modelNew.reactions:
                rxn.id = rxn.id + str(t)
                if rxn.id == "Biomass_tx_" + organTag + str(t):
                    rxn.add_metabolites({BiomassMetLD: 1})
                if rxn.id == "Biomass_tx_" + "LN" + str(t):
                    rxn.add_metabolites({BiomassMetLN: 1})
                self.stoichiometricModel.add_reactions([rxn])
            self.stoichiometricModel.reactions.Total_biomass_tx.add_metabolites({BiomassMetLD: -1, BiomassMetLN: -1})
            self.addLeaf([Organ("leaf" + str(t), position=t, tag = "L"+str(t),connection=self.connections["L"+str(t)])])

            if not self.parameters["singlePhloem"]:
                for met in self.stoichiometricModel.metabolites:
                    if met.compartment == "ph_L1":
                        met2 = self.stoichiometricModel.metabolites.get_by_id(met.id.replace("L"+str(t),self.connections["L"+str(t)]))
                        rxn = Reaction(met.id+"_"+met2.id+"_exchange")
                        rxn.add_metabolites({met:-1,met2:1})
                        rxn.upper_bound = 1000
                        rxn.lower_bound = -1000
                        self.stoichiometricModel.add_reactions([rxn])
                    if met.compartment == "xy_L1":
                        met2 = self.stoichiometricModel.metabolites.get_by_id(met.id.replace("L"+str(t),self.connections["L"+str(t)]))
                        rxn = Reaction(met.id+"_"+met2.id+"_exchange")
                        rxn.add_metabolites({met:-1,met2:1})
                        rxn.upper_bound = 1000
                        rxn.lower_bound = -1000
                        self.stoichiometricModel.add_reactions([rxn])

        if self.parameters["stemPosition"][t-1]:
            self.addStem([Organ("stem" + str(t), position=t, connection=self.connections["ST"+str(t)], tag = "ST"+str(t))])
            if not self.parameters["singlePhloem"]:
                for met in self.stoichiometricModel.metabolites:
                    if met.compartment == "ph_ST1":
                        met2 = self.stoichiometricModel.metabolites.get_by_id(met.id.replace("ST"+str(t),self.connections["ST"+str(t)]))
                        rxn = Reaction(met.id+"_"+met2.id+"_exchange")
                        rxn.add_metabolites({met:-1,met2:1})
                        rxn.upper_bound = 1000
                        rxn.lower_bound = -1000
                        self.stoichiometricModel.add_reactions([rxn])
                    if met.compartment == "xy_ST1":
                        met2 = self.stoichiometricModel.metabolites.get_by_id(met.id.replace("ST"+str(t),self.connections["ST"+str(t)]))
                        rxn = Reaction(met.id+"_"+met2.id+"_exchange")
                        rxn.add_metabolites({met:-1,met2:1})
                        rxn.upper_bound = 1000
                        rxn.lower_bound = -1000
                        self.stoichiometricModel.add_reactions([rxn])

        if self.parameters["shellPosition"][t-1]:
            modelNew = self.shellModel.copy()
            organTag = "SH"
            for met in modelNew.metabolites:
                if not(met.id.endswith("_xy")) and not(met.id.endswith("_ph")):
                    met.id = met.id + str(1)
                    if met.compartment != "f":
                        met.compartment = met.compartment + '1'
            BiomassMet = Metabolite("Biomass_" + organTag + str(1), name = "Biomass_" + organTag + str(1), compartment="b_" + organTag + str(1))
            for rxn in modelNew.reactions:
                rxn.id = rxn.id + str(1)
                if rxn.id == "Biomass_tx_" + organTag + str(1):
                    rxn.add_metabolites({BiomassMet: 1})
                self.stoichiometricModel.add_reactions([rxn])
            self.stoichiometricModel.reactions.Total_biomass_tx.add_metabolites({BiomassMet: -1})

            self.addShell([Organ("shell" + str(1), position=1, tag = "SH"+str(1))])
            if not self.parameters["singlePhloem"]:
                for met in self.stoichiometricModel.metabolites:
                    if met.compartment == "ph_SH"+str(1):
                        met2 = self.stoichiometricModel.metabolites.get_by_id(met.id.replace("SH"+str(1),self.connections["SH"+str(1)]))
                        rxn = Reaction(met.id+"_"+met2.id+"_exchange")
                        rxn.add_metabolites({met:-1,met2:1})
                        rxn.upper_bound = 1000
                        rxn.lower_bound = -1000
                        self.stoichiometricModel.add_reactions([rxn])
                    if met.compartment == "xy_SH"+str(1):
                        met2 = self.stoichiometricModel.metabolites.get_by_id(met.id.replace("SH"+str(1),self.connections["SH"+str(1)]))
                        rxn = Reaction(met.id+"_"+met2.id+"_exchange")
                        rxn.add_metabolites({met:-1,met2:1})
                        rxn.upper_bound = 1000
                        rxn.lower_bound = -1000
                        self.stoichiometricModel.add_reactions([rxn])

        if self.parameters["seedPosition"][t-1]:
            ##DONT ADD MORE THAN 1 SEED###
            if len(self.seeds):
                return
            if not "seedModel" in dir(self):
                self.generateStoichiometricModels(sbmlFile="./PlantCoreMetabolism_v2_0_0.xml",sink_tissues = ["SE",])

            modelNew = self.seedModel.copy()
            organTag = "SE"
            for met in modelNew.metabolites:
                if not(met.id.endswith("_xy")) and not(met.id.endswith("_ph")):
                    met.id = met.id + str(1)
                    if met.compartment != "f":
                        met.compartment = met.compartment + "1"
            BiomassMet = Metabolite("Biomass_" + organTag + str(1), name = "Biomass_" + organTag + str(1), compartment="b_" + organTag + str(1))
            for rxn in modelNew.reactions:
                rxn.id = rxn.id + str(1)
                if rxn.id == "Biomass_tx_" + organTag + str(1):
                    rxn.add_metabolites({BiomassMet: 1})
                self.stoichiometricModel.add_reactions([rxn])
            self.stoichiometricModel.reactions.Total_biomass_tx.add_metabolites({BiomassMet: -1})
            self.seedBiomassWeight = sum([abs(met.formula_weight*self.seedModel.reactions.get_by_id("Biomass_tx_SE").metabolites[met]) for met in self.seedModel.reactions.get_by_id("Biomass_tx_SE").metabolites])/1000

            self.addSeed([Organ("seed" + str(1), position=1, connection=self.connections["S"+str(1)], tag = "SE"+str(1))])
            if not self.parameters["singlePhloem"]:
                for met in self.stoichiometricModel.metabolites:
                    if met.compartment == "ph_SE1":
                        met2 = self.stoichiometricModel.metabolites.get_by_id(met.id.replace("SE"+str(1),self.connections["S"+str(1)]))
                        rxn = Reaction(met.id+"_"+met2.id+"_exchange")
                        rxn.add_metabolites({met:-1,met2:1})
                        rxn.upper_bound = 1000
                        rxn.lower_bound = -1000
                        self.stoichiometricModel.add_reactions([rxn])
                    if met.compartment == "xy_SE1":
                        met2 = self.stoichiometricModel.metabolites.get_by_id(met.id.replace("SE"+str(1),self.connections["S"+str(1)]))
                        rxn = Reaction(met.id+"_"+met2.id+"_exchange")
                        rxn.add_metabolites({met:-1,met2:1})
                        rxn.upper_bound = 1000
                        rxn.lower_bound = -1000
                        self.stoichiometricModel.add_reactions([rxn])
                    if met.compartment == "ph_SH1":
                        met2 = self.stoichiometricModel.metabolites.get_by_id(met.id.replace("SH"+str(1),self.connections["S"+str(1)]))
                        rxn = Reaction(met.id+"_"+met2.id+"_exchange")
                        rxn.add_metabolites({met:-1,met2:1})
                        rxn.upper_bound = 1000
                        rxn.lower_bound = -1000
                        self.stoichiometricModel.add_reactions([rxn])
                    if met.compartment == "xy_SH1":
                        met2 = self.stoichiometricModel.metabolites.get_by_id(met.id.replace("SH"+str(1),self.connections["S"+str(1)]))
                        rxn = Reaction(met.id+"_"+met2.id+"_exchange")
                        rxn.add_metabolites({met:-1,met2:1})
                        rxn.upper_bound = 1000
                        rxn.lower_bound = -1000
                        self.stoichiometricModel.add_reactions([rxn])
    
    def updateConstraints(self,t,verbose=False):

        """
        Function to calculate maximum nitrate uptake rate.

        Parameters
        ----------
        - Nitrate concentration
        - Nitrate uptake Vmax
        - Nitrate uptake Km
        - Root biomass
        - Maintenance respiration requirements


        Returns
        -------
        - Maximum nitrate uptake by root (organ biomass adjusted)
        - Updated maintenance respiration constraints


        Author: Tarit Konuntakiet, Sanu Shameer
        Email: tarit.konuntakiet@seh.ox.ac.uk, sanushameer@gmail.com

        """

        #Root nitrate uptake
        if self.parameters["nitrateConc"] >= 1:
            num = float(self.parameters["NO3UptakeVm"]*self.parameters["nitrateConc"])
            denom = float(self.parameters["NO3UptakeKm"]+self.parameters["nitrateConc"])
        else:
            num = float(self.parameters["NO3UptakeVmLow"]*self.parameters["nitrateConc"])
            denom = float(self.parameters["NO3UptakeKmLow"]+self.parameters["nitrateConc"])

        MichMen = float(num/denom)
        MichMen = float(MichMen*self.roots.root1.mass)
        if "NO3UptakeRate" in self.parameters.keys():
            self.stoichiometricModel.reactions.Nitrate_ec_R1.upper_bound = round(self.parameters["NO3UptakeRate"]*self.roots.root1.mass,5)
        else:
            self.stoichiometricModel.reactions.Nitrate_ec_R1.upper_bound = MichMen
        if verbose:
            print("Nitrate uptake limit =")
            print("self.stoichiometricModel.reactions.Nitrate_ec_R1.upper_bound")


        if self.mode=="daily":
            #Leaf photon uptake rate and RuBisCO carboxylation
            for leaf in self.leaves:
                if leaf.stage == 0: #Leaf is alive
                    self.stoichiometricModel.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p_LD" + str(leaf.position)).upper_bound = self.parameters["VcMax"]*leaf.mass
                    self.stoichiometricModel.reactions.get_by_id("Photon_tx_LD" + str(leaf.position)).upper_bound = self.parameters["Pmax"]*leaf.mass#*self.canopyEffect(0.792, 0.1,t)
                else: #Leaf has senesced
                    self.stoichiometricModel.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p_LD" + str(leaf.position)).upper_bound = 0
                    self.stoichiometricModel.reactions.get_by_id("Photon_tx_LD" + str(leaf.position)).upper_bound = 0

            #Update maintenace respiration requirements
            #Leaves
            for leaf in self.leaves:
                if leaf.stage == 0:
                    self.stoichiometricModel.reactions.get_by_id("ATPase_tx_LD" + str(leaf.position)).upper_bound = self.stoichiometricModel.reactions.get_by_id("ATPase_tx_LD" + str(leaf.position)).lower_bound = (self.parameters["leafMaint"]*leaf.mass/24)*self.parameters["photoperiod"]
                    self.stoichiometricModel.reactions.get_by_id("ATPase_tx_LN" + str(leaf.position)).upper_bound = self.stoichiometricModel.reactions.get_by_id("ATPase_tx_LN" + str(leaf.position)).lower_bound = (self.parameters["leafMaint"]*leaf.mass/24)*(24-self.parameters["photoperiod"])
                else: #Leaf has senseced
                    self.stoichiometricModel.reactions.get_by_id("ATPase_tx_LD" + str(leaf.position)).upper_bound = self.stoichiometricModel.reactions.get_by_id("ATPase_tx_LD" + str(leaf.position)).lower_bound = 0
                    self.stoichiometricModel.reactions.get_by_id("ATPase_tx_LN" + str(leaf.position)).upper_bound = self.stoichiometricModel.reactions.get_by_id("ATPase_tx_LN" + str(leaf.position)).lower_bound = 0
        elif self.mode=="hourly":
            # Ensure total CO2 still includes all CO2 exchange fluxes
            total_co2 = self.stoichiometricModel.problem.Constraint(0,lb=-self.maxFlux,ub=0,name='Total_CO2')
            self.stoichiometricModel.add_cons_vars(total_co2)
            cons_coeffs={}
            for rxn in [x for x in self.stoichiometricModel.reactions if 'CO2_tx' in x.id]:
                cons_coeffs[rxn.forward_variable]=1
                cons_coeffs[rxn.reverse_variable]=-1
            total_co2.set_linear_coefficients(coefficients=cons_coeffs)
            self.stoichiometricModel.solver.update()

            # if round(self.parameters["AssimilationRate"],3)>0:
            #     daytime = True
            # else:
            #     daytime = False
            daytime = self.parameters["AssimilationRate"]>0.052
            print(self.MetaboliteStores)
            print('ASS rate',self.parameters["AssimilationRate"],daytime,self.parameters["AssimilationRate"]<0.052)
            for rxn in self.stoichiometricModel.reactions.query("linker"):
                rxn.lower_bound = -1000
                rxn.upper_bound = 0
                rxn.lower_bound = 0
            for leaf in self.leaves:
                self.stoichiometricModel.reactions.get_by_id("Photon_tx_" + leaf.tag).upper_bound = round(float(self.parameters["Pmax"]*leaf.mass),5)#*self.canopyEffect(0.792, 0.1,t)
                self.stoichiometricModel.constraints.get('Total_CO2').ub = round(float(self.parameters["AssimilationRate"]*leaf.mass),5)#*self.canopyEffect(0.792, 0.1,t)
                self.stoichiometricModel.reactions.get_by_id("ATPase_tx_" + leaf.tag).lower_bound = -1000
                self.stoichiometricModel.reactions.get_by_id("ATPase_tx_" + leaf.tag).upper_bound = float(round((self.parameters["leafMaint"]*leaf.mass),5))
                self.stoichiometricModel.reactions.get_by_id("ATPase_tx_" + leaf.tag).lower_bound = 0#float(round((self.parameters["leafMaint"]*leaf.mass),5))
                self.stoichiometricModel.reactions.get_by_id("NADPHoxp_tx_" + str(leaf.tag)).lower_bound = -1000
                self.stoichiometricModel.reactions.get_by_id("NADPHoxp_tx_" + leaf.tag).upper_bound = round((self.parameters["leafMaint"]*leaf.mass),5)/9
                self.stoichiometricModel.reactions.get_by_id("NADPHoxp_tx_" + str(leaf.tag)).lower_bound = 0#round((self.parameters["leafMaint"]*leaf.mass),5)/9
                self.stoichiometricModel.reactions.get_by_id("NADPHoxc_tx_" + str(leaf.tag)).lower_bound = -1000
                self.stoichiometricModel.reactions.get_by_id("NADPHoxc_tx_" + leaf.tag).upper_bound = round((self.parameters["leafMaint"]*leaf.mass),5)/9
                self.stoichiometricModel.reactions.get_by_id("NADPHoxc_tx_" + str(leaf.tag)).lower_bound = 0#round((self.parameters["leafMaint"]*leaf.mass),5)/9
                self.stoichiometricModel.reactions.get_by_id("NADPHoxm_tx_" + str(leaf.tag)).lower_bound = -1000
                self.stoichiometricModel.reactions.get_by_id("NADPHoxm_tx_" + leaf.tag).upper_bound = round((self.parameters["leafMaint"]*leaf.mass),5)/9
                self.stoichiometricModel.reactions.get_by_id("NADPHoxm_tx_" + str(leaf.tag)).lower_bound = 0#round((self.parameters["leafMaint"]*leaf.mass),5)/9
                # Test new daytime def
                # temp_model = self.stoichiometricModel.copy()
                # if 'Obj_var_biocro' in temp_model.variables:
                #     temp_model.objective = temp_model.problem.Objective(temp_model.variables.Obj_var_biocro,direction='min')
                # else:
                #     print('secondary daytime')
                #     temp_model.objective = {temp_model.reactions.Total_biomass_tx:1}
                #     temp_model.reactions.Starch_storage_p_ST1.lower_bound=0
                #     temp_model.reactions.Starch_storage_p_R1.lower_bound=0
                # daytime = temp_model.slim_optimize()>0 and self.parameters["Pmax"]>0
                # print('D/N:',daytime,temp_model.slim_optimize(), self.parameters["Pmax"])
                if daytime:
                    self.stoichiometricModel.reactions.get_by_id("STARCH_p_linker" + str(leaf.position)).lower_bound = -self.maxFlux
                    self.stoichiometricModel.reactions.get_by_id("STARCH_p_linker" + str(leaf.position)).upper_bound = self.stoichiometricModel.reactions.get_by_id("STARCH_p_linker" + str(leaf.position)).lower_bound = round(float(self.parameters["StarchAccumulationRate"]*leaf.mass),5)
                    #print("CO2 assimilation rate")
                    #print(self.stoichiometricModel.reactions.get_by_id("CO2_tx_L" + str(leaf.position)).upper_bound)
                    #print((24-timeOfDay))
                else:
                    ## Changed!
                    # self.stoichiometricModel.reactions.get_by_id("Photon_tx_" + leaf.tag).upper_bound = round(float(self.parameters["Pmax"]*leaf.mass),5)#0#self.parameters["Pmax"]*leaf.mass*self.canopyEffect(0.792, 0.1,t)
                    # self.stoichiometricModel.reactions.get_by_id("CO2_tx_" + leaf.tag).upper_bound = max(0,round(float(self.parameters["AssimilationRate"]*leaf.mass),5))#0
                    self.stoichiometricModel.reactions.get_by_id("STARCH_p_linker" + str(leaf.position)).lower_bound =-1000
                    self.stoichiometricModel.reactions.get_by_id("STARCH_p_linker" + str(leaf.position)).upper_bound = round(float(-1*(self.MetaboliteStores["STARCH_p_L"+str(leaf.position)])/(self.parameters["nightLength"])),5)
                    self.stoichiometricModel.reactions.get_by_id("STARCH_p_linker" + str(leaf.position)).lower_bound =round(float(-1*(self.MetaboliteStores["STARCH_p_L"+str(leaf.position)])/(self.parameters["nightLength"])),5)
                    #print("Starch degradation rate")
                    #print(-1*(self.MetaboliteStores["STARCH_p_L"+str(leaf.position)])/(24-timeOfDay))
                    #print((24-timeOfDay))
                if verbose:
                    print("---"+str(leaf.id)+"---")
                    print("Photon available ="+str(self.stoichiometricModel.reactions.get_by_id("Photon_tx_L" + str(leaf.position)).upper_bound))
                    print("Assimilation rate constraint ="+str(self.stoichiometricModel.reactions.get_by_id("CO2_tx_L" + str(leaf.position)).upper_bound))
                    print("Starch accumulation rate ="+str(self.stoichiometricModel.reactions.get_by_id("STARCH_p_linker" + str(leaf.position)).upper_bound))
                if self.phloemComposition!="free":
                    if "Phloem_output_rev_tx_L"+str(leaf.position) in [rxn.id for rxn in self.stoichiometricModel.reactions]:
                        self.stoichiometricModel.reactions.get_by_id("Phloem_output_rev_tx_L" + str(leaf.position)).lower_bound = 0
                        self.stoichiometricModel.reactions.get_by_id("Phloem_output_rev_tx_L" + str(leaf.position)).upper_bound = 0
                #self.stoichiometricModel.reactions.get_by_id("Phloem_output_tx_L" + str(leaf.position)).objective_coefficient = 1
        #Stem
        for stem in self.stems:
            self.stoichiometricModel.reactions.get_by_id("ATPase_tx_"+stem.tag).lower_bound=-self.maxFlux
            self.stoichiometricModel.reactions.get_by_id("ATPase_tx_"+stem.tag).upper_bound = self.stoichiometricModel.reactions.get_by_id("ATPase_tx_"+stem.tag).lower_bound = 0#float(self.parameters["stemMaint"]*stem.mass)
            if "STARCH_p_"+stem.tag in self.MetaboliteStores.keys():
                self.stoichiometricModel.reactions.get_by_id("Starch_storage_p_" + stem.tag).lower_bound = min(-self.maxDebt,round(float(-1*(self.MetaboliteStores["STARCH_p_"+stem.tag])),5))
                self.stoichiometricModel.reactions.get_by_id("Starch_storage_p_" + stem.tag).upper_bound = round(float(-1*(self.MetaboliteStores["STARCH_p_"+stem.tag])),5)
            else:
                self.stoichiometricModel.reactions.get_by_id("Starch_storage_p_" + stem.tag).lower_bound = -self.maxDebt

        #Root
        for root in self.roots:
            self.stoichiometricModel.reactions.get_by_id("ATPase_tx_"+root.tag).lower_bound = -self.maxFlux
            self.stoichiometricModel.reactions.get_by_id("ATPase_tx_"+root.tag).upper_bound = self.stoichiometricModel.reactions.get_by_id("ATPase_tx_"+root.tag).lower_bound = 0#float(self.parameters["rootMaint"]*root.mass)
            if "STARCH_p_"+root.tag in self.MetaboliteStores.keys():
                self.stoichiometricModel.reactions.get_by_id("Starch_storage_p_" + root.tag).lower_bound = min(-self.maxDebt,round(float(-1*(self.MetaboliteStores["STARCH_p_"+root.tag])),5))
                self.stoichiometricModel.reactions.get_by_id("Starch_storage_p_" + root.tag).upper_bound = round(float(-1*(self.MetaboliteStores["STARCH_p_"+root.tag])),5)

            else:
                self.stoichiometricModel.reactions.get_by_id("Starch_storage_p_" + root.tag).lower_bound = -self.maxDebt

        # #Fruits
        # for fruit in self.fruits:
        #     self.stoichiometricModel.reactions.get_by_id("ATPase_tx_FR" + str(fruit.position)).upper_bound = self.stoichiometricModel.reactions.get_by_id("ATPase_tx_FR" + str(fruit.position)).lower_bound = self.parameters["fruitMaint"]*fruit.mass

        #Add previous unmetNGAM
        ## NOTE: Think about this...
        for rxn in self.UnmetNGAM.keys():
            self.stoichiometricModel.reactions.get_by_id(rxn).upper_bound = self.stoichiometricModel.reactions.get_by_id(rxn).lower_bound = float(self.stoichiometricModel.reactions.get_by_id(rxn).lower_bound + self.UnmetNGAM[rxn])
        self.UnmetNGAM = dict()
        print('Pre-pickle inspect B')
        # from cobra.io import write_sbml_model
        # write_sbml_model(self.stoichiometricModel,"Model4Testing_2.xml")
        with open('Model4Testing_2'+self.saveTag+'.pkl','wb') as file:
            pkl.dump(self.stoichiometricModel,file)
        print(self.stoichiometricModel.reactions.query("GDWsink"))
        if "Leaf_senesced" in self.parameters.keys():
            self.stoichiometricModel.reactions.get_by_id("GDWsink_L1").lower_bound  = -1*self.parameters["Leaf_senesced"]
            self.stoichiometricModel.reactions.get_by_id("Biomass_tx_L1").lower_bound = -1000
        if "Stem_senesced" in self.parameters.keys():
            self.stoichiometricModel.reactions.get_by_id("GDWsink_ST1").lower_bound  = -1*self.parameters["Stem_senesced"]
            self.stoichiometricModel.reactions.get_by_id("Biomass_tx_ST1").lower_bound = -1000
        if "Root_senesced" in self.parameters.keys():
            self.stoichiometricModel.reactions.get_by_id("GDWsink_R1").lower_bound  = -1*self.parameters["Root_senesced"]
            self.stoichiometricModel.reactions.get_by_id("Biomass_tx_R1").lower_bound = -1000

        #Add TEMPORARY fix to if phloem composition was set incorrectly
        if "Phloem_output_tx_SE1" in self.stoichiometricModel.reactions:
            print('Phloem temp fix...',self.stoichiometricModel.reactions.Phloem_output_tx_SE1.reaction)
            for met in self.stoichiometricModel.reactions.Phloem_output_tx_SE1.reactants:
                if not(met.id.startswith("PROTON")):
                    print(met.id)
                    met_name = met.id[:-6]
                    print(met_name)
                    rxn = Reaction(met_name + "_ph_SE1", name = met_name + "_ph_SE1", subsystem = "Phloem", lower_bound = 0, upper_bound = 1000)
                    new_met = self.stoichiometricModel.metabolites.get_by_id(met.id).copy()
                    new_met.name = met_name + "_ph"
                    new_met.id = met_name + "_ph"
                    print(new_met.id)
                    new_met.compartment = "ph"
                    self.stoichiometricModel.add_reactions([rxn])
                    self.stoichiometricModel.reactions.get_by_id(met_name + "_ph_SE1").add_metabolites({self.stoichiometricModel.metabolites.PROTON_c_SE1: -1,self.stoichiometricModel.metabolites.PROTON_e_SE1: 1,met: 1,new_met: -1})
                    print(self.stoichiometricModel.reactions.get_by_id(met_name + "_ph_SE1").reaction)
            self.stoichiometricModel.reactions.Phloem_output_tx_SE1.remove_from_model()

    def gaussianFunction(self, age, position, organ):
        age = age/(self.parameters[organ + "GrowthDuration"]*self.parameters[organ+"ExpScale"][position-1])
        varStrength = (self.parameters[organ + "VarA"]*math.e)**-(((age-self.parameters[organ + "VarB"])**2)/((2*self.parameters[organ + "VarC"])**(2)))
        return varStrength
    
    def updateTotalBiomass(self,type="default"):
        if type=="default":
            self.stoichiometricModel.reactions.Total_biomass_tx.remove_from_model()
            self.stoichiometricModel.add_reactions([Reaction("Total_biomass_tx", name = "Total_biomass_tx", lower_bound = 0, upper_bound = 1000000)])
            # self.stoichiometricModel.reactions.Total_biomass_tx.objective_coefficient = 1
            if self.mode=="daily":
                for leaf in self.leaves:
                    self.stoichiometricModel.reactions.Total_biomass_tx.add_metabolites({self.stoichiometricModel.metabolites.get_by_id("Biomass_LD" + str(leaf.position)):
                                                                                  -(self.gL*self.leafBiomass*0.75)*(self.gaussianFunction(age=leaf.age, position=leaf.position, organ = "leaf")/self.totSinkStrengthLeaf),
                                                                                 self.stoichiometricModel.metabolites.get_by_id("Biomass_LN" + str(leaf.position)):
                                                                                  -(self.gL*self.leafBiomass*0.25)*(self.gaussianFunction(age=leaf.age, position=leaf.position, organ = "leaf")/self.totSinkStrengthLeaf)})
            elif self.mode=="hourly":
                for leaf in self.leaves:
                    self.stoichiometricModel.reactions.Total_biomass_tx.add_metabolites({self.stoichiometricModel.metabolites.get_by_id("Biomass_L" + str(leaf.position)):
                                                                                          -(self.gL*self.leafBiomass)*(self.gaussianFunction(age=leaf.age, position=leaf.position, organ = "leaf")/self.totSinkStrengthLeaf)})
            self.stoichiometricModel.reactions.Total_biomass_tx.add_metabolites({self.stoichiometricModel.metabolites.Biomass_R1: -(self.gR*self.rootBiomass),
                                                                              self.stoichiometricModel.metabolites.Biomass_ST1: -(self.gS*self.stemBiomass)})
            self.gF = 6.59
            self.gF1 = 0.01
            # if self.reproductiveSwitch == 1:# and round(self.totSinkStrengthFruit,1) > 0:
            #     for fruit in self.fruits:
            #         if fruit.age > self.floweringTime and fruit.age <45:
            #             self.stoichiometricModel.reactions.Total_biomass_tx.add_metabolites(
            #                 {self.stoichiometricModel.metabolites.get_by_id("Biomass_FR" + str(fruit.position)): -(self.gF*fruit.mass)})#*(self.gaussianFunction(age=fruit.age, position=fruit.position, organ = "fruit")/self.totSinkStrengthFruit)})
        elif type=="coeff-from-file":###UPDATE-FIX
            for rxn in self.stoichiometricModel.reactions.query("GDWsink"):
                rxn.lower_bound = 0
                rxn.upper_bound = 0
            if len(self.stoichiometricModel.reactions.query("Total_biomass_tx"))!=0:
                self.stoichiometricModel.reactions.Total_biomass_tx.remove_from_model()
            self.stoichiometricModel.add_reactions([Reaction("Total_biomass_tx", name = "Total_biomass_tx", lower_bound = 0, upper_bound = 1000000)])
            # self.stoichiometricModel.reactions.Total_biomass_tx.objective_coefficient = 1
            for leaf in self.leaves:
                self.stoichiometricModel.reactions.Total_biomass_tx.add_metabolites({self.stoichiometricModel.metabolites.get_by_id("gDW_L" + str(leaf.position)):
                    -10*(self.parameters["Growth_coeff_Leaf"])*(self.gaussianFunction(age=leaf.age, position=leaf.position, organ = "leaf")/self.totSinkStrengthLeaf)})
            self.stoichiometricModel.reactions.Total_biomass_tx.add_metabolites({self.stoichiometricModel.metabolites.gDW_R1: -10*(self.parameters["Growth_coeff_Root"]),
                self.stoichiometricModel.metabolites.gDW_ST1: -10*(self.parameters["Growth_coeff_Stem"])})
            if self.parameters["Growth_coeff_Seed"] > 0.01:
                for seed in self.seeds:
                    MW = 0
                    Rxn = self.stoichiometricModel.reactions.get_by_id("Biomass_tx_"+seed.tag)
                    for met in Rxn.metabolites:
                        MW = MW+(met.formula_weight*Rxn.metabolites[met]*-1)
                    GDW = Metabolite("gDW_"+seed.tag,name="1 gWD",compartment="c_"+seed.tag)
                    Rxn.add_metabolites({GDW:MW/1000})
                    if len(self.stoichiometricModel.metabolites.query("Biomass_SE"))>0:
                        self.stoichiometricModel.metabolites.get_by_id("Biomass_"+seed.tag).remove_from_model()
                    self.stoichiometricModel.reactions.Total_biomass_tx.add_metabolites({self.stoichiometricModel.metabolites.get_by_id("gDW_SE"+str(seed.position)): -10*(self.parameters["Growth_coeff_Seed"])})
                for shell in self.shells:
                    MW = 0
                    Rxn = self.stoichiometricModel.reactions.get_by_id("Biomass_tx_"+shell.tag)
                    for met in Rxn.metabolites:
                        MW = MW+(met.formula_weight*Rxn.metabolites[met]*-1)
                    GDW = Metabolite("gDW_"+shell.tag,name="1 gWD",compartment="c_"+shell.tag)
                    Rxn.add_metabolites({GDW:MW/1000})
                    if len(self.stoichiometricModel.metabolites.query("Biomass_SE"))>0:
                        self.stoichiometricModel.metabolites.get_by_id("Biomass_"+shell.tag).remove_from_model()
                    self.stoichiometricModel.reactions.Total_biomass_tx.add_metabolites({self.stoichiometricModel.metabolites.get_by_id("gDW_SH"+str(shell.position)): -10*(self.parameters["Growth_coeff_Shell"])})
            self.stoichiometricModel.reactions.Total_biomass_tx.upper_bound = 1000
        elif type=="Calloc-from-file":
            for rxn in self.stoichiometricModel.reactions.query("Cmet_sink_"):
                rxn.lower_bound = 0
                rxn.upper_bound = 0
            # for rxn in self.stoichiometricModel.reactions.query("Biomass_tx"):
            #     rxn.objective_coefficient=1
            if len(self.stoichiometricModel.reactions.query("Callocation_ratio_tx"))>0:
                self.stoichiometricModel.reactions.Callocation_ratio_tx.remove_from_model()
            self.stoichiometricModel.add_reactions([Reaction("Callocation_ratio_tx", name = "Callocation_ratio", lower_bound = 0, upper_bound = 1000000)])
            for leaf in self.leaves:
                self.stoichiometricModel.reactions.get_by_id("Callocation_ratio_tx").add_metabolites({self.stoichiometricModel.metabolites.get_by_id("Cmet_L"+str(leaf.position)):
                    -1*(self.parameters["Growth_coeff_Leaf"])*(self.gaussianFunction(age=leaf.age, position=leaf.position, organ = "leaf")/self.totSinkStrengthLeaf)/12})
            self.stoichiometricModel.reactions.get_by_id("Callocation_ratio_tx").add_metabolites({self.stoichiometricModel.metabolites.Cmet_R1:-1*(self.parameters["Growth_coeff_Root"]),
                    self.stoichiometricModel.metabolites.Cmet_ST1: -1*(self.parameters["Growth_coeff_Stem"]),
                    self.stoichiometricModel.metabolites.Cmet_SH1: -1*(self.parameters["Growth_coeff_Shell"])})
            if self.parameters["Growth_coeff_Seed"] > 0.01:
                self.stoichiometricModel.reactions.get_by_id("Callocation_ratio_tx").add_metabolites({self.stoichiometricModel.metabolites.Cmet_SE1: -1*(self.parameters["Growth_coeff_Seed"])})

    def updateBiomassFromSolutionObject(self, sol):
        """Update mass of each organ after a period growth and update the total mass of each organ type

        Parameters
        ----------
        - Optimisation solution
        - Previous organ biomass

        Returns
        -------
        - Mass of organ after a period of growth

        """
        biomassCounter = 0
        for leaf in self.leaves:
            if self.mode=="daily":
                leaf.mass = leaf.mass + self.leafBiomassWeight*(abs(sol["Biomass_tx_LD" + str(leaf.position)]) + abs(sol["Biomass_tx_LN" + str(leaf.position)]))
            else:
                leaf.mass = leaf.mass + self.leafBiomassWeight*(abs(sol["Biomass_tx_L" + str(leaf.position)]))
            biomassCounter = biomassCounter + leaf.mass
        self.leafBiomass = biomassCounter

        for stem in self.stems:
            stem.mass = stem.mass + self.stemBiomassWeight*abs(sol["Biomass_tx_ST1"])*(self.gaussianFunction(age = stem.age, position=stem.position, organ = "stem")/self.totSinkStrengthStem)
        self.stemBiomass = self.stemBiomass + self.stemBiomassWeight*abs(sol["Biomass_tx_ST1"])

        biomassCounter = 0
        for shell in self.shells:
            shell.mass = shell.mass + self.shellBiomassWeight["stage" + str(shell.stage)]*(abs(sol["Biomass_tx_SH" + str(shell.position)]))
            biomassCounter = biomassCounter + shell.mass
        self.shellBiomass = biomassCounter

        self.roots.root1.mass = self.roots.root1.mass + (self.rootBiomassWeight*abs(sol["Biomass_tx_R1"]))
        self.rootBiomass = self.roots.root1.mass
    ## NOTE: Edited to add starch stores and removed leaf for-loop
    def updateMetaboliteStores(self,sol):
        metDict = dict()
        for linker in self.stoichiometricModel.reactions.query("linker"):
                metDict[linker.reactants[0]]=linker

        for met in metDict.keys():
            if met.id in self.MetaboliteStores.keys():
                self.MetaboliteStores[met.id]=self.MetaboliteStores[met.id] + sol[metDict[met].id]
            else:
                self.MetaboliteStores[met.id]=sol[metDict[met].id]
        for stem in self.stems:
            met ='STARCH_p_'+stem.tag
            if met in self.MetaboliteStores.keys():
                self.MetaboliteStores[met]=self.MetaboliteStores[met] + sol[self.stoichiometricModel.reactions.get_by_id("Starch_storage_p_"+stem.tag).id]
            else:
                self.MetaboliteStores[met]=sol[self.stoichiometricModel.reactions.get_by_id("Starch_storage_p_"+stem.tag).id]
        for root in self.roots:
            met ='STARCH_p_'+root.tag
            if met in self.MetaboliteStores.keys():
                self.MetaboliteStores[met]=self.MetaboliteStores[met] + sol[self.stoichiometricModel.reactions.get_by_id("Starch_storage_p_"+root.tag).id]
            else:
                self.MetaboliteStores[met]=sol[self.stoichiometricModel.reactions.get_by_id("Starch_storage_p_"+root.tag).id]

    def simulateHourlyGrowth(self, simulationTime = range(1,85*24), outputCycle = [], outputCycle1 = range(1,85*24),earlyTermination=0,alternate=False):
        leafBiomassOutput = [self.parameters["LbInit"]]
        stemBiomassOutput = [self.parameters["SbInit"]]
        rootBiomassOutput = [self.parameters["RbInit"]]
        # shellBiomassOutput =[self.parameters["RbInit"]]
        fluxOutputs = dict()
        biomassOutputs = dict()
        print("HERE:A")
        self.inspectModel()
        if self.mode == "hourly":
            ageConvertor=24
            self.floweringTime = 8*24
        elif self.mode == "daily":
            ageConvertor=1
            self.floweringTime = 8

        for t in simulationTime:
            self.age=self.age+(t/ageConvertor)
            print("-------")
            print("Age= "+str(self.age))

            #Update age of each organ in days
            for tissue in [self.leaves, self.roots, self.stems, self.fruits]:
                for organ in tissue:
                    organ.age = organ.age + 1./ageConvertor

            #Update organ developmental stage
            for leaf in self.leaves:
                if leaf.age >= self.parameters["leafLifeSpan"]*self.parameters["leafExpScale"][leaf.position-1]:
                    leaf.stage = 1 #Leaf senescence
            
            # #Change in SLA between vegetative and reproductive growth
            # if self.reproductiveSwitch==0:
            #     for fruit in self.fruits:
            #         if fruit.age > self.floweringTime:
            #             self.reproductiveSwitch = 1
            #             self.estimateVcFromPPFD()
            #Check if a complete growth cycle has elapsed
            ####FIXES####
            print("########################")
            print(t)
            print(self.cycleThreshold)
            # if len(self.parameters["leafPosition"])==1581:
            #     self.growthCycle = 936
            #     defList = [0]*4000
            #     self.parameters["leafPosition"] = self.parameters["leafPosition"]+defList
            #     self.parameters["stemPosition"] = self.parameters["stemPosition"]+defList
            #     # self.parameters["shellPosition"] = self.parameters["shellPosition"]+defList
            #     self.parameters["seedPosition"] = self.parameters["seedPosition"]+defList
            #     self.parameters["seedPosition"][937]=1
            print("########################")
            if t >= 1:#self.cycleThreshold:
                self.growthCycle = self.growthCycle + 1
                #if t != simulationTime[-1]:
                #    display.clear_output(wait=True)
                #self.cycleThreshold = self.cycleThreshold + self.parameters["PhyllochronLength"]

                #Check if new organs have appeared and add to stoichiometricModel
                self.createNewOrgansAlongAxis()

            # #Update fruit biomass composition
            # for fruit in self.fruits:
            #     self.stoichiometricModel.reactions.get_by_id("Biomass_tx_FR" + str(fruit.position)).remove_from_model()
            #     self.stoichiometricModel.add_reactions([Reaction("Biomass_tx_FR" + str(fruit.position), name = "Biomass_tx", upper_bound = self.maxFlux, lower_bound = 0)])
            #     self.stoichiometricModel.reactions.get_by_id("Biomass_tx_FR" + str(fruit.position)).add_metabolites({self.stoichiometricModel.metabolites.get_by_id(
            #         "Biomass_FR" + str(fruit.position)): 1})
            #     for met in self.shellBiomassEquations:
            #         if met != "DAE":
            #             self.stoichiometricModel.reactions.get_by_id("Biomass_tx_FR" + str(fruit.position)).add_metabolites({self.stoichiometricModel.metabolites.get_by_id(
            #                 met + "_FR" + str(fruit.position)): -self.shellBiomassEquations[met][1]})

            print("HERE:B")
            self.inspectModel()
            #Update constraints
            print('B: post inspect')
            self.updateConstraints(t)
            print('B: post update')
            #Biomass partitioning
            #print(leaf.age)
            #print(leaf.position)
            self.totSinkStrengthLeaf = sum([self.gaussianFunction(age=leaf.age, position=leaf.position, organ = "leaf") for leaf in self.leaves])
            #print(self.totSinkStrengthLeaf)
            self.totSinkStrengthStem = sum([self.gaussianFunction(age=stem.age, position=stem.position, organ = "stem") for stem in self.stems])
            # self.totSinkStrengthFruit = sum([self.gaussianFunction(age=fruit.age, position=fruit.position, organ = "fruit") for fruit in self.fruits])

            self.updateTotalBiomass(type="coeff-from-file")
            #self.updateTotalBiomass(type="Calloc-from-file")

            #print(self.stoichiometricModel.reactions.Callocation_ratio_tx.reaction)
            #exit()
            ##################################################
            print("HERE:C")
            self.inspectModel()
            with open('Model4Testing'+self.saveTag+'.pkl','wb') as file:
                pkl.dump(self.stoichiometricModel,file)
            #Optimise model with new constraints
            if t in outputCycle:
                signal.signal(signal.SIGALRM,handler)
                signal.alarm(10)
                try:
                    if 'Obj_var_biocro' in self.stoichiometricModel.variables:
                        self.stoichiometricModel.remove_cons_vars(['Obj_var_biocro','Dcons'])
                    # max_biomass = workOutMax(self.stoichiometricModel.copy(),'Total_biomass_tx',5)
                    # print(max_biomass)
                    Dcoeffs = {}
                    sum_var = self.stoichiometricModel.problem.Variable('Obj_var_biocro')
                    self.stoichiometricModel.add_cons_vars(sum_var)
                    for ii in ['Total_biomass_tx']:
                        Dcoeffs[self.stoichiometricModel.reactions.get_by_id(ii).forward_variable]=1
                    for ii in [rxn.id for rxn in self.stoichiometricModel.reactions if 'Starch_storage'in rxn.id]:
                        if sum([self.MetaboliteStores[key] for key in list(self.MetaboliteStores.keys()) if 'STARCH' in key])<0:
                            Dcoeffs[self.stoichiometricModel.reactions.get_by_id(ii).forward_variable]=1
                        else:
                            Dcoeffs[self.stoichiometricModel.reactions.get_by_id(ii).forward_variable]=0.001
                        Dcoeffs[self.stoichiometricModel.reactions.get_by_id(ii).reverse_variable]=-1000
                    Dcoeffs[sum_var]=1
                    Dcons = self.stoichiometricModel.problem.Constraint(0,
                        lb=0,
                        ub=0,
                        name='Dcons')
                    self.stoichiometricModel.add_cons_vars([Dcons])
                    self.stoichiometricModel.solver.update()
                    Dcons.set_linear_coefficients(coefficients=Dcoeffs)
                    max_bio_min_starch=self.stoichiometricModel.problem.Objective(sum_var,
                            direction='min')
                    self.stoichiometricModel.objective = max_bio_min_starch
                    # self.stoichiometricModel.reactions.Total_biomass_tx.lower_bound = max_biomass
                    print('pre-pfba')
                    sol=self.stoichiometricModel.optimize()
                    print('post-pfba')
                except:
                    print('alt solve')
                    self.stoichiometricModel.solver='gurobi'
                    sol=self.stoichiometricModel.slim_optimize()
                    sol=self.stoichiometricModel.optimize()
                    self.stoichiometricModel.solver='glpk'
                    # # print("pFBA returns infeasible solution, so running FBA")
                    # # sol = self.stoichiometricModel.optimize()
                    # if self.parameters["remobilizationPermitted"]:
                    #     for rxn in self.stoichiometricModel.reactions.query("Sucrose_tx"):
                    #         rxn.upper_bound = 1000
                    #         rxn.objective_coefficient = -999
                    #     self.stoichiometricModel.reactions.Total_biomass_tx.upper_bound = 0
                    #     self.stoichiometricModel.reactions.Total_biomass_tx.lower_bound = 0
                    #     sol = cobra.flux_analysis.parsimonious.pfba(self.stoichiometricModel)
                    # else:
                    #     self.stoichiometricModel.reactions.Total_biomass_tx.upper_bound = 0
                    #     self.stoichiometricModel.reactions.Total_biomass_tx.lower_bound = 0
                    #     met2remove = list()
                    #     rxn2remove = list()
                    #     rxn = Reaction("NGAMobjective")
                    #     for organ in [self.leaves+self.stems+self.roots+self.fruits+self.seeds][0]:
                    #         for rtype in ["ATPase","NADPHoxc","NADPHoxp","NADPHoxm"]:
                    #             r = self.stoichiometricModel.reactions.get_by_id(rtype+"_tx_"+organ.tag)
                    #             met = Metabolite("NGAM_"+r.id,compartment="temp")
                    #             met2remove.append(met.id)
                    #             r.add_metabolites({met:1})
                    #             rxn.add_metabolites({met:-1*r.lower_bound})
                    #             r.lower_bound = 0
                    #     self.stoichiometricModel.add_reactions([rxn])
                    #     self.inspectModel()
                    #     # write_sbml_model(self.stoichiometricModel,"Model4Testing.xml")
                    #     with open('Model4Testing.pkl','wb') as file:
                    #         pkl.dump(self.stoichiometricModel,file)
                    #     print(rxn.reaction)
                    #     rxn2remove.append(rxn.id)
                    #     rxn.objective_coefficient = 1
                    #     sol = cobra.flux_analysis.parsimonious.pfba(self.stoichiometricModel)
                    #     print("Fraction of NGAM met =",sol[rxn.id])
                    #     for organ in [self.leaves+self.stems+self.roots+self.fruits+self.seeds][0]:
                    #         for rtype in ["ATPase",]:
                    #             r = self.stoichiometricModel.reactions.get_by_id(rtype+"_tx_"+organ.tag)
                    #             self.UnmetNGAM[r.id] = r.upper_bound - sol[r.id]
                    #     print(self.UnmetNGAM)
                    #     self.stoichiometricModel.reactions.Total_biomass_tx.upper_bound = 1000000
                    #     for met in met2remove:
                    #         self.stoichiometricModel.metabolites.get_by_id(met).remove_from_model()
                    #     for rxn in rxn2remove:
                    #         self.stoichiometricModel.reactions.get_by_id(rxn).remove_from_model()

                    if sol.status != 'optimal':
                        print("solution not optimal")
                        exit()

                fluxOutputs["solution"+str(t)] = sol
                self.updateBiomassFromSolutionObject(sol)
            else:
                print(t)
                if t == earlyTermination:
                    return(self.stoichiometricModel)
                print("HERE:C0")
                if alternate:
                    tmp = self.stoichiometricModel.copy()
                    sol = tmp.optimize()
                    self.stoichiometricModel = tmp.copy()
                else:
                    sol = self.stoichiometricModel.optimize()
                    print("max objective="+str(sol.fluxes["Total_biomass_tx"]))
                print("HERE:C1")
                #print("----------")
                #print("t =" +str(t))
                #print("Biomass =" +str(self.stoichiometricModel.reactions.STARCH_p_linker1.flux))
                #if t>16:
                #    return self
                #    break
                #self.updateBiomassFromCobraModel()
                fluxOutputs["solution"+str(t)] = sol
                self.updateBiomassFromSolutionObject(sol)
                print("HERE:C2")
            signal.alarm(0)
            self.updateMetaboliteStores(sol)
            self.stoichiometricModel.remove_cons_vars('Total_CO2')
            self.recentSolution = sol
            print("HERE:D")
            if t in outputCycle1:
                biomassOutputs["phytomerLeaf" + str(t+12)] = [l.position for l in self.leaves]
                biomassOutputs["massLeaf" + str(t+12)] = [l.mass for l in self.leaves]
                biomassOutputs["phytomerStem" + str(t+12)] = [s.position for s in self.stems]
                biomassOutputs["massStem" + str(t+12)] = [s.mass for s in self.stems]
                # if self.reproductiveSwitch == 1:
                #     biomassOutputs["phytomerFruit" + str(t+12)] = [f.position for f in self.fruits]
                #     biomassOutputs["massFruit" + str(t+12)] = [f.mass for f in self.fruits]

            #Update organ biomass
            leafBiomassOutput.append(self.leafBiomass)
            stemBiomassOutput.append(self.stemBiomass)
            rootBiomassOutput.append(self.rootBiomass)
            # shellBiomassOutput.append(self.shellBiomass)
            print("HERE:E")
            # simulationTime1 = [12]+[i+12 for i in simulationTime]
            #ax1[0].plot(range(len(leafBiomassOutput)),leafBiomassOutput, "k")
            #ax1[1].plot(range(len(stemBiomassOutput)), stemBiomassOutput, "k")
            #ax2[0].plot(range(len(shellBiomassOutput)), shellBiomassOutput, "k")
            #ax2[1].plot(range(len(rootBiomassOutput)), rootBiomassOutput, "k")
            #display.clear_output(wait=True)
            #display.display(pl.gcf())
            #display.clear_output(wait=True)

            #print "leaf biomass=%.5f   root biomass=%.5f   stem biomass=%.5f   fruit biomass=%.5f" %(self.leafBiomass, self.rootBiomass, self.stemBiomass, self.shellBiomass)

            #if t == 103:
            #    tempModel = self.stoichiometricModel

        simulationTime = [12]+[i+12 for i in simulationTime]
        print("HERE:F")
        return fluxOutputs, simulationTime, outputCycle, leafBiomassOutput, stemBiomassOutput, rootBiomassOutput, [0]*len(rootBiomassOutput), biomassOutputs

    def workOutMax(model,rxn_to_opt='Total_biomass_tx',decimal_pts=3):
        feasible= 1
        count=-1
        rxn=model.reactions.get_by_id(rxn_to_opt)
        while feasible:
            count+=1
            rxn.lower_bound = pow(10,-5+count)
            temp=model.slim_optimize()
            if str(temp)=='nan':
                feasible = 0
                count-=1
        magnitude = -5+count
        counts = []
        for dec in range(min(decimal_pts,5+magnitude)):
            feasible=1
            starts_with = sum([pow(10,magnitude-it)*counts[it] for it in range(len(counts))])
            print(starts_with)
            counts+=[0]
            while feasible:
                counts[dec]+=1
                rxn.lower_bound = starts_with+pow(10,magnitude-dec)*counts[dec]
                print(counts,starts_with+pow(10,magnitude-dec)*counts[dec])
                temp=model.slim_optimize()
                if str(temp)=='nan':
                    print('next')
                    feasible = 0
                    counts[dec]-=1
        rxn.lower_bound=0    
        starts_with = sum([pow(10,magnitude-it)*counts[it] for it in range(len(counts))])
        return math.floor(starts_with)
    # def workOutMax(self,rxn_to_opt='Total_biomass_tx',decimal_pts=3):
    #     feasible= 1
    #     count=-1
    #     rxn=self.stoichiometricModel.reactions.get_by_id(rxn_to_opt)
    #     while feasible:
    #         count+=1
    #         rxn.lower_bound = pow(10,-5+count)
    #         temp=self.stoichiometricModel.slim_optimize()
    #         if str(temp)=='nan':
    #             feasible = 0
    #             count-=1
    #     magnitude = -5+count
    #     counts = []
    #     for dec in range(min(decimal_pts,5+magnitude)):
    #         feasible=1
    #         starts_with = sum([pow(10,magnitude-it)*counts[it] for it in range(len(counts))])
    #         print(starts_with)
    #         counts+=[0]
    #         while feasible:
    #             counts[dec]+=1
    #             rxn.lower_bound = starts_with+pow(10,magnitude-dec)*counts[dec]
    #             print(counts,starts_with+pow(10,magnitude-dec)*counts[dec])
    #             temp=self.stoichiometricModel.slim_optimize()
    #             if str(temp)=='nan':
    #                 print('next')
    #                 feasible = 0
    #                 counts[dec]-=1
    #     rxn.lower_bound=0    
    #     starts_with = sum([pow(10,magnitude-it)*counts[it] for it in range(len(counts))])
    #     return math.floor(starts_with)

def handler(signum,frame):
    print('pFBA left hanging')
    # raise Exception
def handler_so(signum,frame):
    print('slim opt left hanging')


class Organ:

    #"""Create new organ objects with age, mass, and position attributes"""

    def __init__(self,id, position, mass = 0, age = 1, stage = 0, connection=None, tag = "tag"):
        self.id = id
        self.age = age
        self.mass = mass
        self.position = position
        self.stage = stage
        self.connection = connection
        self.tag = tag

    def __repr__(self):
        return '%s \nPosition: %d \nAge: %d \nMass: %.4f \nStage: %d' %(self.id, self.position, self.age, self.mass, self.stage)