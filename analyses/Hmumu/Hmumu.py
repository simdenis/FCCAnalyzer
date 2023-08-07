import analysis, functions
import ROOT
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--nThreads", type=int, help="number of threads", default=None)
parser.add_argument("--maxFiles", type=int, help="Max number of files (per dataset)", default=100)
parser.add_argument("--flavor", type=str, choices=["ee", "mumu", "qq"], help="Flavor (ee, mumu, qq)", default="mumu")
parser.add_argument("--jetAlgo", type=str, choices=["kt", "valencia", "genkt"], default="genkt", help="Jet clustering algorithm")
args = parser.parse_args()

functions.set_threads(args)

# define histograms
bins_p_mu = (20000, 0, 200) # 10 MeV bins
bins_m_ll = (3000, 0, 300) # 1 GeV bins
bins_p_ll = (20000, 0, 200) # 10 MeV bins
bins_missingEnergy = (300, 0, 300) # 1 GeV
bins_theta = (500, -5, 5)
bins_phi = (500, -5, 5)
bins_aco = (400, -4, 4)
bins_dR = (1000, 0, 10)

bins_count = (100, 0, 100)
bins_pdgid = (60, -30, 30)
bins_charge = (10, -5, 5)
bins_cat = (10, 0, 10)
bins_resolution = (10000, 0.95, 1.05)

bins_resolution_1 = (20000, 0, 2)
bins_recoil = (20000, 0, 200) # 10 MeV bins 
bins_recoil_fine = (3000, 0, 300) # .1 GeV bins 
bins_cosThetaMiss = (10000, 0, 10000)

jet_energy = (1000, 0, 100) # 100 MeV bins
dijet_m = (2000, 0, 200) # 100 MeV bins
visMass = (2000, 0, 200) # 100 MeV bins
missEnergy  = (2000, 0, 200) # 100 MeV bins

bins_jet = (250, 0, 250)
dijet_m_final = (500, 50, 100) # 100 MeV bins

def build_graph(df, dataset):

    print("build graph", dataset.name)
    results = []

    df = df.Define("weight", "1.0")
    weightsum = df.Sum("weight")
    
    
    df = df.Alias("Particle0", "Particle#0.index")
    df = df.Alias("Particle1", "Particle#1.index")
    df = df.Alias("MCRecoAssociations0", "MCRecoAssociations#0.index")
    df = df.Alias("MCRecoAssociations1", "MCRecoAssociations#1.index")

    df = df.Alias("Muon0", "Muon#0.index")
    df = df.Alias("Electron0", "Electron#0.index")
    
    
    
    # all muons
    df = df.Define("muons_all", "FCCAnalyses::ReconstructedParticle::get(Muon0, ReconstructedParticles)")
    df = df.Define("elec_all", "FCCAnalyses::ReconstructedParticle::get(Electron0, ReconstructedParticles)")
    
    df = df.Define("elec_all_p", "FCCAnalyses::ReconstructedParticle::get_p(elec_all)")
    df = df.Define("muons_all_p", "FCCAnalyses::ReconstructedParticle::get_p(muons_all)")
    df = df.Define("muons_all_theta", "FCCAnalyses::ReconstructedParticle::get_theta(muons_all)")
    df = df.Define("muons_all_phi", "FCCAnalyses::ReconstructedParticle::get_phi(muons_all)")
    df = df.Define("muons_all_q", "FCCAnalyses::ReconstructedParticle::get_charge(muons_all)")
    df = df.Define("muons_all_no", "FCCAnalyses::ReconstructedParticle::get_n(muons_all)")

    
    # cuts on leptons
    #df = df.Define("selected_muons", "FCCAnalyses::excluded_Higgs_decays(muons, MCRecoAssociations0, MCRecoAssociations1, ReconstructedParticles, Particle, Particle0, Particle1)") # was 10
    df = df.Define("muons_sel_p", "FCCAnalyses::ReconstructedParticle::sel_p(25)(muons_all)")
    df = df.Alias("muons", "muons_sel_p") 
    df = df.Define("elec_sel_p", "FCCAnalyses::ReconstructedParticle::sel_p(25)(elec_all)")
    df = df.Alias("elecs", "elec_sel_p") 
    
    df = df.Define("elecs_p", "FCCAnalyses::ReconstructedParticle::get_p(elecs)")
    df = df.Define("muons_p", "FCCAnalyses::ReconstructedParticle::get_p(muons)")
    df = df.Define("muons_theta", "FCCAnalyses::ReconstructedParticle::get_theta(muons)")
    df = df.Define("muons_phi", "FCCAnalyses::ReconstructedParticle::get_phi(muons)")
    df = df.Define("muons_q", "FCCAnalyses::ReconstructedParticle::get_charge(muons)")
    df = df.Define("elec_no", "FCCAnalyses::ReconstructedParticle::get_n(elecs)")
    df = df.Define("muons_no", "FCCAnalyses::ReconstructedParticle::get_n(muons)")


    # lepton kinematic histograms
    results.append(df.Histo1D(("elec_all_p_cut0", "", *bins_p_mu), "elec_all_p"))
    results.append(df.Histo1D(("muons_all_p_cut0", "", *bins_p_mu), "muons_all_p"))
    results.append(df.Histo1D(("muons_all_theta_cut0", "", *bins_theta), "muons_all_theta"))
    results.append(df.Histo1D(("muons_all_phi_cut0", "", *bins_phi), "muons_all_phi"))
    results.append(df.Histo1D(("muons_all_q_cut0", "", *bins_charge), "muons_all_q"))
    #results.append(df.Histo1D(("elec_all_no_cut0", "", *bins_count), "elec_all_no"))
    results.append(df.Histo1D(("muons_all_no_cut0", "", *bins_count), "muons_all_no"))
    
    results.append(df.Histo1D(("elec_p_cut0", "", *bins_p_mu), "elecs_p"))
    results.append(df.Histo1D(("muons_p_cut0", "", *bins_p_mu), "muons_p"))
    results.append(df.Histo1D(("muons_theta_cut0", "", *bins_theta), "muons_theta"))
    results.append(df.Histo1D(("muons_phi_cut0", "", *bins_phi), "muons_phi"))
    results.append(df.Histo1D(("muons_q_cut0", "", *bins_charge), "muons_q"))
    #results.append(df.Histo1D(("elec_no_cut0", "", *bins_count), "elec_no"))
    results.append(df.Histo1D(("muons_no_cut0", "", *bins_count), "muons_no"))
    
  
    
    #########
    ### CUT 0: all events
    #########
    df = df.Define("cut0", "0")
    results.append(df.Histo1D(("cutFlow_cut0", "", *bins_count), "cut0"))
   

    #########
    ### CUT 1: at least a lepton with at least 1 isolated one
    #########
    df = df.Filter("muons_no >= 1").Define("cut1", "1")
    results.append(df.Histo1D(("cutFlow_cut1", "", *bins_count), "cut1"))
    
    
    #########
    ### CUT 2 :at least 2 leptons, and build the resonance
    #########
    df = df.Filter("muons_no >= 2").Define("cut2", "2")
    results.append(df.Histo1D(("cutFlow_cut2", "", *bins_count), "cut2"))
    

    # build the H resonance from the reconstructed particles
    df = df.Define("zbuilder_result", "FCCAnalyses::resonanceBuilder_mass_recoil(125, 91.2, 0.4, 240, false)(muons, MCRecoAssociations0, MCRecoAssociations1, ReconstructedParticles, Particle, Particle0, Particle1)")
    

    df = df.Define("hmm2", "ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>{zbuilder_result[0]}") # the Z
    df = df.Define("hmm_muons2", "ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>{zbuilder_result[1],zbuilder_result[2]}") # the muons
    df = df.Define("hmm_m2", "FCCAnalyses::ReconstructedParticle::get_mass(hmm2)[0]")
    
        

    df = df.Define("hmm_p2", "FCCAnalyses::ReconstructedParticle::get_p(hmm2)[0]")
    df = df.Define("hmm_recoil2", "FCCAnalyses::ReconstructedParticle::recoilBuilder(240)(hmm2)")
    df = df.Define("hmm_recoil_m2", "FCCAnalyses::ReconstructedParticle::get_mass(hmm_recoil2)[0]")

    df = df.Define("hee", "FCCAnalyses::resonanceBuilder(125)(elecs)")
    df = df.Define("hee_elecs", "ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>{hee[1],hee[2]}") # the leptons
    df = df.Define("hee_m", "FCCAnalyses::ReconstructedParticle::get_mass(hee)")
    
    df = df.Define("hmm", "FCCAnalyses::resonanceBuilder(125)(muons)")
    df = df.Define("hmm_muons", "ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>{hmm[1],hmm[2]}") # the leptons
    df = df.Define("hmm_m", "FCCAnalyses::ReconstructedParticle::get_mass(hmm)")
    
    df = df.Define("hmm_no", "FCCAnalyses::ReconstructedParticle::get_n(hmm)")
    df = df.Define("hmm_p", "FCCAnalyses::ReconstructedParticle::get_p(hmm)")
    df = df.Define("hmm_recoil", "FCCAnalyses::ReconstructedParticle::recoilBuilder(240)(hmm)")
    df = df.Define("hmm_recoil_m", "FCCAnalyses::ReconstructedParticle::get_mass(hmm_recoil)[0]")
    df = df.Define("missingEnergy_rp","FCCAnalyses::missingEnergy(240., ReconstructedParticles)")
    df = df.Define("missingEnergy", "missingEnergy_rp[0].energy")
    df = df.Define("hmm_muons_q", "FCCAnalyses::ReconstructedParticle::get_charge(hmm_recoil)")
    df = df.Define("zll_muons_dR", "FCCAnalyses::deltaR(hmm_muons2)")
    df = df.Define("cosTheta_miss", "FCCAnalyses::get_cosTheta_miss(missingEnergy_rp)")
    df = df.Define("acoplanarity", "FCCAnalyses::acoplanarity(muons)")
    df = df.Define("acolinearity", "FCCAnalyses::acolinearity(muons)")
    
    
    results.append(df.Histo1D(("hmm_m_cut2_2", "", *bins_m_ll), "hmm_m2"))
    results.append(df.Histo1D(("hmm_no_cut2", "", *bins_cat), "hmm_no"))
    results.append(df.Histo1D(("hmm_muons_q", "", *bins_charge), "hmm_muons_q"))
    
    df = df.Define("sum_q", "muons_q[0]+muons_q[1]")
    
    results.append(df.Histo1D(("muons_q_cut2", "", *bins_charge), "muons_q"))
    results.append(df.Histo1D(("sum_q", "", *bins_charge), "sum_q"))
    
    #########
    ### CUT 3 :at least 1 resonance (i.e. one opposite sign pair muon)
    #########
    df = df.Filter("hmm_no >= 1").Define("cut3", "3")
    results.append(df.Histo1D(("cutFlow_cut3", "", *bins_count), "cut3"))
    
        #####
    ### CUT 4: cut on momentum
    #####
    results.append(df.Histo1D(("hmm_p_cut4", "", *bins_m_ll), "hmm_p"))
    df = df.Filter("hmm_p[0] > 40 && hmm_p[0] < 72").Define("cut4", "4")
    results.append(df.Histo1D(("cutFlow_cut4", "", *bins_count), "cut4"))
    
   
    #########
    ### CUT 4: recoil cut
    #########  
    results.append(df.Histo1D(("hmm_recoil_m_cut3", "", *(bins_recoil_fine)), "hmm_recoil_m"))
    df = df.Filter("hmm_recoil_m > 86 && hmm_recoil_m < 114").Define("cut6","6")
    results.append(df.Histo1D(("cutFlow_cut6", "", *bins_count), "cut6"))

    


    ####
    ## CUT 5: cut on missing energy
    ####
    results.append(df.Histo1D(("missingEnergy", "", *bins_missingEnergy), "missingEnergy"))
    df = df.Filter(" missingEnergy < 115").Define("cut5", "5")
    results.append(df.Histo1D(("cutFlow_cut5", "", *bins_count), "cut5"))
    
        
    
    df = df.Define("ellipse", "(hmm_recoil_m-100)*(hmm_recoil_m-100)/196 + (hmm_m[0]-114.5)*(hmm_m[0]-114.5)/156 ")
    
    
    
    results.append(df.Histo1D(("ellipse", "", *bins_cosThetaMiss), "ellipse"))
    
    
    
    
    ########
    ### CUT 6: 
    #######
    
    #df = df.Filter("ellipse <=1").Define("cut6","6")
    #results.append(df.Histo1D(("cutFlow_cut6", "", *bins_count), "cut6"))
    
    
        #########
    ### CUT 7 :cut on Higgs mass
    #########
    
    #results.append(df.Histo1D(("hmm_m_cut7", "", *bins_m_ll), "hmm_m"))
    #results.append(df.Histo1D(("hmm_m_cut7_2", "", *bins_m_ll), "hmm_m2"))
    df = df.Filter("hmm_m[0] > 102 && hmm_m[0] < 127").Define("cut7", "7")
    #results.append(df.Histo1D(("cutFlow_cut7", "", *bins_count), "cut7"))
    #"""
    
    df = df.Define("p_leading", "(muons_p[0] > muons_p[1]) ? muons_p[0] : muons_p[1]")
    df = df.Define("p_subleading", "(muons_p[0] < muons_p[1]) ? muons_p[0] : muons_p[1]")
    

    
    ########################
    # Final histograms
    ########################
    results.append(df.Histo1D(("missingEnergy2", "", *bins_missingEnergy), "missingEnergy"))
    results.append(df.Histo1D(("hmm_m_final", "", *bins_m_ll), "hmm_m"))
    results.append(df.Histo1D(("hmm_p_cut4", "", *bins_p_ll), "hmm_p"))
    results.append(df.Histo1D(("muons_p_cut4", "", *bins_p_mu), "muons_p"))
    results.append(df.Histo1D(("hmm_recoil_m", "", *(bins_recoil_fine)), "hmm_recoil_m"))
    results.append(df.Histo1D(("hmm_m_final2", "", *bins_m_ll), "hmm_m2"))
    results.append(df.Histo1D(("cosThetaMiss", "", *bins_cosThetaMiss), "cosTheta_miss"))
    results.append(df.Histo1D(("acoplanarity", "", *bins_aco), "acoplanarity"))
    results.append(df.Histo1D(("acolinearity", "", *bins_aco), "acolinearity"))
    results.append(df.Histo1D(("muons_theta_final", "", *bins_theta), "muons_theta"))
    results.append(df.Histo1D(("muons_phi_final", "", *bins_phi), "muons_phi"))
    results.append(df.Histo1D(("p_leading", "", *bins_p_ll), "p_leading"))
    results.append(df.Histo1D(("p_subleading", "", *bins_p_ll), "p_subleading"))
    results.append(df.Histo1D(("zll_muons_dR_incorrectPairs", "", *bins_dR), "zll_muons_dR"))
    
    results.append(df.Histo1D(("final_muons", "", *bins_cat), "muons_all_no"))
    
    
    
    
    ##### CATEGORIZATION
    
    
    
    
    
    
    df_qq = df.Filter("! ( (muons_no ==2 and elec_no == 0 and missingEnergy>80) or (muons_no ==2 and elec_no == 2 and missingEnergy < 80) or  (muons_no ==4 and elec_no == 0 and missingEnergy < 80) )").Define("CATEGORIZATION", "7")
    #results.append(df.Histo1D(("cutFlow_cut7", "", *bins_count), "CATEGORIZATION"))
    

    df_qq = df_qq.Define("no_muons", "FCCAnalyses::ReconstructedParticle::remove(ReconstructedParticles, muons)")
    df_qq = df_qq.Define("RP_px", "FCCAnalyses::ReconstructedParticle::get_px(no_muons)")
    df_qq = df_qq.Define("RP_py", "FCCAnalyses::ReconstructedParticle::get_py(no_muons)")
    df_qq = df_qq.Define("RP_pz","FCCAnalyses::ReconstructedParticle::get_pz(no_muons)")
    df_qq = df_qq.Define("RP_e", "FCCAnalyses::ReconstructedParticle::get_e(no_muons)")
    df_qq = df_qq.Define("RP_m", "FCCAnalyses::ReconstructedParticle::get_mass(no_muons)")
    df_qq = df_qq.Define("RP_q", "FCCAnalyses::ReconstructedParticle::get_charge(no_muons)")
    
    df_qq = df_qq.Define("pseudo_jets", "FCCAnalyses::JetClusteringUtils::set_pseudoJets(RP_px, RP_py, RP_pz, RP_e)")
    
    df_qq = df_qq.Define("clustered_jets", "JetClustering::clustering_ee_kt(2, 2, 1, 0)(pseudo_jets)")
    df_qq = df_qq.Define("jets", "FCCAnalyses::JetClusteringUtils::get_pseudoJets(clustered_jets)")
    
    df_qq = df_qq.Define("jetconstituents", "FCCAnalyses::JetClusteringUtils::get_constituents(clustered_jets)")
    df_qq = df_qq.Define("jets_e", "FCCAnalyses::JetClusteringUtils::get_e(jets)")
    df_qq = df_qq.Define("jets_px", "FCCAnalyses::JetClusteringUtils::get_px(jets)")
    df_qq = df_qq.Define("jets_py", "FCCAnalyses::JetClusteringUtils::get_py(jets)")
    df_qq = df_qq.Define("jets_pz", "FCCAnalyses::JetClusteringUtils::get_pz(jets)")
    df_qq = df_qq.Define("jets_m", "FCCAnalyses::JetClusteringUtils::get_m(jets)")

    df_qq = df_qq.Define("Tlv_jet", "FCCAnalyses::jetsToTlv(jets_px, jets_py, jets_pz, jets_e)")  # + FCCAnalyses::inv_mass(jetconstituent[0]")
    df_qq = df_qq.Define("inv_mass", "(Tlv_jet[0]+Tlv_jet[1]).M()") 
    
    df_qq = df_qq.Define("jets_e1", "jets_e[0]")
    df_qq = df_qq.Define("jets_e2", "jets_e[1]")
    
    

    
    results.append(df_qq.Histo1D(( "jet_e1", "", *bins_jet), "jets_e1"))
    
    results.append(df_qq.Histo1D(( "jet_e2", "", *bins_jet), "jets_e2"))
    
    results.append(df_qq.Histo1D(("final_mass", "", *bins_m_ll), "inv_mass"))
    

    
    
    
    ####CUT 8
    df_qq = df_qq.Filter(" jets_e2 > 30 and jets_e2 < 50 ").Define("cut8", "8")
    results.append(df_qq.Histo1D(("cutFlow_cut8", "", *bins_count), "cut8"))
    

    ####CUT 9
    df_qq = df_qq.Filter(" jets_e1 < 83 ").Define("cut9", "9")
    results.append(df_qq.Histo1D(("cutFlow_cut9", "", *bins_count), "cut9"))
    
    
    
    results.append(df_qq.Histo1D(("missingEnergy_f", "", *bins_missingEnergy), "missingEnergy"))
    results.append(df_qq.Histo1D(("momentum_f", "", *bins_p_ll), "muons_p"))
    results.append(df_qq.Histo1D(("cosThetaMiss_f", "", *bins_cosThetaMiss), "cosTheta_miss"))
    results.append(df_qq.Histo1D(("acoplanarity_f", "", *bins_aco), "acoplanarity"))
    results.append(df_qq.Histo1D(("acolinearity_f", "", *bins_aco), "acolinearity"))
    results.append(df_qq.Histo1D(("muons_theta_f", "", *bins_theta), "muons_theta"))
    results.append(df_qq.Histo1D(("muons_phi_f", "", *bins_phi), "muons_phi"))
    results.append(df_qq.Histo1D(("p_leading_f", "", *bins_p_ll), "p_leading"))
    results.append(df_qq.Histo1D(("p_subleading_f", "", *bins_p_ll), "p_subleading"))
    results.append(df_qq.Histo1D(("Delta_R", "", *bins_dR), "zll_muons_dR"))
    results.append(df_qq.Histo1D(("mass_f", "", *bins_m_ll), "hmm_m"))
    
    
    
    
    df_nunu = df.Filter("! ( (muons_no ==2 and elec_no == 0 and missingEnergy<80) or (muons_no ==2 and elec_no == 2 and missingEnergy < 80) or  (muons_no ==4 and elec_no == 0 and missingEnergy < 80) )").Define("CATEGORIZATION2", "8")
    

    
    results.append(df_nunu.Histo1D(("mass_nunu", "", *bins_m_ll), "hmm_m"))
    
    
    
    df_mumu = df.Filter("! ( (muons_no ==2 and elec_no == 0 and missingEnergy<80) or (muons_no ==2 and elec_no == 2 and missingEnergy < 80) or  (muons_no == 2 and elec_no == 0 and missingEnergy < 80) )").Define("CATEGORIZATION2", "8")
    
    df_mumu = df_mumu.Define("muons0", "FCCAnalyses::ReconstructedParticle::remove(muons, hmm_muons2)")
    df_mumu = df_mumu.Define("nu_muons", "FCCAnalyses::ReconstructedParticle::get_n(muons0)")
    df_mumu = df_mumu.Filter("nu_muons == 2").Define("muonnumbers", "10")
    
    df_mumu = df_mumu.Define("muons_tlv", "FCCAnalyses::makeLorentzVectors(muons0)")
    df_mumu = df_mumu.Define("mass_mumu", "(muons_tlv[0] + muons_tlv[1]).M()")
    results.append(df_mumu.Histo1D(("nu_muons", "", *bins_cat), "nu_muons"))
    results.append(df_mumu.Histo1D(("mass_mumu", "", *bins_m_ll), "mass_mumu"))
    
    
    
    df_ee = df.Filter("! ( (muons_no ==2 and elec_no == 0 and missingEnergy<80) or (muons_no ==4 and elec_no == 0 and missingEnergy < 80) or  (muons_no == 2 and elec_no == 0 and missingEnergy < 80) )").Define("CATEGORIZATION2", "9")
    
    df_ee = df_ee.Define("elec_nu", "FCCAnalyses::ReconstructedParticle::get_n(elecs)")
    
  
    df_ee = df_ee.Filter("elec_nu == 2").Define("elecnumbers", "10")
    
    df_ee = df_ee.Define("elec_tlv", "FCCAnalyses::makeLorentzVectors(elecs)")
    df_ee = df_ee.Define("massee", "(elec_tlv[0] + elec_tlv[1]).M()")
    results.append(df_ee.Histo1D(("elec_nu", "", *bins_cat), "elec_nu"))
    results.append(df_ee.Histo1D(("mass_ee", "", *bins_m_ll), "massee"))
    
    
    return results, weightsum
    
    
    
   

if __name__ == "__main__":

    baseDir = "/data/submit/cms/store/fccee/winter2023_official/IDEA/"
    wzp6_ee_nunuH_Hmumu_ecm240 = {"name": "wzp6_ee_nunuH_Hmumu_ecm240", "datadir": f"{baseDir}wzp6_ee_nunuH_Hmumu_ecm240/",  "xsec": 1.005e-05}
    wzp6_ee_eeH_Hmumu_ecm240 = {"name": "wzp6_ee_eeH_Hmumu_ecm240", "datadir": f"{baseDir}wzp6_ee_eeH_Hmumu_ecm240/",  "xsec": 1.558e-06}
    wzp6_ee_ccH_Hmumu_ecm240 = {"name": "wzp6_ee_ccH_Hmumu_ecm240", "datadir": f"{baseDir}wzp6_ee_ccH_Hmumu_ecm240/",  "xsec": 5.079e-06}
    wzp6_ee_tautauH_Hmumu_ecm240 = {"name": "wzp6_ee_tautauH_Hmumu_ecm240", "datadir": f"{baseDir}wzp6_ee_tautauH_Hmumu_ecm240/",  "xsec": 1.469e-06}
    wzp6_ee_bbH_Hmumu_ecm240 = {"name": "wzp6_ee_bbH_Hmumu_ecm240", "datadir": f"{baseDir}wzp6_ee_bbH_Hmumu_ecm240/",  "xsec": 6.521e-06}
    wzp6_ee_qqH_Hmumu_ecm240 = {"name": "wzp6_ee_qqH_Hmumu_ecm240", "datadir": f"{baseDir}wzp6_ee_qqH_Hmumu_ecm240/",  "xsec": 1.161e-05}
    wzp6_ee_ssH_Hmumu_ecm240 = {"name": "wzp6_ee_ssH_Hmumu_ecm240", "datadir": f"{baseDir}wzp6_ee_ssH_Hmumu_ecm240/",  "xsec": 6.519e-06}
    wzp6_ee_mumuH_Hmumu_ecm240 = {"name": "wzp6_ee_mumuH_Hmumu_ecm240", "datadir": f"{baseDir}wzp6_ee_mumuH_Hmumu_ecm240/",  "xsec": 1.472e-06}
    #p8_ee_WW_ecm240 = {"name" : "p8_ee_WW_ecm240", "datadir": f"{baseDir}p8_ee_WW_ecm240/", "xsec": 16.4385}
    p8_ee_ZZ_ecm240 = {"name" : "p8_ee_ZZ_ecm240", "datadir": f"{baseDir}p8_ee_ZZ_ecm240/", "xsec": 1.35899}
    p8_ee_WW_ecm240 = {"name" : "p8_ee_WW_ecm240", "datadir": f"{baseDir}p8_ee_WW_mumu_ecm240/", "xsec":  0.25792}
    datasets = [wzp6_ee_nunuH_Hmumu_ecm240, wzp6_ee_eeH_Hmumu_ecm240 , wzp6_ee_ccH_Hmumu_ecm240, wzp6_ee_tautauH_Hmumu_ecm240, wzp6_ee_bbH_Hmumu_ecm240, wzp6_ee_qqH_Hmumu_ecm240, wzp6_ee_ssH_Hmumu_ecm240, wzp6_ee_mumuH_Hmumu_ecm240, p8_ee_ZZ_ecm240, p8_ee_WW_ecm240]
    
    result = functions.build_and_run(datasets, build_graph, "tmp/output_mumu.root", maxFiles=args.maxFiles, norm=True, lumi=7200000)


