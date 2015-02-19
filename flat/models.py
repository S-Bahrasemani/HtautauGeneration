from ROOT import TLorentzVector

from rootpy.tree import TreeModel, FloatCol, IntCol, BoolCol, CharCol
from rootpy.vector import LorentzRotation, LorentzVector, Vector3, Vector2
from rootpy.extern.hep import pdg
from rootpy import log
ignore_warning = log['/ROOT.TVector3.PseudoRapidity'].ignore(
    '.*transvers momentum.*')

from .taudecay import TauDecay

class EventModel(TreeModel):
    runnumber = IntCol()
    evtnumber = IntCol()
    weight = FloatCol()
    hadhad = IntCol() # 1 or 0
    lephad = IntCol() # 1 or 0
    leplep = IntCol() # 1 or 0
    
class FourMomentum(TreeModel):
    pt = FloatCol()
    p = FloatCol()
    et = FloatCol()
    e = FloatCol()
    eta = FloatCol(default=-1111)
    phi = FloatCol(default=-1111)
    m = FloatCol()

    @classmethod
    def set(cls, this, other):
        if isinstance(other, TLorentzVector):
            vect = other
        else:
            vect = other.fourvect
        this.pt = vect.Pt()
        this.p = vect.P()
        this.et = vect.Et()
        this.e = vect.E()
        this.m = vect.M()
        with ignore_warning:
            this.phi = vect.Phi()
            this.eta = vect.Eta()

class TrueMet(FourMomentum):
    @classmethod
    def set(cls, this, miss1, miss2):
        FourMomentum.set(this, miss1 + miss2)

class TrueTau(FourMomentum + FourMomentum.prefix('vis_')):
    nProng = IntCol(default=-1111)
    nPi0 = IntCol(default=-1111)
    charge = IntCol()
    flavor = CharCol()
    pdgId = IntCol(default=-1111)
    index = IntCol()
    @classmethod
    def set_vis(cls, this, other):
        if isinstance(other, TLorentzVector):
            vect = other
        else:
            vect = other.fourvect
        this.vis_pt = vect.Pt()
        this.vis_p = vect.P()
        this.vis_et = vect.Et()
        this.vis_e = vect.E()
        this.vis_m = vect.M()
        with ignore_warning:
            this.vis_phi = vect.Phi()
            this.vis_eta = vect.Eta()

class TrueTauBlock(TrueTau.prefix('tau1_') + TrueTau.prefix('tau2_') + TrueMet.prefix('met_')):
    
    dR_taus = FloatCol()
    dEta_taus = FloatCol()
    dPhi_taus = FloatCol()
    
    dPhi_taus_met = FloatCol()
    dPhi_tau1_met = FloatCol()
    dPhi_tau2_met = FloatCol()

    pt_sum_taus_met = FloatCol()
    pt_tot_taus_met = FloatCol()

    pt_sum_tau1_tau2 = FloatCol() # TO BE SET
    pt_tot_tau1_tau2 = FloatCol() # TO BE SET
    pt_ratio_tau1_tau2 = FloatCol()

    transverse_mass_tau1_tau2 = FloatCol() # TO BE SET
    transverse_mass_tau1_met = FloatCol() # TO BE SET
    transverse_mass_tau2_met = FloatCol() # TO BE SET

    @classmethod 
    def set(cls, tree, tau1, tau2):


        TrueTau.set(tree.tau1, tau1.fourvect)
        TrueTau.set(tree.tau2, tau2.fourvect)

        TrueTau.set_vis(tree.tau1, tau1.decay.fourvect_vis)
        TrueTau.set_vis(tree.tau2, tau2.decay.fourvect_vis)


        tree.tau1.index = tau1.index
        tree.tau1.charge = tau1.charge
        tree.tau1.flavor = 'l' if tau1.decay.leptonic else 'h'
        if tree.tau1.flavor == 'l':
            tree.tau1.pdgId = pdg.mu if tau1.decay.leptonic_muon else pdg.e
        else:
            tree.tau1.nProng = tau1.decay.nprong
            tree.tau1.nPi0s = tau1.decay.nneutrals
        
        tree.tau2.index = tau2.index
        tree.tau2.charge = tau2.charge
        tree.tau2.charge = tau2.charge
        tree.tau2.flavor = 'l' if tau2.decay.leptonic else 'h'
        if tree.tau2.flavor == 'l':
            tree.tau2.pdgId = pdg.mu if tau2.decay.leptonic_muon else pdg.e
        else:
            tree.tau2.nProng = tau2.decay.nprong
            tree.tau2.nPi0s = tau2.decay.nneutrals

        MET = tau1.decay.fourvect_missing + tau2.decay.fourvect_missing
        TrueMet.set(tree.met, tau1.decay.fourvect_missing, tau2.decay.fourvect_missing)

        vis_tau1 = tau1.decay.fourvect_vis
        vis_tau2 = tau2.decay.fourvect_vis
        tree.dR_taus = vis_tau1.DeltaR(vis_tau2)
        tree.dEta_taus = abs(vis_tau1.Eta() - vis_tau2.Eta())
        tree.dPhi_taus = abs(vis_tau1.DeltaPhi(vis_tau2))
        
        vis_taus = vis_tau1 + vis_tau2

        tree.dPhi_taus_met = abs(vis_taus.DeltaPhi(MET))
        tree.dPhi_tau1_met = abs(vis_tau1.DeltaPhi(MET))
        tree.dPhi_tau2_met = abs(vis_tau2.DeltaPhi(MET))
        
        tree.pt_sum_taus_met = (vis_taus + MET).Pt()
        tree.pt_tot_taus_met = 0 #( tlv_tau1 + tlv_lep2 + tlv_met ).Pt() / ( lep1_pt + lep2_pt + met_et ) # TO BE SET

        tree.pt_sum_tau1_tau2 = vis_taus.Pt()
        tree.pt_tot_tau1_tau2 = 0 # TO BE SET

        tree.transverse_mass_tau1_tau2 = vis_taus.Mt() 
        tree.transverse_mass_tau1_met = (vis_tau1 + MET).Mt()
        tree.transverse_mass_tau2_met = (vis_tau2 + MET).Mt()


        if vis_tau2.Pt() != 0:
            tree.pt_ratio_tau1_tau2 = vis_tau1.Pt() / vis_tau2.Pt()
        else:
            tree.pt_ratio_tau1_tau2 = 0



class TrueJet(FourMomentum):
    index = IntCol()



class TrueJetBlock(TrueJet.prefix('jet1_') + TrueJet.prefix('jet2_')):
    
    dEta_jet1_jet2 = FloatCol()
    mass_je1_jet2 = FloatCol()
    sum_PT_jet1_jet2_vistau1_vistau2_ETmiss = FloatCol()
    
    product_jet1_eta_jet2_eta = FloatCol()
    
    @classmethod 
    def set(cls, tree, jet1, jet2):
        if jet1 is not None:
            tree.jet1.index = jet1.index
            FourMomentum.set(tree.jet1, jet1)
        if jet2 is not None:
            tree.jet2.index = jet2.index
            FourMomentum.set(tree.jet2, jet2)
            tree.dEta_jet1_jet2 = abs(jet1.eta - jet2.eta)
        

def get_model():
    model = EventModel + TrueTauBlock + FourMomentum.prefix('higgs_') + TrueJetBlock
    return model
