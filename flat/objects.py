from .mixins import *
from .taudecay import TauDecay
from . import log; log = log[__name__]

def define_objects(tree):

    tree.define_collection(
        name="taus",
        prefix="tau_",
        size="tau_n",
        mix=MCTauFourMomentum)
    tree.define_collection(
        name="mc",
        prefix="mc_",
        size="mc_n",
        mix=MCParticle)
    tree.define_collection(
        name="jets",
        prefix="jet_AntiKt4TruthJets_",
        size="jet_AntiKt4TruthJets_n")
