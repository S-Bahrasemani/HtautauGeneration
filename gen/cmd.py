import re

JO_PATTERN = re.compile('(MC12.(?P<id>\d+).PowhegPythia8_AU2CT10_(?P<mode>gg|VBF)H(?P<mass>\d+)_tautau.py)?$')
DIR_PATTERN = re.compile('(MC12.(?P<id>\d+).PowhegPythia8_AU2CT10_(?P<mode>gg|VBF)H(?P<mass>\d+)_tautau_run(?P<seed>\d+))?$')

def get_dir_mode(d):
    match = re.match(DIR_PATTERN, d)
    if match:
        return match.group('mode')
    else:
        return None

def get_dir_mass(d):
    match = re.match(DIR_PATTERN, d)
    if match:
        return int(match.group('mass'))
    else:
        return None


def generate_cmd(jo, seed=1234, nevents=5000):
    """
    EVNT GENERATION COMMAND LINE
    """
    match = re.match(JO_PATTERN, jo)
    if match:
        datasetid = match.group('id')
        run_dir = jo.replace('.py', '_run{0}'.format(seed))
        
        cmd_args = [
            "Generate_trf.py",
            "ecmEnergy=13000",
            "runNumber={0}".format(datasetid),
            "firstEvent=0",
            "maxEvents={0}".format(nevents), 
            "randomSeed={0}".format(seed), 
            "jobConfig={0}".format(jo),
            "outputEVNTFile=EVNT_{0}_seed{1}.root".format(datasetid, seed),
            "evgenJobOpts=MC12JobOpts-00-14-70_v8.tar.gz"
            ]
        cmd = ' '.join(cmd_args)
        cmd = 'mkdir prod/{0} && mv joboptions/{1} prod/{0}/ && cd prod/{0} && '.format(run_dir, jo) + cmd
        return cmd
    else:
        return None


def d3pd_cmd(run_dir, input_root):
    """
    EVNT -> NTUP_TRUTH COMMAND LINE
    """
    cmd_args = [
        'Reco_trf.py',
        'preExec="rec.doApplyAODFix.set_Value_and_Lock(False);from D3PDMakerConfig.D3PDMakerFlags import D3PDMakerFlags;D3PDMakerFlags.TruthWriteEverything=True"',
        'inputEVNTFile={0}'.format(input_root),
        'outputNTUP_TRUTHFile=ntup.{0}'.format(input_root),
        ]
    cmd = ' '.join(cmd_args)
    cmd = 'cd {0} && '.format(run_dir) + cmd
    return cmd


def flat_cmd():
    """
    TRUTH -> FLAT TREE COMMAND LINE
    """

    pass
