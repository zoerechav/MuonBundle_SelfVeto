#!/usr/bin/env python
import numpy as np
import time
import argparse

# veto and MCEq things
from nuVeto.nuveto import passing
from nuVeto.utils import Units
import crflux.models as pm
from nuVeto.mu import interp

start_time = time.time()


def modified_fsigmoid(emu, a, b, x1, k1, c):
    """p_light parametrisation: sigmoid and gaussian peak

    Parameters
    ----------
    emu : float
        muon bundle energy
    a : float
        position of sigmoid
    b : float
        slope of sigmoid
    x1 : float
        position of peak
    k1 : float
        width of peak
    c : float
        height of peak

    Returns
    -------
    float
        p_light for the given energy
    """
    x = emu
    return 1 / (1 + np.exp(-b * (np.log10(x) - np.log10(a)))) + c * np.exp(
        -(np.log10(x) - np.log10(x1))**2 / k1)


def get_pl_params(flavor, depth, cos_zen, shift):
    """get the parameters to parametrise p_light for a given bin

    Parameters
    ----------
    flavor : str
        neutrino flavor (mu or e, antineutrinos will be treated as the
        corresponding neutrinos)
    depth : float
        depth in the detector in km
    cos_zen : float
        cos of zenith of the incoming event
    shift : str
        either base or m2,m1,p1,p2 for the pls shifted by statistical errors

    Returns
    -------
    floats
        five parameters to parametrise pl in the bin
    """
    if "_bar" in flavor:
        flavor.strip("_bar")
    # read file containing all pl params
    pars = np.genfromtxt("/home/jhellrung/calc_PF/PL_params.txt",
                         delimiter=",",
                         dtype=None,
                         names=True,
                         encoding=None)
    # available values
    flavors = np.unique(pars["flavor"])
    depths = np.unique(pars["depth"])
    cos_zens = np.unique(pars["cos_zen"])

    # check if flavor is available
    if flavor not in flavors:
        raise ValueError(f"flavor {flavor} is not available, only {flavors}")

    # choose available values closest to arguments
    if cos_zen not in cos_zens:
        print(f"cos_zen {cos_zen} not in available values")
        cos_zen = cos_zens[(np.abs(cos_zens - cos_zen)).argmin()]
        print(f"chosing nearest value {cos_zen}")

    if depth not in depths:
        print(f"depth {depth} not in available values")
        depth = depths[(np.abs(depths - depth)).argmin()]
        print(f"chosing nearest value {depth}")

    # select pl-parameters for arguments
    pars = pars[pars["shift"] == shift]
    pars = pars[pars["flavor"] == flavor]
    pars = pars[pars["depth"] == depth]
    pars = pars[pars["cos_zen"] == cos_zen]
    if len(pars) != 1:
        raise KeyError(
            f"There are {len(pars)} sets of pl-params for the given variables"
        )
    a = pars[0]["a"]
    b = pars[0]["b"]
    x1 = pars[0]["x1"]
    k1 = pars[0]["k1"]
    c = pars[0]["c"]
    return a, b, x1, k1, c


# dicts with available models
hadmodel_dict = {
    'SIBYLL': 'SIBYLL2.3c',
    'EPOSLHC': 'EPOS-LHC',
    'QGSJET': 'QGSJet-II-04',
    'DPMJET': 'DPMJET-III-19.1',
}

pmodel_dict = {
    'GST3': (pm.GaisserStanevTilav, "3-gen"),
    'GST4': (pm.GaisserStanevTilav, "4-gen"),
    'GHandHGH3a': (pm.CombinedGHandHG, "H3a"),
    'GHandHGH4a': (pm.CombinedGHandHG, "H4a"),
    'GaisserH3a': (pm.HillasGaisser2012, "H3a"),
    'GaisserH4a': (pm.HillasGaisser2012, "H4a"),
    'PolyGonato': (pm.PolyGonato, "poly-gonato"),
    'Thunman': (pm.Thunman, "TIG"),  # broken
    'ZatsepinSokolskaya': (pm.ZatsepinSokolskaya, "default"),
    'ZatsepinSokolskayaPamela': (pm.ZatsepinSokolskaya, "pamela"),
    'GaisserHonda': (pm.GaisserHonda, "GH"),
    'GlobalSplineFitBeta': (pm.GlobalSplineFitBeta, "GSF spl"),
}

density_dict = {
    'January': ('PL_SouthPole', 'January'),
    'June': ('SouthPole', 'June'),
    'August': ('PL_SouthPole', 'August'),
    'December': ('SouthPole', 'December'),
}

parser = argparse.ArgumentParser(
    description='Inputs for passing fraction calculator')
parser.add_argument("-cr",
                    dest="cr_model",
                    type=str,
                    default='GaisserH4a',
                    help="Mass Composition Model")
parser.add_argument("-had",
                    dest="had_model",
                    type=str,
                    default='SIBYLL',
                    help="Hadronic Interaction type")
parser.add_argument("-flux",
                    dest="flux",
                    type=str,
                    default="conv",
                    help="neutrino flux type (conv or pr)")
parser.add_argument("-flavor",
                    dest="flavor",
                    type=str,
                    default="mu",
                    help="neutrino flavor (mu, mu_bar, e or e_bar)")
parser.add_argument("-depth",
                    dest="depth",
                    type=float,
                    default=1.4,
                    help="depth in km")
parser.add_argument("-density",
                    dest="density",
                    type=str,
                    default="January",
                    help="month for the density model")
parser.add_argument("-preach",
                    dest="preach",
                    type=str,
                    default='/home/jhellrung/nuVeto/nuVeto/resources/mu/mmc/ice_allm97.pklz',
                    help="preach file to be used")
parser.add_argument("-outdir",
                    dest="outdir",
                    type=str,
                    default='/home/jhellrung/calc_PF/PFs/',
                    help="folder to write output")
args = parser.parse_args()

outdir = args.outdir
if outdir[-1] != "/":
    outdir += "/"

# setup our parameters of interest
preach_file = args.preach
flavor = args.flavor
flux = args.flux
depth = args.depth
hadr = hadmodel_dict[args.had_model]
cr_model = args.cr_model

# flavors = ["mu", "mu_bar", "e",  "e_bar"]
# fluxes = ["conv", "pr"]
# depths = [1.4, 2.0, 2.1]
shifts = ["base", "p1", "p2", "m1", "m2"]

energies = np.logspace(1, 7, 51)  # 10GeV-1PeV
# round to have nice predictable behavior when selecting PL
cos_zeniths = np.round(np.linspace(0, 1, 11), 3)

# put all settings in array to save later
settings = np.array(
    [preach_file, flavor, flux, depth, hadr, cr_model, args.density])

print()
print("Preach:", preach_file)
print("flavor:", flavor)
print("flux:", flux)
print("depth:", depth)
print("HI model:", hadr)
print("CR model:", cr_model)

neut = flux + " " + flavor

# calculate PF table for all shifted pls
for shift in shifts:
    outfile_name = f"PF_{shift}_neut_type_{neut}_at_depth_{depth}km"
    outfile = outdir + outfile_name

    PFs = np.empty((len(energies), len(cos_zeniths)))
    for i in range(len(cos_zeniths)):
        # get pl params and define the correct pl
        a, b, x1, k1, c = get_pl_params(flavor, depth, cos_zeniths[i], shift)

        def pl(emu):
            return modified_fsigmoid(emu, a, b, x1, k1, c)

        # convolve preach and plight to prpl
        prpl = interp(preach_file, pl)

        # calculate PF for all energies
        for j in range(len(energies)):
            PFs[j, i] = passing(
                energies[j],
                cos_zeniths[i],
                kind=f"{flux} nu_{flavor}",
                pmodel=pmodel_dict[args.cr_model],
                hadr=hadr,
                depth=depth * Units.km,
                density=('CORSIKA', density_dict[args.density]),
                prpl=prpl)
    # save the PF table (including the bins and the settings)
    np.savez(outfile,
             PFs=PFs,
             energies=energies,
             cos_zeniths=cos_zeniths,
             settings=settings)
    print("saving numpy array: " + outfile)
print(f"took {(time.time()-start_time)/60:.2f} mins")
