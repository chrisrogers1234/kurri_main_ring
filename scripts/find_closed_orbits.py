import sys
import json
sys.path.insert(1, "scripts")
from opal_tracking import OpalTracking
import xboa.Common as common
from xboa.Hit import Hit
from xboa.algorithms.closed_orbit import EllipseClosedOrbitFinder

import ROOT

def reference(energy):
    hit_dict = {}
    hit_dict["pid"] = 2212
    hit_dict["mass"] = common.pdg_pid_to_mass[2212]
    hit_dict["charge"] = 1
    hit_dict["x"] = 4600.
    hit_dict["kinetic_energy"] = energy
    return Hit.new_from_dict(hit_dict, "pz")

def plot_iteration(i, iteration, energy):
    canvas, hist, graph, fit = iteration.plot_ellipse("x", "px", "mm", "MeV/c")
    hist.SetTitle('KE='+str(energy)+' iter='+str(i))
    canvas.Update()
    canvas.Print("plots/closed_orbit-i_"+str(i)+"-ke_"+str(energy)+".png")

def find_closed_orbit(energy, seed):
    print "ENERGY", energy
    tmp_dir = "tmp/find_closed_orbits/"
    common.substitute('lattices/KurriMainRingTuneComparison/KurriMainRingTuneComparison.in', tmp_dir+'/Kurri_ADS_Ring.tmp', {'__energy__':energy})
    ref_hit = reference(energy)
    tracking = OpalTracking(tmp_dir+'/Kurri_ADS_Ring.tmp', tmp_dir+'/disttest.dat', ref_hit, 'PROBE1.loss', "/home/cr67/OPAL/source/opal_dev/src/trunk/src/opal", tmp_dir+"/log")
    mass = common.pdg_pid_to_mass[2212]
    seed_hit = ref_hit.deepcopy()
    seed_hit["x"] = seed[0]
    finder = EllipseClosedOrbitFinder(tracking, seed_hit)
    generator = finder.find_closed_orbit_generator(["x", "px"], 5)
    for i, iteration in enumerate(generator):
        print iteration.points, iteration.centre
        if i == 0:
            plot_iteration(i, iteration, energy)
        if i >= 10 or abs(iteration.points[0][0] - iteration.centre[0]) < 0.1 and abs(iteration.points[0][1]) < 0.1:
            break
    plot_iteration(i, iteration, energy)

    return tracking.last[0]

if __name__ == "__main__":
      next_seed = [4602., 0.0]
      fout = open('find_closed_orbit.out', 'w')
      energy_list = [11]+range(15, 151, 2)
      for i, energy in enumerate(energy_list):
          is_batch = len(energy_list) > 5 and i%(len(energy_list)/5) != 0
          print "Batch", is_batch, len(energy_list), i%(len(energy_list)/5)
          ROOT.gROOT.SetBatch(is_batch)
          hit_list = find_closed_orbit(energy, next_seed)
          next_seed = [hit_list[0]["x"], 0.]
          output = [energy]+[[hit["x"], hit["t"]] for hit in hit_list]
          print >> fout, json.dumps(output)
          fout.flush()
      if len(energy_list) < 5:
          raw_input()

