import glob
import json
import sys
import os
import shutil
import ROOT
from opal_tracking import OpalTracking
import xboa.common as common
from xboa.hit import Hit
from xboa.algorithms.tune import FFTTuneFinder

class Tune(object):
    def __init__(self, file_name):
        """
        Find the tune. If file_name is specified, use that file_name, otherwise
        generate a new one by tracking
        """
        self.closed_orbits_cached = None
        self.tmp_dir = None
        self.unique_id = 0
        self.just_plot = file_name != None
        if self.just_plot:
            self.opal = "/bin/echo"
        else:
            self.opal = "/home/cr67/OPAL/source/opal_dev/src/trunk/src/opal"
        self.lattice_src = "lattices/KurriMainRingTuneComparison/KurriMainRingTuneComparison.in"
        self.lattice = "/Tune.tmp"
        self.beam_file = "/disttest.dat"
        self.log_file = "/log"
        if self.just_plot:
            self.output_filename = file_name
        else:
            self.output_filename = "PROBE1.loss"
        self.nturns = 100.1
        self._load_closed_orbits("closed_orbits.ref")
        self.output = "tunes.out"

    def find_tune(self):
        fout = open(self.output, "w")
        for energy, position in sorted(self.closed_orbits_cached.iteritems())[10::5]:
            print "Finding tune at", energy, "MeV and closed orbit x=", position, "mm"
            hit = self._reference(energy)
            hit["x"] = position
            tune_info = {
                "energy":energy,
            }
            for axis, delta in [("x", 10.), ("y", 1.)]:
                self._temp_dir()
                tracking = OpalTracking(self.tmp_dir+self.lattice,
                                        self.tmp_dir+self.beam_file,
                                        self._reference(energy),
                                        self.output_filename,
                                        self.opal,
                                        self.tmp_dir+self.log_file)
                finder = FFTTuneFinder()
                if self.just_plot:
                    tracking.do_tracking = False
                    finder.run_tracking(axis, delta, hit, tracking) # DUMMY
                else:
                    common.substitute(
                        self.lattice_src, 
                        self.tmp_dir+self.lattice, {
                            "__energy__":energy,
                            "__nturns__":self.nturns,
                            "__beamfile__":self.tmp_dir+self.beam_file,
                            "__stepsize__":0.1,
                    })
                    finder.do_tracking = True
                    finder.run_tracking(axis, delta, hit, tracking)
                    os.rename(self.output_filename, self.tmp_dir+self.output_filename)
                canvas, hist, graph = finder.plot_signal(axis+" signal - KE: "+str(energy)+" MeV")
                self._print_canvas(canvas, axis, "signal", energy)
                try:
                    tune = finder.get_tune(self.nturns/10.)
                    for i in range(len(finder.k_mag_x)):
                        k_x = finder.k_mag_x[i]
                        k_y = finder.k_mag_y[i]
                        k_t = finder._sft(k_x)
                        #print str(round(k_x, 3)).ljust(10), str(round(k_y, 3)).ljust(10), str(round(k_t, 3)).ljust(10)
                    tune = finder.get_tune(1./self.nturns/100.)
                except:
                    sys.excepthook(*sys.exc_info())
                canvas, hist, graph = finder.plot_fft(axis+" FFT - KE: "+str(energy)+" MeV")
                graph.SetMarkerStyle(6)
                graph.Draw()
                self._print_canvas(canvas, axis, "fft", energy)
                tune_info[axis+"_tune"] = tune
                tune_info[axis+"_signal"] = finder.u
            for key in sorted(tune_info.keys()):
                if "signal" not in key:
                    print "   ", key, tune_info[key]
            print >> fout, json.dumps(tune_info)
            fout.flush()

    def _temp_dir(self):
        self.tmp_dir = "tmp/tune/"+str(self.unique_id)+"/"
        if self.just_plot:
            return
        if self.unique_id == 0:
            for tmp_dir in glob.glob("tmp/tune/*"):
                shutil.rmtree(tmp_dir)
        os.makedirs(self.tmp_dir)
        self.unique_id += 1

    def _load_closed_orbits(self, filename):
        fin = open(filename)
        closed_orbits = [json.loads(line) for line in fin.readlines()]
        closed_orbits_energy = [orbit[0] for orbit in closed_orbits]
        closed_orbits_x = [orbit[1:][0][0] for orbit in closed_orbits]
        closed_orbits_dict = dict(zip(closed_orbits_energy, closed_orbits_x))
        self.closed_orbits_cached = closed_orbits_dict

    def _reference(self, energy):
        hit_dict = {}
        hit_dict["pid"] = 2212
        hit_dict["mass"] = common.pdg_pid_to_mass[2212]
        hit_dict["charge"] = 1
        hit_dict["x"] = 4600.
        hit_dict["kinetic_energy"] = energy
        return Hit.new_from_dict(hit_dict, "pz")

    def _print_canvas(self, canvas, axis, name, energy):
        name = "plots/"+axis+"_"+name+"_energy="+str(energy)
        for format in ["png", "root"]:
            canvas.Print(name+"."+format)

def main():
    tune = Tune(file_name = "tmp/tune_step=0.1_nturns=100_dx=1./PROBE1.loss")
    tune.find_tune()
    raw_input()

if __name__ == "__main__":
    main()
    print "Finished"
    #raw_input()










