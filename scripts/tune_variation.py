import math
import numpy
import copy
import glob
import json
from xboa import Common

def parse_file(file_name, heading, types):
    if len(heading) != len(types):
        raise KeyError("Heading mismatched to types in parse_file "+file_name)
    fin = open(file_name)
    line = fin.readline()
    data = {}
    for item in heading:
        data[item] = []
    line = fin.readline()[:-1]
    while line != "":
        words = line.split()
        if len(words) != len(heading):
            print "Line\n  "+line+"\nmismatched to heading\n  "+str(heading)+"\nin parse_file "+file_name
        else:
            words = [types[i](x) for i, x in enumerate(words)]
            for i, item in enumerate(heading):
                data[item].append(words[i])
        line = fin.readline()[:-1]
    print "Got data from file "+file_name
    return data

def parse_track_file(directory = "."):
    file_name = directory+"/OpalRingTest-trackOrbit.dat"
    heading = ["id", "x", "px", "y", "py", "z", "pz"]
    types = [str]+[float]*(len(heading)-1)
    data = parse_file(file_name, heading, types)
    data = r_phi_track_file(data)
    return data

def r_phi_track_file(data):
    data = copy.deepcopy(data)
    data["r"] = range(len(data["x"]))
    data["phi"] = range(len(data["x"]))
    for i in range(len(data["r"])):
        data["r"][i] = (data["x"][i]**2+data["y"][i]**2.)**0.5
        data["phi"][i] = math.degrees(math.atan2(data["x"][i], data["y"][i]))
    return data

def parse_probe_file(directory = "."):
    file_name = directory+"/PROBE1.loss"
    # Element x (mm),  y (mm),  z (mm),  px ( ),  py ( ),  pz ( ), id,  turn,  time (ns) 
    heading = ["element", "x", "y", "z", "px", "py", "pz", "id", "turn", "time"]
    types = [str]+[float]*6+[int, int, float]
    return parse_file(file_name, heading, types)   
  
def make_xy_plot(data, canvas = None):
    new_canvas = canvas == None
    if new_canvas:
        canvas = Common.make_root_canvas("x vs y")
    hist, graph = Common.make_root_graph("x_vs_y", data["x"], "x [mm]", data["y"], "y [mm]", sort=False, xmin=-3000., xmax=3000., ymin=-3000., ymax=3000.)
    if new_canvas:
        hist.Draw()
    graph.Draw('p')
    canvas.Update()
    print "made plot"

def make_rphi_plot(data, canvas = None):
    new_canvas = canvas == None
    if new_canvas:
        canvas = Common.make_root_canvas("r vs phi")
    hist, graph = Common.make_root_graph("r_vs_phi", data["phi"], "phi [^{o}]", data["r"], "r [mm]", sort=False, xmin=-45., xmax=45., ymin=2100., ymax=2600.)
    if new_canvas:
        hist.Draw()
    graph.Draw()
    canvas.Update()
    print "made plot"

def make_phibz_plot(data, canvas = None):
    new_canvas = canvas == None
    bz = [bz/10. for bz in data["bz"]]
    if new_canvas:
        canvas = Common.make_root_canvas("by vs phi")
    hist, graph = Common.make_root_graph("by_vs_phi", data["phi"], "phi [^{o}]", bz, "bz [T]", sort=False, xmin=-45., xmax=45., ymin=-1., ymax=1.)
    if new_canvas:
        hist.Draw()
    graph.Draw()
    canvas.Update()
    print "made plot"


def make_horizontal_phase_space_plot(data, canvas = None):
    data = copy.deepcopy(data)
    new_canvas = canvas == None
    if new_canvas:
        pz = (11.**2+938**2)**0.5
        canvas = Common.make_root_canvas("x vs px")
        data["x"] = [x-2350. for x in data["x"]]
        data["px"] = [px*pz for px in data["px"]]
    hist, graph = Common.make_root_graph("x_vs_px", data["x"], "x [mm]", data["px"], "px [MeV/c]", sort=False, xmin=-250, xmax=+250, ymin=-25., ymax=+25.)
    if new_canvas:
        hist.Draw()
    graph.Draw('p')
    graph.SetMarkerStyle(6)
    canvas.cd()
    points = [(data["x"][i], data["px"][i]) for i in range(len(data["z"]))]
    mean, cov = Common.fit_ellipse(points, 100.)
    try:
        fit_function = Common.make_root_ellipse_function(mean, cov, [1., 10., 100., 400.], -250., +250., -25., +25.)
        fit_function.Draw("SAME")
    except numpy.linalg.linalg.LinAlgError:
        pass
    canvas.Update()
    canvas.Print("dynamic_aperture_x.png")
    print "made plot"

def make_vertical_phase_space_plot(data, canvas = None):
    data = copy.deepcopy(data)
    new_canvas = canvas == None
    if new_canvas:
        canvas = Common.make_root_canvas("y vs py")
    hist, graph = Common.make_root_graph("y_vs_py", data["z"], "z [mm]", data["pz"], "z'", sort=False, ymin=-25., ymax=+25., xmin=-250, xmax=+250)
    if new_canvas:
        hist.Draw()
    graph.Draw('p')
    graph.SetMarkerStyle(6)
    canvas.Update()
    print "made plot"

def _norm_to_1(array):
    max_array = max([abs(x) for x in array])
    return [x/max_array for x in array]

def _find_tune(f_spec):
    #return [float(i)/len(f_spec) for i, x in enumerate(f_spec) if x > 0.5]
    return [f_spec.index(max(f_spec[1:len(f_spec)/2]))/float(len(f_spec))]

def get_tune(a1, a2, make_plots=True):
    print 'get_tune'
    w1 = numpy.fft.fft(numpy.array(a1))
    re_w1 = [z.real for z in w1]
    im_w1 = [z.imag for z in w1]
    mag_w1 = [(z*z.conjugate()).real**0.5 for z in w1]
    w2 = numpy.fft.fft(numpy.array(a2))
    re_w2 = _norm_to_1([z.real for z in w2])
    im_w2 = _norm_to_1([z.imag for z in w2])
    mag_w2 = _norm_to_1([(z*z.conjugate()).real**0.5 for z in w2])
    tune = (_find_tune(mag_w1), _find_tune(mag_w2))
    if make_plots:
        a1 = _norm_to_1(a1)
        canvas = Common.make_root_canvas('fft')
        hist, graph_mag1 = Common.make_root_graph('horizontal fft', [float(x)/len(mag_w1) for x in range(len(mag_w1))], 'frequency', mag_w1, 'magnitude', ymin=-max(mag_w1)*1.1, ymax=max(mag_w1)*1.1)
        hist.Draw()       
        graph_mag1.SetLineColor(1)
        graph_mag1.Draw()
        hist, graph_x1 = Common.make_root_graph('x', [float(x)/len(a1) for x in range(len(a1))], 'turn', a1, 'var', ymin=-1.1, ymax=1.1)
        graph_x1.Draw('p')
        hist, graph_re1 = Common.make_root_graph('re(z)', [float(x)/len(re_w1) for x in range(len(re_w1))], 'tune', re_w1, 'var')
        graph_re1.SetLineColor(2)
        graph_re1.Draw()
        hist, graph_im1 = Common.make_root_graph('im(z)', [float(x)/len(im_w1) for x in range(len(im_w1))], 'tune', im_w1, 'var')
        graph_im1.SetLineColor(4)
        graph_im1.Draw()
        """
        hist, graph2 = Common.make_root_graph('vertical fft', [float(x)/len(mag_w2) for x in range(len(mag_w2))], 'frequency', mag_w2, 'magnitude')
        graph2.SetLineColor(4)
        graph2.Draw()
        """
        Common.make_root_legend(canvas, [graph_mag1, graph_re1, graph_im1])
        canvas.Update()
    return tune

def plot(data):
    ke, tune_x, tune_y = [], [], []
    for directory in sorted(data.keys()): 
        my_data = data[directory]
        ke.append(my_data['kinetic_energy'])
        try:
            tune_x.append(2.-my_data['x_tune'][1.][0])
            tune_y.append(2.+my_data['y_tune'][1.][0])
        except KeyError:
            continue
        print ke[-1], tune_x[-1], tune_y[-1]
    canvas = Common.make_root_canvas('tune_variation_with_energy')
    hist, graph1 = Common.make_root_graph('horizontal tune', ke, 'Energy [MeV]', tune_x, 'fractional tune', ymin=1.0, ymax=3.0)
    hist.Draw()
    graph1.SetLineColor(2)
    graph1.Draw()
    hist, graph2 = Common.make_root_graph('vertical tune', ke, 'Energy [MeV]', tune_y, 'fractional tune', ymin=1.0, ymax=3.0)
    graph2.SetLineColor(4)
    graph2.Draw()
    Common.make_root_legend(canvas, [graph1, graph2])
    canvas.Update()

def main():
    tune_data = {}
    for directory in sorted(glob.glob("runs_nturns=100_step=0.01-ns/OpalRingTest_ke=*")):
        #track_data = parse_track_file(directory)
        try:
            probe_data = parse_probe_file(directory)
        except IOError:
            continue
        tune_data[directory] = {}
        print probe_data['x'][0]
        x_tunes, y_tunes = {}, {}
        tune_data[directory]['particles'] = {}
        tune_data[directory]['closed_orbit'] = probe_data['x'][0]
        particle_number = 1
        x1 = [x-probe_data['x'][0] for i, x in enumerate(probe_data['x']) if probe_data['id'][i] == particle_number]
        particle_number = 2
        y1 = [y-probe_data['y'][0] for i, y in enumerate(probe_data['z']) if probe_data['id'][i] == particle_number]
        tune_data[directory]['particles'][particle_number] = len(x1)
        if len(x1) > 4 and abs(x1[0]) > 0.1 and len(y1) > 4 and abs(y1[0]) > 0.1:
            print  directory[-5:]
            x_tunes[y1[0]], y_tunes[y1[0]] = get_tune(x1, y1, directory[-5:] == "ke=11")
        tune_data[directory]['kinetic_energy'] = float(directory[directory.find('ke=')+3:])
        tune_data[directory]['x_tune'] = x_tunes
        tune_data[directory]['y_tune'] = y_tunes
    for key_1 in sorted(tune_data.keys()):
        print key_1
        for key_2 in sorted(tune_data[key_1].keys()):
            print '   ', key_2, tune_data[key_1][key_2]
    plot(tune_data)

    #make_horizontal_phase_space_plot(probe_data)
    #make_vertical_phase_space_plot(probe_data)

if __name__ == "__main__":
    main()
    print "press CR to continue"
    raw_input()

