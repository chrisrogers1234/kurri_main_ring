import math
import copy
import glob
import ROOT
from xboa import Common

class RootObjects:
    histograms = []
    canvases = []
    graphs = []

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
    file_name = glob.glob(directory+"/*-trackOrbit.dat")[0]
    heading = ["id", "x", "px", "y", "py", "z", "pz", "bx", "by", "bz"]
    types = [str]+[float]*(len(heading)-1)
    data = parse_file(file_name, heading, types)
    data = r_phi_track_file(data)
    return data

def r_phi_track_file(data):
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

def make_phi_plot(data, y_var="bz", canvas = None):
    new_canvas = canvas == None
    bz = [bz/10. for bz in data[y_var]]
    phi = [f for f in data["phi"]]
    for i, f in enumerate(phi):
        while f > 30.:
            f -= 30. 
        while f < 0.:
            f += 30. 
        phi[i] = f
    if new_canvas:
        canvas = Common.make_root_canvas(y_var+" vs phi")
    hist, graph = Common.make_root_graph(y_var+"_vs_phi", phi, "phi [^{o}]", bz, y_var+" [T]", sort=False, xmin=-5., xmax=35., ymin=-1., ymax=1.)
    if new_canvas:
        hist.Draw()
    graph.Draw('p')
    canvas.Update()
    print "made plot"


def make_horizontal_phase_space_plot(data, canvas = None):
    data = copy.deepcopy(data)
    new_canvas = canvas == None
    if new_canvas:
        pz = (11.**2+938**2)**0.5
        canvas = Common.make_root_canvas("x vs px")
        data["x"] = [x for x in data["x"]]
        data["px"] = [px for px in data["px"]]
    hist, graph = Common.make_root_graph("x_vs_px", data["x"], "x [mm]", data["px"], "x'", sort=False)#, xmin=-250, xmax=+250, ymin=-25., ymax=+25.)
    if new_canvas:
        hist.Draw()
    graph.Draw('p')
    graph.SetMarkerStyle(6)
    canvas.cd()
    points = [(data["x"][i], data["px"][i]) for i in range(len(data["z"]))]
    mean, cov = Common.fit_ellipse(points, 100.)
    fit_function = Common.make_root_ellipse_function(mean, cov, [1., 10., 100., 400.], -250., +250., -1., +1.)
    fit_function.Draw("SAME")
    canvas.Update()
    canvas.Print("dynamic_aperture_x.png")
    print "made plot"

def make_vertical_phase_space_plot(data, canvas = None):
    data = copy.deepcopy(data)
    new_canvas = canvas == None
    if new_canvas:
        canvas = Common.make_root_canvas("y vs py")
    hist, graph = Common.make_root_graph("y_vs_py", data["z"], "z [mm]", data["pz"], "z'", sort=False) #, ymin=-1., ymax=+1.e-2, xmin=-250, xmax=+250)
    if new_canvas:
        hist.Draw()
    graph.Draw('p')
    graph.SetMarkerStyle(6)
    canvas.Update()
    print "made plot"


def plot_x_y_projection(step_list):
    canvas = ROOT.TCanvas("x_y_projection", "x_y_projection")
    axes = ROOT.TH2D("x_y_projection_axes", ";x [mm];y [mm]",
                     1000, -5500., 5500.,
                     1000, -5500., 5500.)
    axes.SetStats(False)
    graph = ROOT.TGraph(len(step_list))
    canvas.Draw()
    axes.Draw()
    for i in range(len(step_list["x"])):
        graph.SetPoint(i, step_list["x"][i], step_list["y"][i])
    graph.Draw("l")
    canvas.Update()
    RootObjects.histograms.append(axes)
    RootObjects.canvases.append(canvas)
    RootObjects.graphs.append(graph)
    plot_beam_pipe(3900., 5500., canvas)
    return canvas, axes, graph

def plot_beam_pipe(inner_radius, outer_radius, canvas=None):
    n_steps = 361 # number of azimuthal steps
    n_periods = 12

    if canvas == None:
        canvas = ROOT.TCanvas("beam_pipe", "beam_pipe")
        canvas.Draw()
        axes = ROOT.TH2D("beam_pipe_axes", ";x [mm];y [mm]",
                         1000, -5500., 5500.,
                         1000, -5500., 5500.)
        axes.Draw()
        RootObjects.histograms.append(axes)
        RootObjects.canvases.append(canvas)
    canvas.cd()
    graph_inner = ROOT.TGraph(n_steps)
    graph_outer = ROOT.TGraph(n_steps)
    for i in range(n_steps):
        graph_inner.SetPoint(i,
                             inner_radius*math.sin(i/float(n_steps-1)*2.*math.pi),
                             inner_radius*math.cos(i/float(n_steps-1)*2.*math.pi))
        graph_outer.SetPoint(i,
                             outer_radius*math.sin(i/float(n_steps-1)*2.*math.pi),
                             outer_radius*math.cos(i/float(n_steps-1)*2.*math.pi))
    for i in range(n_periods):
        graph = ROOT.TGraph(2)
        graph.SetPoint(0,
                       inner_radius*math.sin(i/float(n_periods)*2.*math.pi),
                       inner_radius*math.cos(i/float(n_periods)*2.*math.pi))
        graph.SetPoint(1,
                       outer_radius*math.sin(i/float(n_periods)*2.*math.pi),
                       outer_radius*math.cos(i/float(n_periods)*2.*math.pi))
        graph.Draw("l")
        RootObjects.graphs.append(graph)
    graph_inner.Draw("l")
    graph_outer.Draw("l")
    canvas.Update()
    RootObjects.graphs.append(graph_inner)
    RootObjects.graphs.append(graph_outer)

def main():

    directory = "tmp/tune/0/"
    #track_data = parse_track_file(directory)
    #plot_x_y_projection(track_data)
    #for field in "bx", "by", "bz":
    #    make_phi_plot(track_data, field)
    probe_data = parse_probe_file(directory)
    make_horizontal_phase_space_plot(probe_data)
    directory = "tmp/tune/1/"
    probe_data = parse_probe_file(directory)
    make_vertical_phase_space_plot(probe_data)

if __name__ == "__main__":
    main()
    print "press CR to continue"
    raw_input()

